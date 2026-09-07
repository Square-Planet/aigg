/// @file CxDefs.h Low-level header, mostly defining preprocessor macros, and
/// some low-level types (e.g., `cstr_buf_t`, `cx_runtime_error_t`) and functions (e.g.,
/// related to raising/communicating error conditions).
///
/// Most `cx` library components include this header (directly or indirectly).
#ifndef CX_DEFS_H
#define CX_DEFS_H

#include <cstddef> // for size_t

/**
 * @defgroup CX_CONFIG_MACROS C++ preprocessor macros influencing the behavior of `CxDefs.h`, `CxUnicode.h`, and other cx library components
 *
 */

#ifndef CX_NO_STDEXCEPT
   #include <exception> // for std::exception in cx_assert_fail_t and cx_alloc_fail_t
   namespace ct { typedef std::exception cx_exception_t; }
#endif


// Define macros for restricted pointer extensions. By defining a pointer as
// "restricted", we promise the compiler that the pointer is not aliased in
// the current scope. E.g., that in constructs like
//
//   void fn(double *CX_RESTRICT pOut, double const *CX_RESTRICT pIn) {
//      pOut[0] = pIn[0];
//      pOut[1] = pIn[1];
//   },
//
// the two assignments can be scheduled out of order, because no pOut[x] and
// pIn[y] cannot possibly point to the same memory location.
//
// Note:
//   - CX_RESTRICT now univerally goes behind the star ("double *CX_RESTRICT p;").
//     That is, the *pointer* is restricted, not the data it points to.
//     (in Olden Thymes(tm) this was handled differently in different compilers)
#ifdef CX_RESTRICT
   // CX_RESTRICT already pre-defined --> leave it alone
#elif defined (RESTRICT)
   // there is a generic RESTRICT macro from host program --> use that one
   #define CX_RESTRICT RESTRICT
#elif defined(__GNUC__)  || defined(__clang__) // g++ or anything which pretends to be g++ (e.g., clang, intel c++ on linux)
   #define CX_RESTRICT __restrict__
#elif defined(_MSC_VER) // microsoft c++, or anything which pretends to be msvc (e.g., intel c++ on windows)
   #define CX_RESTRICT __restrict
#else
   #define CX_RESTRICT
#endif

#define IR_RP CX_RESTRICT

#ifdef INCLUDE_ABANDONED
// #ifndef RESTRICT
//    #define RESTRICT IR_RP
// #endif
#endif // INCLUDE_ABANDONED

// In some cases we need to inline driver functions as a prerequisite
// for inlining virtual functions. If the driver functions are large,
// compilers might be inclined to ignore our "inline" directive.
// Set up some macros to tell them that we really mean it.
//
// (yes, in the few places where this is used it *does* matter. Yes, I did
// measure. E.g., IR's 3ix shell driver abstractions use this virtual/force-
// inline mechanism (in lieu of a complex template setup) to treat various
// combinations of base/1st/2nd derivative integrals and no-contract/ab- inline-
// contract/c-inline contract with common code; without this, the abstractions
// lead to significant overhead in total integral timings in some cases)
#ifdef __GNUC__ // g++
   #define IR_FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
   #define IR_FORCE_INLINE __forceinline
#else
   #define IR_FORCE_INLINE inline
#endif


#ifdef __GNUC__
   #define IR_NO_INLINE __attribute__ ((noinline))
#elif defined(_MSC_VER)
   #define IR_NO_INLINE __declspec(noinline)
#else
   #define IR_NO_INLINE
#endif


// Some other assorted defines...


// Data prefetch for reading and writing.
#ifdef __GNUC__
   #define IR_PREFETCH_R(p) __builtin_prefetch(p, 0, 1)
   #define IR_PREFETCH_W(p) __builtin_prefetch(p, 1, 1)
#elif defined(_MSC_VER)
   #include <xmmintrin.h>
   #define IR_PREFETCH_R(p) _mm_prefetch((char*)p, 1)
   #define IR_PREFETCH_W(p) _mm_prefetch((char*)p, 1)
#else
   #define IR_PREFETCH_R(p) static_cast<void>(0)
   #define IR_PREFETCH_W(p) static_cast<void>(0)
#endif
// ^-- Note: All explicit prefetches should be subjected to reviews.
// - In places where they are used, they definitely did make a significant
//   positive difference at the time of writing the respective codes and on the
//   respective platforms.
// - However, explicit prefetching is dirty business, and in general at least as
//   likely to make code slower than faster. If it actually helps depends *a lot*
//   on the concrete hardware and on what the compiler concretely translates
//   the surrounding code into
// - So prefetch-optimized code does not age well, and the respective places
//   should be reviewed from time to time.


// Make some macros for C++ standard version... in case we have optional
// support for some features beyond C++98/C++11.
//
// Background:
// - Until ~2021, I ensured that all my code runs in C++98-mode.
//   (rant: Mostly to make it accessible to "High-performance-computing"
//   clusters, which for some reason beyond my imagination tend to be set up
//   with horribly outdated compilers and operating systems. In my experience,
//   seeing 10 year old compilers in action is very common, especially on
//   University or science-lab maintained clusters. I can't even recount how
//   many compilers I have seen which did not even support the instruction sets
//   of the CPUs they are running on...)
//
// - However, in the meantime 10 years have passed since C++11 came around. So
//   most of my codes will now assume C++11 to be present (in particular,
//   wmme/MicroScf, IboView,  csf-fci, ITF2, migrid, permlib, mq1c++).
//   Special-purpose programs may also assume C++17 to be available (e.g., aigg)
//
// - The exception to this are the following codes, which will remain fully
//   compatible with C++98:
//   + The IR integral core
//   + The base-level codes in cx (e.g., CxDefs, CxMemoryStack, CxPodArray,
//     CxIntrusivePtr)
#if __cplusplus >= 201103L
   #define IR_HAVE_CXX11
   #define CX_HAVE_CXX11
#endif
#if __cplusplus >= 201402L
   #define IR_HAVE_CXX14
   #define CX_HAVE_CXX14
#endif
#if __cplusplus >= 201703L
   #define IR_HAVE_CXX17
   #define CX_HAVE_CXX17
#endif
// ^-- The concrete values for the C++ versions were recycled from:
// https://stackoverflow.com/questions/2324658/how-to-determine-the-version-of-the-c-standard-used-by-the-compiler

// Set custom IR_DEBUG/CX_DEBUG macro depending on whether NDEBUG and/or _DEBUG
// are set. _DEBUG is a standard macro in the VC world, whereas standard C
// uses NDEBUG (formally) only to control assertions. My own programs also
// always set _DEBUG in debug builds and NDEBUG in release builds.
#if defined(IR_DEBUG) && !defined(CX_DEBUG)
   #define CX_DEBUG IR_DEBUG
#endif
#if defined(CX_DEBUG) && !defined(IR_DEBUG)
   #define IR_DEBUG CX_DEBUG
#endif
#ifndef IR_DEBUG
   #if (defined(_DEBUG) && !defined(NDEBUG))
      // _DEBUG explicitly set, but no NDEBUG --> definitely debug mode
      #define IR_DEBUG 2
   #elif (!defined(_DEBUG) && defined(NDEBUG))
      // NDEBUG explicitly set, but no _DEBUG --> definitely release mode
      #define IR_DEBUG 0
   #elif (!defined(_DEBUG) && !defined(NDEBUG))
      // neither _DEBUG nor NDEBUG set. Probably included in a non-cgk program
      // being compiled in debug mode. Turn on simple debug code (including assertions)
      #define IR_DEBUG 1
   #elif (defined(_DEBUG) && defined(NDEBUG))
      #warning "CxDefs.h is confused... both _DEBUG and NDEBUG macros are set. Turning off CX/IR debug code."
      #define IR_DEBUG 0
   #endif
   #define CX_DEBUG IR_DEBUG
#endif // IR_DEBUG already set from the outside



#if __cplusplus >= 201103L
   /// `#define`s to `constexpr` in >= C++11 mode and to `inline` in C++98 mode
   #define CX_CONSTEXPR_FN constexpr
   /// `#define`s to `constexpr` in >= C++11 mode and to `const` in C++98 mode
   #define CX_CONSTEXPR_VAR constexpr
   /// Macro CX_CONSTEXPR_IF_NDEBUG `#define`s to `constexpr` if compiled in
   /// NDEBUG mode (that's the standard macro for disabling assertions) and to
   /// `inline` if not.
   ///
   /// The point of this is that constexpr functions in C++11-mode are *very* limited as to what they
   /// can contain, and regular assertions are *NOT* one of those things. So, normally, we either
   /// have to either remove lots of potentially very useful debug checking code (e.g., range
   /// checking code in arrays), or remove constexpr declarations on many functions which otherwise
   /// very well could be constexpr'd. This macro works around the issue by #defining to constexpr
   /// only if assertions are disabled, so the debug code can still be added, as long as it either
   /// consists of `assert()`s only, or is otherwise wrapped in a #ifndef NDEBUG ... #endif block.
   #if NDEBUG
      #define CX_CONSTEXPR_IF_NDEBUG constexpr
   #else
      #define CX_CONSTEXPR_IF_NDEBUG inline
   #endif
   /// `#define`s to `static constexpr` in >= C++11 mode and to `static const` in C++98 mode
   #define CX_STATIC_CONSTEXPR_VAR static CX_CONSTEXPR_VAR
   /// `#define`s to `[[no_unique_address]]` in >= C++11 mode and to nothing in C++98 mode
   /// (note, however, that attribute `[[no_unique_address]]` is only required by the standard
   /// since >= C++20. However, since C++11 we can still set it, but it may be ignored.)
   #define CX_NO_UNIQUE_ADDRESS [[no_unique_address]]
   // ^-- it's only defined as necessarily there in >= C++20 (yes, it took the committee a full 20
   //     years to finally understand that enforcing their idiotic "every object need to have a
   //     unique address (unless ....)" is indeed a problem). However, starting in C++11, we have
   //     the attribute syntax, and if the implementation indeed does not know it, it is technically
   //     just supposed to be ignored.
   #define CX_NOEXCEPT noexcept
   #define CX_NOTHROW_OR_NOEXCEPT noexcept
#else
   // No constexpr functions in C++98. The next best thing is inline functions.
   #define CX_CONSTEXPR_FN inline
   #define CX_CONSTEXPR_VAR const
   #define CX_CONSTEXPR_IF_NDEBUG inline
   #define CX_STATIC_CONSTEXPR_VAR static CX_CONSTEXPR_VAR
   #define CX_NO_UNIQUE_ADDRESS
   #define CX_NOEXCEPT
   #define CX_NOTHROW_OR_NOEXCEPT throw()
   // ^-- there are two different macros because specifying 'throw()' on a function
   //     is generally a terrible idea---the only reason to ever use it is when it
   //     is strictly required. E.g., because a base class declaration uses it,
   //     and we wish to override a function of it. Here: for std::exception.
   //
   //     (^-- Background rant: the exception specification feature of C++98 is
   //     severely mis-designed, and almost(?) never makes any sense to use.
   //     In particular, it can only make both the compiled code and runtime
   //     behavior worse, never better:
   //       That is the case because 'throw()' does not ACTUALLY mean that no
   //     exception will be thrown... the C++ standard requires that *if* an
   //     exception is thrown despite 'throw()', rather than terminating the
   //     program (or just declaring this to result in 'undefined behavior'),
   //     the compiler is instead supposed to generate code to catch that
   //     exception and instead raise a highly unhelpful std::unexpected
   //     exception. So it requires extra runtime-code for handling, rather than
   //     less! And not only that: if the contract actually ends up being
   //     violated, then by implicitly discarding the original exception, the
   //     exception translation also obfuscates the reason for the program
   //     failure. Wait... program failure? Yes, of course the program is going
   //     to abort anyway: no one catches std::unexpected; and even if there was
   //     a catch(), there is no meaningful way to handle it! After all, its very
   //     occurrence establishes that there is a programming(!) error, and
   //     incorrect code cannot be handled in any sensible manner at runtime.
   //       Basically, the only non-negative aspect of a throw specification I can
   //     see is that it communicates *to the programmer* which exceptions they
   //     should consider handling for dealing with failure paths of a certain piece
   //     of code. ...and even in that aspect the throw(...) specification offers
   //     nothing of value over a comment "throws: ...".
   //     Well, just another case in C++ of "...what were they thinking?!")
#endif

// Function attribute promising that a function does not return (e.g., because
// it definitely raises an exception, terminates the program, contains an
// infinite loop, or similar)
#ifndef CX_NORETURN
   #if __cplusplus >= 201103L
      // since C++11 there's a standard attribute for that.
      #define CX_NORETURN [[noreturn]]
   #else
      #if (defined  __GNUC__)
         // tell the compiler that this function will never return
         #define CX_NORETURN __attribute__((noreturn))
      #elif defined(_MSC_VER)  // microsoft c++ (other win32 compilers might emulate this syntax)
         #define CX_NORETURN __declspec(noreturn)
      #else
         #define CX_NORETURN
      #endif
   #endif
#endif

#if (defined  __GNUC__)
   // 'insert software breakpoint here'.
   #define DEBUG_BREAK __asm__("int $0x03");
#elif defined(_MSC_VER)  // microsoft c++ (other win32 compilers might emulate this syntax)
   #define DEBUG_BREAK __asm{ int 3 };
#else
   #define DEBUG_BREAK
#endif



/// Ask the compiler to not warn about a given symbol if it is not actually used
/// (for example, because this is only a partial IR core and the symbol is used
/// in some code not supplied with the partial version. Or some symbol is used
/// only if the code is compiled in IR_DEBUG mode)
///
/// TODO:
/// - unfortunately, on VC++ this style of suppressing warnings (void-cast)
///   currently causes some other warnings about the statement having no effect...
/// - In C++11 mode we could change it to `[](...){}(x);` (as in the va-args macro
///   `CX_SUPPRESS_UNUSED_WARNINGS(...)` below).
///   This does work fine on my up-to-date g++ and clang, but I did not yet
///   check what this does in VC and other compilers, including the ancient
///   g++'s we find so frequently on HPC clusters...
#define IR_SUPPRESS_UNUSED_WARNING(x) (void)x
#define IR_SUPPRESS_UNUSED_WARNING2(x,y) {(void)x; (void)y;}
#define IR_SUPPRESS_UNUSED_WARNING3(x,y,z) {(void)x; (void)y; (void)z;}
#define IR_SUPPRESS_UNUSED_WARNING4(x,y,z,w) {(void)x; (void)y; (void)z; (void)w;}
#if __cplusplus >= 201103L
   #define CX_SUPPRESS_UNUSED_WARNINGS(...) [](...){}(__VA_ARGS__)
   // ^-- pretty, isn't it? Hooray for C++ \o/
   template <typename... Args>
   inline void cx_suppress_unused(Args&&...) noexcept {}
#endif


#define IR_TOKEN_PASTE2(x, y) x ## y
/// Macro-combination for concatenating two C preprocessor tokens, in such a way
/// that they are actually expanded.
#define IR_TOKEN_PASTE(x, y) IR_TOKEN_PASTE2(x, y)
/// Macro for making an anonymous local variable name (it just concats a generic fixed name-prefix with `__LINE__`)
#define IR_UNIQUE_NAME IR_TOKEN_PASTE(_dummy, __LINE__)


#define IR_TOKEN_QUOTE2(x) #x
/// Macro-combination for quoting a number-valued preprocessor macro as a
/// c-string, in such a way that it is actually expanded.
/// (e.g., `IR_TOKEN_QUOTE(__LINE__)` would expand to a `char const *` literal containing
/// the line number of the code which invoked the macro)
#define IR_TOKEN_QUOTE(x) IR_TOKEN_QUOTE2(x)


#ifdef CX_ASSERT
   // If we're asked to, re-define assert(). The main reason for doing this is
   // that the standard version has a habit of hard-crashing MPI executables
   // (e.g., Molpro), and in a way which makes it near-impossible to debug them.
   // Additionally, in some compilers `assert()`s do *not* go away, even if you **do**
   // `#define NDEBUG`!
   //
   // (...and with a custom one, one can put in a debug_break into the assert
   // fail function, which helps a lot if using a visual debugger like VC)
   #undef assert
   #undef assert_rt
   #undef cx_assert_rt
   #undef cx_assert
#endif // CX_ASSERT


#ifdef CX_NO_STDEXCEPT
namespace ct {
   struct cx_exception_t {
   #if __cplusplus >= 201103L
      cx_exception_t() CX_NOEXCEPT {}
      cx_exception_t(cx_exception_t& const) = default;
      cx_exception_t& operator=(cx_exception_t& const) = default;
      cx_exception_t(cx_exception_t&&) = default;
      cx_exception_t& operator=(cx_exception_t&&) = default;
   #endif
      virtual ~cx_exception_t() {};
      virtual char const *what() const CX_NOTHROW_OR_NOEXCEPT = 0;
   };
}
#endif

#ifdef IR_HAVE_CXX11
   // If possible (requires _Pragma form of #pragma, which is only available in
   // >= C++11), set up a pair of macros to suppress warnings about violations
   // of exception specifications in a specific region of code (e.g., a throw()
   // in a function marked as noexcept).
   //
   // - We use this as intentional behavior in a number of situation to
   //   "communicate" hard error conditions which would result in a program
   //   termination in any case---for example an assertion failure, or as an
   //   error-case in a C++11 constexpr function (the latter can only consist of
   //   static_asserts (which cannot involve function paramters!) and return
   //   statements, so `(cond)? value : throw(...)` is the only viable way of
   //   communicating errors related to function arguments).
   // - We do this because is the only standard-compatible way of creating a
   //   termination which more likely than not results in an error message being
   //   shown in one form or another. There are some issues, however: e.g., in
   //   Win32 GUI environments, the VC runtime does not translate the messages
   //   of uncaught exceptions leading to program termination into a human-
   //   readable form, iirc---but neither would any output to `cerr` be shown in
   //   this case, so this is at least not a degression.
   // - Note that raising exceptions in noexcept situations is semi-well-defined
   //   behavior accoding to the standard---it results in the termination of the
   //   program via a call to std::terminate (the implementation-defined part is
   //   whether stack unwinding happens before that).
   // - If it happens during the evaluation of a constexpr expression in a
   //   consteval context, it is a compile error (which is good!)
   #if defined(__GNUC__) && !defined(__clang__)
      #define CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN \
         _Pragma("GCC diagnostic push") \
         _Pragma("GCC diagnostic ignored \"-Wterminate\"")
      #define CX_SUPPRESS_NOEXCEPT_WARNING_END \
         _Pragma("GCC diagnostic pop")
   #elif defined(__clang__)
      #define CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN \
         _Pragma("clang diagnostic push") \
         _Pragma("clang diagnostic ignored \"-Wexceptions\"")
      #define CX_SUPPRESS_NOEXCEPT_WARNING_END \
         _Pragma("clang diagnostic pop")
   #elif defined(_MSC_VER)
      #define CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN \
         _Pragma("warning(push)") \
         _Pragma("warning(disable: 4297)")
      #define CX_SUPPRESS_NOEXCEPT_WARNING_END \
         _Pragma("warning(pop)")
   #else
      #define CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN
      #define CX_SUPPRESS_NOEXCEPT_WARNING_END
   #endif
#else
   #define CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN
   #define CX_SUPPRESS_NOEXCEPT_WARNING_END
#endif
// #define CX_SUPPRESS_NOEXCEPT_WARNING_FOR_EXPR(x) CX_SUPPRESS_NOEXCEPT_WARNING_BEGIN x CX_SUPPRESS_NOEXCEPT_WARNING_END
// ^-- hmpf. Apparently gcc does not allow putting pragmas into the middle of expressions, unlike just about any other compiler framework


#if defined(CX_ASSERT_EXCEPT)
   // If CX_ASSERT_EXCEPT is set, failed assertions will emit an exception,
   // rather than calling CxAssertFail() to abort the program.
   //
   // Advantages of this mode:
   // - The CxDefs header is sufficient to get working assertions. No separate
   //   CxAssertFail() function needs to be implemented (neither our default one
   //   in CxAssertFail.cpp, no any program-specific replacement).
   // - This is easier to track in a debugger.
   // Disadvantages:
   // - This makes it impossible to use assertions inside destructors and some
   //   other code places which are supposed to not throw.
   // - It can confuse the compiler, as it creates more complex potential
   //   control flows (but since this is for debug mode only, with a few
   //   exceptions, it is probably not that big of a deal)

   // An exception class which may be thrown by own assertion macros, in case
   // they are used.
   namespace ct {
      struct cx_assert_fail_t : public cx_exception_t {
         explicit cx_assert_fail_t(char const *msg) : m_msg(msg) {}
         char const *what() const CX_NOTHROW_OR_NOEXCEPT { return m_msg; };
         char const *m_msg;
      };
   }
   // ^-- Note that you might require something to the degree of SetUnhandledExceptionFilter() or
   // std::set_terminate to actually display the message of the unhandled exceptions retained in
   // release builds, expecially if using this in GUI programs or other non-console programs.
   // See post of Logan Capaldo at https://stackoverflow.com/questions/3774316/c-unhandled-exceptions
   // and https://en.cppreference.com/w/cpp/error/set_terminate.
   //
   // Idea is to create a construct like this:
   // ```
   //     #include <exception>
   //     extern "C" { struct FILE; int fflush(FILE *stream); int printf(char const* format, ...); }
   //     static void (*s_pLastTerminateHandler)() = 0;
   //     static void TerminateHandlerWithExceptionPrint()
   //     {
   //        try { throw; }
   //        catch(cx_assert_fail_t const &e) { printf("%s%s\n", "Terminating after runtime-assertion failure:\n", e.what()); fflush(0); }
   //        catch(std::exception const &e) { printf("%s%s\n", "Terminating after unhandled exception: ", e.what()); fflush(0); }
   //        catch(...) { return s_pLastTerminateHandler(); }
   //     }
   //
   //     void CxActivateCustomTerminateHandler() {
   //        s_pLastTerminateHandler = std::set_terminate(TerminateHandlerWithExceptionPrint);
   //     }
   // ```
   // And inside the handler, do whatever needs to be done to display the message
   // of the error condition in such a way that users are likely to actually see it
#endif // CX_ASSERT_EXCEPT

#ifndef cx_assert_rt
   // cx_assert_rt ("_rt" = runtime): these assertions stays active even in release builds.
   #ifdef CX_ASSERT_EXCEPT
      // the exception variant is easier for debugging and maintaining (CxAssertFail not strictly needed).
      // However, it means that one can't put assertions in certain places without being swamped
      // with warnings (e.g., in destructors or below other nowthrow functions) and the throw() may
      // cause the compiler to misbehave regarding optimizations
      #define cx_assert_rt(x) (x)? static_cast<void>(0) : throw ct::cx_assert_fail_t(__FILE__ ":" IR_TOKEN_QUOTE(__LINE__)  ": Assertion failed: '" #x "'"))
   #else
//       void CX_NORETURN CxAssertFail(char const *pExpr, char const *pFile, int iLine) CX_NOEXCEPT;
      void CxAssertFail(char const *pExpr, char const *pFile, int iLine);
      #define cx_assert_rt(x) ((x)? static_cast<void>(0) : CxAssertFail(#x,__FILE__,__LINE__))
      #pragma xref add(CxAssertFail.cpp)
      // ^-- note: if you don't like that particular implementation, you can
      // just make your own and skip CxAssertFail.cpp entirely.
   #endif
#endif

// Already an assert() macro around? (e.g., by including <cassert>)
#ifdef assert
   // use that one for our own assertions
   #define cx_assert(x) assert(x)
#else
   // no assert() macro defined. Make a new cx_assert macro (unless already pre-
   // defined), and if CX_ASSERT preprocessor directive is set, set it as global
   // assert macro, too.
   #ifndef cx_assert
      #if (IR_DEBUG >= 1)
         #define cx_assert(x) cx_assert_rt(x)
      #else
         #define cx_assert(x) static_cast<void>(0)
      #endif
   #endif // cx_assert
   #ifdef CX_ASSERT
      #define assert(x) cx_assert(x)
      #define assert_rt(x) cx_assert_rt(x)
   #endif
#endif

#undef _CxAssert
#undef _CxAssertRt
#undef _IrAssert
#undef _IrAssertRt
#define _CxAssert(x) cx_assert(x)
#define _CxAssertRt(x) cx_assert_rt(x)
#define _IrAssert(x) cx_assert(x)
#define _IrAssertRt(x) cx_assert_rt(x)


#if (defined  __GNUC__ && defined __linux__ && IR_DEBUG >= 2)
   // to allow putting "mcheck_check_all();" (glibc's explicit heap consistency check) anywhere.
   #include <mcheck.h>
#endif



#ifndef IR_CURRENT_FUNCTION_NAME
   #ifdef IR_HAVE_CXX11
      // it's a standard predef since C++11, although with the highly non-useful
      // definition of "an implementation-defined string" (C++11 §8.4.1)
      #define IR_CURRENT_FUNCTION_NAME __func__
   #else
      #if defined(__GNUC__)
         #define IR_CURRENT_FUNCTION_NAME __PRETTY_FUNCTION__
      #elif !defined(__PRETTY_FUNCTION__) && defined(_MSC_VER)
         #define IR_CURRENT_FUNCTION_NAME __FUNCSIG__
      #else
         #define IR_CURRENT_FUNCTION_NAME __FUNCTION__
      #endif
   #endif
#endif


// some global functions for indicating/raising error conditions.
namespace ct {
   /// Holds a copy of a character string of explicitly length, storing it in a
   /// raw buffer. This object is copyable and movable; it can be used to make
   /// copies of 0-terminated c-strings obtained as literals or temporaries.
   /// Primary use by `cx` itself is for representing error/exception messages.
   struct cstr_buf_t {
      struct TakeOverTag{};

      /// Create string buffer storing a heap-hosted copy of 0-terminated c-string `str`
      cstr_buf_t(char const *str) CX_NOEXCEPT;
      /// Create string buffer storing a heap-hosted copy of input c-string [str, str+len].
      /// This function allocates a buffer of length `len + 1` and adds a 0-terminator automatically.
      cstr_buf_t(char const *str, std::size_t len) CX_NOEXCEPT;

      /// Create string buffer, attaching it to an externally created heap-hosted 0-terminated c-string `str`.
      /// `str` must be 0-terminated.
      /// `str` must freeable by std::free(str) (i.e., allocated with std::malloc, std::callc or similar)
      cstr_buf_t(TakeOverTag, char const *str) CX_NOEXCEPT;
      /// Create string buffer, attaching it to an externally created heap-hosted c-string `str` of length `len`.
      /// `str[len-1]` MUST be 0! (i.e., the final string subsection must be 0-terminated. Intermediate 0-values are okay, though.).
      /// `str` must freeable by std::free(str) (i.e., allocated with std::malloc, std::callc or similar)
      cstr_buf_t(TakeOverTag, char const *str, std::size_t len) CX_NOEXCEPT;
      cstr_buf_t(cstr_buf_t const &) CX_NOEXCEPT;
      cstr_buf_t& operator=(cstr_buf_t const &) CX_NOEXCEPT;
   #if __cplusplus >= 201103L
      cstr_buf_t(cstr_buf_t &&) CX_NOEXCEPT;
      cstr_buf_t& operator=(cstr_buf_t &&) CX_NOEXCEPT;
   #endif
      ~cstr_buf_t();
      operator char const *() const CX_NOTHROW_OR_NOEXCEPT { return m_str; }
      char const *data() const CX_NOTHROW_OR_NOEXCEPT { return m_str; }
      char *data() CX_NOTHROW_OR_NOEXCEPT { return m_str; }
      size_t size() const CX_NOTHROW_OR_NOEXCEPT { return m_len; }
      friend void swap(cstr_buf_t &a, cstr_buf_t &b) CX_NOEXCEPT;
   private:
      inline void _Release();
      inline void _Init();
      void _Reset(char const *str, std::size_t len);
      inline void _Swap(cstr_buf_t &);
      inline void _Copy(cstr_buf_t const &);
      char *m_str;
      std::size_t m_len;
   };

   /// Basic runtime error exception class for cx library components
   /// (to allow distinguishing them from host program exceptions, if
   /// needed)
   struct cx_runtime_error_t : public cx_exception_t {
      cx_runtime_error_t(char const *msg) CX_NOEXCEPT : m_msg(msg) {}
      char const *what() const CX_NOTHROW_OR_NOEXCEPT { return m_msg; }
   private:
      cstr_buf_t m_msg;
   };
#ifdef INCLUDE_ABANDONED
//    struct cx_runtime_error_t : public cx_exception_t {
//       cx_runtime_error_t(char const *msg) CX_NOEXCEPT;
//       cx_runtime_error_t(cx_runtime_error_t const &);
//       cx_runtime_error_t& operator=(cx_runtime_error_t const &);
//    #if __cplusplus >= 201103L
//       cx_runtime_error_t(cx_runtime_error_t&&);
//       cx_runtime_error_t& operator=(cx_runtime_error_t&&);
//    #endif
//       ~cx_runtime_error_t() CX_NOTHROW_OR_NOEXCEPT;
//       char const *what() const CX_NOTHROW_OR_NOEXCEPT;
//    private:
//       char const *m_msg;
//    };
#endif // INCLUDE_ABANDONED
   // Why not std::bad_alloc, you ask? Well, std::bad_alloc does not allow specifying
   // a failure message... however, in our case, this condition is more likely to occur
   // upon a failure to allocate a 50 GB array, rather than in a condition where
   // not enough memory is left to do the string formatting to make one.
   struct cx_alloc_fail_t : public cx_runtime_error_t {
      explicit cx_alloc_fail_t(char const *msg) : cx_runtime_error_t(msg) {}
   };
   double _fSizeMb(std::size_t nBytes);
   // Throws cx_alloc_fail_t.
   void _RaiseMemoryAllocError(char const *pFnName, std::size_t ValueTypeSize, std::size_t NewReserveCount);
   void _RaiseMemoryAllocError(char const *pFnName, std::size_t ValueTypeSize, std::size_t NewReserveCount, std::size_t PrevReservedCount);
}


// Next: some trivial functions to avoid having to `#include <algorithm>` just about everywhere.
// (...it drags in a lot of other stuff, and what exactly that is is
// implementation defined. This creates a portability problem in terms of making
// it hard to see if `#include`s of std lib stuff are formally missing.)
namespace ct {
   /// swap two values of a primitive type (pointers, integers, etc.).
   /// WARNING: copies by value! No move semantics.
   template <class FType>
   void _CxSwap1v(FType &A, FType &B){
      FType t = A; A = B; B = t;
   }

   /// Return larger of two values of a primitive type (pointers, integers, etc.).
   /// WARNING: args and return type passed by value!
   template <class FType>
   CX_CONSTEXPR_FN FType _CxMaxQv(FType a, FType b) {
      return (a < b)? b : a; // <-- we often define the < operator only.
   }

   /// Return smaller of two values of a primitive type (pointers, integers, etc.).
   /// WARNING: args and return type passed by value!
   template <class FType>
   CX_CONSTEXPR_FN FType _CxMinQv(FType a, FType b) {
      return (b < a)? b : a; // <-- we often define the < operator only.
   }

   namespace predicates {
   #ifdef INCLUDE_ABANDONED
   //    template<typename _Arg1, typename _Arg2, typename _Result>
   //    struct binary_function_t {
   //       typedef _Arg1 first_argument_type;
   //       typedef _Arg2 second_argument_type;
   //       typedef _Result result_type;
   //    };
   //
   //    template<class T> struct less_t : binary_function_t<T, T, bool> {
   //       bool CX_CONSTEXPR_FN operator()(T const &a, T const &b) const { return a < b; }
   //    };
   //    template<class T> struct equal_t { bool CX_CONSTEXPR_FN operator()(T const &a, T const &b) const { return a == b; } };
   //    template<class T> struct not_equal_t { bool CX_CONSTEXPR_FN operator()(T const &a, T const &b) const { return a != b; } };
   #endif // INCLUDE_ABANDONED
      template<class T> struct less_t { bool CX_CONSTEXPR_FN operator()(T const &a, T const &b) const { return a < b; } };
      template<class T> struct greater_t { bool CX_CONSTEXPR_FN operator()(T const &a, T const &b) const { return b < a; } };
   }
} // namespace ct



/// Given operators a == b and a < b, defines the remaining comparison operators a != b,  a <= b, a > b, a >= b.
/// Use "friend static" in _TemplateOrFriendPrefix to define this inside class scope.
/// Use empty token (>= C++11) for global/namespace scope for concrete types
/// Use "template <...>" in _TemplateOrFriendPrefix to define class template parameters for _Type.
#define CX_ADD_DERIVED_COMPARE_OPS(_TemplateOrFriendPrefix, _Type) \
   _TemplateOrFriendPrefix bool CX_CONSTEXPR_FN operator !=(_Type const &a, _Type const &b) noexcept { return !(a == b); } \
   _TemplateOrFriendPrefix bool CX_CONSTEXPR_FN operator  >(_Type const &a, _Type const &b) { return b < a; } \
   _TemplateOrFriendPrefix bool CX_CONSTEXPR_FN operator <=(_Type const &a, _Type const &b) { return !(b < a); } \
   _TemplateOrFriendPrefix bool CX_CONSTEXPR_FN operator >=(_Type const &a, _Type const &b) { return !(a < b); }


/// Defines a minimal iterator interface to allow use in range-based for loops
/// and generic container printing mechanism. This one is for objects which have
/// a working [] operator to access elements, but not necessarily much else.
/// Examples are raw-pointer based arrays (with non-unit stride, otherwise one
/// can just use the pointers directly).
/// Usage: Add:
///
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,const_iterator,const)
///     CX_IMPLEMENT_ITERATOR_IndexInBracket(TypeOfClassHere,iterator,)
///
/// to class declaration to make both const and non-const iterators which wrap
/// calls to operator [].
#define CX_IMPLEMENT_ITERATOR_IndexInBracket(HostType,IteratorType,ConstOrNothing) \
   struct IteratorType { \
      HostType ConstOrNothing \
         *m_pArrayObj; \
      size_t \
         m_Index; \
      IteratorType &operator ++() { ++m_Index; return *this; } \
      IteratorType &operator --() { --m_Index; return *this; } \
      bool operator != (IteratorType const& other) const { \
         return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      bool operator == (IteratorType const& other) const { \
         return this->m_pArrayObj == other.m_pArrayObj || this->m_Index != other.m_Index; \
      } \
      auto operator *() -> decltype((*m_pArrayObj)[m_Index]) { return (*m_pArrayObj)[m_Index]; } \
   }; \
   IteratorType begin() ConstOrNothing { return IteratorType{this, 0}; } \
   IteratorType end() ConstOrNothing { return IteratorType{this, size()}; }
   // ^--- note: == operator not strictly required.



#endif // CX_DEFS_H
