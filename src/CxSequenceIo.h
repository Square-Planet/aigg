/// @file CxSequenceIo.h Contains implementations of ostream-based (and maybe later
/// `fmt::BasicWriter` based) generic implementations of formatted output operators (<<) for generic
/// sequences and a simple interface to make additional ones for other containers (`std::vector`,
/// `std::set`, `std::map`, etc.).
///
/// # Example (see `sequence_io_example.cpp`):
/// ```
///     #include <vector>
///     #include <list>
///     #include <set>
///     #include <map>
///
///     // tell CxSequenceIo to #include <array> and <tuple> and define <<
///     // operators for them (if SEQUENCE_IO_PREDEFINE_LEVEL is unset, it
///     // #includes only what is necessary for io ops; see below for options)
///     #define SEQUENCE_IO_PREDEFINE_LEVEL 10
///     #include "CxSequenceIo.h"
///
///     // import CxSequenceIo.h's default `operator <<`s for std::tuple
///     // and std::array into the current namespace
///     using ct::sequence_io::operator <<;
///     using ct::sequence_io::open_sep_close_t;
///
///     // define output operators for generic std::vectors, lists, and sets.
///     template<class T, class A>
///     std::ostream &operator << (std::ostream &out, std::vector<T,A> const &seq) {
///        return ct::sequence_io::TSequenceWriter({"[",", ","]"}).WriteList(out, seq);
///        // ^-- first argument is a ct::sequence_io::open_sep_close_t object.
///        // This is equivalent to:
///        //    return ct::sequence_io::TSequenceWriter(
///        //       ct::sequence_io::open_sep_close_t{"[",", ","]"}
///        //    ).WriteList(out, seq);
///        // The WriteList function just returns the ostream object (by reference).
///     }
///
///     template<class T, class A>
///     std::ostream &operator << (std::ostream &out, std::list<T,A> const &seq) {
///        return ct::sequence_io::TSequenceWriter({"[",", ","]"}).WriteList(out, seq);
///     }
///
///     template<class T, class A>
///     std::ostream &operator << (std::ostream &out, std::set<T,A> const &seq) {
///       return ct::sequence_io::TSequenceWriter({"{",", ","}"}).WriteList(out, seq);
///     }
///
///     // If SEQUENCE_IO_PREDEFINE_LEVEL is undefined, or defined to >= 1,
///     // the WriteTuple method is available, which can be used to emit entries
///     // of generic types implementing get/std::tuple_size/std::tuple_element.
///     // We here use it to implement an operator for std::pair, as an example.
///     template<class F, class S>
///     std::ostream &operator << (std::ostream &out, std::pair<F,S> const &tpl) {
///       return ct::sequence_io::TSequenceWriter({"",": ",""}).WriteTuple(out, tpl);
///     }
///
///
///     #include <iostream>
///
///     int main() {
///        std::cout << "-- " __FILE__ << std::endl;
///        std::tuple<int, char const*, char> some_tuple{10, "whee", 'X'};
///        std::array<int,4>  some_array{20,1,30,23};
///        std::set<int>      some_set{20,1,30,23};
///        std::vector<int>   some_vec{20,1,30,23};
///        std::list<int>     some_list{20,1,30,23};
///        std::pair<std::string,int> some_pair{"whooo", 123};
///        std::cout << "Output of `some_tuple`:  " << some_tuple << std::endl;
///        std::cout << "Output of `some_array`:  " << some_array << std::endl;
///        std::cout << "Output of `some_set`:    " << some_set << std::endl;
///        std::cout << "Output of `some_vec`:    " << some_vec << std::endl;
///        std::cout << "Output of `some_list`:   " << some_list << std::endl;
///        std::cout << "Output of `some_pair`:   " << some_pair << std::endl;
///     }
///
///     // Prints:
///     // ```
///     //     -- sequence_io_example.cpp
///     //     Output of `some_tuple`:  (10, "whee", 'X')
///     //     Output of `some_array`:  [20, 1, 30, 23]
///     //     Output of `some_set`:    {1, 20, 23, 30}
///     //     Output of `some_vec`:    [20, 1, 30, 23]
///     //     Output of `some_list`:   [20, 1, 30, 23]
///     //     Output of `some_pair`:   "whooo": 123
///     // ```
/// ```
///
/// # Adding predefined generic operators via `SEQUENCE_IO_PREDEFINE_LEVEL` #define
///
/// By default, `CxSequenceIo.h` only `#include`s support code for making your
/// own operators. However, you may `#define` the preproessor macro
/// `SEQUENCE_IO_PREDEFINE_LEVEL` to an integer level before including this to
/// make `CxSequenceIo.h` (or set it up in your project makefile) also pre-
/// define generic output operators for various standard containers
/// (this mechanism is used to allow control the level of #include pollution).
///
/// Concretely: predefine `SEQUENCE_IO_PREDEFINE_LEVEL` to an integer as follows:
///
/// + `#if (SEQUENCE_IO_PREDEFINE_LEVEL <= 0)`:
///
///    - At level 0, `CxSequenceIo.h` will `#include` only what is necessary to
///      define ostream-based formatting routines for various kinds of
///      containers or other container-like sequences.
///    - ...but it WILL NOT actually declare & define actual generic realizations of
///      those for any standard containers, nor for `std::tuple` or `std::array`
///      (which, strictly speaking, are not containers)
///    - Note: each container is three lines (see below), so there may not
///      actually be a terribly good reason to define them like this (apart
///      from generic std::tuple output or for debugging or development work)
///    - At level 0, it will also not define the generic `TSequenceWriter::WriteTuple`
///      function, which normally can be used to emit generic C++ types implementing
///      the structured binding interface (`std::tuple_element`, `std::tuple_size`, `std::get`).
///    - At this moment there is no guarantee that it will `#include` any concrete
///      std headers apart from `<ostream>`
///
/// +  if `SEQUENCE_IO_PREDEFINE_LEVEL` undefined or `>= 1` (i.e., default behavior)
///
///    - As level 0, but also `#include "CxMetaLight.h"` and define a generic
///      `TSequenceWriter::WriteTuple`, which can be used to emit generic C++
///      types implementing the structured binding interface
///      (`std::tuple_element`, `std::tuple_size`, `std::get`).
///      E.g., an operator for generic tuples can then be implemented as:
///      ```
///          template<class... Args>
///          std::ostream &operator << (std::ostream &out, std::tuple<Args...> const &A) {
///             return TSequenceWriter({"(",", ",")"}).WriteTuple(out, A);
///          }
///      ```
///      This is the same operator whichis predefined at level >= 10.
///
/// + `#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 10)`:
///
///    - At level 10, add `#include`s for `<tuple>` and `<array>` and define
///      generic `std::ostream` output operators for those
///
/// + `#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 20)`:
///
///    - At level 20, add `#include`s for the most widely used standard
///      containers (namely, `<vector>`, `<map>`, `<set>`, and `<list>`), and
///      define generic `operator << (ostream &, ...)` for them. Concretely, the
///      generic operators are defined for `std::vector`, `std::map`,
///      `std::multimap`, `std::set`, `std::multiset`, and `std::list`.
///    - These still sit in namespace `ct::sequence_io`, so they should not by
///      themselves interfere with other code. To actually make them available,
///      you may want to add a
///      ````
///          using ct::sequence_io::operator <<;
///      ````
///      to a suitable local scope.
///
/// + `#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 100)`:
///
///    - At level 30, this file fill attempt to define working `operator <<`s
///      (still in namespace `ct::sequence_io`) for all standard containers,
///      including rather rarely used ones (`std::deque`(*),
///      `std::unordered_multiset`, `std::forward_list`, etc...)
///
/// [*] (cgk's mystery: ...`std::deque` is actually a cool data structure, from
/// an algorithms & data structures point of view. But somehow, at least in
/// my past, it just happened soooo exceedingly rarely that I would actualy
/// want to use one in a practical problem. Is this bias? Am I missing
/// something? Or is it just really not all that common to actually want
/// both FIFO storage and random access at the same time? And, just for the
/// record: I *have* used plenty of `std::multiset`s and `std::multimap`s.
/// In fact, I have probably used more `std::rope`s (back when SGI STL and
/// STLport and stuff still were a thing) than I have used `std::deque`s...)


#ifndef CX_SEQUENCE_IO_H
#define CX_SEQUENCE_IO_H

#include <type_traits> // for remove_const / remove_reference
#include <utility> // for std::pair
#include <exception>
// #include <iostream> // FIXME: remove this
#ifdef INCLUDE_ABANDONED
// #include <algorithm> // for std::max
// #include <functional> // for std::function output wrapper to allow container element overrides.
#endif // INCLUDE_ABANDONED
#include <ostream>
#include <string>
#include <sstream> // for support of alignment/padding for non-basic types
#include <memory> // for std::unique_ptr

#if (defined(CX_UNICODE) && CX_UNICODE >= 1)
   #include "CxUnicode.h" // for UTF-8 support, in particular regarding column width measurement of strings.
   namespace ct { namespace sequence_io {
      using ct::unicode::ColumnWidthQ;
   } }
#else
   namespace ct { namespace sequence_io {
      std::size_t ColumnWidthQ(char const *, std::size_t n) { return n; }
      std::size_t ColumnWidthQ(char const *cstr) { return std::char_traits<char>::length(cstr); }
   } }

#endif

#ifndef SEQUENCE_IO_PREDEFINE_LEVEL
   #define SEQUENCE_IO_PREDEFINE_LEVEL 1
#endif

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 1)
   #include "CxMetaLight.h" // for template index lists (index_list_t, index_list_for_t); used to implement TSequenceWriter::WriteTuple().
   // ^-- note that this does *not* actually #include <tuple> or <array>! (on
   // this machine, both are 25kloc headers with ~650kb source code implementing
   // a huge template mess. I still do not get how they managed that for
   // std::array...). CxMetaLight only brings in <utiltiy> and <type_traits>, which
   // are about 3kloc together and anyway included in just about any standard
   // container header.
#endif
#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 10)
   #include <array>
   #include <tuple>
#endif

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 20)
   // if we're trashing the include space so throughly anyway, I guess it makes
   // little difference whether one brings in all standard containers, too.
   #include <vector>
   #include <map>
   #include <set>
   #include <list>
#endif

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 30)
   #include <deque>
   #include <forward_list>
   #include <unordered_set>
   #include <unordered_map>
#endif



namespace ct {
namespace sequence_io {

struct FFormatError : public std::exception {
   char const *m_msg;
   explicit FFormatError(char const *pmsg) : m_msg(pmsg) {}
   virtual char const *what() const noexcept { return m_msg; }
};



/// Encapsulates the set of io-flags (entry width/precision,ios flags) of a `std::ostream`',
/// and allowed retrieving them from / activating them on a stream.
struct ostream_flags_t {
   int m_prec, m_width;
   std::ios::fmtflags m_flags;
   char m_fill;

   void copy_from(std::ostream *pstr) { m_prec = pstr->precision(); m_width = pstr->width(); m_flags = pstr->flags(); m_fill = pstr->fill(); }
   void enact(std::ostream *pstr) { pstr->precision(m_prec); pstr->width(m_width); pstr->flags(m_flags); pstr->fill(m_fill); }
   explicit ostream_flags_t(std::ostream *pstr) { copy_from(pstr); }
   ostream_flags_t() : m_prec(0), m_width(0), m_flags(std::ios::fmtflags{}), m_fill(' ') {};

   ostream_flags_t &flags(std::ios::fmtflags flags_) { m_flags = flags_; return *this; }
   ostream_flags_t &fill(char f) { m_fill = f; return *this; }
   ostream_flags_t &precision(int p) { m_prec = p; return *this; }
   ostream_flags_t &width(int w) { m_width = w; return *this; }
   ostream_flags_t &setf(std::ios::fmtflags on, std::ios::fmtflags mask) { m_flags = (m_flags & (~mask)) | on; return *this; }
   ostream_flags_t &unsetf(std::ios::fmtflags mask) { m_flags = (m_flags & (~mask)); return *this; }
   char fill() const { return m_fill; }
   // (todo: maybe unset std::ios::unitbuf before emitting arrays?)
};


/// An object which makes a copy of an `std::ostream`'s flags & co (e.g., number
/// formats) on construction, and restores them when it goes out of scope.
struct ostream_flag_guard_t {
   std::ostream *pout;
   ostream_flags_t flags;
   explicit ostream_flag_guard_t(std::ostream *pout_) : pout(pout_), flags(pout) {}
   ~ostream_flag_guard_t() { flags.enact(pout); }
};


// // Defines a minimal iterator interface to allow use in range-based for loops
// // and generic container printing mechanism. This one is for objects which have
// // a working [] operator to access elements, but not necessarily much else.
// // Examples are raw-pointer based arrays (with non-unit stride, otherwise one
// // can just use the pointers directly).
// //
// // Example to put in class:
// //    typedef TMinimalIterator_Bracket<FThisType const> const_iterator;
// //    const_iterator begin() const { return const_iterator{this, 0}; }
// //    const_iterator end() const { return const_iterator{this, size()}; }
// template<class FArrayLike, class FElement = typename FArrayLike::value_type>
// struct TMinimalIterator_Bracket {
//    TMinimalIterator_Bracket &operator ++() { ++m_Index; return *this; }
//    bool operator != (TMinimalIterator_Bracket const& other) const {  // <-- only != is strictly required (== optional).
//       return this->m_pArrayObj != other.m_pArrayObj || this->m_Index != other.m_Index;
//    }
//    FElement operator *() { return (*m_pArrayObj)[m_Index]; }
// protected:
//    FArrayLike
//       *m_pArrayObj;
//    size_t
//       m_Index;
// };

// ^-- UPDATE: that doesn't work because the class is not complete at that point.
//     Doing this requires a macro. I made one in CxIterTools.h. See:
//     CX_IMPLEMENT_ITERATOR_IndexInBracket

/// Wrapper for giving raw-pointer based arrays an interface which is sufficient
/// for WriteList (implements what is needed for the range-based for).
template<class FElement>
struct TArrayViewStrided {
   FElement const
      *m_pData;
   size_t
      // number of array elements
      m_Size;
   ptrdiff_t
      // in-memory distance between two adjacent elements we should consider
      // as adjacent. That is, actual element #i is at m_pData[i * m_Stride]
      m_Stride = 1;
   size_t size() const { return m_Size; }
   FElement const &operator[] (size_t i) { return m_pData[m_Stride * i]; }
//    typedef TMinimalIterator_Bracket<TArrayViewStrided const> const_iterator;
//    const_iterator begin() const { return const_iterator{this, 0}; }
//    const_iterator end() const { return const_iterator{this, size()}; }
#ifdef CX_IMPLEMENT_ITERATOR_IndexInBracket
   CX_IMPLEMENT_ITERATOR_IndexInBracket(TArrayViewStrided,const_iterator,const)
#endif
};


/// Emit a string (given in terms of iterator sequence), quoting it and (by
/// default) performing some more common escape character replacements (e.g.,
/// actual newlines by the character sequence '\n' and the QuoteLeft/QuoteRight
/// characters by backslash-ed versions thereof)
template <class FStringIt>
void WriteQuoted(std::ostream &out, FStringIt const &first, FStringIt const &last, char QuoteLeft, char QuoteRight, bool EscapeQuotesOnly=false) {
   out << QuoteLeft;
   auto const &try_emit_escape = [&](char c, char from, char const *to, char to2) -> bool {
      if (c == from) {
         out << to;
         if (to2 != 0) out << to2;
         return true;
      }
      return false;
   };
   for (FStringIt it = first; it != last; ++ it) {
      auto c = *it;
      if (!EscapeQuotesOnly) {
         if (try_emit_escape(c, '\n', "\\n", 0)) continue;
         if (try_emit_escape(c, '\t', "\\t", 0)) continue;
         if (try_emit_escape(c, '\0', "\\0", 0)) continue;
         if (try_emit_escape(c, '\\', "\\\\", 0)) continue;
      }
      if (try_emit_escape(c, QuoteLeft, "\\", QuoteLeft)) continue;
      if (try_emit_escape(c, QuoteRight, "\\", QuoteRight)) continue;
      // still here? That means it's not one of the to-escape ones above.
      out << c;
      // btw: this is not supposed to be safe at all. Just to deal with some
      // more common escapes to help with debugging. DIY escapes are highly NOT
      // recommended in situations where you might encounter malicious data...
   }
   out << QuoteRight;
}

/// Emit the template element verbatim using their own `operator <<`; this
/// includes `std::string`, `char`, `char const *`, which this function
/// **does not** quote or escape.
///
/// Note: this function is declared `inline` instead of `static` only to suppress
/// annoying compiler warnings.
template<class FElement>
inline void DefaultWriteElement_Unquoted(std::ostream &out, size_t iArg, FElement const &Arg) {
   out << Arg; (void)iArg;
}


/// Emit `std::string`, `char`, and `char const *` using the `WriteQuoted`
/// function, and all other types as verbatim using their own `operator <<
/// (ostream&, ...)` functions.
///
/// Note: this function is declared `inline` instead of `static` only to suppress
/// annoying compiler warnings.
template<class FElement>
inline void DefaultWriteElement(std::ostream &out, size_t iArg, FElement const &Arg) {
   out << Arg; (void)iArg;
}

template<>
inline void DefaultWriteElement<std::string>(std::ostream &out, size_t iArg, std::string const &Arg) {
   WriteQuoted(out, Arg.begin(), Arg.end(), '\"', '\"'); (void)iArg;
}

template<>
inline void DefaultWriteElement<char const *>(std::ostream &out, size_t iArg, char const * const &pArg) {
   // assume 0-terminated string. count how long it is.
//    size_t n = 0; for (char const *p = pArg; *p;) { n += 1; p += 1; }
   size_t n = std::char_traits<char>::length(pArg);
   WriteQuoted(out, pArg, pArg+n, '\"', '\"'); (void)iArg;
}

template<>
inline void DefaultWriteElement<char>(std::ostream &out, size_t iArg, char const &cArg) {
   // assume 0-terminated string. count how long it is.
   WriteQuoted(out, &cArg, &cArg+1u, '\'', '\''); (void)iArg;
}


/// Meta-type wrapping `DefaultWriteElement()` as the `.write(ostream&, size_t iArg, T const &Arg)` static member function.
template<class FElement>
struct TDefaultElementWriter {
   static std::ostream &write(std::ostream &out, size_t iArg, FElement const &Arg) { DefaultWriteElement(out, iArg, Arg); return out; }
};


/// Meta-type wrapping `DefaultWriteElement_Unquoted()` as the `.write(ostream&, size_t iArg, T const &Arg)` static member function.
template<class FElement>
struct TDefaultElementWriter_Unquoted {
   static std::ostream &write(std::ostream &out, size_t iArg, FElement const &Arg) { DefaultWriteElement_Unquoted(out, iArg, Arg); return out; }
};


// // suppress some unused-warnings. I think it is not particularly concerning if
// // certain of these functions are not, in fact, used in every single compilation unit.
// // update: took this out again and made DefaultWriteElement non-static instead. For i/o
// // stuff it doesn't much matter, I guess, although it does mean that we bring in those
// // functions even if never used anywhere.
// inline void please_disregard_this() {
//    (void)DefaultWriteElement<char const *>;
//    (void)DefaultWriteElement<std::string>;
// }


// DefaultWriteElement is specialized for the base types (but takes
// const refs as arguments to the elements). To get the right
// function we therefore need to strip off the reference and const-ness.

/// C++11 equivalent of C++14's `std::remove_reference_t`
template<class FDeclType> using without_ref_t = typename std::remove_reference<FDeclType>::type;
/// C++11 equivalent of C++20's `std::remove_cvref`
template<class FDeclType> using without_cvref_t = typename std::remove_cv<typename std::remove_reference<FDeclType>::type>::type;



/// Helper class to identify whether a sequence entry is some sort of elementary type
/// for which using iostream.width/align directly is likely to do a somewhat sensible
/// thing.
template<class T> struct _IsElementaryStreamType { enum { value = false }; };
template<> struct _IsElementaryStreamType<bool>               { enum { value = true }; };
template<> struct _IsElementaryStreamType<short>              { enum { value = true }; };
template<> struct _IsElementaryStreamType<unsigned short>     { enum { value = true }; };
template<> struct _IsElementaryStreamType<int>                { enum { value = true }; };
template<> struct _IsElementaryStreamType<unsigned int>       { enum { value = true }; };
template<> struct _IsElementaryStreamType<long>               { enum { value = true }; };
template<> struct _IsElementaryStreamType<unsigned long>      { enum { value = true }; };
template<> struct _IsElementaryStreamType<long long>          { enum { value = true }; };
template<> struct _IsElementaryStreamType<unsigned long long> { enum { value = true }; };
template<> struct _IsElementaryStreamType<float>              { enum { value = true }; };
template<> struct _IsElementaryStreamType<double>             { enum { value = true }; };
template<> struct _IsElementaryStreamType<long double>        { enum { value = true }; };
template<> struct _IsElementaryStreamType<void*>              { enum { value = true }; };
template<> struct _IsElementaryStreamType<char const*>        { enum { value = true }; };
template<> struct _IsElementaryStreamType<std::string>        { enum { value = true }; };

/// Encapsulation for field formatting (alignment, width, precision, etc) associated with elements of a type
struct FFieldFormat {
   int width; //< -1 or element width (using padding)
   int prec; //< entry 'precision' (number of post-decimal digits for floats for 'f' and 'e' types).
   char type; //< one of: 0 (auto), f (fixed), e (scientific), g (general), x (hex)
   char align; //< one of: 0 (auto), < (left), > (right), ^ (center)   (note: f/e/g use right/internal-align by default)
   char const *fill; //< 0 (for auto/space) or char sequence to use for padding.
   FFieldFormat(int width_ = -1, int prec_ = -1, char type_ = 0, char align_ = 0, char const *fill_ = 0)
   : width(width_), prec(prec_), type(type_), align(align_), fill(fill_) {}
};


/// Abstraction class for factoring out template-type-specific code for emitting a single element of a sequence.
struct _FEmitElementProxy {
   virtual void operator ()(std::ostream &out) const = 0;
};

template<class T, class FWriteElementFn>
struct _TEmitElementProxy : public _FEmitElementProxy {
   size_t i;
   T const *ei;
   FWriteElementFn const *pEmitElement;
   _TEmitElementProxy(size_t i_, T const &ei_, FWriteElementFn const &pEmitElement_) : i(i_), ei(&ei_), pEmitElement(&pEmitElement_) {}

   virtual void operator ()(std::ostream &out) const /* override */ {
      (*pEmitElement)(out, i, *ei);
   }
};

template<class T, class FWriteElementFn>
auto _BindElementEmitter(size_t i, T const &ei, FWriteElementFn const &pEmitElement) -> _TEmitElementProxy<T, FWriteElementFn> { return {i, ei, pEmitElement}; }


struct FFieldFormatContext {
#ifdef INCLUDE_ABANDONED
//    template<class OuterT>
//    explicit FFieldFormatContext(FFieldFormat const &ffInner, FFieldFormat const &ffOuter, std::ostream &out)
//       : m_InnerFormat(ffInner), m_OuterFormat(ffOuter) { Init(out, _IsElementaryStreamType<without_cvref_t<OuterT> >::value); }
//    template<class OuterT>
#endif // INCLUDE_ABANDONED
   explicit FFieldFormatContext(FFieldFormat const &ffInner, FFieldFormat const &ffOuter, std::ostream &out)
      : m_InnerFormat(ffInner), m_OuterFormat(ffOuter) { _Init(out, false); }


   template<class OuterT, class FWriteElementFn>
   std::ostream &Emit(std::ostream &out, size_t i, OuterT const &ei, FWriteElementFn pEmitElement) {
      return _EmitElement(out, _BindElementEmitter(i, ei, pEmitElement));
   }

   std::ostream &_EmitElement(std::ostream &out, _FEmitElementProxy const &_pEmitBoundElement) {
      using std::get;
      if (!m_ProxyPad) {
         ostream_flag_guard_t _flag_restore_toplevel(&out);
         // no special operations required --- enact prepared flags and output directly to target stream `out`.
         m_EntryFlags.enact(&out);
         _pEmitBoundElement(out);
      } else {
         // need to perform some sort of non-trivial outer-element formatting (e.g., aligning entire vectors,
         // not just the individual vector elements). With iostreams that can only be done by
         // rendering the data into a temporary buffer (...and it has to be something to the degree
         // of a std::stringstream or equivalent) and then from there out to the target buffer.

         if (0) {
            ostream_flag_guard_t _flag_restore_toplevel(&out);
            // reset previous buffer content
            m_sstr->str({});
            // ^- fun fact: that may or may not reallocate its previous buffer memory---it is implementation
            //    defined (apparently the VC RTL reallocates, while the g++ RTL does not). As far as I can
            //    see, there is no standard-specified way to reset a stringstream *without* reallocation of its buffer.
            //    So this fully conforms to the primary iostream design priciple: "what were they thinking?!"
            m_EntryFlags.enact(&*m_sstr);
            _pEmitBoundElement(*m_sstr);
            out.unsetf(std::ios::adjustfield);
            out.width(0);
            out << m_sstr->rdbuf();
         } else {
            m_sstr->clear(std::ios::goodbit); // clear error state flags (incl. possible eof from last run)
            m_sstr->seekp(0, std::ios::beg); // just reset 'put' position to start
            // ^-- ... note that this does not neither reset the buffer length, nor clear out the
            //     previous data in the stream. So m_sstr.str() will, in general, produce garbage after
            //     just doing this.
            m_EntryFlags.enact(&*m_sstr);
            _pEmitBoundElement(*m_sstr);
            size_t n = m_sstr->tellp();
            size_t
               padleft = 0,
               padright = 0,
               padtotal = 0;
            if (m_OuterFormat.width > ptrdiff_t(n)) {
               padtotal = size_t(m_OuterFormat.width) - n;
            }
            if (m_OuterFormat.align == '<')
               padright = padtotal;
            else if (m_OuterFormat.align == '>')
               padright = 0;
            else if (m_OuterFormat.align == '^')
               padright = padtotal / 2;
            padleft = padtotal - padright;
            for (size_t i = 0; i != padleft / m_FillWidth; ++ i)
               out.write(m_OuterFormat.fill, m_FillSize);

            m_sstr->seekg(0, std::ios::beg);
//             std::array<char,256> buf1;
            size_t const buf1_size = 256;
            char buf1[buf1_size];
            while (n != 0 && out.good() && m_sstr->good()) {
               size_t n1 = (n < buf1_size)? n : buf1_size;
               m_sstr->read(&buf1[0], n1);
               if (m_sstr->fail()) throw FFormatError("FFieldFormatContext::Emit(): m_sstr->read() failed unexpectedly.");
               out.write(&buf1[0], n1);
               n -= n1;
            }
            for (size_t i = 0; i != padright / m_FillWidth; ++ i)
               out.write(m_OuterFormat.fill, m_FillSize);
            // FIXME: this is not measuring column width of inner string correctly!
         }
      }
      return out;
   }
#ifdef INCLUDE_ABANDONED
//    template<class OuterT, class FWriteElementFn>
//    std::ostream &Emit(std::ostream &out, size_t i, OuterT const &ei, FWriteElementFn pEmitElement) {
//       using std::get;
//       if (!m_ProxyPad) {
//          ostream_flag_guard_t _flag_restore_toplevel(&out);
//          // no special operations required --- enact prepared flags and output directly to target stream `out`.
//          m_EntryFlags.enact(&out);
//          pEmitElement(out, i, ei);
//       } else {
//          // need to perform some sort of non-trivial outer-element formatting (e.g., aligning entire vectors,
//          // not just the individual vector elements). With iostreams that can only be done by
//          // rendering the data into a temporary buffer (...and it has to be something to the degree
//          // of a std::stringstream or equivalent) and then from there out to the target buffer.
//
//          if (0) {
//             ostream_flag_guard_t _flag_restore_toplevel(&out);
//             // reset previous buffer content
//             m_sstr->str({});
//             // ^- fun fact: that may or may not reallocate its previous buffer memory---it is implementation
//             //    defined (apparently the VC RTL reallocates, while the g++ RTL does not). As far as I can
//             //    see, there is no standard-specified way to reset a stringstream *without* reallocation of its buffer.
//             //    So this fully conforms to the primary iostream design priciple: "what were they thinking?!"
//             m_EntryFlags.enact(&*m_sstr);
//             pEmitElement(*m_sstr, i, ei);
//             out.unsetf(std::ios::adjustfield);
//             out.width(0);
//             out << m_sstr->rdbuf();
//          } else {
//             m_sstr->clear(std::ios::goodbit); // clear error state flags (incl. possible eof from last run)
//             m_sstr->seekp(0, std::ios::beg); // just reset 'put' position to start
//             // ^-- ... note that this does not neither reset the buffer length, nor clear out the
//             //     previous data in the stream. So m_sstr.str() will, in general, produce garbage after
//             //     just doing this.
//             m_EntryFlags.enact(&*m_sstr);
//             pEmitElement(*m_sstr, i, ei);
//             size_t n = m_sstr->tellp();
//             size_t
//                padleft = 0,
//                padright = 0,
//                padtotal = 0;
//             if (m_OuterFormat.width > ptrdiff_t(n)) {
//                padtotal = size_t(m_OuterFormat.width) - n;
//             }
//             if (m_OuterFormat.align == '<')
//                padright = padtotal;
//             else if (m_OuterFormat.align == '>')
//                padright = 0;
//             else if (m_OuterFormat.align == '^')
//                padright = padtotal / 2;
//             padleft = padtotal - padright;
//             for (size_t i = 0; i != padleft / m_FillWidth; ++ i)
//                out.write(m_OuterFormat.fill, m_FillSize);
//
//             m_sstr->seekg(0, std::ios::beg);
//             std::array<char,256> buf1;
//             while (n != 0 && out.good() && m_sstr->good()) {
//                size_t n1 = (n < buf1.size())? n : buf1.size();
//                m_sstr->read(&buf1[0], n1);
//                if (m_sstr->fail()) throw FFormatError("FFieldFormatContext::Emit(): m_sstr->read() failed unexpectedly.");
//                out.write(&buf1[0], n1);
//                n -= n1;
//             }
//             for (size_t i = 0; i != padright / m_FillWidth; ++ i)
//                out.write(m_OuterFormat.fill, m_FillSize);
//             // FIXME: this is not measuring column width of inner string correctly!
//          }
//       }
//       return out;
//    }
#endif // INCLUDE_ABANDONED

protected:
   // format description for inner-most elements (e.g., one might pass a "{:12.6f}" to a sequence of
   // numerical vectors here). Must be representable by standard iostream formats.
   FFieldFormat    m_InnerFormat;
   FFieldFormat    m_OuterFormat;  // format description for elements of the current-level sequence itself
   ostream_flags_t m_EntryFlags;   // settings to enact on target iostream before emitting an output element.
   bool            m_Numeric;      // set if entries have a numeric format specifier set (e.g., 'x', 'f', 'e')
   std::size_t     m_FillWidth;    // column width of m_OuterFormat.fill fill-character sequence
   std::size_t     m_FillSize;     // strlen of m_OuterFormat.fill fill-character sequence
   bool            m_ProxyPad;     // indicates that the given combination of format flags cannot be natively performed by iostream

   using buffer_t = std::stringstream;
   std::unique_ptr<buffer_t> m_sstr;

   void _Init(std::ostream &out, bool _OuterElementTypeIsBasic) {
      using std::get;
//       m_EntryFlags = {0};  // (has a default c'tor now)
      m_Numeric = false;
      m_FillWidth = 1;
      m_FillSize = 1;
      m_ProxyPad = false;
//       if (!_OuterElementTypeIsBasic)
//          m_ProxyPad = true;

      int w = m_InnerFormat.width;
      int p = m_InnerFormat.prec;
      m_EntryFlags.width((w >= 0)? w : int(out.width()));
      m_EntryFlags.precision((p >= 0)? p : int(out.precision()));

      char const *fi = m_InnerFormat.fill;
      if (fi == 0) {
         m_EntryFlags.fill(out.fill());
      } else {
         if (1 != std::char_traits<char>::length(fi)) throw FFormatError("FFieldFormatContext::_Init(): non-unit-length "
            "fill characters not allowed for inner format (not iostream-representable).");
         m_EntryFlags.fill(fi[0]);
      }

      m_EntryFlags.flags(out.flags());
      switch (m_InnerFormat.type) {
         // switch mode to "general" (typically omits trailing zeros, unlike 'f')
         case 'g': { m_EntryFlags.unsetf(std::ios::floatfield); m_Numeric = true; break; };
         case 'f': { m_EntryFlags.setf(std::ios::fixed, std::ios::floatfield); m_Numeric = true; break;}
         case 'e': { m_EntryFlags.setf(std::ios::scientific, std::ios::floatfield); m_Numeric = true; break; }
         case 'x': { m_EntryFlags.setf(std::ios::hex, std::ios::basefield); m_Numeric = true; break; }
         case 'd': { m_EntryFlags.setf(std::ios::dec, std::ios::basefield); m_Numeric = true; break; }
         case  0 : { break; }
         case 's': { break; }
         default: throw FFormatError("FFieldFormatContext::_Init(): unrecognied inner type format.");
      }
      char a = m_InnerFormat.align;
      if (a == 0 && m_Numeric)
         // for numeric format types, right-align by default.
         a = '>';
      if (w >= 0 && a != 0) {
         // for numeric types with '0' filling, replace right-align by 'internal' align
         // (fill padding beween sign and rest of number, not before sign)
         if (m_Numeric && m_EntryFlags.fill() == '0' && a == '>')
            a = '=';
         switch (a) {
            case '>': m_EntryFlags.setf(std::ios::right, std::ios::adjustfield); break;
            case '=': m_EntryFlags.setf(std::ios::internal, std::ios::adjustfield); break;
            case '<': m_EntryFlags.setf(std::ios::left, std::ios::adjustfield); break;
            case '^': throw FFormatError("FFieldFormatContext::_Init(): center-align not allowed for inner format (not iostream-representable)");
//          case '^': m_EntryFlags.unsetf(std::ios::adjustfield); m_ProxyPad = true; break;  // center align: not supported by iostream.
            default: throw FFormatError("FFieldFormatContext::_Init(): unrecognied alignment format.");
         }
      }
      char const *fo = m_OuterFormat.fill;
      if (fo != 0 && m_OuterFormat.width > 0) {
         m_FillSize = std::char_traits<char>::length(fo);
         m_FillWidth = ColumnWidthQ(fo);
         if (m_FillWidth == 0) throw FFormatError("FFieldFormatContext::_Init(): cannot pad output with a zero-width fill string.");
         m_ProxyPad = true;
      }
      if (m_OuterFormat.type != 0 || m_OuterFormat.prec != -1)
         throw FFormatError("FFieldFormatContext::_Init(): non-default type/prec formats not allowed for outer type");

      if (m_ProxyPad) {
         // some combination of options cannot be processed by iostream formats -- need manual hacking to work-around it.
         m_sstr = std::unique_ptr<buffer_t>(new buffer_t);
      }
      [](...){}(_OuterElementTypeIsBasic); // suppress unused warning.
   }

};


/// Main interface for defining output operators. This class encapsulates and
/// contains generic members (most of which are optional) for controlling the
/// formatted output of sequences.
struct TSequenceWriter {
   struct FSeparators {
      char const *open;
      char const *sep;
      char const *close;
      template<std::size_t i> friend constexpr char const *get(FSeparators const &osc) {
         return (i == 0)? osc.open : ((i == 1)? osc.sep : osc.close);
      }
   };
//    typedef std::array<char const *, 3> FSeparators;
   /// Separators between the actual direct elements of the container
   FSeparators m_SequenceSeptors = {"[",",","]"}; // LeftMidRight
   /// Separators for parts of an individual element; not normally used,
   /// except for cases like associative containers (`std::map` and co),
   /// in which we'd separate `(key : value)` pairs instead of just printing
   /// them as 2-tuples.
   FSeparators m_ElementSeptors = {"",": ",""};
   FFieldFormat m_FieldFormatInner = {-1, -1, 0, 0, 0};
   FFieldFormat m_FieldFormatOuter = {-1, -1, 0, 0, 0};

   TSequenceWriter() {}
   explicit TSequenceWriter(FSeparators LeftMidRight_) : m_SequenceSeptors(LeftMidRight_) {}
   explicit TSequenceWriter(FSeparators LeftMidRight_, FFieldFormat FieldFormatInner_, FFieldFormat FieldFormatOuter_ = {-1, -1, 0, 0, 0}) : m_SequenceSeptors(LeftMidRight_), m_FieldFormatInner(FieldFormatInner_), m_FieldFormatOuter(FieldFormatOuter_) {}
   explicit TSequenceWriter(FSeparators Sequence_LeftMidRight_, FSeparators Element_LeftMidRight_) : m_SequenceSeptors(Sequence_LeftMidRight_), m_ElementSeptors(Element_LeftMidRight_) {}

#ifdef INCLUDE_ABANDONED
//    // arguments: target stream, index of element in sequence, element
//    template<class FElement>
//    using TWriteElementFn = std::function<void (std::ostream &, size_t, FElement const &)>;

//    // change to iterator-based. Would allow sets and lists, but no longer the packed
//    // index.... unless I add some sort of fake iterator class...
//    template<class FArrayLike, class FWriteElementFn>
//    std::ostream &WriteList(std::ostream &out, FArrayLike const &A, FWriteElementFn &&pEmitElement) const {
//       using std::get;
//       ostream_flag_guard_t _flag_restore_outer(&out);
//       out << get<0>(m_SequenceSeptors);
//
//       ostream_flags_t entry_flags(&out);
//       {
//          ostream_flag_guard_t _flag_restore_inner(&out);
//          int w = get<0>(m_FieldFormat);
//          if (w >= 0)
//             out.width(w);
//          if (get<1>(m_FieldFormat) >= 0)
//             out.precision(get<1>(m_FieldFormat));
//          char c = get<2>(m_FieldFormat);
//          if (c == 'g')
//             out.unsetf(std::ios::floatfield); // switch mode to "general" (typically omits trailing zeros, unlike 'f')
//          else if (c == 'f')
//             out.setf(std::ios::fixed, std::ios::floatfield);
//          else if (c == 'e')
//             out.setf(std::ios::scientific, std::ios::floatfield);
//          char a = get<3>(m_FieldFormat);
//          if (a == 0 && (c == 'g' || c == 'f' || c == 'e'))
//             a = '>'; // for numeric format types use right-align by default
//          if (w >= 0 && a != 0) {
//             if (a == '>')
//                out.setf(std::ios::right, std::ios::adjustfield);
//             else if (a == '<')
//                out.setf(std::ios::left, std::ios::adjustfield);
//          }
//          entry_flags.copy_from(&out);
//       }
//
//       size_t i = 0;
//       for (auto const &Ai : A) {
//          if (i != 0)
//             out << get<1>(m_SequenceSeptors);
//          {
//             entry_flags.enact(&out);
//             pEmitElement(out, i, Ai);
// //          out << Ai;
//             out.width(0);
//          }
//          i += 1;
//       }
//       out << get<2>(m_SequenceSeptors);
//       return out;
//    }
#endif // INCLUDE_ABANDONED


   // change to iterator-based. Would allow sets and lists, but no longer the packed
   // index.... unless I add some sort of fake iterator class...
   template<class FArrayLike, class FWriteElementFn>
   std::ostream &WriteList(std::ostream &out, FArrayLike const &A, FWriteElementFn &&pEmitElement) const {
      using std::get;
      ostream_flag_guard_t _flag_restore_outer(&out);
      out.unsetf(std::ios::unitbuf); // <- disable flush after each io operation.
      auto _EmitRaw = [&](std::ostream &out, char const *p) -> std::ostream & {
         out.write(p, std::char_traits<char>::length(p));
         return out;
      };
      _EmitRaw(out, get<0>(m_SequenceSeptors));

      FFieldFormatContext ffctx(m_FieldFormatInner, m_FieldFormatOuter, out);

      size_t i = 0;
      for (auto const &Ai : A) {
         if (i != 0)
            _EmitRaw(out, get<1>(m_SequenceSeptors));
         ffctx.Emit(out, i, Ai, pEmitElement);
         i += 1;
      }
      _EmitRaw(out, get<2>(m_SequenceSeptors));
      return out;
   }



   template<class FArrayLike>
   std::ostream &WriteList(std::ostream &out, FArrayLike const &A) const {
      return WriteList(out, A, &DefaultWriteElement<without_cvref_t<decltype(*A.begin())> >);
   }

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 1)
   template<template<class E> class TElementWriter, cx::index_t... Is, class FTupleLike>
   void WriteTupleImpl(std::ostream &out, FTupleLike const &A, cx::index_list_t<Is...>) const {
      using std::get;
      // that's apparently the recommended way before C++17 fold expressions
      // to perform statements on parameter packs...  The C-array
      // is there to ask the compiler to execute the statements in original sequence.
      // To make it an int array, we use the comma operator (..., 0).
      int dummy[] = {
         (
            TElementWriter<cx::without_cvref_t<decltype(get<Is>(A))> >::write(
               (out << ((Is == 0)? "" : get<1>(m_SequenceSeptors))),
               Is,
               get<Is>(A)
            ),
            0
         )...
      };
      [](...){}(dummy); // suppress unused warning
   }

#ifdef INCLUDE_ABANDONED
//    template<cx::index_t... Is, class FTupleLike>
//    void WriteTupleImpl(std::ostream &out, FTupleLike const &A, cx::index_list_t<Is...>) const {
//       using std::get;
//       int dummy[] = {
//          // that's apparently the recommended way before C++17 fold expressions
//          // to perform statements on parameter packs...  The C-array
//          // is there to ask the compiler to execute the statements in original sequence.
//          // To make it an int array, we use the comma operator (..., 0).
//          (
//             (out << ((Is == 0)? "" : get<1>(m_SequenceSeptors)) << get<Is>(A)),
//             0
//          )...
//       };
//       (void)dummy; // suppress unused warning
//    }
//    template<class... Args>
//    std::ostream &WriteTuple(std::ostream &out, std::tuple<Args...> const &A) const {
//       ostream_flag_guard_t _flag_restore(&out);
//       using std::get;
//       out << get<0>(m_SequenceSeptors);
//       WriteTupleImpl(out, A, cx::index_list_for_t<Args...>());
//       out << get<2>(m_SequenceSeptors);
//       return out;
//    }
#endif // INCLUDE_ABANDONED
   template<class FTupleLike, template<class> class TElementWriter = TDefaultElementWriter>
   std::ostream &WriteTuple(std::ostream &out, FTupleLike const &A) const {
      ostream_flag_guard_t _flag_restore(&out);
      using std::get;
      out << get<0>(m_SequenceSeptors);
      size_t constexpr N = std::tuple_size<cx::without_cvref_t<FTupleLike> >::value;
      WriteTupleImpl<TElementWriter>(out, A, cx::index_list_of_size_t<N>());
      out << get<2>(m_SequenceSeptors);
      return out;
   }
#endif // #if (SEQUENCE_IO_PREDEFINE_LEVEL >= 1)

};


   template<class K, class V>
   void DefaultWriteMapElement(std::ostream &out, size_t i, std::pair<K,V> const &ei) {
      // map's value_type is std::pair<Key, Val>. Print them separated by `:`, similar to a python dict literal.
   //    out << ei.first << ": " << ei.second;
      DefaultWriteElement<without_cvref_t<K> >(out, i, ei.first);
      out << ": ";
      DefaultWriteElement<without_cvref_t<V> >(out, i, ei.second);
   }


   typedef TSequenceWriter::FSeparators open_sep_close_t;

} // namespace sequence_io
} // namespace ct


#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 10)
namespace ct {
namespace sequence_io {

   template<class T, size_t N> // <-- WARNING: this will not match any std::array<...> if you have another integral type than size_t there! (e.g., unsigned).
   std::ostream &operator << (std::ostream &out, std::array<T,N> const &A) {
      return TSequenceWriter({"[",", ","]"}).WriteList(out, A);
   }


   template<class... Args>
   std::ostream &operator << (std::ostream &out, std::tuple<Args...> const &A) {
      return TSequenceWriter({"(",", ",")"}).WriteTuple(out, A);
   }

} // namespace sequence_io
} // namespace ct
#endif // #if (SEQUENCE_IO_PREDEFINE_LEVEL >= 10)


namespace ct {
namespace sequence_io {

#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 20)
// define fully generic << operators for some more common standard containers: vector, list, set, map.

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::vector<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::list<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::set<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::multiset<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::map<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::multimap<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

#endif // if SEQUENCE_IO_PREDEFINE_LEVEL >= 20


#if (SEQUENCE_IO_PREDEFINE_LEVEL >= 30)
   // define fully generic << operators for rarely encountered containers

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::deque<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::forward_list<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::unordered_set<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq);
   }

   template<class T, class A> std::ostream &operator << (std::ostream &out, std::unordered_multiset<T,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::unordered_map<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

   template<class K, class V, class A> std::ostream &operator << (std::ostream &out, std::unordered_multimap<K,V,A> const &seq) {
      return ct::sequence_io::TSequenceWriter(open_sep_close_t{"{|",", ","|}"}).WriteList(out, seq, &DefaultWriteMapElement<K,V>);
   }

#endif // if SEQUENCE_IO_PREDEFINE_LEVEL >= 20



} // namespace sequence_io
} // namespace ct



// need to either define these before including format.h or put them into namespace std to make
// them accessible via argument-dependent-lookup in fmt::format/fmt::print...
// at least formally, and clang++ actually enforces it:
// /home/cgk/dev/cx/format.h:2349:6: error: call to function 'operator<<' that is neither visible in the template definition nor found by argument-dependent lookup
// namespace std { using ct::sequence_io::operator <<; }

// template<class T, size_t N>
// std::ostream &operator << (std::ostream &out, std::array<T,N> const &A) {
//    ct::TSequenceWriter({"[",", ","]"}).WriteList(out, A);
//    return out;
// }

#endif // CX_SEQUENCE_IO_H
