#ifndef CX_OS_INT
#define CX_OS_INT

#include <string>
#include <vector>

// Some operating system specific stuff, mainly related to files and paths.

// Set some OS-specific macros. These are TARGET_OS_* and CX_PATH_SEPARATOR (the
// default path separator character, as a preprocessor macro string-literal)
//
// Distinction here is based on:
// https://stackoverflow.com/questions/5919996/how-to-detect-reliably-mac-os-x-ios-linux-windows-in-c-preprocessor
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
   // some sort of windows. Note: _WIN32 is also there for x64 compilation.
   #define CX_TARGET_OS_WINDOWS
   #define CX_PATH_SEPARATOR "\\" // backslash
   #define CX_PATH_ALL_SEPARATORS "\\/" // backslash or front slash (windows allows using both for relative paths)
   #define CX_PATH_NUM_SEPARATORS 2
#elif __APPLE__
    #include <TargetConditionals.h>
   // this sets stuff like TARGET_OS_MAC.
   #ifdef TARGET_OS_MAC
      // MacOS has POSIX interfaces.
      #define CX_TARGET_OS_MAC
      #define CX_TARGET_OS_POSIX
   #endif
   #define CX_PATH_SEPARATOR "/" // front slash
#elif __linux__
   // some sort of linux
   #define CX_TARGET_OS_LINUX
   #define CX_TARGET_OS_POSIX
   #define CX_PATH_SEPARATOR "/" // front slash
#elif __unix__ // all unices not caught above
   // other unix
   #define CX_TARGET_OS_UNIX
   #define CX_TARGET_OS_POSIX
   #define CX_PATH_SEPARATOR "/" // front slash
#elif defined(_POSIX_VERSION)
   // non-unix with POSIX interface
   #define CX_TARGET_OS_POSIX
   #define CX_PATH_SEPARATOR "/" // front slash
#endif

#ifndef CX_PATH_ALL_SEPARATORS
   #define CX_PATH_ALL_SEPARATORS CX_PATH_SEPARATOR
   #define CX_PATH_NUM_SEPARATORS 1
#endif



namespace ct {
   typedef std::string FFileName;
   typedef std::pair<FFileName, FFileName> FFileName_PathAndBase;
   typedef std::pair<FFileName, FFileName> FFileName_StemAndExt;
   typedef std::vector<FFileName> FFileNames;

#ifdef INCLUDE_ABANDONED
#if 0
   namespace path {
      typedef std::string::iterator
         str_it;
      typedef std::string::const_iterator
         cstr_it;

      struct is_path_separator_c
      {
         inline bool operator() (char c) const
         {
         #ifdef _WIN32
            return c == '/' || c == '\\';
         #else
            return c == '/';
         #endif
         }
      };


      struct is_equal_c
      {
         is_equal_c(char c_)  : m_TargetChar(c_) {}
         
         inline bool operator() (char c) const {
            return c == m_TargetChar;
         }
      protected:
         char
            m_TargetChar;
      };


      static std::string const s_EmptyDummyStrToMakeDefinedIterators;


      struct string_slice
      {
         cstr_it first, last;

         string_slice() { first = s_EmptyDummyStrToMakeDefinedIterators.end(); last = first; };
         string_slice(cstr_it itFirst_, cstr_it itLast_) : first(itFirst_), last(itLast_) { _CxAssert(first <= last); }
         string_slice(std::string const &s) : first(s.begin()), last(s.end()) {}

         // make this object reference the external std::string 's'. Beware that 's'
         // needs to stay alive for the lifetime of *this!
         void assign_ref(std::string const &s) { first = s.begin(); last = s.end(); }

         cstr_it begin() const { return first; }
         cstr_it end() const { return last; }
         size_t size() const { return last - first; }
         bool empty() const { return first == last; }
         void clear() { first = last; }

         std::string to_str(unsigned Flags = 0) const;

         cstr_it find(char c, size_t offset = 0, unsigned flags=0) const;
         cstr_it find(char const *c, size_t offset = 0, unsigned flags=0) const;
         template<class FCharPred>
         cstr_it find_c(FCharPred const &pred, size_t offset = 0, unsigned flags=0) const;

         cstr_it rfind(char c) const;
         cstr_it rfind(char const *c) const;
         template<class FCharPred>
         cstr_it rfind_c(FCharPred const &pred) const;

         string_slice substr(size_t iFirst, size_t nCount) const {
            _CxAssert(iFirst + nCount <= size());
            return string_slice(first + iFirst, first + iFirst + nCount);
         }

         string_slice substr(cstr_it iFirst, size_t nCount) const {
            _CxAssert(iFirst + nCount <= end());
            return string_slice(iFirst, iFirst + nCount);
         }

         bool operator == (char const *p) const;
         bool operator != (char const *p) const { return !this->operator == (p); }
         
         // similar to startswith(), but returns an iterator to the end of the
         // matching sequence in *this in case of a successful match.
         // Returns this->first unless the sequence in `search` matches completely.
         // NOTE:
         // - this returns `first` on failure, not `last`! `last` is a valid output
         //   if *this is a complete match with `search`.
         cstr_it match_at_start(char const *search) const;
         cstr_it match_at_start(string_slice const &search) const;
         // similar to endsswith(), but returns an iterator to the start of the
         // matching sequence in *this in case of a successful match.
         // Returns this->last unless the sequence in `search` matches completely.
         cstr_it match_at_end(char const *search) const;
         cstr_it match_at_end(string_slice const &search) const;
         
         bool startswith(char const *p) const;
         bool startswith(string_slice const &search) const;
         bool endswith(char const *pStr) const;
         bool endswith(string_slice const &search) const;
         bool startswith(char c) const;
         bool endswith(char c) const;

         // checks if *this starts with pSearch; if not, returns `false`; if yes,
         // increments this->first to the first character behind the matching string,
         // and returns `true`. This is a convenience wrapper around try_skip_literal(),
         // which can also handle some more complex arrangements.
         //
         // Example:
         //   string_slice sl("df-rhf");
         //   bool df = false;
         //   if (sl.remove_if_at_start("df-")) {
         //      df = true;
         //   }
         //
         //   At the end, `sl` will be covering the string "rhf" only, and, `df`
         //   will be `true`. If instead, `sl` is initialized to "rhf" only at
         //   the start, it will still cover only "rhf" at the end, but `df`
         //   will be false.
         bool remove_if_at_start(char const *pSearch);
         // checks if *this ends with pSearch; if not, returns `false`; if yes,
         // decrements this->last to the first character of the matching string,
         // and returns `true`.
         bool remove_if_at_end(char const *pSearch);
         
         void resize(size_t n) {
            last = first + n;
         };

         // skip all characters contained in pCharsToSkip (0-terminated) at begin and end.
         // If pCharsToSkip is 0, this skips whitespace.
         void trim(char const *pCharsToSkip=0);
         // skip any characters in pCharsToSkip (if pCharsToSkip != 0) or whitespace (if pCharsToSkip==0) at end only
         void trim_right(char const *pCharsToSkip=0);
         // skip any characters in pCharsToSkip (if pCharsToSkip != 0) or whitespace (if pCharsToSkip==0) at start only
         void trim_left(char const *pCharsToSkip=0);

         FPair rsplit1(char Delim, size_t Flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;
         FPair split1(char Delim, size_t Flags = SPLIT_SkipWhitespace | SPLIT_SkipParenthesis) const;

         pchar_range to_pchar_range() const {
            if (empty())
               return pchar_range(0,0);
            else
               return pchar_range(&*first, &*first + size());
               // ^- why &*first + size() instead of &*last? the latter would
               // formally dereference a potentially invalid iterator, and some
               // debug CRTs might take offence to that.
         }
         
         operator std::string () const { return to_str();  };
         operator pchar_range () const { return this->to_pchar_range();  };

         char const &operator [] (size_t i) const { _CxAssert(i < size()); return first[i]; };
      protected:
         // advance iterator iPos, to next position.
         // Flags may specify SPLIT_SkipParenthesis and/or SPLIT_SkipWhitespace
         // to skip over blocks between matching parenthesis and/or whitespace
         inline void increment_it(cstr_it &iPos, unsigned Flags) const;
         void skip_over_flagged(cstr_it &iPos, unsigned Flags) const;
      };

      std::ostream &operator << (std::ostream &out, string_slice const &s);
   }
#endif
#endif // INCLUDE_ABANDONED
   
   namespace path {
      /// Split file name into path component (including separator) and base name
      /// (e.g., `ct::path::split("/tmp/whee.xyz")` --> `{"/tmp/", "whee.xyz"}`)
      FFileName_PathAndBase split(FFileName const &FileName);
      
      /// Split file name into base (including path) and file extension, including '.'
      /// (e.g., `ct::path::split_ext("/tmp/whoo/whee.xyz")` --> `{"/tmp/whoo/whee", ".xyz"}`).
      ///
      /// Notes:
      /// - This assumes that the file extension is last part of the file name after the last '.'.
      /// - If the file does not actually have a logical extension (a unfortunate practice in many unixes),
      ///   but still contains a '.' with some other meaning, this will yield erratic results
      ///   (e.g., for a file "read.me.first" it would return ("read.me", ".first"))
      FFileName_StemAndExt split_ext(FFileName const &FileName);
      
      /// Split file name into individual path components (excluding inner connecting separators) and base name.
      /// @param MaxSplit If specified, only up to this many sub-paths are split off, starting at the deepest level 
      ///
      /// Examples:
      /// ```
      ///   >> ct::path::split_path("/tmp/whoo/whee.xyz")
      ///   {"/tmp", "whoo", "whee.xyz"}
      ///
      ///   >> ct::path::split_path("/opt/programs/iboview2021/bin/iboview", 2)
      ///   {"/opt/programs/iboview2021", "bin", "iboview"}
      /// ```
      FFileNames split_path(FFileName const &FileName, size_t MaxSplit = std::string::npos);
      
      /// Join path `a` and path fragment or file name `b` using the OS' path separator
      FFileName join(FFileName const &a, FFileName const &b);

      /// Join path fragments and/or base file name in L using the OS' path separator
      FFileName join(FFileNames const &L);

      /// Removes redundant path separators, and resolves '.' and '..' references, if possible
      FFileName normalize_path(FFileName const &FileName);
      /// Examples:
      /// ```
      ///   >> ct::path::normpath("/opt/programs/iboview2021/bin/../bases")
      ///   {"/tmp", "whoo", "whee.xyz"}
      ///
      ///   >> ct::path::split_path("/opt/programs/iboview2021/bin/iboview", 2)
      ///   {"/opt/programs/iboview2021", "bin", "iboview"}
      /// ```

      /// Checks whether a file (given by name) could be opened for reading (if it can't this typically means it does not exist)
      bool file_readable_q(char const *pFileName);

      /// Checks whether a file (given by name) could be opened for reading (if it can't this typically means it does not exist)
      bool file_readable_q(std::string const &FileName);

      /// Returns `PathName` with an added trailing path separator, unless PathName already ends in one
      /// (e.g., `with_trailing_path_separator("/tmp") --> "/tmp/"`, but `"/tmp/"` returned unchanged)
      FFileName with_trailing_path_separator(FFileName &&PathName);
      /// Returns `PathName` with an added trailing path separator, unless PathName already ends in one
      /// (e.g., `with_trailing_path_separator("/tmp") --> "/tmp/"`, but `"/tmp/"` returned unchanged)
      FFileName with_trailing_path_separator(FFileName const &PathName);
      /// Returns `PathName` with trailing path separator removed, provided it ends in one,
      /// (e.g., `without_trailing_path_separator("/tmp/") --> "/tmp"`, but `"/tmp"` returned unchanged)
      /// and provided the trailing path separator is not part of a root-directory or network-share
      /// declaration or similar (e.g., does not turn "/" into "" or leading double-backslash (part of
      /// UNC/SMB network share name) into a single backslash)
      FFileName without_trailing_path_separator(FFileName &&PathName);
      /// Checks whether PathName ends with a path separator, and if yes, removes it
      FFileName without_trailing_path_separator(FFileName const &PathName);
      
   }
   
   FFileName GetExePath();
   std::string FmtLastError();
}


#endif // CX_OS_INT
