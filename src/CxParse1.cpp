#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <cstdlib> // for strtol

// #include <cstring> // for strlen
#include <cctype> // for isalnum, isalpha, etc.
#include "format.h"

//#include "fmt.h"
#include "CxParse1.h"

namespace ct {

char const
   *g_pDefaultPropertyList_OpenDelimCloseDef = "{;}:";

static unsigned _SanitizedPropertyListFlags(unsigned Flags_);

// TOSTR_QuotesRequired itself implies TrimQuotes; that also means that we cannot check if
// QuotesRequired was set by just using bool(Flags & TOSTR_QuotesRequired). Therefore
// this one here is created.
unsigned const TOSTR_QuotesRequiredFlagOnly = string_slice::TOSTR_QuotesRequired & ~string_slice::TOSTR_TrimQuotes;

struct FParenthesisPair { char Open; char Close; };
// ^- note: that's an aggregate type; so it can be { ... }-initialized
// even in C++98 (same with FQuotePair and FEscapePair below)
static FParenthesisPair const s_ParenthesisPairs[] = { {'(',')'}, {'{', '}'}, {'<','>'}, {'[',']'} };
static size_t const s_nParenthesisTypes = sizeof(s_ParenthesisPairs)/sizeof(s_ParenthesisPairs[0]);

struct FQuotePair { char Open; char Close; };
static FQuotePair const s_QuotePairs[] = { {'"','"'}, {'\'', '\''} };
static size_t const s_nQuoteTypes = sizeof(s_QuotePairs)/sizeof(s_QuotePairs[0]);

struct FEscapePair { char Placeholder; char Replacement; };
static FEscapePair const s_BasicEscapeSequences[] = { {'\\', '\\'}, {'\'', '\''}, {'\"', '\"'}, {'n', '\n'}, {'t','\t'} };
static size_t const s_nBasicEscapeSequences = sizeof(s_BasicEscapeSequences)/sizeof(s_BasicEscapeSequences[0]);



// Used in place of `char` to represent short immutable multi-byte sequences of
// raw chars (e.g., many utf-8 encoded Unicode code-points). Employs some dirty
// tricks for interacting with normal-looking code (e.g. the operator == with
// a 'char const &' which actually takes the object address and looks beyond that...)
//
// Note: we might want to change this into something with a fixed maximum size.
// currently it is a bit of a hybrid of static/dynamic storage. Or even to templatize
// this thing over the storage data type.
struct short_str_t {
   // --- UPCOMING ---
   // c/p-style interface which simply forwards array operations to member
   // object m_data.
   using value_type = char;
   using pointer = value_type*;
   using iterator = pointer;
   using const_pointer = value_type const*;
   using const_iterator = const_pointer;

   explicit short_str_t(value_type const *pSeq) : m_data(pSeq, pSeq + std::char_traits<char>::length(pSeq)) {}
   explicit short_str_t(std::initializer_list<value_type> init) : m_data{init} {}
   template<class _Iterator>
   explicit short_str_t(_Iterator first, _Iterator last) : m_data{first, last} {}
   explicit short_str_t(string_slice const &sl) : m_data{sl.begin(), sl.end()} {}

   bool empty() const { return m_data.empty(); }
   size_t size() const { return m_data.size(); }
   value_type &operator [] (size_t i) { return m_data[i]; }
   value_type const &operator [] (size_t i) const { return m_data[i]; }
   value_type *data() { return m_data.data(); }
   value_type const *data() const { return m_data.data(); }

   iterator begin() { return &*m_data.begin(); }
   iterator end() { return &*m_data.end(); }
   const_iterator begin() const { return &*m_data.begin(); }
   const_iterator end() const { return &*m_data.end(); }

   bool operator == (short_str_t const &other) const { return m_data == other.m_data; }
   bool operator != (short_str_t const &other) const { return !(*this == other); }
//    void operator = (short_str_t const &other) { m_data = other.m_data; }
//    void operator = (short_str_t &&other) { m_data = std::move(other.m_data); }
   // ^-- default implementations of those and copy/move ctor should do the trick

   // --- UPCOMING ---
   // The core functionality this object. In particular, this concerns
   // interactions with char const * pointers and dereferenced cstr_it
   // objects. These are meant to allow formulating string algorithms
   // very similarly to how they would be formulated if short_str_t would
   // be replaced by an actual elementary char object.

   // returns whether ∀i ∈ {0,…,A.size()-1} : (&MemRef)[i] == A[i].
   // ...so: WARNING: (a) this accesses memory beyond MemRef itself,
   // (b) making sure that there *are* at least A.size() chars of meaningful
   // data at &MemRef is the caller's burden.
   // Also, this is meant for *short* strings, so it intentionally does
   // not early-out on the memory comparison (at least formally).
   friend bool operator == (char const &MemRef, short_str_t const &A) {
      char const
         *pMemRef = &MemRef,
         *pA = A.data();
      size_t N = A.size();
      bool bRes = true;
      for (size_t i = 0; i < N; ++ i)
         bRes = bRes && (pA[i] == pMemRef[i]);
      return bRes;
      // ^-- this weird construction is meant to communicate something to the compiler:
      //
      //     """yes, you can definitely assume that the memory at (&MemRef)[i] {i = 0,…,(N-1)} can
      //     be dereferenced, regardless of any conditionals involving whatever data (&MemRef)[...]
      //     and A.data()[...] may hold"""
   }
   friend bool operator == (short_str_t const &A, char const &MemRef) { return MemRef == A; }
   friend bool operator != (short_str_t const &A, char const &MemRef) { return !(MemRef == A); }
   friend bool operator != (char const &MemRef, short_str_t const &A) { return !(MemRef == A); }

   // returns whether the sequence [first,last) is long enough to contain
   // at least one copy of *this.
   template<class _Iterator>
   bool would_fit(_Iterator first, _Iterator last) const { return first + size() <= last; }
   // fits_in?

   // returns whether the sequence [first,last) is too short to contain
   // at least one copy of *this.
   template<class _Iterator>
   bool longer_than(_Iterator first, _Iterator last) const { return first + size() > last; }
   // ^--   longer_than_slice? longer_than_range? too_long_for?
protected:
   enum {
      // max amount of chars we can store in TArray's local static data buffer
      // without dynamic allocations. Note: that also means these have to be
      // copied around if these guys here are passed around by value.
      _StaticSize = 16
   };
   using array_t = TArray<char, _StaticSize>;
   array_t m_data;
   // TODO: maybe look into this --> https://en.cppreference.com/w/cpp/types/aligned_storage
};



std::string ReadFile(std::string FileName)
{
   std::ifstream
      File(FileName.c_str());
   File.seekg(0, std::ifstream::end);
   for (; File.good();) {
      std::string
         Result(File.tellg(), '\0');
      File.seekg(0, std::ifstream::beg);
      File.read(&Result[0], Result.size());
      if (!File.good())
         break;
      return Result;
   }
   throw std::runtime_error("Failed to read file '" + FileName + "'.");
}

// cstr_it FindIdentifierEnd( cstr_it iPos, cstr_it iEnd )
// {
//    while( iPos != iEnd && ( std::isalnum(*iPos) || '_' == *iPos )  )
//       ++ iPos;
//    return iPos;
// };


bool FIdentifierDecl::IsCharAllowedAsFirst(char c) const {
   return std::isalpha(static_cast<unsigned char>(c)) ||
      (Enabled(AllowUnderscores) && !Enabled(ForbidUnderscoreAsFirst) && '_' == c);
}


bool FIdentifierDecl::IsCharAllowed(char c, char cLast) const {
   if (std::isalpha(static_cast<unsigned char>(c)))
      return true;
   if (Enabled(AllowDigits) && std::isdigit(static_cast<unsigned char>(c)))
      return true;
   if (Enabled(AllowUnderscores) && c == '_')
      return true;
   if (Enabled(AllowMinus | AllowColon | AllowCaret)) {
      if (std::isalnum(static_cast<unsigned char>(cLast))) {
         // these ones only allowed if last character was alpha or digit
         if (Enabled(AllowMinus) && c == '-') return true;
         if (Enabled(AllowCaret) && c == '^') return true;
      }
      if (Enabled(AllowColon) && c == ':') return true;
   }
   return false;
}


cstr_it FindIdentifierEnd(cstr_it iPos, cstr_it iEnd, FIdentifierDecl const &Decl)
{
   if (iPos == iEnd || !Decl.IsCharAllowedAsFirst(*iPos))
      return iPos;
   for ( ; ; ) {
      char cLast = *iPos;
      ++ iPos;
      if (iPos == iEnd || !Decl.IsCharAllowed(*iPos, cLast))
         break;
   }
   return iPos;
}

bool IsIdentifier(cstr_it iPos, cstr_it iEnd, FIdentifierDecl const &Decl) {
   return FindIdentifierEnd(iPos, iEnd, Decl) == iEnd;
}

bool IsIdentifier(string_slice sl, FIdentifierDecl const &Decl) {
   return FindIdentifierEnd(sl.first, sl.last, Decl) == sl.last;
}


// cstr_it FindIdentifierEnd(cstr_it iPos, cstr_it iEnd, FIdentifierDecl const &Decl)
// {
//    if (iPos == iEnd || !Decl.IsCharAllowedAsFirst(*iPos))
//       return iPos;
//    ++ iPos;
//    while (iPos != iEnd && Decl.IsCharAllowed(*iPos))
//       ++ iPos;
//    return iPos;
// }

#ifdef INCLUDE_ABANDONED
// cstr_it FindIdentifierEnd(cstr_it iPos, cstr_it iEnd, FIdentifierDecl const &Decl)
// {
//    if (iPos == iEnd)
//       return iPos;
//    char
//       c = *iPos, cLast;
//    if (!Decl.IsCharAllowedAsFirst(c))
//       return iPos;
//    for ( ; ; ) {
//       cLast = c;
//       ++ iPos;
//       c = *iPos;
//       if (iPos == iEnd && !Decl.IsCharAllowed(c, cLast))
//          break;
//    }
//    return iPos;
// }
#endif // INCLUDE_ABANDONED


// assumes that *pFirst == ChOpen. Returns pLast if no closing braket found.
// TODO: add options to skip over quote pairs?
char const *FindClosingParenthesis(char Open, char Close, char const *pFirst, char const *pLast)
{
   _CxAssert(pFirst < pLast && *pFirst == Open);
   _CxAssert(Open != Close);
   int
      // number of brackets still open.
      iLevel = 0;
   for (; pFirst != pLast; ++pFirst) {
      if (*pFirst == Open)
         iLevel += 1;
      if (*pFirst == Close)
         iLevel -= 1;
      if (iLevel == 0)
         break;
   }
   return pFirst;
   IR_SUPPRESS_UNUSED_WARNING(Open);
}


// Assumes that *pFirst == ChOpen. Returns pLast if no closing braket found.
// Allows multi-byte sequences for pOpen and pClose (such as UTF-8-encoded
// unicode parenthesis or separators. Hooray for ⟨a|b⟩!), but they do need
// to fit into a `short_str_t`.
char const *FindClosingParenthesis(char const *pOpen, char const *pClose, char const *pFirst, char const *pLast)
{
   short_str_t Open(pOpen), Close(pClose);
   _CxAssert(pFirst < pLast && *pFirst == Open);
   _CxAssert(Open != Close);
   int
      // number of brackets still open.
      iLevel = 0;
   for (; pFirst <= pLast - Close.size(); ++pFirst) {
//       if (pLast - pFirst >= ptrdiff_t(Open.size()) && *pFirst == Open)
      if (Open.would_fit(pFirst, pLast) && Open == *pFirst)
         iLevel += 1;
      if (*pFirst == Close)
         // ^-- min-remaining-size requirement fixed by adjusting
         //     upper loop limit. No extra check needed in this one.
         iLevel -= 1;
      if (iLevel == 0)
         return pFirst;
   }
   return pLast;
   // ^-- note subtle change of loop iLevel == 0 case and return
   //     value here (compared to straight char version)
   IR_SUPPRESS_UNUSED_WARNING(Open);
}



// assumes that *pFirst == ChOpen. Returns pLast if no closing quote found.
// Unlike in FindClosingParenthesis, the symbols for opening and closing quotes
// (Open, Close) may be identical.
//
// Notes:
// - This function will ignore inner closing quote symbols escaped with a
//   backslash
// - However, it will *NOT* remove the escape symbols from the inner quotes, if
//   present; note that doing so would require modification of the input string,
//   which string_slice functions cannot do in the actual analysis stage,
//   because they operate on the input string, and it is constant.
// - To interpret the inner escaped characters, use
//   string_slice.to_str(string_slice::TOSTR_QuotesRequired | string_slice::TOSTR_InterpretEscapes))
//   As this function *does* create a new string as return value, it can
//   perform the escape conversions there.
char const *FindClosingQuote(char Open, char Close, char const *pFirst, char const *pLast)
{
   _CxAssert(pFirst < pLast && *pFirst == Open);
   // skip over the opening quote at *pFirst.
   ++ pFirst;
   int
      iEscapeLevel = 0;
   for (; pFirst != pLast; ++pFirst) {
      if (*pFirst == '\\') {
         iEscapeLevel += 1;
         continue;
      }
      if (*pFirst == Close && iEscapeLevel == 0)
         return pFirst;
      iEscapeLevel = 0;
   }
   _CxAssert(pFirst == pLast);
   return pLast;
   IR_SUPPRESS_UNUSED_WARNING(Open);
}



std::ostream &operator << (std::ostream &out, string_slice const &s)
{
   out << s.to_str();
   return out;
}

bool is_whitespace_cxp1(char c) {
   return c == ' ' || c == '\t';
}


static bool is_contained_in_or_whitespace(char c, char const *pCharsToTest) {
   if (pCharsToTest) {
      for (char const *pt = pCharsToTest; *pt != 0; ++ pt)
         if (*pt == c)
            return true;
      return false;
   } else {
      return is_whitespace_cxp1(c);
   }
}


cstr_it SkipWhiteSpace( cstr_it iPos, cstr_it iEnd ){
   while (iPos != iEnd && is_whitespace_cxp1(*iPos))
      ++ iPos;
   return iPos;
}


// string_slice::string_slice(char const *p)
//    : first(p, p + std::strlen(p))
// {
// }


void string_slice::trim_left(char const *pCharsToSkip)
{
   while (first != last && is_contained_in_or_whitespace(*first, pCharsToSkip))
      ++ first;
}


void string_slice::trim_right(char const *pCharsToSkip)
{
   if (first == last)
      return;
   while (first != (last-1) && is_contained_in_or_whitespace(*(last-1), pCharsToSkip))
      -- last;
}


void string_slice::trim(char const *pCharsToSkip)
{
   _CxAssert(first <= last);
   if (pCharsToSkip == 0) {
      while (first != last && is_whitespace_cxp1(*first))
         ++ first;
      if (first == last)
         return;
      while (first != (last-1) && is_whitespace_cxp1(*(last-1)))
         -- last;
   } else {
      trim_left(pCharsToSkip);
      trim_right(pCharsToSkip);
   }
}


char string_slice::trim_quotes(char const *pQuoteChars, bool QuotesRequired)
{
   // remove leading and trailing whitespace in any case---BEFORE quote removal.
   // (i.e., whitespace *inside* the quotes will stay)
   trim();

   // TODO: possibly synchronize this with new global QuotePairs mechanism.
   // in skip_over_flagged.
   if (pQuoteChars == 0)
      pQuoteChars = "''\"\""; // these come in matching pairs for front/back, and the str is null-terminated.
   size_t
      nQuoteChars = 0;
   for (char const *p = pQuoteChars; *p != 0; ++p)
      nQuoteChars += 1;
   if ((nQuoteChars & 1) != 0)
      throw std::runtime_error("string_slice::trim_quotes: set of input quote character pairs must have even length, but I got <" + std::string(pQuoteChars) + ">.");
   if (size() >= 2) {
      for (size_t iQuotePair = 0; iQuotePair != nQuoteChars/2; ++ iQuotePair) {
         char
            Open = pQuoteChars[2*iQuotePair],
            Close = pQuoteChars[2*iQuotePair+1];
         if (*first == Open && *(last-1) == Close) {
            ++ first;
            -- last;
            return Open;
         }
      }
      // if we arrived here, the quote match above has not occurred for any
      // quote pair.
      if (QuotesRequired)
         throw std::runtime_error("string_slice::trim_quotes: quotes are set as mandatory, but input string <" + to_str() + "> is not enclosed in quotes.");
   }
   if (QuotesRequired) {
      throw std::runtime_error("string_slice::trim_quotes: expected input str to be quoted, but got <" + to_str() + ">");
   }
   return 0;
}


char const *string_slice::try_skip_literal(char const *pSequence, unsigned Flags)
{
   _CxAssert(pSequence != 0);
   if (Flags & SKIP_SkipWhitespaceBefore) {
      string_slice copy_of_this = *this;
      copy_of_this.trim_left();
      return try_skip_literal(pSequence, Flags & (~SKIP_SkipWhitespaceBefore));
   }

   if (Flags == 0) {
      // no special things to consider. Take a shortcut.
      cstr_it it = this->match_at_start(pSequence);
      if (*pSequence == 0 || it != first) {
         // ^-- *pSequence == 0 test: skipping an empty string always succeeds.
         first = it;
         return pSequence;
      }
      return 0;
   }

   // this is a case with some non-trivial options. Proceed like in
   // string_slice::startswith, but with consideration of some additional
   // optional details.
   char const
      *p = pSequence;
   cstr_it it = first;
   for (; *p && it != last; ++it, ++p) {
      bool is_equal;
      if (bool(Flags & SKIP_IgnoreCase)) {
         is_equal = (std::tolower(static_cast<unsigned char>(*p)) == std::tolower(static_cast<unsigned char>(*it)));
      } else {
         is_equal = (*p == *it);
      }
      if (!is_equal)
         return 0;
   }
   if (!(*p == 0))
      return 0;
   if (it == last) {
      // *this is equal to pSequence --- nothing following.
      // That's an successful match, but we need not process anything more,
      // and the string is empty after removing the sequence.
      first = it;
      return pSequence;
   }

   char
      // character directly following the end of the matched sequence.
      next = *it;
   bool
      match_accepted = true;
   if (bool(Flags & SKIP_TokenUnderscore) && next == '_')
      match_accepted = false;
   if (bool(Flags & SKIP_TokenMinus) && next == '-')
      match_accepted = false;
   if (bool(Flags & SKIP_TokenAlpha) && std::isalpha(next))
      match_accepted = false;
   if (bool(Flags & SKIP_TokenDigit) && std::isdigit(next))
      match_accepted = false;
   if (!match_accepted) {
      // *this did start with the target sequence, but we processed it
      // in full-token mode, and the next character in *this would have
      // been part of a corresponding token to match
      return 0;
   }

   first = it;
   if (bool(Flags & SKIP_SkipWhitespaceAfter))
      trim_left();
   return pSequence;
}


cstr_it string_slice::match_at_start(char const *p) const
{
   cstr_it it = first;
   for (; *p && it != last; ++it, ++p)
      if (*p != *it)
         return first;
   if (*p != 0)
      return first;
   return it;
}


cstr_it string_slice::match_at_start(string_slice const &search) const
{
   size_t n = search.size();
   if (this->size() < n)
      return first;
   size_t i = 0;
   for (; i < n; ++ i)
      if (this->first[i] != search.first[i])
         return first;
   return first + i;
}


bool string_slice::startswith(char const *p, unsigned Flags) const
{
   if (Flags == 0) {
      return this->startswith(p);
   } else {
      string_slice copy_of_this = *this;
      return bool(copy_of_this.try_skip_literal(p, Flags));
   }
}


bool string_slice::startswith(char c) const
{
   return (!empty() && (*first) == c);
}


bool string_slice::endswith(char c) const
{
   return (!empty() && (*(last-1)) == c);
}


bool string_slice::startswith(char const *search) const
{
   _CxAssert(search != 0);
   return *search == 0 || first != match_at_start(search);
}


bool string_slice::startswith(string_slice const &search) const
{
//    _CxAssert(!search.empty());
   return search.empty() || first != match_at_start(search);
}

// bool string_slice::startswith(char const *p) const
// {
//    cstr_it it = first;
//    for (; *p && it != last; ++it, ++p)
//       if (*p != *it)
//          return false;
//    return *p == 0;
// }
//
//
// bool string_slice::startswith(string_slice const &search) const
// {
//    size_t n = search.size();
//    if (this->size() < n)
//       return false;
//    for (size_t i = 0; i < n; ++ i)
//       if (this->first[i] != search->first[i])
//          return false;
//    return true;
// }

// bool string_slice::startswith(char c) const
// {
//    char const buf[2] = {c, 0};
//    return this->startswith(&buf[0]);
// }
//
// bool string_slice::endswith(char c) const
// {
//    char const buf[2] = {c, 0};
//    return this->endswith(&buf[0]);
// }

cstr_it string_slice::match_at_end(char const *pStr) const
{
   if (empty())
      return last;
   char const *p = pStr;
   while (*p != 0)
      ++p;
   string_slice
      slEndOfThis(last - (p - pStr), last);
   if (slEndOfThis == pStr)
      return slEndOfThis.first;
   else
      return last;
}


cstr_it string_slice::match_at_end(string_slice const &search) const
{
   size_t n = search.size();
   if (this->size() < n)
      return last;
   cstr_it match_start = this->last - n;
   for (size_t i = 0; i < n; ++ i)
      if (match_start[i] != search.first[i])
         return last;
   return match_start;
}

bool string_slice::endswith(char const *search) const {
   _CxAssert(search != 0);
   return *search == 0 || last != match_at_end(search);
}

bool string_slice::endswith(string_slice const &search) const {
   return search.empty() || last != match_at_end(search);
}

// bool string_slice::endswith(char const *pStr) const
// {
//    if (empty())
//       return false;
//    char const *p = pStr;
//    while (*p != 0)
//       ++p;
//    return string_slice(last - (p - pStr), last) == pStr;
// }
//
//
// bool string_slice::endswith(string_slice const &search) const
// {
//    size_t n = search.size();
//    if (this->size() < n)
//       return false;
//    cstr_it start = this->last - n;
//    for (size_t i = 0; i < n; ++ i)
//       if (start[i] != search->first[i])
//          return false;
//    return true;
// }


bool string_slice::startswith(string_slice const &search, unsigned Flags) const
{
   if (Flags == 0)
      return this->startswith(search);
   else {
      // TODO: this variant is not properly implemented --- would require update of
      // try_skip_literal // try_skip_literals first. That's a bit much for now.
      // So just defer to the pchar version. But, least in debug mode, make sure
      // there are no 0-characters inside the sequence which defines the target
      // string. These are used as 0-terminators in the pchar case...
      _CxAssert(find('\0') == last);
      return this->startswith(&*search.first, Flags);
   }
}


bool looks_like_property_list(string_slice const &sl, unsigned PropertyListFlags_, char const *pOpenDelimCloseDef, FIdentifierDecl const &IdentifierDecl)
{
   unsigned
      PropertyListFlags = _SanitizedPropertyListFlags(PropertyListFlags_);
   cstr_it
      itNameEnd = FindIdentifierEnd(sl.first, sl.last, IdentifierDecl);
   if (bool(PropertyListFlags & PROPLIST_UnnamedList) && itNameEnd != sl.first)
      // this one is supposed to have a *no name*, but has one
      return false;
   if (bool(PropertyListFlags & PROPLIST_NameRequired) && itNameEnd == sl.first)
      // this one is supposed to have a name, but does not
      return false;
   cstr_it
      itPropBeg = itNameEnd;
   char
      Open = pOpenDelimCloseDef[0], Close = pOpenDelimCloseDef[2];
   if ((itPropBeg == sl.last) || !(*itPropBeg == pOpenDelimCloseDef[0]))
      return false;
   char const
      *pPropBeg = &*itPropBeg,
      *pPropEnd = FindClosingParenthesis(Open, Close, pPropBeg, &*sl.last);
   cstr_it
      itPropEnd = itPropBeg + (pPropEnd - pPropBeg);
   if (itPropEnd >= sl.last)
      return false;
   _CxAssert(*pPropEnd == Close);
   ++ itPropEnd;
   return itPropEnd == sl.last;
}


bool looks_like_property_list(string_slice const &sl, unsigned PropertyListFlags, char const *pOpenDelimCloseDef)
{
   FIdentifierDecl
      IdentifierDecl(FIdentifierDecl::AllowMinus | FIdentifierDecl::AllowDigits);
   return looks_like_property_list(sl, PropertyListFlags, pOpenDelimCloseDef, IdentifierDecl);
}


bool string_slice::remove_if_at_start(char const *pSearch)
{
   _CxAssert(pSearch != 0);
   return 0 != try_skip_literal(pSearch, 0);
}


bool string_slice::remove_if_at_end(char const *pSearch)
{
   _CxAssert(pSearch != 0);
   if (*pSearch == 0)
      return true; // removing an empty string always succeeds.
   cstr_it
      it = match_at_end(pSearch);
   if (it != last) {
      last = it;
      return true;
   } else
      return false;
}


template<class TSliceList, class TDelimList>
size_t string_slice::split_any(TSliceList &out, TDelimList *pOutDelims, char const *pDelims, size_t nDelims, size_t Flags) const
{
//    bool const Print = false;
//    if ( Print ) xout << "splitting '" << *this << "' at '" << std::string(pDelims, pDelims + nDelims) << "' into: [";

   out.clear();
   out.push_back(string_slice(first,first));

   if (pOutDelims)
      pOutDelims->clear();

   cstr_it
      iPos = first;

//    while (out.size() < nSplitStaticStorage) {
   while (out.size() < out.max_size()) {
      char Delim = 0;
      for ( ; ; ) {
//       while ( iPos != last && *iPos != Delim ){
         if (iPos == last)
            break;
         for (size_t i = 0; i < nDelims; ++ i)
            if (*iPos == pDelims[i])
               Delim = pDelims[i];
         if ((0 != (Flags & SPLIT_SkipWhitespace)) && is_whitespace_cxp1(Delim)) {
            // skip additional whitespace until we reach something which is not a whitespace.
            while (iPos != last && is_whitespace_cxp1(*iPos))
               iPos += 1;
            iPos -= 1;
            break;
         }
         if (Delim != 0)
            break;
         // ^- FIXME cgk 2020-04-27: reconsider logic in this again. The two breaks look a bit shady.
         //    The iPos -= 1 without check for iPos==last is fine, however, because if iPos cannot
         //    start out as last here due to the 'if (iPos == last) break' just above.

//          if (0 != (Flags & SPLIT_SkipParenthesis)) {
//             for (size_t iParensType = 0; iParensType < s_nParenthesisTypes; ++ iParensType) {
//                char
//                   Open = s_ParenthesisPairs[iParensType].Open,
//                   Close = s_ParenthesisPairs[iParensType].Close;
//                if (*iPos == Open) {
//                   char const
//                      *pFirst = &*iPos,
//                      *pLast = &*(last-1)+1,
//                      *p = FindClosingParenthesis(Open, Close, pFirst, pLast);
//                   if (p == pLast) {
//                      std::stringstream str;
//                      str << "imbalanced parenthesis '" << Open << "'/'" << Close << "' in '" << *this << "'.";
//                      throw std::runtime_error(str.str());
//                   }
//                   iPos += p - pFirst;
//                   -- iPos;
//                }
//             }
//          }
//          ++ iPos;
         if (0 != (Flags & SPLIT_SkipQuotes))
            skip_over_flagged(iPos, SPLIT_SkipQuotes);
         if (0 != (Flags & SPLIT_SkipParenthesis))
            skip_over_flagged(iPos, SPLIT_SkipParenthesis);
            // ^- note: this is not 100% functionally equivalent to the code above; it will
            // end up on the closing parenthesis itself, not one char before it ("-- iPos").
            // but I currently am of the impression that this is the right thing to do.
         ++ iPos;
         out.back().last = iPos;
      }
      if (iPos == last) {
         if (out.back().last == out.back().first)
            out.pop_back();
         if (0 != (Flags & SPLIT_SkipWhitespace))
            for (size_t i = 0; i < out.size(); ++ i)
               out[i].trim();
//          if (Print) {
//             for (size_t i = 0; i < out.size(); ++ i) {
//                if (i != 0)
//                   xout << ", ";
//                xout << "'" << out[i] << "'";
//             }
//             xout << "]"<< std::endl;
//          }
         return out.size();
      }
      if (pOutDelims)
         pOutDelims->push_back(Delim);
      ++ iPos;
      out.push_back(string_slice(iPos,iPos));
   }

   std::stringstream str;
   str << "parsing error: split count overflow in '" << *this << "' while scanning for '"  << std::string(pDelims, pDelims + nDelims) << "'.";
   throw std::runtime_error(str.str());
}

template<class TSliceList>
size_t string_slice::split_list(TSliceList &out, char ListStart, char ListEnd, char ListDelim) const
{
   string_slice
      sub = *this;
   sub.trim();
   _CxAssert((ListStart == 0) == (ListEnd == 0));
   if (ListStart != 0 && ListEnd != 0) {
      // cut out list start and list end indicators (normaly [...]).
      if (sub.size() < 2 || sub[0] != ListStart || sub[sub.size()-1] != ListEnd)
         throw std::runtime_error(fmt::format("Expected argument list ({}..{}..{}), but got '{}'", ListStart, ListDelim, ListEnd, sub.to_str()));
      sub.first += 1;
      sub.last -= 1;
   }
   // split rest normally.
   return sub.split(out, ListDelim, SPLIT_SkipWhitespace | SPLIT_SkipParenthesis);
}




template<class TSliceList>
size_t string_slice::split_lines(TSliceList &out, size_t flags) const
{
   long_split_delim_list
      *pDelims = 0;
   split_any(out, pDelims, "\n", 1, flags & (~SPLIT_SkipWhitespace));
   // fix up line endings in case they were made by "\n\r", "\r\n", or similar.
   for (typename TSliceList::iterator it = out.begin(); it != out.end(); ++ it) {
//       if (it->startswith('\r'))
//          it->first += 1;
//       if (it->endswith('\r'))
//          it->last -= 1;
      while (it->startswith('\r') || it->startswith('\0'))
         it->first += 1;
      while (it->endswith('\r') || it->endswith('\0'))
         it->last -= 1;
      if (0 != (flags & SPLIT_SkipWhitespace))
         it->trim();
   }
   return out.size();
}


void string_slice::remove_trailing_path_separators()
{
   // what is the point of this?
   // When dealing with path names,
   // in some situations we need to keep the last separator
   // to obtain the correct absolute path.
   //
   // E.g., in '/home/cgk/stuff/whee.xyz', removing the trailing '/' is fine.
   // in '/whee.xyz' it is not (because with the separator it will denote the root directory,
   // without it, it will denote a directory relative to the local one).
   //
   // Similarly, both in
   //    in lib/whee.xyz
   // and in
   //    f:\whee.xyz
   // we need to keep the directory, despite the fact that there is only one.

   is_path_separator_c ips;
   for (;first != last;) {
      // return if string is empty or last character is not a path separator.
      if (empty() || !ips(*(last-1)))
         break;
      bool
         must_be_kept_anyway = false;
      // so the last character *is* some sort of path separator.
      if (size() == 1) {
         // this means the string *starts* with this.
         // means it's relative to some root directory. We must keep it.
         must_be_kept_anyway = true;
   #ifdef _WIN32
      } else if (size() >= 2 && *(last-2) == ':') {
         // logical drive designation followed by a path separator.
         must_be_kept_anyway = true;
   #endif
      } else {
      }
      if (must_be_kept_anyway)
         break;
      else {
         last -= 1;
      }
   }
}


string_slice::FPair string_slice::split_path() const
{
   // find the last '/' separator within the string content
   // (directory separator)
   cstr_it
      itDirSep = rfind('/');
#ifdef _WIN32
   // Windows is okay with both front-slash and back-slash as directory
   // separators. So look for a back-slash, too, and if it occurs farther
   // out than the last front slash (or if none was found),
   // take that one.
   cstr_it
      itBackSlash = rfind('\\');
   if (itBackSlash != last && (itDirSep == last || itBackSlash > itDirSep))
      itDirSep = itBackSlash;
#endif
   if (itDirSep != last) {
      FPair
         r(string_slice(first,itDirSep+1), string_slice(itDirSep+1,last));
      r.first.remove_trailing_path_separators();
      return r;
   } else {
      // no directory separator. Return original file name.
      return FPair(string_slice(cstr_it(),cstr_it()), *this);
   }
}


// explicitly instantiate the core split function for the combinations of sequence types we will use.
template size_t string_slice::split_any<short_split_result,short_split_delim_list>(short_split_result &, short_split_delim_list *, char const *, size_t, size_t) const;
template size_t string_slice::split_any<long_split_result,long_split_delim_list>(long_split_result &, long_split_delim_list *, char const *, size_t, size_t) const;
template size_t string_slice::split_any<long_split_result,short_split_delim_list>(long_split_result &, short_split_delim_list *, char const *, size_t, size_t) const;

template size_t string_slice::split_list<short_split_result>(short_split_result &out, char ListStart, char ListEnd, char ListDelim) const;
template size_t string_slice::split_list<long_split_result>(long_split_result &out, char ListStart, char ListEnd, char ListDelim) const;
template size_t string_slice::split_lines<long_split_result>(long_split_result &out, size_t) const;






// void string_slice::split( split_result &out, char Delim, bool SkipWhitespace, bool SkipBraces ) const
// {
//    bool const Print = false;
//    if ( Print ) xout << "splitting '" << *this << "' at '" << Delim << "' into: [";
//
//    out.clear();
//    out.push_back( string_slice(first,first) );
//
//    cstr_it
//       iPos = first;
//
//    while ( out.size() < nSplitStaticStorage ) {
//       while ( iPos != last && *iPos != Delim ){
//          if ( SkipBraces ) {
//             for ( size_t iParensType = 0; iParensType < s_nParenthesisTypes; ++ iParensType ) {
//                char
//                   Open = s_ParenthesisPairs[iParensType].Open,
//                   Close = s_ParenthesisPairs[iParensType].Close;
//                if ( *iPos == Open ) {
//                   char const
//                      *pFirst = &*iPos,
//                      *pLast = &*(last-1)+1,
//                      *p = FindClosingParenthesis( Open, Close, pFirst, pLast );
//                   if ( p == pLast ) {
//                      std::stringstream str;
//                      str << "imbalanced parenthesis '" << Open << "'/'" << Close << "' in '" << *this << "'.";
//                      throw std::runtime_error(str.str());
//                   }
//                   iPos += p - pFirst;
//                   -- iPos;
//                }
//             }
//          }
//          ++ iPos;
//          out.back().last = iPos;
//       }
//       if ( iPos == last ) {
//          if ( SkipWhitespace )
//             for ( size_t i = 0; i < out.size(); ++ i )
//                out[i].trim();
//          if ( Print ) {
//             for ( size_t i = 0; i < out.size(); ++ i ) {
//                if ( i != 0 )
//                   xout << ", ";
//                xout << "'" << out[i] << "'";
//             }
//             xout << "]"<< std::endl;
//          }
//          return;
//       }
//       ++ iPos;
//       out.push_back( string_slice(iPos,iPos) );
//    }
//
//    std::stringstream str;
//    str << "parsing error: split count overflow in '" << *this << "' while scanning for '"  << Delim << "'.";
//    throw std::runtime_error(str.str());
// };

static std::string interpret_escapes(cstr_it first, cstr_it last, FEscapePair const *pEscapePairs, size_t nEscapePairs)
{
   std::string out;
   _CxAssert(last >= first);
   char const EscapeChar = '\\';

   out.reserve(last - first);
   bool
      NextCharIsEscaped = false;
   for (cstr_it it = first; it != last; ++it) {
      if (NextCharIsEscaped) {
         // last it set this one to an escape. see if we have the current char
         // as an escape sequence.
         bool Handled = false;
         for (size_t iEscape = 0; iEscape != nEscapePairs; ++ iEscape) {
            FEscapePair const
               &EscapePair = pEscapePairs[iEscape];
            if (*it == EscapePair.Placeholder) {
               // ...yes. replace escape character itself and current
               // character by target expression
               out.push_back(EscapePair.Replacement);
               Handled = true;
               break;
            }
         }
         if (!Handled) {
            // ...no. replicate both the '\' and the current character
            // to output string.
            out.push_back(EscapeChar);
            out.push_back(*it);
         }
         NextCharIsEscaped = false;
      } else {
         if (*it == EscapeChar) {
            NextCharIsEscaped = true;
         } else {
            NextCharIsEscaped = false;
            out.push_back(*it);
         }
      }
   }
   // is there still an escape character left which was not processed?
   // Might happen at the end of a string.
   if (NextCharIsEscaped) {
      // if yes, replicate it
      out.push_back(EscapeChar);
   }
   return out;
}


std::string string_slice::to_str(unsigned Flags) const
{
   if (bool(Flags & TOSTR_TrimQuotes)) {
      unsigned
         // to_str() flags other than TrimQuotes and QuotesRequired (which
         // are handled here) which still need processing afterwards.
         RemainingFlags = Flags & ~(TOSTR_QuotesRequired | TOSTR_TrimQuotes);
      string_slice copy_of_this = *this;
      char
         TrimmedOpenQuote = copy_of_this.trim_quotes(0, bool(Flags & TOSTR_QuotesRequiredFlagOnly));
      if (bool(Flags & TOSTR_QuotesRequiredFlagOnly) && !bool(Flags & TOSTR_InterpretEscapes)) {
         // if QuotesRequired is set, we have promised to evaluate escaped inner
         // quotes of the same type as the trimmed outer quote.
         // If TOSTR_InterpretEscapes is also set, this will happen in the next
         // step automatically. Otherwise, we should do it separately here.
         FEscapePair
            TrimmedQuoteEscape = { TrimmedOpenQuote, TrimmedOpenQuote };
         std::string
            out1 = interpret_escapes(first, last, &TrimmedQuoteEscape, 1);

         // was this the only flag set?
         if (RemainingFlags == 0) {
            // ...yes, nothing else to do, return directly after interpreting the
            // inner quotes here.
            return out1;
         } else {
            // ...no, other conversion flags must still be processed. So make a
            // new string_slice object for out1, and recurse to that one's
            // to_str() to handle the remaining flags.
            return string_slice(out1).to_str();
         }
      } else {
         return copy_of_this.to_str(RemainingFlags);
      }
   } else if (bool(Flags & TOSTR_InterpretEscapes)) {
      return interpret_escapes(first, last, &s_BasicEscapeSequences[0], s_nBasicEscapeSequences);
   } else {
      return std::string(first,last);
   }
}


bool string_slice::try_convert_to_int(ptrdiff_t *pr) const
{
   char const
      *pEnd = 0;
   long
      r_ = strtol(&*first, const_cast<char**>(&pEnd), 0);
   if (pEnd != &*(last-1) + 1)
      return false;
   if (pr)
      *pr = static_cast<ptrdiff_t>(r_);
   return true;
}


bool string_slice::try_convert_to_size(size_t *pr) const
{
   char const
      *pEnd = 0;
   long
      r_ = strtol(&*first, const_cast<char**>(&pEnd), 0);
   if (pEnd != &*(last-1) + 1)
      return false;
   if (r_ < 0)
      return false;
   size_t
      r = size_t(r_);
   if (pr)
      *pr = static_cast<ptrdiff_t>(r);
   return true;
}


bool string_slice::try_convert_to_float(double *pr) const
{
   char const
      *pEnd = 0;
   double
      r = strtod(&*first, const_cast<char**>(&pEnd));
   if (pEnd != &*(last-1) + 1)
      return false;
   if (pr)
      *pr = r;
   return true;
}


// bool string_slice::try_convert_to_bool(bool *pr) const
// {
//    std::string ss(first,last);
//    for (size_t i = 0; i < ss.size(); ++ i)
//       ss[i] = (char)std::tolower(static_cast<unsigned char>(ss[i]));
//    if (ss == "on" || ss == "true" || ss == "1" || ss == "yes") {
//       if (pr) *pr = true;
//       return true;
//    }
//    if (ss == "off" || ss == "false" || ss == "0" || ss == "no") {
//       if (pr) *pr = false;
//       return true;
//    }
//    return false;
// }

bool string_slice::try_convert_to_tribool(ETriBool *pr) const
{
   std::string ss(first,last);
   for (size_t i = 0; i < ss.size(); ++ i)
      ss[i] = (char)std::tolower(static_cast<unsigned char>(ss[i]));
   if (ss == "on" || ss == "true" || ss == "1" || ss == "yes") {
      if (pr) *pr = TRIBOOL_True;
      return true;
   }
   if (ss == "off" || ss == "false" || ss == "0" || ss == "no") {
      if (pr) *pr = TRIBOOL_False;
      return true;
   }
   if (ss == "auto" || ss == "undef" || ss == "none" || ss == "-1" || ss == "maybe") {
      if (pr) *pr = TRIBOOL_NotDefined;
      return true;
   }
   return false;
}


bool string_slice::try_convert_to_bool(bool *pr) const
{
   ETriBool tb;
   if (!try_convert_to_tribool(&tb))
      return false;
   if (!(tb == TRIBOOL_True || tb == TRIBOOL_False))
      return false;
   *pr = (tb == TRIBOOL_True)? true : false;
   return true;
}


ptrdiff_t string_slice::to_int() const
{
   ptrdiff_t r(0);
   if (!this->try_convert_to_int(&r))
      throw std::runtime_error("expected base10 integer, but got '" + to_str() + "'." );
   return r;
}

size_t string_slice::to_size() const
{
   size_t r(0);
   if (!this->try_convert_to_size(&r))
      throw std::runtime_error("expected non-negative base10 integer, but got '" + to_str() + "'." );
   return r;
}


double string_slice::to_float() const
{
   double r(0);
   if (!this->try_convert_to_float(&r))
      throw std::runtime_error("expected floating point number, but got '" + to_str() + "'." );
   return r;
}


bool string_slice::to_bool() const
{
   bool r;
   if (!this->try_convert_to_bool(&r))
      throw std::runtime_error("expected boolean argument, but got '" + to_str() + "'." );
   return r;
}

ETriBool string_slice::to_tribool() const
{
   ETriBool r;
   if (!this->try_convert_to_tribool(&r))
//       throw std::runtime_error("expected tribool argument (i.e., true/on/yes/1, false/off/no/0, or auto/none/undef/maybe/-1), but got '" + to_str() + "'." );
      throw std::runtime_error("expected tribool argument (i.e., true/on/yes/1, false/off/no/0, or auto/none/undef/maybe/-1), but got '" + to_str() + "'." );
   return r;
}


void string_slice::skip_over_flagged(cstr_it &iPos, unsigned Flags) const
{
   _CxAssert(iPos < last);
   if (0 != (Flags & SPLIT_SkipWhitespace)) {
      // skip additional whitespace until we reach something which is not a
      // whitespace (or the end of the string slice)
      while (iPos != last && is_whitespace_cxp1(*iPos))
         iPos += 1;
   }

   if (0 != (Flags & SPLIT_SkipQuotes)) {
      // If we ended up on a opening quote, skip to the corresponding closing quote.
      //
      // The difference to the (cleaner) parenthesis mechanism below is that the
      // symbols for opening and closing quotes corresponding to each other will
      // typically identical (e.g., ".."), and therefore cannot be cleanly nested.
      //
      // These is some limited support of including quotes inside a string by
      // means of python-style quote escape sequences: A quote symbol preceeded
      // by a backslash will not be be considered as a viable closing quote, and
      // therefore skipped over. This allows specifying strings such as:
      //
      //   "The drama's name is \"Escaped quotes are hard to read by humans\"."
      //
      // (...of course, getting these escapes and quotes into the program input
      // without interpretation may itself be a bit of a challenge if processed
      // by bash or Python or similar)
      //
      // TODO: synchronize this with previous trim_quotes mechanism.
      for (size_t iQuoteType = 0; iQuoteType < s_nQuoteTypes; ++ iQuoteType) {
         char
            Open = s_QuotePairs[iQuoteType].Open,
            Close = s_QuotePairs[iQuoteType].Close;
         if (iPos != last && *iPos == Open) {
            char const
               *pFirst = &*iPos,
               *pLast = &*(last-1)+1,
               *p = FindClosingQuote(Open, Close, pFirst, pLast);
            if (p == pLast) {
               std::stringstream str;
               str << "imbalanced quote " << Open << "..." << Close << " in '" << *this << "'.";
               throw std::runtime_error(str.str());
            }
            _CxAssert(*p == Close);
            iPos += (p - pFirst); // should be on the closing quote now
            _CxAssert(*iPos == Close);
            _CxAssert(iPos < last);
         }
      }
   }

   if (0 != (Flags & SPLIT_SkipParenthesis)) {
      // if we ended up on a opening parenthesis, skip to the corresponding
      // closing parenthesis
      for (size_t iParensType = 0; iParensType < s_nParenthesisTypes; ++ iParensType) {
         char
            Open = s_ParenthesisPairs[iParensType].Open,
            Close = s_ParenthesisPairs[iParensType].Close;
         if (iPos != last && *iPos == Open) {
            char const
               *pFirst = &*iPos,
               *pLast = &*(last-1)+1,
               *p = FindClosingParenthesis(Open, Close, pFirst, pLast);
            if (p == pLast) {
               std::stringstream str;
               str << "imbalanced parenthesis '" << Open << "'/'" << Close << "' in '" << *this << "'.";
               throw std::runtime_error(str.str());
            }
            _CxAssert(*p == Close);
            iPos += (p - pFirst); // should be on the closing parens now
            _CxAssert(*iPos == Close);
            _CxAssert(iPos < last);
         }
      }
   }
}

#ifdef INCLUDE_ABANDONED
// cstr_it string_slice::find(char const *c, std::size_t offset) const
// {
//    for (cstr_it it = first + offset; it != last; ++it)
//       if (string_slice(it, last).startswith(c))
//          return it;
//    return last;
// }
//
//
// cstr_it string_slice::find(char c, std::size_t offset) const
// {
//    return std::find(first + offset, last, c);
// }
#endif // INCLUDE_ABANDONED

cstr_it string_slice::find(char const *c, std::size_t offset, unsigned flags) const
{
   for (cstr_it it = first + offset; it < last; increment_it(it, flags))
      if (string_slice(it, last).startswith(c))
         return it;
   return last;
}


cstr_it string_slice::find(char c, std::size_t offset, unsigned flags) const
{
   _CxAssert(first + offset <= last);
   for (cstr_it it = first + offset; it < last; increment_it(it, flags))
      if (*it == c)
         return it;
   return last;
}


cstr_it string_slice::rfind(char c) const
{
   return rfind_c(is_equal_c(c));
}


cstr_it string_slice::rfind(char const *c) const
{
   cstr_it it(last - 1);
   for ( ; it >= first; --it) {
      if (string_slice(it, last).startswith(c))
         return it;
   }
   return last;
}


cstr_it string_slice::find_closing_parenthesis(char open, char close) const
{
   // Can't find anything in an empty string slice --> fail.
   if (this->empty()) return this->first;

   if (open == 0)
      open = *first;
   else if (*first != open)
      // An explicit opening parenthesis was given, but there
      // isn't one of it's type at the start of *this --> fail
      return this->first;

   // If no `close` was given, attempt to auto-identify it based on
   // the standard pairs and what we have at *open.
   for (size_t ipar = 0; close == 0 && ipar < s_nParenthesisTypes; ++ ipar) {
      if (s_ParenthesisPairs[ipar].Open == open) {
         // found it!
         close = s_ParenthesisPairs[ipar].Close;
         break;
      }
   }

   if (close == 0)
      // Still no closing parenthesis type? That means we failed to
      // identify what to search for and cannot proceed --> fail.
      return this->first;
   _CxAssert(close != open);

   // okay, we got what we need. We now just pass control to the free
   // pchar version of FindClosingParenthesis(), and convert back whatever it
   // produces into `cstr_it`s.
   {
      char const
         *pFirst = &*first,
         *pLast = &*last,
         *pClosePar = FindClosingParenthesis(open, close, pFirst, pLast);
      // We need some conversion beyond the types:
      // (a) ...the pchar routine returns pLast if not no closing parenthesis is
      //     found, rather than pFirst as we do here
      // (b) ...if it does find the closing parenthesis, it returns an iterator
      //     to the closing parenthesis *itself*, rather than the character after
      _CxAssert(pClosePar == pLast || *pClosePar == close);
      if (pClosePar == pLast || *pClosePar != close)
         return first; // no close parens --> fail
      // convert pLast to cstr_it, via iterator differences.
      cstr_it itClosePar = first + (pLast - pFirst);
      _CxAssert(*itClosePar == close);
      return itClosePar;
   }
}


cstr_it string_slice::find_closing_parenthesis(char const *pOpen, char const *pClose) const {
//    char_traits::length
   char const
      *pFirst = &*first,
      *pLast = &*last,
      *pClosePar = FindClosingParenthesis(pOpen, pClose, pFirst, pLast);
   // We need some conversion beyond the types:
   // (a) ...the pchar routine returns pLast if not no closing parenthesis is
   //     found, rather than pFirst as we do here
   // (b) ...if it does find the closing parenthesis, it returns an iterator
   //     to the closing parenthesis *itself*, rather than the character after
   if (pClosePar == pLast)
      return first; // no close parens --> fail
   // convert pLast to cstr_it, via iterator differences.
   cstr_it itClosePar = first + (pLast - pFirst);
   _CxAssert(string_slice(itClosePar,last).startswith(pClose));
   return itClosePar;
}



bool string_slice::operator==(char const *p) const
{
   for (cstr_it it = first; it != last; ++it, ++p) {
      if (*p == 0)
         // reference string ended before string covered in *this.
         return false;
      if (*p != *it)
         return false;
   }
   return *p == 0;
}


bool string_slice::is_equal_except_case(char const *p) const
{
   for (cstr_it it = first; it != last; ++it, ++p) {
      if (*p == 0)
         // reference string ended before string covered in *this.
         return false;
      if (std::tolower(*p) != std::tolower(*it))
         return false;
   }
   return *p == 0;
}



string_slice::FPair string_slice::rsplit1(char Delim, size_t Flags) const
{
   long_split_result
      sr;
   this->split(sr, Delim, Flags);
   string_slice
      nothing(last,last);
   if (sr.empty())
      return FPair(nothing,nothing);
   else if (sr.size() == 1)
      return FPair(sr[0], nothing);
   else {
      _CxAssert(sr.size() >= 2);
      return FPair(string_slice(sr.front().first, sr[sr.size()-2].last), sr.back());
   }
}


string_slice::FPair string_slice::split1(char Delim, size_t Flags) const
{
   long_split_result
      sr;
   this->split(sr, Delim, Flags);
   string_slice
      nothing(last,last);
   if (sr.empty())
      return FPair(nothing,nothing);
   else if (sr.size() == 1)
      return FPair(sr[0], nothing);
   else {
      _CxAssert(sr.size() >= 2);
      return FPair(sr.front(), string_slice(sr[1].first, sr.back().last));
   }
}






void SplitLinesAndRemoveComments(FLineList &out, std::string const &Text)
{
   size_t
      nLines1 = 0;
   for ( cstr_it it = Text.begin(); it != Text.end(); ++ it )
      if ( *it == '\n' ) nLines1 += 1;
   out.reserve(nLines1);

   size_t
      iLineBegin,
      iLineEnd = 0;
   size_t
      iLine = 0;
   string_slice
      Line;
   for ( ; iLineEnd < Text.size() && iLineEnd != std::string::npos; ) {
      ++ iLine;
      iLineBegin = iLineEnd;
      while (iLineBegin < Text.size() && Text[iLineBegin] == '\n')
         ++ iLineBegin;
      iLineEnd = Text.find('\n', iLineBegin);
      if (iLineEnd != std::string::npos)
         Line = string_slice(Text.begin() + iLineBegin, Text.begin() + iLineEnd);
      else
         Line = string_slice(Text.begin() + iLineBegin, Text.end());
      cstr_it
         iLineContentEnd = Line.find( "//" );
      if (iLineContentEnd != Line.last)
         Line.last = iLineContentEnd;
      while (Line.last > Line.first && is_whitespace_cxp1(*(Line.last-1)))
         -- Line.last;
      if (Line.empty())
         continue;
      out.push_back(FLine(Line,FCodeLoc(iLine)));
   }
}

std::ostream &operator << (std::ostream &out, FLine const &Line){
   out << "\nin line " << Line.Loc.iLine << ": '" << Line.Text << "'.";
   return out;
}


static bool _PropertyListFlagsValidQ(unsigned Flags) {
   if (bool(Flags & PROPLIST_IdentifierWithOptionalArgs) && !bool(Flags & PROPLIST_NameRequired))
      return false;
   if ((bool(Flags & PROPLIST_NameRequired) && bool(Flags & PROPLIST_UnnamedList)))
      return false;
   return true;
}


static void _AssertPropertyListFlagsAreValid(unsigned Flags) {
   if (!_PropertyListFlagsValidQ(Flags))
//       throw std::logic_error(fmt::format("CxParse1: invalid combination of property list flags (Flags = {:04x}).", Flags));
      throw std::runtime_error(fmt::format("CxParse1: invalid combination of property list flags (Flags = {:04x}).", Flags));
}


// return Flags_ extended by other flags which are implied by deeper ones, if necessary,
// and perform some basic consistency checks (e.g., not specifying PROPLIST_UnnamedList and
// PROPLIST_NameRequired at the same time, etc.)
static unsigned _SanitizedPropertyListFlags(unsigned Flags_) {
   unsigned Flags = Flags_;
   if (bool(Flags & PROPLIST_IdentifierWithOptionalArgs))
      Flags |= PROPLIST_NameRequired;
   _AssertPropertyListFlagsAreValid(Flags);
   return Flags;
}



FPropertyListStr::FPropertyListStr(string_slice const &Input, char const *pAssumedName, unsigned Flags, char const *pOpenDelimCloseDef)
{
//    _CxAssert(!(bool(Flags & PROPLIST_IdentifierWithOptionalArgs) && bool(Flags & PROPLIST_UnnamedList)));
//    _CxAssert(!(bool(Flags & PROPLIST_NameRequired) && bool(Flags && PROPLIST_UnnamedList)));
   Flags = _SanitizedPropertyListFlags(Flags);
   static std::string EmptyString = "";
   char
      cOpen = '{',
      cDelim = ';',
      cClose = '}',
      cDef = ':';
   if (pOpenDelimCloseDef != 0) {
      cOpen = pOpenDelimCloseDef[0];
      cDelim = pOpenDelimCloseDef[1];
      cClose = pOpenDelimCloseDef[2];
      if (pOpenDelimCloseDef[3] != 0)
         cDef = pOpenDelimCloseDef[3];
      _CxAssert(cDelim != cDef);
      _CxAssert(cOpen != cClose); // <- otherwise the error checks below must be adjusted (iPos+1, Input.last-1 etc.)
   }
   bool
      AllowOneEntryWithoutName = (0 != (Flags & PROPLIST_AllowOneEntryWithoutName)),
      AllowEntriesWithoutName = AllowOneEntryWithoutName || (0 != (Flags & PROPLIST_AllowEntriesWithoutName));

   string_slice
      Data;
   if (!bool(Flags & PROPLIST_UnnamedList)) {
      // this is a (potentially) named list. Find name and opening parenthesis
      // to designate start of data.
      if (bool(Flags & PROPLIST_IdentifierWithOptionalArgs)) {
         // pursue option to make a named "list" without { } provided the
         // name is itself a well-defined identifier. Useful if adding
         // actual options is optional. See also looks_like_property_list() function.
         FIdentifierDecl
            IdentifierDecl(FIdentifierDecl::AllowDigits | FIdentifierDecl::AllowMinus);
         cstr_it
            itNameEnd = FindIdentifierEnd(Input.first, Input.last, IdentifierDecl);
         Name = string_slice(Input.first, itNameEnd);
         Data = string_slice(itNameEnd, Input.last);
         Name.trim();
         if (Name.empty()) {
            // empty list name. should we do something about that?
         }
         if (Data.empty()) {
            // empty property list. no braces, either. And that's okay in this mode.
         } else {
            if (!(Data.size() >= 2 && Data.startswith(cOpen) && Data.endswith(cClose)))
               throw std::runtime_error(fmt::format("expected 'name{}..{},..{}' property list, but found '{}'.", cOpen, cDelim, cClose, Input.to_str()));
            // skip open and close braces.
            Data.first += 1;
            Data.last -= 1;
         }
      } else {
         cstr_it iPos = Input.find(cOpen);
         if (!Input.endswith(cClose) || iPos == Input.last)
            throw std::runtime_error(fmt::format("expected 'name{}..{},..{}' property list, but found '{}'.", cOpen, cDelim, cClose, Input.to_str()));
         Name = string_slice(Input.first, iPos);
         Name.trim();

         Data = string_slice(iPos+1,Input.last-1);
      }
      if (bool(Flags & PROPLIST_NameRequired) && Name.empty()) {
         throw std::runtime_error(fmt::format("expected 'name{}..{},..{}' property list, with mandatory name, but found '{}'.", cOpen, cDelim, cClose, Input.to_str()));
      }
   } else {
      // this is a definitely unnamed list.
      if (Input.startswith(cOpen) && Input.endswith(cClose)) {
         // unnamed list starts/ends with cOpen/cClose. So entry data
         // is in the middle---skip parenthesis.
         Data = string_slice(Input.first+1, Input.last-1);
      } else {
         // unnamed list does NOT start/end with cOpen/cClose.
         // Were we told that is okay?
         if (0 != (Flags & PROPLIST_AllowOmitOpenClose)) {
            // yes. so assume the entire input consists of straight-forward data
            // specifications only
            Data = Input;
            // (^- note that this *can* lead to problems if the
            // individual list entries would be allowed to start/end with
            // cOpen/cClose, too. so do not use the flag in this case!)
         } else {
            // omission of outer {...} is not allowed. So this is an input error.
            throw std::runtime_error(fmt::format("expected '{}..{},..{}' unnamed property list, but found '{}'.", cOpen, cDelim, cClose, Input.to_str()));
         }
      }
   }

   unsigned
      FindFlags = 0;
   if (bool(Flags & PROPLIST_AllowComplexKeys)) {
      FindFlags |= string_slice::SPLIT_SkipParenthesis;
   }

   long_split_result
      sr;
   Data.split(sr, cDelim, string_slice::SPLIT_SkipWhitespace | string_slice::SPLIT_SkipParenthesis | string_slice::SPLIT_SkipQuotes);
   for (size_t i = 0; i < sr.size(); ++ i) {
      if (!sr[i].empty()) { // <- empty property list is allowed.
         push_back(FPropertyStr());
         FPropertyStr
            &out = back();
         cstr_it
            iPos = sr[i].find(cDef, 0, FindFlags);
         if (iPos == sr[i].end()) {
            // store as a name-less entry...
            if (!AllowEntriesWithoutName)
               throw std::runtime_error(fmt::format("expected '{}' in property declaration of '{}'.", cDef, sr[i].to_str()));
            out.Name = string_slice(EmptyString);
            out.Content = sr[i];
            if (AllowOneEntryWithoutName)
               // e.g., where entries without name would set the default value.
               // in this case we would allow at most one of those in a given specification.
               AllowEntriesWithoutName = false;
         } else {
            out.Name = string_slice(sr[i].first, iPos);
            out.Name.trim();
            out.Content = string_slice(iPos+1, sr[i].last);
            out.Content.trim();
         }
      }
   }

   if (pAssumedName && Name != pAssumedName) {
      std::stringstream str;
      str << "expected '" << pAssumedName << "' but found '" << Name << "'.";
      throw std::runtime_error(str.str());
   }
}

} // namespace ct


// c/p'd from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <string.h> // for strlen

namespace ct {

// // trim from left
// void TrimLeft(std::string &s) {
//    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
// }

// // trim from right
// void TrimRight(std::string &s) {
//    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
// }
void TrimLeft(std::string &s) {
   s.erase(
      s.begin(),
      std::find_if(s.begin(), s.end(),
         [](unsigned char c) { return !std::isspace(c); })
   );
}

void TrimRight(std::string &s) {
   s.erase(
      std::find_if(s.rbegin(), s.rend(),
         [](unsigned char c) { return !std::isspace(c); }
      ).base(),
      s.end()
   );
}


// trim from both ends
void Trim(std::string &s) {
   TrimLeft(s);
   TrimRight(s);
}

bool StartsWith(std::string const &s, std::string const &prefix) {
//    return string_slice(s).startswith(string_slice(prefix));
   return 0 == s.compare(0, prefix.size(), prefix);
}

bool StartsWith(std::string const &s, char const *prefix) {
   return string_slice(s).startswith(prefix);
//    return 0 == s.compare(0, strlen(prefix), prefix);
}

bool StartsWith(std::string const &s, char prefix) {
   if (s.empty()) return false;
   return s[0] == prefix;
}

bool EndsWith(std::string const &s, std::string const &suffix)
{
   return string_slice(s).endswith(string_slice(suffix));
//    if (s.empty() || suffix.size() > s.size()) return false;
//    return 0 == s.compare(s.size() - suffix.size(), suffix.size(), suffix);

}

bool EndsWith(std::string const &s, char const *suffix)
{
   return string_slice(s).endswith(suffix);
//    size_t suffixlen = strlen(suffix);
//    if (s.empty() || suffixlen > s.size()) return false;
//    return 0 == s.compare(s.size() - suffixlen, suffixlen, suffix);
}


void StripLineComment(std::string &s, char const *prefix)
{
   size_t iPos = s.find(prefix);
   if (iPos == std::string::npos)
      return; // not in there.
   // strip off everything after the comment indicator, and then whitespace to the right.
   s.resize(iPos);
   TrimRight(s);
}


// makes a string lower case. Will break with localized multibyte characters
std::string tolower(std::string const &in)
{
   std::string
      Result(in);
   std::string::iterator
      it;
   for (it = Result.begin(); it != Result.end(); ++ it)
      *it = std::tolower(*it);
   return Result;
}


std::string toupper(std::string const &in)
{
   std::string
      Result(in);
   std::string::iterator
      it;
   for (it = Result.begin(); it != Result.end(); ++ it)
      *it = std::toupper(*it);
   return Result;
}


// returns string 'in' with all spaces (but not tabs or newlines or other whitespace characters) removed.
std::string stripwhitespace(std::string const &in)
{
   std::string
      Result;
   Result.reserve(in.size());
   for (size_t i = 0; i < in.size(); ++ i)
      if (' ' != in[i])
         Result.push_back(in[i]);
   return Result;
}


bool IsEqual_CaseInsensitive(std::string const &a, std::string const &b)
{
   if (a.size() != b.size())
      return false;
   for (size_t i = 0; i < a.size(); ++i)
      if (std::tolower(a[i]) != std::tolower(b[i]))
         return false;
   return true;
}

bool IsEqual_CaseInsensitive(std::string const &a, char const *b)
{
   for (size_t i = 0; i < a.size(); ++i) {
      if (b[i] == 0)
         // that means b terminated before a did. so strs have different length,
         // and cannot be equal.
         return false;
      if (std::tolower(a[i]) != std::tolower(b[i]))
         return false;
   }
   return true;
}

bool IsEqual_CaseInsensitive(char const *a, std::string const &b)
{
   return IsEqual_CaseInsensitive(b, a);
}


// FRuntimeError::FRuntimeError(std::string const &Reason)
//    : FBase(Reason)
// {}
//
// FUnsupportedError::FUnsupportedError(std::string const &Reason)
//    : FBase("Unsupported combination of options/methods: " + Reason)
// {}



} // namespace ct
