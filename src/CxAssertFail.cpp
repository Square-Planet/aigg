// This is a reference implementation of the "CxAssertFail" function used in
// CxDefs.h's assertion functions. It is only used to make the distribution
// self-contained: a function with signature
//
//    CxAssertFail(char const *pExpr, char const *pFile, int iLine)
//
// needs to be *somewhere* in the program, but it does not have to be this one.
// If you have a different way of handling assertion failures, you need
// not use/compile in this file, but can instead just substitute your own
// implementation.
//
// UPDATE:
// - Some other error raising function implementations also made their way here...

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <algorithm> // for std::swap
#include "CxDefs.h"
#ifdef INCLUDE_ABANDONED
// #include "format.h"
// // ^-- hm... format.h was not required without the extra _RaiseMemoryAllocError
// // stuff. On the other hand, I probably really don't want to make a program
// // without this anyway.
// // hmpf... whatever. Let's remove it.
#endif // INCLUDE_ABANDONED
using std::size_t;


void /*CX_NORETURN*/ CxAssertFail(char const *pExpr, char const *pFile, int iLine)
{
   std::stringstream
      str;
   if (pFile) {
      str << pFile << ":" << iLine << ": ";
   }
   str << "Assertion Failed: '" << pExpr << "'";
   std::cerr << str.str() << std::endl;
#ifdef _DEBUG
   DEBUG_BREAK
#endif // _DEBUG
   throw std::runtime_error(str.str());
}

#ifdef INCLUDE_ABANDONED
//          std::size_t len = 0;
//          for (; msg[len] != 0; ++ len) {
//          }
//          len = std::char_traits<char>::length(string);
#endif // INCLUDE_ABANDONED
namespace ct {
   cstr_buf_t::cstr_buf_t(char const *msg) CX_NOEXCEPT {
      _Init(); if (msg) { _Reset(msg, std::char_traits<char>::length(msg)); }
   }
   cstr_buf_t::cstr_buf_t(char const *msg, std::size_t len) CX_NOEXCEPT { _Init(); _Reset(msg, len); }
   cstr_buf_t::cstr_buf_t(cstr_buf_t const &other) CX_NOEXCEPT : m_str(0), m_len(0) { _Init(); _Copy(other); }
   cstr_buf_t &cstr_buf_t::operator=(cstr_buf_t const &other) CX_NOEXCEPT { _Init(); _Copy(other); return *this; }
   cstr_buf_t::cstr_buf_t(cstr_buf_t::TakeOverTag, char const *msg) CX_NOEXCEPT {
      _Init(); if (msg) { m_str = const_cast<char *>(msg); m_len = std::char_traits<char>::length(msg); };
   }
   cstr_buf_t::cstr_buf_t(TakeOverTag, char const *msg, std::size_t len) CX_NOEXCEPT {
      _Init(); if (bool(msg) && len > 0) { m_str = const_cast<char *>(msg); m_len = len; /*_CxAssert(len >= 1 && m_str[len] == 0);*/ }
   }
   #if __cplusplus >= 201103L
   cstr_buf_t::cstr_buf_t(cstr_buf_t &&other) CX_NOEXCEPT { _Init(); _Swap(other); }
   cstr_buf_t &cstr_buf_t::operator=(cstr_buf_t &&other) CX_NOEXCEPT { _Init(); _Swap(other); return *this; }
   #endif
   cstr_buf_t::~cstr_buf_t() { _Release(); }
   void swap(cstr_buf_t &a, cstr_buf_t &b) CX_NOEXCEPT { a._Swap(b); }

   void cstr_buf_t::_Init() {
      m_str = 0;
      m_len = 0;
   }
   void cstr_buf_t::_Release() {
      std::free(m_str);
      m_str = 0;
      m_len = 0;
   }
   void cstr_buf_t::_Reset(char const *str, std::size_t len) {
      _Release();
      m_len = len;
      m_str = static_cast<char*>(std::malloc(m_len+1));
      for (std::size_t i = 0; i < m_len; ++i)
         m_str[i] = str[i];
      m_str[m_len] = 0; // add 0-terminator.
   }
   void cstr_buf_t::_Swap(cstr_buf_t &other) {
      std::swap(m_str, other.m_str);
      std::swap(m_len, other.m_len);
   }
   void cstr_buf_t::_Copy(cstr_buf_t const &other) {
      _Reset(other.m_str, other.m_len);
   }
}



namespace ct {
   double _fSizeMb(size_t nBytes) {
      return double(nBytes)/double(size_t(1u) << size_t(20u));
   }

#ifdef INCLUDE_ABANDONED
//    void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount) {
//       throw cx_alloc_fail_t(fmt::format(
//             "{}: memory allocation failed ({:.2f} MB requested; {} objects of size {}).",
//             pFnName, _fSizeMb(ValueTypeSize * NewReserveCount), NewReserveCount, ValueTypeSize
//          ).c_str()
//       );
//    }
//
//
//    void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount, size_t PrevReservedCount) {
//       throw cx_alloc_fail_t(fmt::format(
//             "{}: memory allocation failed ({:.2f} MB requested; {} objects of size {}; reallocation: {}; prev size: {:.2f} MB).",
//             pFnName, _fSizeMb(ValueTypeSize * NewReserveCount), NewReserveCount, ValueTypeSize,
//             (PrevReservedCount != 0)? "yes" : "no", _fSizeMb(ValueTypeSize * PrevReservedCount)
//          ).c_str()
//       );
//    }
#endif // INCLUDE_ABANDONED
   static const size_t _MagicAbsence = static_cast<size_t>(0xbadc0de);

   void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount, size_t PrevReservedCount) {
      std::stringstream ss;
      ss << pFnName << ": " << "memory allocation failed ("
         << _fSizeMb(ValueTypeSize * NewReserveCount) << " MB requested"
         << "; " << NewReserveCount << " objects of size " << ValueTypeSize;
      if (PrevReservedCount != _MagicAbsence) { ss
         << "; " << "reallocation: " << ((PrevReservedCount != 0)? "yes" : "no")
         << "; " << "prev size: " << _fSizeMb(ValueTypeSize * PrevReservedCount) << " MB";
      }
      ss << ").";
      throw cx_alloc_fail_t(ss.str().c_str());
   }

   void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount) {
      _RaiseMemoryAllocError(pFnName, ValueTypeSize, NewReserveCount, _MagicAbsence);
   }
}



#include <algorithm> // for std::stable_sort

// well... this one *most definitely* should not be here. It is declared in
// CxTypes.h, but that one does not have an implementation file, and it doesn't
// really fit anywhere else, either.
namespace ct {

   struct FPredArgSort {
      FPredArgSort(double const *pVals_, size_t nValSt_, bool Reverse_) : pVals(pVals_), nValSt(nValSt_), Reverse(Reverse_) {};
      bool operator () (size_t i, size_t j) const { return Reverse? (iVal(i) > iVal(j)) : (iVal(i) < iVal(j));  }
   private:
      double iVal(size_t i) const { return pVals[nValSt * i]; }
      double const *pVals;
      size_t nValSt;
      bool Reverse;
   };

   // find permutation which sorts values in pVals[nValSt*i].
   void ArgSort1(size_t *pOrd, double const *pVals, size_t nValSt, size_t nVals, bool Reverse) {
      for (size_t i = 0; i < nVals; ++ i)
         pOrd[i] = i;
      std::stable_sort(pOrd, pOrd+nVals, FPredArgSort(pVals, nValSt, Reverse));
   }

} // namespace ct
