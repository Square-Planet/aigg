#ifndef CX_PACKED_INDEX_H
#define CX_PACKED_INDEX_H

#include <type_traits>
#include <array>
#include <utility>
// #include <functional> // for std::plus
// #include <algorithm> // for std::max
#include <cstdint>

#include "CxDefs.h" // for assert
#include "CxIterTools.h" // for template index lists (index_sequence_t, index_sequence_for_t)


namespace ct {


// helper class to select a native unsigned integer type with a given minimum number of bits
template<size_t nBits> struct TSelectType_FastUnsignedInt_MinNumBits {
   static_assert(nBits <= 64, "sorry, don't know an unsigned integer capable of handling this many bits. there may be implementation-specific versions, though, or some type emulation in boost (e.g., <boost/multiprecision/cpp_bin_float.hpp>, although I do not know how fast that is).");
   typedef
      typename std::conditional<nBits <= 8, std::uint_fast8_t,
      typename std::conditional<nBits <= 16, std::uint_fast16_t,
      typename std::conditional<nBits <= 32, std::uint_fast32_t,
      typename std::conditional<nBits <= 64, std::uint_fast64_t,
      void >::type>::type>::type>::type
      // ^--- there is also `std::conditional_t`, which would allow omitting the
      // extra `typename's and `::type`s, but it was only introduced in C++17.
      type;
};



// e.g., <unsigned short, 4>
template<class FEntryType, unsigned nEntriesT, unsigned nBitsPerEntryT = sizeof(FEntryType)*8>
struct TPackedIndexN {
   static_assert(!std::is_signed<FEntryType>::value, "TPackedIndexN<T,N,W> is only for packed storage of unsigned integers.");

   enum {
      nEntries = nEntriesT,
      nBits = nBitsPerEntryT
   };
   typedef FEntryType
      value_type;
   typedef typename TSelectType_FastUnsignedInt_MinNumBits<nBits * nEntries>::type
      // select a native unsigned integer capable of representing at
      // least `nGen*nBits` bits total
      FPackedType;
   typedef std::array<value_type, nEntries>
      FIndexN;

   TPackedIndexN() : m_PackedData(0) {}
   // construct from index array
   template<class T>
   TPackedIndexN(std::array<T,nEntries> const &IndexN) : m_PackedData(_PackIndexN(IndexN)) {}

   // construct from brace-enclosed initializer list
   template<class T>
   TPackedIndexN(std::initializer_list<T> vl) { _InitFromInitList(vl); }

   // construct from function arguments; one for each index
   template<class... Args>
   explicit TPackedIndexN(Args ...args) : m_PackedData(_PackIndex(args...)) {
      static_assert(sizeof...(args) <= nEntries, "cannot use more indices in initialization than fixed for template");
   }

   FIndexN Unpacked() const { return _UnpackIndex<FIndexN>(m_PackedData); }
//    operator FIndexN () const { return _UnpackIndex(m_PackedData); }
   // …and the same thing as an implicit conversion operator to `std::array<T,nEntriesT>`, too.
   template<class T>
   operator std::array<T,nEntriesT> () const { return _UnpackIndex<std::array<T,nEntriesT> >(m_PackedData); }

   // compute and return the i`th element from the packed data
   value_type operator [] (unsigned i) const {
      _CxAssert(i <= size());
      return static_cast<value_type>((m_PackedData >> (i*nBits)) & EntryMask);
   }

   bool operator == (TPackedIndexN const &other) const { return this->m_PackedData == other.m_PackedData; }
   bool operator != (TPackedIndexN const &other) const { return this->m_PackedData != other.m_PackedData; }
   bool operator < (TPackedIndexN const &other) const { return this->m_PackedData < other.m_PackedData; }

   // element-wise addition of sub-indices.
   TPackedIndexN operator + (TPackedIndexN const &other) const {
      for (size_t i = 0; i < size(); ++ i)
         _CxAssert(_index1_in_range((*this)[i] + other[i]));
      return TPackedIndexN(this->m_PackedData + other.m_PackedData);
   }

   size_t constexpr size() const { return nEntries; }

   // add a minimal iterator interface to allow use in range-based for and generic container printing mechanism.
   struct const_iterator {
      TPackedIndexN const *pOwner;
      unsigned iArg;
      const_iterator &operator ++() { ++iArg; return *this; }
      bool operator != (const_iterator const& other) const { return this->pOwner != other.pOwner || this->iArg != other.iArg; }
      TPackedIndexN::value_type operator *() const { return (*pOwner)[iArg]; }
   };
   const_iterator begin() const { return const_iterator{this, 0}; }
   const_iterator end() const { return const_iterator{this, nEntries}; }

public:
   // this one is for compatibility with Eigen in a special application.
   static TPackedIndexN Zero() { return TPackedIndexN(FPackedType(0)); }

   FEntryType sum() const { return this->FoldLeft(FEntryType(0), [](FEntryType a, FEntryType b) {return a + b;}); }
   FEntryType min() const { return this->FoldLeft(FEntryType(-1), [](FEntryType a, FEntryType b) {return a <= b ? a : b;}); }
   FEntryType max() const { return this->FoldLeft(FEntryType(0), [](FEntryType a, FEntryType b) {return a >= b ? a : b;}); }

   // compute  Op(…Op(Op(StartValue, m[0]), m[1]), m[2], …)
   // Note: the two arguments of the operator need not be the same type.
   // e.g., `StartValue` as bool, which is then accumulated, is fine in principle.
   template<class FResultType, class FBinaryOp>
   FResultType FoldLeft(FResultType StartValue, FBinaryOp const &BinaryOp) const {
      for (unsigned i = 0; i < nEntries; ++ i)
         StartValue = BinaryOp(StartValue, (*this)[i]);
      return StartValue;
   }

   FPackedType _PackedData() const { return m_PackedData; }
protected:
   static constexpr FPackedType
      EntryMask = (FPackedType(1) << nBits) - 1;
   FPackedType
      m_PackedData;
   explicit TPackedIndexN(FPackedType PackedData_) : m_PackedData(PackedData_) {}

   static bool _index1_in_range(value_type index) { return size_t(index) < (size_t(1)<<nBits); }

   // implementation for initialization from brace-enclosed initializer list
   template<class T>
   void _InitFromInitList(std::initializer_list<T> InitList) {
      _CxAssert(InitList.size() == nEntries);
      m_PackedData = 0;
      size_t iArg = 0;
      for (T const &Arg : InitList) {
         if (iArg < nEntries) {
            m_PackedData += FPackedType(Arg) << (nBits * iArg);
            iArg += 1;
         }
      }
   }

   // implementation for initialization from function argument pack
   template<unsigned iArg = 0>
   static FPackedType _PackIndex(value_type Arg) {
      static_assert(iArg < nEntries, "exceeded fixed number of sub-indices");
      _CxAssert(_index1_in_range(Arg));
      return FPackedType(Arg) << (nBits * iArg);
   }

   template<unsigned iArg = 0, class... FOtherArgs>
   static FPackedType _PackIndex(value_type Arg, FOtherArgs ...OtherArgs) {
      return _PackIndex<iArg>(Arg) | _PackIndex<iArg+1>(OtherArgs...);
   }


   // implementation for initialization from static-size array/tuple.
#if __cplusplus >= 201703L
   template<unsigned iArg = nEntries-1>
   FPackedType _PackIndexN(FIndexN const &ArgN) {
      using std::get;
      FPackedType tmp = _PackIndex<iArg>(get<iArg>(ArgN));
      if constexpr (iArg > 0)  // <-- if constexpr requires >= C++17. Certainly makes these kinds of things less cumbersome.
         tmp += _PackIndexN<iArg-1>(ArgN);
      return tmp;
   }
#else
   template<unsigned iArg> struct _TPackIndexN_Helper {
      struct Return{}; struct Recurse{}; // <-- no partial template specialization outside namespaces
      static FPackedType ReturnOrRecurse(Recurse, FPackedType Acc, FIndexN const &ArgN) { return Acc | _TPackIndexN_Helper<iArg-1>::Eval(ArgN); }
      static FPackedType ReturnOrRecurse(Return, FPackedType Acc, FIndexN const &) { return Acc; }
      static FPackedType Eval(FIndexN const &ArgN) {
         using std::get;
         typedef typename std::conditional<iArg == 0, Return, Recurse>::type FDecider;
         return ReturnOrRecurse(FDecider(), _PackIndex<iArg>(get<iArg>(ArgN)), ArgN);
      }
   };
   static FPackedType _PackIndexN(FIndexN const &ArgN) {
      return _TPackIndexN_Helper<nEntries-1>::Eval(ArgN);
   }
#endif

   // implementation into unpacking into array-like objects.
   template<class FResultType, cx::index_t iArg>
   static FResultType _UnpackIndex1(FPackedType pi) {
      return (pi >> (nBits*iArg)) & EntryMask;
      // ^-- one could leave off the mask for the last entry and the iArg for the first.
   }

   template<class FArrayLike, cx::index_t... Is>
   static FArrayLike _UnpackIndexImpl(FPackedType pi, cx::index_list_t<Is...>) {
      return {_UnpackIndex1<typename FArrayLike::value_type, Is>(pi)...};
   }

   template<class FArrayLike>
   static FArrayLike _UnpackIndex(FPackedType pi) {
      return _UnpackIndexImpl<FArrayLike>(pi, cx::make_index_list<nEntries>());
   }
};

// #pragma pack(pop)



} // namespace ct

// #if __cplusplus >= 201703L
// define get<n>(e), tuple_size<>, tuple_element<> to provide support for C++17
// 'structured binding'. I guess it doesn't hurt to define the tuple-interface
// stuff even if there is no structured binding. I think some other algorithms
// may also use them.
namespace std {
   template<class T, unsigned N, unsigned W>
   struct tuple_size<ct::TPackedIndexN<T,N,W> >{ enum { value = N }; };

   template<size_t iElemT, class T, unsigned N, unsigned W>
   struct tuple_element<iElemT, ct::TPackedIndexN<T,N,W> >{ typedef typename ct::TPackedIndexN<T,N,W>::value_type type; };
}

namespace ct {
   template<size_t iElemT, class T, unsigned N, unsigned W>
   typename TPackedIndexN<T,N,W>::value_type get(TPackedIndexN<T,N,W> const &pi) {
      return pi[iElemT];
   }
}
// #endif


#ifdef CX_SEQUENCE_IO_H
namespace ct {
   template<class T, unsigned N, unsigned W>
   std::ostream &operator << (std::ostream &out, ct::TPackedIndexN<T,N,W> const &I) {
      return ct::sequence_io::TSequenceWriter({"[",",","]"}).WriteList(out, I);
   }
} // namespace ct
#endif // CX_SEQUENCE_IO_H

#endif // CX_PACKED_INDEX_H

