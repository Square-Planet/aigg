#ifndef CX_PODARRAY_H
#define CX_PODARRAY_H

#if __cplusplus >= 201103L
   #define CX_ARRAY_HAVE_CXX11
#endif

#include <new> // for placement new
#include <string.h> // for memcpy
#include <stdlib.h> // for malloc/free
#ifdef INCLUDE_ABANDONED
// #include <functional> // for std::less
// #include <algorithm> // for swap
#endif // INCLUDE_ABANDONED
#ifdef CX_ARRAY_HAVE_CXX11
   #include <utility> // for std::forward
   #include <initializer_list>
#endif
#include "CxDefs.h"



namespace ct {

#ifdef INCLUDE_ABANDONED
// #ifdef __GNUC__ // g++ or anything which pretends to be g++ (e.g., clang, intel c++ on linux)
//    // default number of static elements in TArray. Should be 0, but this is apparently
//    // of disputed legality in the C++ standard. VC++ doesn't like that ("illegal size 0 array")
//    #define CX_DEFAULT_ARRAY_STATIC_STORAGE 0
// #else
//    #define CX_DEFAULT_ARRAY_STATIC_STORAGE 8
// #endif
// ^- UPDATE: fixed now with the TStaticStorage workaround below.
#endif // INCLUDE_ABANDONED

#ifndef CX_DEFAULT_ARRAY_STATIC_STORAGE
   #define CX_DEFAULT_ARRAY_STATIC_STORAGE 0
#endif // CX_DEFAULT_ARRAY_STATIC_STORAGE



namespace detail {
   // workaround for this:
   //
   //    CxPodArray.h:285:7: warning: ISO C++ forbids zero-size array [-Wpedantic]
   //
   // VC apparently does not like this at all, either.
   // So instead of adding the zero-length array, we wrap the static storage
   // into a separate template class (TStaticStorage), and then do a partial
   // specialization for the case _StaticStorage == 0 does not have the zero-
   // size array member at all. I think this one should be legal C++ without any
   // caveats.

   template<class _ValueType, std::size_t _StaticSize, class _IndexType = std::size_t>
   struct TStaticStorage {
      _ValueType m_StaticData[_StaticSize];
      _ValueType *static_data() noexcept { return &this->m_StaticData[0]; }
      constexpr _ValueType const *static_data() const noexcept { return &this->m_StaticData[0]; }
      constexpr _IndexType static_capacity() const noexcept { return static_cast<_IndexType>(_StaticSize); }
   };

   template<class _ValueType, class _IndexType>
   struct TStaticStorage<_ValueType, std::size_t(0), _IndexType> {
      _ValueType *static_data() noexcept { return 0; }
      constexpr _ValueType const *static_data() const noexcept { return 0; }
      constexpr _IndexType static_capacity() const noexcept { return 0; }
   };
}

// raises std::bad_alloc
void _RaiseMemoryAllocError(char const *pFnName, size_t ValueTypeSize, size_t NewReserveCount, size_t PrevReservedCount);
double _fSizeMb(size_t nBytes);


/// A dynamic array of POD (``plain old data'') types that can be copied via a
/// memcpy and cleared with memset(...,0,..), and require neither constructors
/// no destructors to work.
///
/// - Main point for this is that std::vector does not allow containing C-style
///   arrays (e.g., double [3]), because C-style arrays are not assignable.
/// - Additionally, std::vector can be *very* slow when allocating large amounts
///   of data, because the data is set to zero on resize. This class explicitly
///   DOES NOT DO THAT: It has RAII semantics, but non-explicitly touched data
///   is just random. This behavior is ITENTIONAL! It makes a large difference
///   for large data blocks (in particular, it means that they need not be
///   committed by the OS before being written to).
/// - Also, this class neither calls constructors nor destructors for its data
///   elements, so this should most definitely NOT be used on classes which
///   would be broken due to this.
///
/// Template arguments:
/// - _StaticStorage: If '_StaticStorage' is zero, the class is
///   effectively a 'buffer-ptr + size' pair. Otherwise, the class contains
///   a local buffer for '_StaticStorage' elements which will be used
///   for allocation-less storage of the first few elements until a dynamic
///   allocation is required when the static storage size is exceeded.
///
/// TODO:
/// - The class has become complex enough to warrant splitting into a dumb
///   data implementation class (with char* ptrs) and an actual typed
///   class, I think. But it is a bit messy and dangerous, due to the
///   strict aliasing rules. Would need some actual time to think through.
///   Reconsider if I get really bored...
template<class _ValueType, std::size_t _StaticStorage = CX_DEFAULT_ARRAY_STATIC_STORAGE, class _IndexType = std::size_t>
struct TArray : protected detail::TStaticStorage<_ValueType, _StaticStorage, _IndexType>
{
   using value_type = _ValueType;
   using pointer = _ValueType*;
   using reference = _ValueType&;
   using iterator = _ValueType*;
   using const_pointer = _ValueType const*;
   using const_reference = _ValueType const&;
   using const_iterator = _ValueType const*;
   using index_type = _IndexType;
   using difference_type = std::ptrdiff_t;

   using value_t = _ValueType;
   using index_t = _IndexType;


   iterator begin() noexcept { return m_pData; }
   iterator end() noexcept { return m_pDataEnd; }
   const_iterator begin() const noexcept { return m_pData; }
   const_iterator end() const noexcept { return m_pDataEnd; }

   value_t &front() noexcept { return *m_pData; }
   value_t &back() noexcept { return *(m_pDataEnd-1); }
   value_t const &front() const noexcept { return *m_pData; }
   value_t const &back() const noexcept { return *(m_pDataEnd-1); }

   value_t *data() noexcept { return m_pData; }
   value_t const *data() const noexcept { return m_pData; }

   value_t &operator[] (index_t i) noexcept { return m_pData[i]; }
   value_t const &operator[] (index_t i) const noexcept { return m_pData[i]; }

   index_t size() const noexcept { return m_pDataEnd - m_pData; }
//    index_t capacity() const noexcept { return _CxMaxQv(index_t(m_pDataLimit - m_pData), index_t(_StaticStorage)); }
   index_t capacity() const noexcept { return static_cast<index_t>(m_pDataLimit - m_pData); }
   index_t max_size() const noexcept { return index_t(-1)/sizeof(value_t); }
   bool empty() const noexcept { return m_pData == m_pDataEnd; }

   /// WARNING: contrary to std::vector, this function *DOES NOT*
   /// initialize (or touch, for that matter) the newly created data.
   /// If you want that behavior, use the other resize function!
   /// It also does not call any destructors in case the number of
   /// element decreases.
   void resize(index_t n) { reserve(n); m_pDataEnd = m_pData + n; }

   void resize(index_t n, value_t t) {
      index_t old_size = size();
      resize(n);
      for (index_t i = old_size; i < n; ++ i)
         m_pData[i] = t;
   }

   // change the number of controlled elements to zero, but **do not**
   // release the array storage space they used
   void clear() noexcept { resize(0); }

   // memset the entire array data to 0.
   void clear_data() noexcept {
      memset(m_pData, 0, sizeof(value_t)*std::size_t(size()));
   }

   void resize_and_clear(index_t n) { resize(n); clear_data(); }

   // Clear array data and release all memory.
   //
   // Note: This frees all controlled data, but does not invalidate the array
   // object itself. Calling release() and then refilling an array is fine.
   void release() noexcept { _ReleaseBuffer(); }

   void push_back(value_t const &t) {
      _ReserveSpaceForPushBack();
      *m_pDataEnd = t;
      ++m_pDataEnd;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   template<class... Args>
   void emplace_back(Args&&... args) {
      _ReserveSpaceForPushBack();
      new(m_pDataEnd) value_t{std::forward<Args>(args)...};
      ++m_pDataEnd;
   }
#endif

   void pop_back() noexcept {
      _CxAssert(!empty());
      --m_pDataEnd;
   }

#ifdef INCLUDE_ABANDONED
//    value_t pop_back() noexcept {
//       _CxAssert(!empty());
//       m_pDataEnd -= 1;
//       return *m_pDataEnd;
//    }
//    ^-- doesn't work with arrays.
//    enum {
//       // if given as flag to reserve(), reserve may reduce the array's capacity;
//       // otherwise reserve(n) with n < capacity() are ignored.
//       // Capacity will never be reduced below the actual number of stored items.
//       RESERVE_AllowShrink = 0x0001
//    };
//
//    void reserve(index_t NewReservedSize, unsigned Flags = 0) {
//       if (bool(Flags & RESERVE_AllowShrink)) {
//          if (NewReservedSize < size())
//             NewReservedSize = size();
//       } else {
//          if (NewReservedSize == 0)
//             return;
//       }
#endif // INCLUDE_ABANDONED

   void reserve(index_t nCells) {
      if (m_pData == 0 && _StaticStorage != 0) {
         // No supposed to happen, I think. But I there might be cases where
         // `*this` was part of a structure which got summarily initialized via
         //
         //     memset(…, 0, size)...
         //
         // Now, allowing this is not necessarily a terrible idea. But if done,
         // this definitely requires extra care and some other checks
         // (e.g., if the object is properly destroyed, or also just `::free`d).
         //
         // In any case: assuming this *was* a `memset`, we can just re-init the
         // structure here, to establish the expected structure of internal data members.
         _CxAssert(m_pDataEnd == 0 && m_pDataLimit == 0);
         _CxAssert(!"TArray::reserve(): pData == 0, despite non-zero static storage size. "
               "Structure memset()-ed over? If yes, check other usage patterns...");
         _InitDataMembers();
      }

      // space already sufficient as-is?
      if (nCells <= capacity()) return;
      index_t nElem = size();
      _CxAssert(nElem <= capacity() && nElem < nCells);
      value_t *pNewData;
      if (_HeapArrayQ() && capacity() < 2*nElem) {
         // ^-- note: this test implicitly uses that `(capacity() < 2*nElem)`
         //     automatically implies `nElem > 0`

         // Attempt a reallocation of the present memory block if...
         // ...there are already data elements stored (nElem > 0),
         // ...current data is on dynamic memory,
         // ...and the currently reserved memory block size does not
         //    significantly exceed the memory size of the actual data.
         // The last condition is for the case in which the reallocation has to
         // resort to an alloc & copy -- in that case we'd be better off
         // memcpy'ing the only actual amount of data ourselves, rather than
         // risking realloc copying a large amount of unnecessary random data
         // (whatever the empty space [m_pDataEnd … m_pDataLimit] contains).
         pNewData = static_cast<value_t*>(::realloc(m_pData, sizeof(value_t) * nCells));
      } else {
         pNewData = static_cast<value_t*>(::malloc(sizeof(value_t) * nCells));
         if (pNewData != 0) {
            if (nElem != 0)
               ::memcpy(pNewData, m_pData, sizeof(value_t) * nElem);
            _ReleaseBuffer();
         }
      }
      if (pNewData == 0)
         return _RaiseMemoryAllocError(__func__, sizeof(value_t), nCells, std::size_t(m_pDataLimit - m_pData));
#ifdef INCLUDE_ABANDONED
//             throw std::runtime_error(fmt::format("TArray<{}>: memory allocation failed ({:.2f} MB requested; reallocation: {}; prev size: {:.2f} MB).",
//                sizeof(value_t), _fSizeMb(NewReservedSize), (m_pData != 0)? "yes" : "no", _fSizeMb(_nReserved())));
#endif // INCLUDE_ABANDONED
      m_pData = pNewData;
      m_pDataEnd = m_pData + nElem;
      m_pDataLimit = m_pData + nCells;
   }


   void shrink_to_fit(index_t nCells) {
      if (_StaticArrayQ()) return; // static buffer cannot be shrunken.
      index_t nElem = size();
      if (nCells < nElem)
         // Ensure capacity remains sufficient to store existing elements
         nCells = nElem;
      if (nElem < this->static_capacity()) {
         // The data is not on the static buffer currently, but with the current
         // level of filling, it could be. So move it there and completely
         // discard the dynamic allocation.
         _CxAssert(!_StaticArrayQ());
         ::memcpy(this->static_data(), m_pData, sizeof(value_t) * nElem);
         _ReleaseBuffer(); // <-- re-attaches data array to static buffer
         _CxAssert(m_pData == this->static_data() && m_pDataEnd == m_pData && m_pDataLimit == m_pData + this->static_capacity());
         m_pDataEnd  = m_pData + nElem;
      } else if (capacity() >= (nElem + (3u + nElem)/4u) || (capacity() - nElem) * sizeof(value_t) >= (size_t(16) << 10u)) {
         // fulfill request (via reallocation, if possible) only if either filling
         // level is < 75% of capacity or the shrink_to_fit would free up at least
         // 16 kiB
         _CxAssert(!_StaticArrayQ());
         value_t
            *pNewData = static_cast<value_t*>(::realloc(m_pData, sizeof(value_t) * nCells));
         if (pNewData == 0)
            return _RaiseMemoryAllocError(__func__, sizeof(value_t), nCells, std::size_t(m_pDataLimit - m_pData));
         m_pData = pNewData;
         m_pDataEnd = m_pData + nElem;
         m_pDataLimit = m_pData + nCells;
      }
   }


   TArray() noexcept {
      _InitDataMembers();
   }

   TArray(TArray const &other) {
      _InitDataMembers();
      *this = other;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray(TArray &&other) noexcept {
      _InitDataMembers();
      _MoveArrayData(other);
   }
#endif // CX_ARRAY_HAVE_CXX11

   // WARNING: this function **does not** initialize the array content! (with intention!)
   explicit TArray(index_t n) {
      _InitDataMembers();
      resize(n);
   }

   // initializes the the array to size `n`, with all `n` entries set to `t`.
   TArray(index_t n, value_t t) {
      _InitDataMembers();
      assign(n, t);
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray(std::initializer_list<value_t> init) noexcept {
      // @noexcept: since this type is meant for pods, aggregates, trivially
      // constructible types/etc, the only likely exception for this is a memory
      // allocation failure. And if we get the data from a initializer list,
      // there is little chance of it, say, requiring several GB of memory which
      // the machine does not have. Additionally, we'd near certainly have to
      // terminate the program in any case, unless this is one of the very few
      // placed where the code actually (1) expects large allocations, (2)
      // has viable alternatives, and (3) actually handles the allocation failure
      // in by using them.
      _InitDataMembers();
      assign(init.begin(), init.end());
   }
#endif // CX_ARRAY_HAVE_CXX11

   template<class _Iterator>
   TArray(_Iterator first, _Iterator last) {
      _InitDataMembers();
      assign(first, last);
   }

   ~TArray() noexcept {
      _ReleaseBuffer();
   }

   TArray const &operator = (TArray const &other) {
      if (this == &other) return *this;
      resize(other.size());
      memcpy(m_pData, other.m_pData, sizeof(value_t) * size());
      return *this;
   }

#ifdef CX_ARRAY_HAVE_CXX11
   TArray &operator = (TArray &&other) noexcept {
      _CxAssert(this != &other);
      // ^-- probably not supposed to happen. But opinions are somewhat divided:
      //     See answer of Howard Hinnant at:
      //     https://stackoverflow.com/questions/9322174/move-assignment-operator-and-if-this-rhs
      //     In any case, we do provide a formally correct implementation for that case just next.
      if (this == &other) {
         // this kills *this === other, but does leave it in a "valid but unspecified"
         // state. This is a valid outcome for an operation like `x = std::move(x)`.
         // Note that even with killing *this in this variant, the `swap(x,x)` operation
         // will still work out fine and yield the correct result of `x` staying as it
         // was in input.
         _ReleaseBuffer();
      } else {
         _InitDataMembers();
         _MoveArrayData(other);
      }
      return *this;
   }
#endif // CX_ARRAY_HAVE_CXX11

   void swap(TArray &other) noexcept {
      if (this == &other || (this->empty() && other.empty()))
         return;
      if (other._HeapArrayQ() && this->_HeapArrayQ()) {
         // both arrays are hosted on existing, heap-allocated dynamic buffers
         // --> some pointer swappery is sufficient.
         _CxSwap1v(m_pData, other.m_pData);
         _CxSwap1v(m_pDataEnd, other.m_pDataEnd);
         _CxSwap1v(m_pDataLimit, other.m_pDataLimit);
      } if (this->empty() && other._HeapArrayQ()) {
         // take over other arrays heap buffer (note: _MoveArrayData also clears out the source array)
         _MoveArrayData(other);
      } if (other.empty() && this->_HeapArrayQ()) {
         // other array is empty, but we have heap data. Move it over there,
         // (note: _MoveArrayData also clears out the source, here *this).
         other._MoveArrayData(*this);
      } else {
         // at least one of the array data sets is on static buffers. This makes
         // this complicated. We actually need to physically copy stuff around
         // in this case.
         TArray tmp(*this);
         *this = other;
         other = tmp;
      }
   }


   template<class _RandomIt>
   void assign(_RandomIt begin, _RandomIt end) {
      resize(end - begin);
      for (index_t i = 0; i < size(); ++ i)
         m_pData[i] = value_t(begin[i]);
   }

   // replaces content of *this with `n` copies of value `t`.
   void assign(index_t n, value_t t) {
      resize(n);
      for (index_t i = 0; i < n; ++ i)
         m_pData[i] = t;
   }


#ifdef CX_ARRAY_HAVE_CXX11
   void assign(std::initializer_list<value_t> init) noexcept {
      assign(init.begin(), init.end());
   }
#endif // CX_ARRAY_HAVE_CXX11


   template<class _RandomIt>
   void insert(iterator pos, _RandomIt begin, _RandomIt end) {
      // insert position must be within range of existing data values.
      _CxAssert(m_pData <= pos && pos <= m_pDataEnd);
      // convert everything from pointers to integers -- we will resize the
      // array, and this will invalidate the pointers.
      index_t
         n = end - begin,                  // insert-seq size
         N = this->end() - this->begin(),  // own array size
         ipos = pos - this->begin();       // own array index to place insert-seq at
      if (n == 0) return; // nothing to insert.
      resize(n + N);
      // new subsequence elements take space at `data[ipos:(ipos+n)]`
      // --> move all data elements currently at `data[ipos:N]` backwards
      //     by `n` data cells, ⟶ into data[(ipos+n):(N+n)].
      for (index_t i = N + n; i != ipos + n; ) { -- i;
         this->m_pData[i] = this->m_pData[i-n];
      }
      for (index_t i = 0; i < n; ++ i)
         this->m_pData[ipos+i] = value_t(begin[i]);
   }

   // first compare array size; then lexicographical comparison of array elements for
   // same-size arrays. Employs member type's operator <.
   bool operator < (TArray const &other) const noexcept {
      if (this->size() < other.size()) return true;
      if (other.size() < this->size()) return false;
      _CxAssert(this->size() == other.size());
      index_t N = size();
      for (index_t i = 0; i != N; ++ i) {
         value_t const
            &ai = this->m_pData[i],
            &bi = other.m_pData[i];
         if (ai < bi) return true;
         if (bi < ai) return false;
      }
      return false;
   }

   bool operator != (TArray const &other) const noexcept {
      if (other.size() != this->size()) return true;
      for (index_t i = 0; i != size(); ++ i)
         if (this->m_pData[i] != other.m_pData[i]) return true;
      return false;
   }

   bool operator == (TArray const &other) const noexcept {
      return !(*this != other);
   }
private:
   value_t
      // start of controlled array. may be 0 if no data is contained.
      *m_pData,
      // end of actual data
      *m_pDataEnd,
      // end of reserved space (i.e., of the space *this owns)
      *m_pDataLimit;

   // If `true`, the back-end array is (still) self-hosted within `m_StaticData`
   // (i.e., no external memory allocations so far).
   // If `false`, `m_pData` points to an external memory block, acquired by
   // heap allocation, and `*this` is responsible for freeing it.
   constexpr bool _StaticArrayQ() noexcept { return _StaticStorage != 0 && m_pData == this->static_data(); }
   // return `true` the back-end array exists and has been heap-allocated (as opposed to using internal static storage)
   constexpr bool _HeapArrayQ() noexcept { return m_pData != 0 && m_pData != this->static_data(); }

#ifdef INCLUDE_ABANDONED
//    value_t
//       // if the total number of elements is smaller than _StaticStorage, then
//       // this buffer may be used to store them.
//       m_StaticStorage[_StaticStorage];
#endif // INCLUDE_ABANDONED

   double _fSizeMb(size_t n) {
      return _fSizeMb(sizeof(value_t) * n);
   }

   void _ReleaseBuffer() {
      if (_HeapArrayQ()) ::free(m_pData);
      _InitDataMembers(); // <-- re-attach array pointers to static storage array (if we have one)
   }

   void _InitDataMembers() {
      // On init, attach our data pointers into the local static storage array
      // (so if there is a static storage, we'll start out with valid pointers
      // and a non-zero reservation size to begin with)
      m_pData = this->static_data();
      m_pDataEnd = m_pData;
      m_pDataLimit = m_pData + this->static_capacity();
   }

   void _MoveArrayData(TArray &other) {
      if (other._HeapArrayQ()) {
         // other array's data is on dynamic memory, not on static buffers. just
         // move around the pointers.
         _ReleaseBuffer();
         m_pData = other.m_pData;
         m_pDataEnd = other.m_pDataEnd;
         m_pDataLimit = other.m_pDataLimit;
         other._InitDataMembers();
      } else {
         // other array's data is on its static memory buffer---which we obviously
         // cannot take control of in here. So there is no way to avoid physically copying
         // that other array's data.
         assign(other.begin(), other.end());
         other.clear();
      }
   }


   // exponentially increases array size if space is not enough, for average logarithmic
   // time complexity of push_back/emplace_back
   inline void _ReserveSpaceForPushBack() {
      _CxAssert(m_pDataLimit >= m_pDataEnd && m_pDataEnd >= m_pData);
      if (m_pDataEnd == m_pDataLimit) {
         index_t n = (3u*this->size())/2u;
         // hm.. maybe should come up with more reasonable starting size and
         // groth factor. Note: *this itself requires at least the three 8-byte
         // pointers, so we probably should not make a heap array with less than
         // 32 bytes in length.
         // The construction in the next line is triggered in the _StaticSize == 0
         // case, and then chooses the smallest number of elements for which
         // the target allocation size is at least 32 bytes.
         if (n == 0) n = index_t((32u + sizeof(value_type) - 1u)/sizeof(value_type));
         _CxAssert(n >= index_t(1));
         reserve(n);
      }
   }
};


// Note to self: @swap:
// --------------------
// The recommended usage patter for a free template swap<T> function is this:
//
// 1. When declaring swappable classes `T`, make sure that a free-standing
//    swap(T&, T&) function is declared *in the same namespace*. This way in
//    a call of `swap(a,b)` where `a` or `b` are of type `T`, this function
//    will automatically be considered as a candidate (even without explicit
//    namespace qualification) by means of argument-dependent lookup (ADL).
//
// 2. In functions using swap on unknown/template types, *do not* use qualified
//    namespaces to refer to the swap function (e.g., "std::swap"), but rather
//    bring in viable candidates into the local namespace:
//
//       template<class T>
//       void some_function(T &a, T &b) { ...
//           using std::swap;
//           swap(A, B);
//       }
//
// 3. The actual swap() resolved in there then may or may not be std::swap,
//    depending on actual #includes and, most importantly, the concrete type
//    of `T` and functions handling it.
//      In effect, the "using std::swap;" just enables *considering* all the std
//    namespace specializations as candidates in the function resolution
//    ---however, if a `swap` is defined in the same namespace as `T`, it will
//    *also* be considered as candidate (via ADL), even if `T`s namespace has
//    not been explicitly imported.
template<class value_t>
void swap(TArray<value_t> &A, TArray<value_t> &B) {
   A.swap(B);
}


} // namespace ct



namespace ct {

/// Array class with fixed maximum size which stores its elements
/// in-place (i.e., no allocations).
///
/// For technical reasons, all _StaticSize array elements are default-constructed on
/// construction and only destroyed when TArrayFix is (i.e., technically,
/// all array elements are alive all the time even if size() < max_size()).
template<class _ValueType, unsigned int _StaticSize, class _IndexType = std::size_t>
struct TArrayFix {
   typedef _ValueType value_t;
   typedef _IndexType index_t;

   typedef value_t value_type;
   typedef value_t *iterator;
   typedef value_t const *const_iterator;
   typedef index_t size_t;

   // compiler-generated copy-ctor and assignment op should work.
   TArrayFix() { _InitDataMembers(); };

   explicit TArrayFix(index_t nEntries) { _InitDataMembers(); resize(nEntries); };
   TArrayFix(index_t nEntries, value_t const &Scalar) { _InitDataMembers(); assign(nEntries, Scalar); };

   template<class _InputIt>
   TArrayFix(_InputIt first, _InputIt last) { _InitDataMembers(); assign(first, last); }
#ifdef CX_ARRAY_HAVE_CXX11
   TArrayFix(std::initializer_list<value_t> init) { _InitDataMembers(); assign(init); };
   // ^-- intentionally not "explicit"
#endif // CX_ARRAY_HAVE_CXX11


   value_t &operator[](index_t i)             { _CxAssert(i < m_Size); return m_Data[i]; }
   value_t const &operator[](index_t i) const { _CxAssert(i < m_Size); return m_Data[i]; }

   value_t &back()              { _CxAssert(m_Size != 0); return m_Data[m_Size - 1]; }
   value_t const &back() const  { _CxAssert(m_Size != 0); return m_Data[m_Size - 1]; }
   value_t &front()             { _CxAssert(m_Size != 0); return m_Data[0]; }
   value_t const &front() const { _CxAssert(m_Size != 0); return m_Data[0]; }

   bool operator == (TArrayFix const &other) const {
      if (size() != other.size())
         return false;
      for (index_t i = 0; i != size(); ++i)
         if ((*this)[i] != other[i])
            return false;
      return true;
   }

   bool operator!=(TArrayFix const &other) const { return !this->operator==(other); }

   index_t size() const     { return m_Size; }
   bool empty() const       { return m_Size == 0; }
   index_t capacity() const { return _StaticSize; }
   index_t max_size() const { return _StaticSize; }
   void clear()             { resize(0); }

   void push_back(value_t const &t) { _CxAssert(m_Size < _StaticSize); m_Data[m_Size] = t; ++m_Size; }
   void pop_back()                  { _CxAssert(m_Size > 0); --m_Size; }

   void resize(index_t NewSize)     { _CxAssert(NewSize <= _StaticSize); m_Size = NewSize; }

   void resize(index_t NewSize, value_t const &value) {
      _CxAssert(NewSize <= _StaticSize);
      for (index_t i = m_Size; i < NewSize; ++i)
         m_Data[i] = value;
      m_Size = NewSize;
   }

   iterator erase(iterator first, iterator last) {
      _CxAssert(first >= begin() && first <= last && last <= end());
      m_Size -= last - first;
      for (iterator it = first; it < end(); ++it, ++last)
         *it = *last;
      return first;
   }

   iterator erase(iterator itWhere) { return erase(itWhere, itWhere + 1); };

   // assign scalar value to every element of *this
   void operator = (value_t const &Value) {
      for (iterator it = begin(); it != end(); ++it)
         *it = Value;
   }

   template<class _RandomIt>
   void assign(_RandomIt first, _RandomIt last) {
      resize(last - first);
      for (index_t i = 0; i < m_Size; ++i)
         m_Data[i] = first[i];
   }
   void assign(size_t nCount, value_t const &Value) { resize(nCount); *this = Value; }
#ifdef CX_ARRAY_HAVE_CXX11
   void assign(std::initializer_list<value_t> init) { assign(init.begin(), init.end()); }
#endif
   value_t *data()              { return &m_Data[0]; };
   value_t const *data() const  { return &m_Data[0]; };

   iterator begin()             { return &m_Data[0]; }
   iterator end()               { return &m_Data[m_Size]; }
   const_iterator begin() const { return &m_Data[0]; }
   const_iterator end() const   { return &m_Data[m_Size]; }

protected:
   value_t
      m_Data[_StaticSize];
   index_t
      m_Size; // actual number of elements (<= _StaticSize).

   void _InitDataMembers() {
      // that's it already.
      m_Size = 0;
   }
};

// lexicographically compare two arrays. -1: A < B; 0: A == B; +1: A > B.
template<class value_t, unsigned int MaxN, class FPred>
int Compare(TArrayFix<value_t, MaxN> const &A, TArrayFix<value_t, MaxN> const &B, FPred less = predicates::less_t<value_t>())
{
   typedef TArrayFix<value_t, MaxN> FArray;
   if (A.size() < B.size())
      return -1;
   if (B.size() < A.size())
      return +1;
   for (typename FArray::index_t i = 0; i < A.size(); ++i) {
      if (less(A[i], B[i]))
         return -1;
      if (less(B[i], A[i]))
         return +1;
   }
   return 0;
}

} // namespace ct


#ifdef INCLUDE_ABANDONED
//    // TODO: rewire this for empty-base-class optimization, deriving TArray from
//    // TStaticStorage. This requires:
//    // - moving m_StaticStorage and m_DataIsOnStaticBuffer into the base-class object
//    // - possibly creating a `std::size_t const _static_storage_size = _StaticStorage`
//    //   in the base-class object
//    // - accessing either of them only indirectly.
//    // We can probably remove the operators from TStaticStorage, because the
//    // resulting object is more or less equivalent to std::array (≥C++11), so
//    // there is little reason to retain it as a separate object.
//    typedef value_t *pointer;
//    typedef value_t &reference;
//    typedef value_t *iterator;
//    typedef value_t const *const_pointer;
//    typedef value_t const &const_reference;
//    typedef value_t const *const_iterator;
//    typedef std::size_t index_t;
//    typedef std::ptrdiff_t difference_type;
//    typedef value_t value_type;
#endif // INCLUDE_ABANDONED


#endif // CX_PODARRAY_H
