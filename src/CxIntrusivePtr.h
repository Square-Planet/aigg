#ifndef CX_INTRUSIVE_PTR
#define CX_INTRUSIVE_PTR

#include <stddef.h> // for ptrdiff_T
#include <algorithm> // for std::swap.

// A re-implementation of the boost intrusive ptr. reason is that due to
// recent changes in gcc, some problematic template order- of-definition
// rules are now enforced, which leads to the situation that
// add_ref/release *must* be declared before intrusive_ptr.hpp is
// included... something we can of course not guarantee in general!
//
// Additionally, this allows us to drop one more dependency on boost.
//
// Actual pointer relies on free-standing:
//      void intrusive_ptr_add_ref(T *p);
//      void intrusive_ptr_release(T *p);
// functions, just as the boost ptr.
//
// Objects derived from class FIntrusivePtrDest1 support those methods
// automatically, but deriving from this class is not mandatory.
// (as long as the corresponding add_ref/release functions are declared
// before this header is included).

namespace ct {

struct FIntrusivePtrDest1
{
   FIntrusivePtrDest1() : m_RefCount(0) {}
   inline virtual ~FIntrusivePtrDest1() = 0;

   ptrdiff_t mutable
      m_RefCount;
//    friend void ::intrusive_ptr_add_ref(FIntrusivePtrDest1 const *p);
//    friend void ::intrusive_ptr_release(FIntrusivePtrDest1 const *p);
};

inline FIntrusivePtrDest1::~FIntrusivePtrDest1()
{
}

} // namespace ct

inline void intrusive_ptr_add_ref(ct::FIntrusivePtrDest1 const *p) {
   p->m_RefCount += 1;
}

inline void intrusive_ptr_release(ct::FIntrusivePtrDest1 const *p) {
#ifdef assert
   assert(p->m_RefCount != 0);
#endif
   p->m_RefCount -= 1;
   if (p->m_RefCount == 0)
      delete p;
}

namespace ct {

template <class T>
struct TIntrusivePtr
{
   typedef T
      *pointer;
   typedef T
      element_type;
   typedef TIntrusivePtr
      this_type;
   // we do it just like boost: this is a pointer-to-member, which can be
   // implicitly converted to bool, but not to, say, other integral types or
   // floating point numbers. returning this instead of actual bools thus prevents
   // some unwanted implicit conversions, e.g.:
   //      TIntrusivePtr<something> p;
   //      double f = p;
   // is not legal with 'pointer to member', but would be legal with 'bool'.
   typedef pointer
      this_type::*unspecified_bool_type;

//
// constructors, destructors, and assignment operators
//
   TIntrusivePtr() : m_ptr(0) {}

   // construct from raw pointer of same type
   TIntrusivePtr(pointer p, bool add_ref = true) : m_ptr(p)  {
      if (m_ptr != 0 && add_ref)
         intrusive_ptr_add_ref(m_ptr);
   }

   TIntrusivePtr(TIntrusivePtr const &p) : m_ptr(p.get())  {
      if (m_ptr != 0)
         intrusive_ptr_add_ref(m_ptr);
   }

   // construct from related intrusive pointer
   template <class U>
   TIntrusivePtr(TIntrusivePtr<U> const &p) : m_ptr(p.get())  {
      if (m_ptr != 0)
         intrusive_ptr_add_ref(m_ptr);
   }

   ~TIntrusivePtr() {
      if (m_ptr != 0)
         intrusive_ptr_release(m_ptr);
   }

   // assign via copy constructor: x = bla -> x.swap(type(bla))
   TIntrusivePtr &operator = (pointer other) {
      this_type(other).swap(*this);
      return *this;
   };

   TIntrusivePtr &operator = (TIntrusivePtr const &other) {
      this_type(other).swap(*this);
      return *this;
   }

   template <class U>
   TIntrusivePtr &operator = (TIntrusivePtr<U> const &other) {
      this_type(other).swap(*this);
      return *this;
   }

//
// accessors
//
   pointer &get() { return m_ptr; }
   pointer const &get() const { return m_ptr; }

   T &operator*() const {  return *m_ptr; }

   pointer &operator->() { return m_ptr; }
   pointer const &operator->() const { return m_ptr; }

   operator unspecified_bool_type () const {
      return m_ptr == 0? 0: &this_type::m_ptr;
   }

   bool operator !() const { return m_ptr == 0; }

   void swap(TIntrusivePtr &other) {
      std::swap(m_ptr, other.m_ptr);
   }
private:
   T *m_ptr;
};

template<class U, class V>
inline void swap(TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   a.swap(b);
}



// comparison operators: ==, !=, and < between the intrusive pointers
// and raw pointers (in all permutations).
template<class U, class V>
inline bool operator == (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() == b.get();
}

template<class U, class V>
inline bool operator != (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() != b.get();
}

template<class U, class V>
inline bool operator < (TIntrusivePtr<U> const &a, TIntrusivePtr<V> const &b) {
   return a.get() < b.get();
}

template<class U, class V>
inline bool operator == (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() == b;
}

template<class U, class V>
inline bool operator != (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() != b;
}
template<class U, class V>
inline bool operator < (TIntrusivePtr<U> const &a, V const *b) {
   return a.get() < b;
}

template<class U, class V>
inline bool operator == (U const *a, TIntrusivePtr<V> const &b) {
   return a == b.get();
}

template<class U, class V>
inline bool operator != (U const *a, TIntrusivePtr<V> const &b) {
   return a != b.get();
}

template<class U, class V>
inline bool operator < (U const *a, TIntrusivePtr<V> const &b) {
   return a < b.get();
}


} // namespace ct

#endif // CX_INTRUSIVE_PTR
