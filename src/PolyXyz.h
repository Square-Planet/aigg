#ifndef POLY_XYZ_H
#define POLY_XYZ_H

#include "Aigg.h"
#include <iosfwd>
#include <map>
#include <tuple>
#include <cstdint>
#include "CxTypes.h"
// TODO:
// - remove the eigen typedef for the monomial (as alternative, a
//   straight std::array with 3 unsigned shorts would probably do the trick).
// - replace implementation by the one from CxPackedIndex. It's more general
//   and better cross-checked.
// - change the integer sizes. Apparently on 64bit linux, *literally all*
//   the uint_fastXX_t types are typedef'ed to a 64bit unsigned long,
//   even including uint_fast8_t. Also, there now are "...at least 64bit"
//   uint types in standard C++, unlike in C++98.

// #include "CxVec3.h"


#define POLYN_PACKED_INDEX
// #define POLYN_USE_HASH_MAP


#ifdef POLYN_PACKED_INDEX
// todo: replace by CxPackedIndex version, and extend to 16bits per component.
struct FMonomialN {
   enum {
      nGen = AIG_SPACE_AXES,
      nBits = 10,
      Mask = (size_t(1) << nBits) - 1
   };
   typedef unsigned value_type;
   
   FMonomialN() : m_PackedPow(0) {}
   
   FMonomialN(value_type ix, value_type iy, value_type iz)
      : m_PackedPow(FPackedType(ix) | (FPackedType(iy) << nBits) | (FPackedType(iz) << (2*nBits)))
   {
      assert((ix < size_t(1)<<nBits) && (iy < size_t(1)<<nBits) && (iz < size_t(1)<<nBits));
   }
   
   template<class FUnsignedInt>
   explicit FMonomialN(std::array<FUnsignedInt,nGen> Exp) : FMonomialN(Exp[0], Exp[1], Exp[2]) {
      static_assert(!std::is_signed<FUnsignedInt>::value, "Monomials must have non-negative powers. No signed integers, please.");
   }
   
   value_type operator [] (unsigned i) const {
      assert(i <= 2);
      return (m_PackedPow >> (i*nBits)) & Mask;
   }

   bool operator == (FMonomialN const &other) const { return this->m_PackedPow == other.m_PackedPow; }
   bool operator != (FMonomialN const &other) const { return this->m_PackedPow != other.m_PackedPow; }
   bool operator < (FMonomialN const &other) const { return this->m_PackedPow < other.m_PackedPow; }
   
   FMonomialN operator + (FMonomialN const &other) const {
      return FMonomialN(this->m_PackedPow + other.m_PackedPow);
   }
   constexpr size_t size() const { return nGen; }
   
   static FMonomialN Zero() { return FMonomialN(FPackedType(0)); }
   
//    unsigned Degree() const { return (*this)[0] + (*this)[1] + (*this)[2]; }
protected:
   typedef size_t
      FPackedType;
   FPackedType
      m_PackedPow;
   explicit FMonomialN(FPackedType PackedPow_) : m_PackedPow(PackedPow_) {}
}; 

typedef std::less<FMonomialN> FMonomialCompareLess;

#else

typedef Eigen::Matrix<int, 3, 1>
   FMonomialN;

// less predicate ordering the monomials so that we can use them in a map.
struct FMonomialCompareLess {
   bool operator() (FMonomialN const &a, FMonomialN const &b) const
   {
#if __cplusplus > 199711L
      // C++11 or later?
      return std::make_tuple(a[0], a[1], a[2]) < std::make_tuple(b[0], b[1], b[2]);
#else
      for (unsigned i = a.size()-1; i < a.size(); --i)
         if (a[i] < b[i])
            return true;
         else if (a[i] > b[i])
            return false;
      return false; // never reached, just to please the compiler.
#endif
   }
};
   
#endif // POLYN_PACKED_INDEX



#ifdef POLYN_USE_HASH_MAP
#include <unordered_map>

   
static size_t const FPolynomialN_MaxOrder = 1024; // <- with this one can fit three orders into even a 32bit int.

namespace std {
   template <>
   class hash<FMonomialN> {
   public:
      size_t operator() (const FMonomialN &m) const {
#ifdef POLYN_PACKED_INDEX
         return size_t(m[0]) + FPolynomialN_MaxOrder*(size_t(m[1]) + FPolynomialN_MaxOrder*(size_t(m[2])));
#else
         return size_t(m.m_PackedPow);
#endif
      }
   };
}

#endif // POLYN_USE_HASH_MAP


// #if __cplusplus >= 201703L
// define get<n>(e) and tuple_size<> and tuple_element<> for FMonomialN to
// provide support for C++17 'structured binding'. Since all elements have the
// same type (unsigned int) and are anyway tiny intermediates, this is nowhere
// near as painful as normally.

namespace std {
   template<>
   struct tuple_size<FMonomialN>{ enum { value = FMonomialN::nGen }; };

   template<size_t iElemT>
   struct tuple_element<iElemT, FMonomialN >{ typedef unsigned type; };
}

template<size_t iElemT>
unsigned get(FMonomialN const &m) {
   return m[iElemT];
}

// #endif




extern FScalar g_ThrDefaulPolyCoeffNeglect;

struct FPolynomialN : public ct::FIntrusivePtrDest {
   // represents term x^p[0] * y^p[1] * z^p[2].
   enum {
      // number of generators (the N in FPolynomialN)
      nGen = 3
   };
   typedef FMonomialN
      FMonomial;
   
   typedef double
      FEstimatedNormType;
   typedef std::function<FEstimatedNormType(FMonomialN const&)>
      // a function taking a FMonomialN object and returning an estimate (therefore
      // double, instead of FScalar) of its squared norm in a application-specific manner.
      // Used (optionally) for term screening, by allowing to take account of the
      // importance of the term being screened itself (e.g, a term with a small coefficient
      // may still need to be retained if the attached monomial has a high weight)
      FEstimateSquaredNormFn;
   struct FCullCriteria {
      FScalar
         m_ThrCoeffNeglectAbs,
         // controls weight (coeff * term norm) with relation to heighest weight in the polynomial.
         // 0 -> don't use; -1 -> use same as absolute threshold.
         // Note: relative-threshold purge() is expensive and not normally done automatically.
         m_ThrCoeffNeglectRel;
      FEstimateSquaredNormFn
         m_pEstimateSquaredNormFn;

      explicit FCullCriteria(FScalar const &ThrCullAbsolute = g_ThrDefaulPolyCoeffNeglect, FScalar const &ThrCullRelative = 0, FEstimateSquaredNormFn const &pEstimateSquaredNorm = FEstimateSquaredNormFn())
      : m_ThrCoeffNeglectAbs(ThrCullAbsolute), m_ThrCoeffNeglectRel(ThrCullRelative), m_pEstimateSquaredNormFn(pEstimateSquaredNorm) {};

      // return whether to retain the monomial `m` with coefficient `c` according to the criteria of *this
      bool RetainQ(FMonomialN const &m, FScalar const &c) const;
      bool RetainAllQ() const { return m_ThrCoeffNeglectAbs == 0; }
      // replace unset/default thresholds in *this by data setup in `other`; otherwise
      // update *this to harder thresholds.
      // (note: allows free-standing functions and default initialized intermediates
      // to get cull criteria propagated from computations with arguments which have
      // explicitly set them)
      void Update(FCullCriteria const &other);

      FScalar const &GetThrCoeffNeglect() const { return m_ThrCoeffNeglectAbs; }
      void SetThrCoeffNeglect(FScalar const &thr) { m_ThrCoeffNeglectAbs = thr; }
      
      void SetNormEstimateFn(FEstimateSquaredNormFn const &pFn) { m_pEstimateSquaredNormFn = pFn; }
   public:
      FCullCriteria(FCullCriteria const &other) = default;
      FCullCriteria(FCullCriteria &&other) = default;
      FCullCriteria& operator = (FCullCriteria const &other) = default;
   };
   
//    explicit FPolynomialN(FScalar const &ThrCullAbsolute)
//    : m_Cull(ThrCullAbsolute) {}
   // no terms -> zero polynomial.
   explicit FPolynomialN(FCullCriteria const &Cull_ = FCullCriteria())
   : m_Cull(Cull_) {}
#ifdef INCLUDE_ABANDONED
//    : m_ThrCoeffNeglect(ThrZero), m_pEstimateSquaredNormFn(pEstimateSquaredNorm) {};
//    FPolynomialN(FScalar Coeff, FMonomial const &Monomial, FScalar const &ThrZero = g_ThrDefaulPolyCoeffNeglect);
#endif // INCLUDE_ABANDONED
   
   // creates a polynomial that represents a monomial, e.g. 1.2 * x^2 y or
   // something.
   FPolynomialN(FScalar const &Coeff, FMonomial const &Monomial, FCullCriteria const &Cull = FCullCriteria());

   // adds Factor * `other` to the current polynomial
   void Add(FScalar const &Factor, FPolynomialN const &other);
   // adds Factor * `other` to the current polynomial
   void Add(FScalar const &Factor, FMonomial const &Monomial);
   // multiplies this polynomial by Coeff * x^M[0] * y^M[1] * z^M[2];
   void Multiply(FScalar const &Coeff, FMonomial const &Monomial);

   void operator += (FPolynomialN const &other) {
      Add(FScalar(1), other);
   }
   FPolynomialN operator + (FPolynomialN const &other) const {
      FPolynomialN Result(*this);
      Result += other;
      return Result;
   }
   void operator -= (FPolynomialN const &other) {
      Add(FScalar(-1), other);
   }
   FPolynomialN operator - (FPolynomialN const &other) const {
      FPolynomialN Result(*this);
      Result -= other;
      return Result;
   }

   FPolynomialN operator * (FPolynomialN const &other) const;

   bool operator == (FPolynomialN const &other) const;
   
   bool operator != (FPolynomialN const &other) const {
      return !this->operator == (other);
   }
   
   void operator *= (FScalar f) {
      Multiply(f, FMonomial::Zero());
   }

   void operator *= (FMonomial const &m) {
      Multiply(FScalar(1), m);
   }

   void Print(std::ostream &out) const;
   
   bool IsZero() const { return m_Terms.empty(); }
   // ^- everything == 0 should have been Purged() already. So if anything is left,
   //    the result is non-zero.
   
   // delete all terms x^i y^j z^k for which any of i,j, or k is uneven.
   void EraseMonomialsWithOddPowers();
   
   FPolynomialN Purged() const {
      FPolynomialN r(*this);
      r.Purge();
      return r;
   }
#ifdef INCLUDE_ABANDONED
//    // return a copy of *this purged with a target threshold of `thr`.
//    // After the Purge(), the target threshold will be reset to the
//    // initial value of this->GetThrCoeffNeglect()
//    FPolynomialN Purged(FScalar thr) const {
//       if (thr <= 0) thr = this->GetThrCoeffNeglect();
//       FPolynomialN r(*this);
//       r.SetThrCoeffNeglect(thr);
//       r.Purge();
//       r.SetThrCoeffNeglect(this->GetThrCoeffNeglect());
//       return r;
//    }
#endif // INCLUDE_ABANDONED
   // return a copy of *this purged with a cull parameters of `Cull`.
   // After the Purge(), the cull parameters will be reset to the
   // initial value of this->GetCullCriteria()()
   FPolynomialN Purged(FCullCriteria const &Cull) const {
      FPolynomialN r(*this);
      r.SetCullCriteria(Cull, true); // true: purge now/
      r.SetCullCriteria(this->GetCullCriteria(), false);
      return r;
   }
   
   // removes terms with near-zero prefactors (most likely resulting
   // from numerical noise.)
   void Purge();

#ifdef INCLUDE_ABANDONED
//    // Set Threshold m_ThrCoeffNeglect below which coefficients will
//    // be deleted in Purge. Notes:
//    // - defaults to g_ThrAlmostZero
//    // - can be changed via SetThrCoeffNeglect()
//    //   (to afford higher accuracy in intermediates of calculations)
//    // - with m_ThrCoeffNeglect == 0, Purge will become a non-op and all
//    //   terms will be retained. With SetThrCoeffNeglect(-1), m_ThrCoeffNeglect
//    //   will be reset to its default value.
//    // - will be typicall inherited as min(a,b) for operations relating a and b.
//    // - if a smaller threshold is changes to an (effective) larger threshold,
//    //   Purge() will be called by this function and terms over the old
//    //   threshold, but below the new one, will be deleted.
//    void SetThrCoeffNeglect(FScalar ThrZero, bool AutoPurge=true);
//    FScalar GetThrCoeffNeglect() const { return m_ThrCoeffNeglect; }
#endif // INCLUDE_ABANDONED
   void SetCullCriteria(FCullCriteria const &Cull, bool AutoPurge=true);
   FCullCriteria const &GetCullCriteria() const { return m_Cull; }
   
   FScalar const &GetThrCoeffNeglect() const { return m_Cull.GetThrCoeffNeglect(); }
protected:
//    FScalar
//       m_ThrCoeffNeglect;
//    FEstimateSquaredNormFn
//       m_pEstimateSquaredNormFn;
   FCullCriteria
      m_Cull;
   
#ifdef POLYN_USE_HASH_MAP
   typedef std::unordered_map<FMonomial, FScalar>
      FTermList;
#else
   typedef std::map<FMonomial, FScalar, FMonomialCompareLess>
      FTermList;
#endif
   FTermList
       // maps monomials to their prefactors. The entire polynomial is the sum
       // of all terms contained herein.
       m_Terms;
public:
   void swap(FPolynomialN &other) {
      using std::swap;
      swap(m_Cull, other.m_Cull);
      m_Terms.swap(other.m_Terms);
   }
   
   typedef FTermList::value_type value_type;
   typedef FTermList::iterator iterator;
   typedef FTermList::const_iterator const_iterator;
   const_iterator begin() const { return m_Terms.begin(); }
   const_iterator end() const { return m_Terms.end(); }

   size_t size() const { return m_Terms.size(); }
   // these are technically valid polynomials---which evaluate to zero.
   bool empty() const { return m_Terms.empty(); }
   // delete all controlled terms and set the polynomial to zero.
   void clear() { m_Terms.clear(); }

public:
   FPolynomialN(FPolynomialN &&other) = default;
   FPolynomialN(FPolynomialN const &other) = default;
   FPolynomialN& operator = (FPolynomialN const &other) = default;
   FPolynomialN& operator = (FPolynomialN &&other) = default;
};
typedef ct::TIntrusivePtr<FPolynomialN>
   FPolynomialNPtr;
typedef ct::TIntrusivePtr<FPolynomialN const>
   FPolynomialNCptr;


inline void swap(FPolynomialN &a, FPolynomialN &b) {
   a.swap(b);
}


inline FPolynomialN operator * (FScalar f, FPolynomialN const &other)
{
   FPolynomialN Result(other);
   Result.Multiply(f, FPolynomialN::FMonomial::Zero());
   return Result;
}


inline std::ostream &operator << (std::ostream &out, FPolynomialN const &p) {
   p.Print(out);
   return out;
}

#ifdef INCLUDE_ABANDONED
// return A*B*C, but checks in which order to multiply.
// inline FPolynomialN Multiply(FPolynomialN const &A, FPolynomialN const &B, FPolynomialN const &C)
// {
//    // we want to combine the two items with the *largest* number of terms first.
//    // UPDATE: ...I think that would only be the right choice if we assume that
//    // in this case the most terms cancel in the first product. Not clear if that
//    // is the case. Maybe we should actually compute it?
//    if (C.size() <= A.size() && C.size() <= B.size())
//       return (A*B)*C;
//    if (B.size() <= A.size() && B.size() <= C.size())
//       return (A*C)*B;
//    return (B*C)*A;
// }
#endif // INCLUDE_ABANDONED

// return A*B*C, but checks in which order to multiply.
FPolynomialN Multiply(FPolynomialN const &A, FPolynomialN const &B, FPolynomialN const &C);

// return A^na * B^nb * C^nc, forming suitable intermediate results
FPolynomialN MultiplyPow(FPolynomialN const &A, size_t na, FPolynomialN const &B, size_t nb, FPolynomialN const &C, size_t nc);


// returns A*A*...*A  (n-times). This function bisects across 'n' and caches
// intermediate products (so they exist)
FPolynomialN pow(FPolynomialN const &A, size_t n);


struct FMonomialN_ExportOrderPred { bool operator() (FMonomialN const &A, FMonomialN const &B) const; };   
struct FMonomialN_SmallIntegralsFirstOrder { bool operator() (FMonomialN const &A, FMonomialN const &B) const; };   
struct FMonomialN_BasicOrder { bool operator() (FMonomialN const &A, FMonomialN const &B) const; }; // creates an explicit tuple of (ix,iy,iz) as sort key
// struct FMonomialN_FastOrder { bool operator() (FMonomialN const &A, FMonomialN const &B) const; }; // uses FMonomialN::operator < () const.
typedef std::less<FMonomialN> FMonomialN_FastOrder;



#endif // POLY_XYZ_H
