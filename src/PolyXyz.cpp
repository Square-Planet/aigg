#include "Aigg.h"
#include "PolyXyz.h"
#include <ostream>
#include <math.h>

static bool const FPolynomialN_AutoPurge = true;
FScalar g_ThrDefaulPolyCoeffNeglect = g_ThrAlmostZero;

// static bool g_DebugUseFactorScreening = true;
// ^-- yields strange symmetrization errors in some cases. And is not very clean in any case.
static bool g_DebugUseFactorScreening = false;

#ifdef _for_each
   #undef _for_each
#endif
#define _for_each(it,cont) for (it = cont.begin(); it != cont.end(); ++it)


void FPolynomialN::Print(std::ostream &out) const
{
//    bool FlatPows = true;
//    FlatPows = false;
   int PowFormat = 2;

   FPolynomialN::FTermList::const_iterator itTerm;
   bool EmittedSomethingAtAll = false; // set to true if any term is printed.
   _for_each(itTerm, m_Terms)
   {
      // attempt to move unary minus sign operator of coefficient into binary
      // operator joining terms, replacing the default addition with a subtraction.
//       int iSign = (itTerm->second < 0) ? -1 : +1;
      int iSign = 1;
      if (itTerm != m_Terms.begin()) {
         if (itTerm->second < 0) {
            out << " - ";
            iSign *= -1; // if present, will cancel negative sign of coeff
         } else {
            out << " + ";
         }
      }
      char const *pScalarOp = 0;
      bool EmittedSomethingForTerm = false;
      {
         FScalar fp = (iSign*itTerm->second);
         if (fp < 0) { // can happen if first coeff is negative
            out << "-"; fp = -fp;
         }
         if (!IsAlmostZero(fp-1)) {
            out << fp;
//             pScalarOp = "*";
            pScalarOp = " ";
            EmittedSomethingForTerm = true;
         }
      }

      if (itTerm->first != FMonomial::Zero() || PowFormat == 2) {
         if (pScalarOp) {
            out << pScalarOp;
            pScalarOp = 0;
         }
         char Components[] = {'x', 'y', 'z', 'w'};
         assert(sizeof(Components)/sizeof(Components[0]) >= FPolynomialN::nGen);
         char const *pCompJoinOp = 0;
         for (unsigned i = 0; i < FPolynomialN::nGen; ++i) {
            if (itTerm->first[i] != 0 || PowFormat == 2) {
               if (PowFormat == 0) {
                  // emit monomials like "x^2 y^4 z"
                  if (pCompJoinOp) { out << pCompJoinOp; pCompJoinOp = 0; }
                  out << Components[i];
                  EmittedSomethingForTerm = true;
                  if (itTerm->first[i] != 1)
                     out << "^" << itTerm->first[i];
                  pCompJoinOp = " ";
               } if (PowFormat == 2) {
                  // emit monomials like "[4,2,6]"
                  if (i == 0) out << "[";
//                   out << Components[i];
                  if (pCompJoinOp) { out << pCompJoinOp; pCompJoinOp = 0; }
                  EmittedSomethingForTerm = true;
                  out << itTerm->first[i];
                  pCompJoinOp = ",";
                  if (i == FPolynomialN::nGen-1) out << "]";

               } else {
                  // emit monomials like "xx yyyy z"
                  unsigned np = itTerm->first[i];
                  if (np != 0) {
                     if (pCompJoinOp) { out << pCompJoinOp; pCompJoinOp = 0; }
                     EmittedSomethingForTerm = true;
                     pCompJoinOp = " ";
                  }
                  for (size_t ip = 0; ip != np; ++ ip) {
                     out << Components[i];
                  }
               }
            }
         }
      }
      
      if (!EmittedSomethingForTerm) {
         // Normally we omit print 1 coefficients and a (0,0,0) monomials.
         // However, here we have both at the same time, and then we can't
         // just leave everything away. So we need to print a 1 now.
         if (pScalarOp) {
            out << pScalarOp;
            pScalarOp = 0;
         }
         out << 1;
         EmittedSomethingForTerm = true;
      }
      if (EmittedSomethingForTerm)
         EmittedSomethingAtAll = true;
   }
   if (!EmittedSomethingAtAll) {
      out << "0";
   }
}


bool FPolynomialN::FCullCriteria::RetainQ(FMonomialN const &m, FScalar const &c) const
{
   using std::abs;
   using std::sqrt;
   assert(m_ThrCoeffNeglectAbs >= 0);
   if (RetainAllQ())
      return true;
   if (!m_pEstimateSquaredNormFn)
      return abs(c) >= m_ThrCoeffNeglectAbs;
   else
      return sqr(double(c)) * double(m_pEstimateSquaredNormFn(m)) >= sqr(double(m_ThrCoeffNeglectAbs)); // <-- note: right one should have target norm of object absorbed. so there are norms and coeffs on both sides.
//       return abs(double(c)) * sqrt(double(m_pEstimateSquaredNormFn(m))) >= abs(double(m_ThrCoeffNeglectAbs));
}

void FPolynomialN::FCullCriteria::Update(FCullCriteria const &other)
{
   this->m_ThrCoeffNeglectAbs = std::min(this->m_ThrCoeffNeglectAbs, other.m_ThrCoeffNeglectAbs);
   if (!bool(this->m_pEstimateSquaredNormFn) && bool(other.m_pEstimateSquaredNormFn))
      this->m_pEstimateSquaredNormFn = other.m_pEstimateSquaredNormFn;
//    if (other.m_ThrCoeffNeglectRel > 0) {
//       // `other` has a relative threshold explicitly set. copy/update it.
//       if (this->m_ThrCoeffNeglectRel > 0)
//          // both relative thresholds explicitly set. retain the harder one.
//          this->m_ThrCoeffNeglectRel = std::min(this->m_ThrCoeffNeglectRel, other.m_ThrCoeffNeglectRel);
//       else
//          // *this thr is <= 0 (set to default or disable). Copy from other object.
//          this->m_ThrCoeffNeglectRel = other.m_ThrCoeffNeglectRel;
//    }
}


void FPolynomialN::Add(FScalar const &Factor, FPolynomialN::FMonomial const &Monomial)
{
   using std::abs;
   if (Factor == 0 || (g_DebugUseFactorScreening && (abs(Factor) < m_Cull.m_ThrCoeffNeglectAbs))) {
      return;
   } else {
      FTermList::iterator
         it = this->m_Terms.find(Monomial);
      // have such a term already? 
      if (it != this->m_Terms.end()) {
         // yes. add the new coefficient
         it->second += Factor;
         // does this result in a almost zero-coefficient term?
         if (!m_Cull.RetainQ(it->first, it->second))
            // yes. remove the monomial completely.
            this->m_Terms.erase(it);
#ifdef INCLUDE_ABANDONED
//          if (abs(it->second) < m_Cull.m_ThrCoeffNeglectAbs)
//             // yes. remove the monomial completely.
//             this->m_Terms.erase(it);
#endif // INCLUDE_ABANDONED
      } else {
         // it's a new monomial. Just insert it to the term list with its target factor.
         if (m_Cull.RetainQ(Monomial, Factor))
            this->m_Terms.insert(FTermList::value_type(Monomial, Factor));
      }
   }
   // note: no need for call to Purge(). The above logic should cover all
   // cases where the Add() leads to vanishing terms.
}


void FPolynomialN::Add(FScalar const &Factor, FPolynomialN const &other)
{
   using std::abs;
   if (g_DebugUseFactorScreening && (abs(Factor) < std::min(m_Cull.m_ThrCoeffNeglectAbs, other.m_Cull.m_ThrCoeffNeglectAbs))) {
      return;
   } else {
      FTermList::const_iterator itOtherTerm;
      _for_each(itOtherTerm, other.m_Terms)
      {
         FTermList::iterator itThisTerm;
         // if we have such a term already, just add the coefficient. Otherwise
         // copy the term from the other polynomial.
         itThisTerm = this->m_Terms.find(itOtherTerm->first);
         if (itThisTerm != this->m_Terms.end()) {
            itThisTerm->second += Factor*itOtherTerm->second;
         } else {
            this->m_Terms[itOtherTerm->first] = Factor*itOtherTerm->second;
         }
      }
      m_Cull.Update(other.m_Cull);
      if (FPolynomialN_AutoPurge)
         Purge();
   }
}


// creates a polynomial that represents a monomial, e.g. 1.2 * x^2 y or
// something.
FPolynomialN::FPolynomialN(FScalar const &Coeff, FMonomial const &Monomial, FPolynomialN::FCullCriteria const &Cull)
   : m_Cull(Cull)
{
   m_Terms[Monomial] = Coeff;
   Purge();
}


void FPolynomialN::Multiply(FScalar const &Coeff, FMonomial const &Monomial)
{
   FTermList NewTerms;

   FTermList::iterator itTerm;
   _for_each(itTerm, this->m_Terms)
   {
      // add powers, multiply coefficients.
      NewTerms[itTerm->first + Monomial] = itTerm->second * Coeff;
   }
   m_Terms.swap(NewTerms); // <- can't do the replacement inline due to the map
                           // structure.
   if (FPolynomialN_AutoPurge)
      Purge();
}


// FPolynomialN FPolynomialN::operator * (FPolynomialN const &other) const
// {
//    FPolynomialN Result;
//    FTermList::const_iterator itOtherTerm;
//    _for_each(itOtherTerm, other.m_Terms)
//    {
//       // sequentially multiply monomials of other polynomial by terms of
//       // *this and add the results.
//       FPolynomialN ThisMulOtherTerm(*this);
//       ThisMulOtherTerm.Multiply(itOtherTerm->second, itOtherTerm->first);
//       Result += ThisMulOtherTerm;
//    }
//    Result.Purge();
//    return Result;
// }


FPolynomialN FPolynomialN::operator * (FPolynomialN const &other) const
{
   FPolynomialN
      Result(m_Cull);
   Result.m_Cull.Update(other.m_Cull);

   FTermList::const_iterator
      itThisTerm, itOtherTerm;
   _for_each(itThisTerm, this->m_Terms)
      _for_each(itOtherTerm, other.m_Terms)
         Result.m_Terms[itThisTerm->first + itOtherTerm->first] += itThisTerm->second * itOtherTerm->second;
   if (FPolynomialN_AutoPurge)
      Result.Purge();
//    Result.m_Cull.m_ThrCoeffNeglectAbs = std::min(m_Cull.m_ThrCoeffNeglectAbs, other.m_Cull.m_ThrCoeffNeglectAbs);
   return Result;
}

#ifdef INCLUDE_ABANDONED
// void FPolynomialN::SetThrCoeffNeglect(FScalar ThrZero, bool AutoPurge)
// {
//    if (ThrZero < 0)
//       ThrZero = g_ThrDefaulPolyCoeffNeglect;
//    bool
//       // change to larger threshold -> some terms might have fallen below the
//       // new one. Delete those.
//       NeedsPurge = AutoPurge && (ThrZero > m_Cull.m_ThrCoeffNeglectAbs);
//    m_Cull.m_ThrCoeffNeglectAbs = ThrZero;
//    
//    if (NeedsPurge)
//       Purge();
// }
// 
// void FPolynomialN::Purge()
// {
//    assert(m_Cull.m_ThrCoeffNeglectAbs >= 0);
//    if (m_Cull.m_ThrCoeffNeglectAbs == 0)
//       return;
//    FTermList::iterator itTerm, itNextTerm;
//    for (itTerm = m_Terms.begin(); itTerm != m_Terms.end();) {
//       itNextTerm = itTerm;
//       ++itNextTerm;
// 
//       if (IsAlmostZero(itTerm->second, m_Cull.m_ThrCoeffNeglectAbs))
//          m_Terms.erase(itTerm, itNextTerm); // itTerm now invalidated.
// 
//       itTerm = itNextTerm;
//    }
// }
#endif // INCLUDE_ABANDONED

void FPolynomialN::SetCullCriteria(FCullCriteria const &Cull, bool AutoPurge)
{
   m_Cull = Cull;
   if (AutoPurge)
      Purge();
}

void FPolynomialN::Purge()
{
   if (m_Cull.RetainAllQ())
      return;
   FTermList::iterator itTerm, itNextTerm;
   for (itTerm = m_Terms.begin(); itTerm != m_Terms.end();) {
      itNextTerm = itTerm;
      ++itNextTerm;

      if (!m_Cull.RetainQ(itTerm->first, itTerm->second))
         m_Terms.erase(itTerm, itNextTerm); // itTerm now invalidated.

      itTerm = itNextTerm;
   }
}

#ifdef INCLUDE_ABANDONED
// #ifdef POLYN_USE_HASH_MAP
// 
// bool FPolynomialN::operator == (FPolynomialN const &other) const
// {
//    if (this->m_Terms.size() != other.m_Terms.size())
//       return false;
//    FScalar
//       ThrCmp = std::min(m_Cull.m_ThrCoeffNeglectAbs, other.m_Cull.m_ThrCoeffNeglectAbs);
//    FTermList::const_iterator
//       itThis = this->m_Terms.begin(),
//       itOther;
//    for ( ; itThis != this->m_Terms.end(); ++ itThis) {
//       itOther = other.m_Terms.find(itThis->first);
//       if (itOther == other.m_Terms.end())
//          return false;
// 
//       if (!IsAlmostZero(itThis->second - itOther->second, ThrCmp))
//          return false;
//    }
//    return true;
// }
// 
// #else
// bool FPolynomialN::operator == (FPolynomialN const &other) const
// {
//    if (this->m_Terms.size() != other.m_Terms.size())
//       return false;
//    FScalar
//       ThrCmp = std::min(m_Cull.m_ThrCoeffNeglectAbs, other.m_Cull.m_ThrCoeffNeglectAbs);
//    FTermList::const_iterator
//       itThis = this->m_Terms.begin(),
//       itOther = other.m_Terms.begin();
//    for ( ; itThis != this->m_Terms.end(); ) {
//       if (itThis->first != itOther->first)
//          return false;
//       if (!IsAlmostZero(itThis->second - itOther->second, ThrCmp))
//          return false;
//       
//       ++ itThis;
//       ++ itOther;
//    }
//    return true;
// }
// #endif // POLYN_USE_HASH_MAP
#endif // INCLUDE_ABANDONED

#ifdef POLYN_USE_HASH_MAP

bool FPolynomialN::operator == (FPolynomialN const &other) const
{
   if (this->m_Terms.size() != other.m_Terms.size())
      return false;
   FCullCriteria
      TestCrtiera(m_Cull);
   TestCrtiera.Update(other.m_Cull);
   
   FTermList::const_iterator
      itThis = this->m_Terms.begin(),
      itOther;
   for ( ; itThis != this->m_Terms.end(); ++ itThis) {
      itOther = other.m_Terms.find(itThis->first);
      if (itOther == other.m_Terms.end())
         return false;

      // According to the merged culling criteria, would we retain the monomial
      // if it had the *this vs `other` coefficient difference as prefactor?
      if (TestCrtiera.RetainQ(itThis->first, itThis->second - itOther->second))
         // yes, difference term would be retained. So it is not small enough to pretend it is zero.
         return false;
   }
   return true;
}

#else

bool FPolynomialN::operator == (FPolynomialN const &other) const
{
   if (this->m_Terms.size() != other.m_Terms.size())
      return false;

   FCullCriteria
      TestCrtiera(m_Cull);
   TestCrtiera.Update(other.m_Cull);

   FTermList::const_iterator
      itThis = this->m_Terms.begin(),
      itOther = other.m_Terms.begin();
   for ( ; itThis != this->m_Terms.end(); ) {
      if (itThis->first != itOther->first)
         return false;

      if (TestCrtiera.RetainQ(itThis->first, itThis->second - itOther->second))
         return false;
      
      ++ itThis;
      ++ itOther;
   }
   return true;
}

#endif // POLYN_USE_HASH_MAP


void FPolynomialN::EraseMonomialsWithOddPowers()
{
   FTermList::iterator itTerm, itNextTerm;
   for (itTerm = m_Terms.begin(); itTerm != m_Terms.end();) {
      itNextTerm = itTerm;
      ++itNextTerm;
      
      FMonomial const
         &m = itTerm->first;

      if (m[0] % 2 != 0 || m[1] % 2 != 0 || m[2] % 2 != 0)
         m_Terms.erase(itTerm, itNextTerm); // itTerm now invalidated.

      itTerm = itNextTerm;
   }
}


// typedef std::map<size_t, FPolynomialNCptr>
//    FPolynomialNPowerMap;
typedef std::map<size_t, FPolynomialNPtr>
   FPolynomialNPowerMap;

FPolynomialNCptr ComputePowN(size_t n, FPolynomialNPowerMap &Int)
{
   // is this a result we already have?
   FPolynomialNPowerMap::const_iterator
      itInt = Int.find(n);
   if (itInt != Int.end())
      return itInt->second;
   // nope. So make it now.
   assert(n >= 2); // should already have 0 and 1 in the map to start with.
   size_t
      n1 = n/2,
      n2 = n - n1;   
   // compute A^n1 and A^n2
   FPolynomialN const
      &Left = *ComputePowN(n1, Int),
      &Right = *ComputePowN(n2, Int);
   // and now A^n = A^n1 * A^n2 itself. And remember the result.
   FPolynomialNPtr
      res = FPolynomialNPtr(new FPolynomialN(Left*Right));
   Int[n] = res;
   return res;
}

   
// returns A*A*...*A  (p-times). This function bisects across 'p' and caches
// intermediate products (so they exist)
FPolynomialN pow(FPolynomialN const &A, size_t n)
{
   if (n == 0)
      // that's 1.
      return FPolynomialN(1, FPolynomialN::FMonomial(0,0,0), A.GetCullCriteria());
   if (n == 1)
      // A^1 -> return A itself.
      return A;
   if (n == 2)
      return A*A;
   if (0) {
      size_t n1 = n/2;
      return pow(A, n1) * pow(A, n-n1);
   } else {
      // something larger. Make cache for intermediate results,
      // and add '1' and a copy of A to it.
      FPolynomialNPowerMap
         Int; // <- cache for intermediate results
      Int[0] = FPolynomialNPtr(new FPolynomialN(1, FPolynomialN::FMonomial(0,0,0)));
      Int[1] = FPolynomialNPtr(new FPolynomialN(A));
      Int[0]->SetCullCriteria(A.GetCullCriteria());
      // ^- should propagate to rest automagically
      return *ComputePowN(n, Int); 
   }
}


// typedef Eigen::Array<int, 3, 1> FPow3Index;
// 
// // less predicate ordering the monomials so that we can use them in a map.
// struct FPow3IndexCompareLess {
//    bool operator() (FPow3Index const &a, FPow3Index const &b) const
//    {
//       for (unsigned i = a.size()-1; i < a.size(); --i)
//          if (a[i] < b[i])
//             return true;
//          else if (a[i] > b[i])
//             return false;
//       return false; // never reached, just to please the compiler.
//    }
// };
// 
// 
// typedef std::map<FPow3Index, FPolynomialNCptr, FPow3IndexCompareLess>
//    FPolynomialNPowerM_Map;
// 
// FPolynomialNCptr ComputePowMN(FPow3Index nabc, FPolynomialNPowerM_Map &Int)
// {
//    typedef FPolynomialNPowerM_Map
//       FMap;
//    // is this a result we already have?
//    FMap::const_iterator
//       itInt = Int.find(nabc);
//    if (itInt != Int.end())
//       return itInt->second;
//    // nope. So make it now.
//    unsigned
//       na = nabc[0],
//       nb = nabc[1],
//       nc = nabc[2];
//    // reduce largest exponent we have. First find which one that is.
//    unsigned
//       iDirMax = 0,
//       nMax = na;
//    if (nb >= nMax) { nMax = nb; iDirMax = 1; }
//    if (nc >= nMax) { nMax = nc; iDirMax = 2; }
//    
//    FPow3Index
//       // power indices of the two intermediate terms we will multiply together
//       nabc_Left,
//       nabc_Right;
//    
//    if (nMax == 1) {
//       // all exponents are 1 or 0. reduce highest one which is not yet zero.
//       cx_assert_rt(iDirMax >= 1); // <- shouldn't be here in case na<=1 nb=0 nc=0.
//       nabc_Left = nabc;
//       nabc_Left[iDirMax] -= 1;
//       nabc_Right = FPow3Index(0,0,0);
//       nabc_Right[iDirMax] += 1;
//    } else {
//       nabc_Left = nabc / 2;
//       nabc_Right = nabc - nabc_Left;
//    }
//    
//    // compute A^nabc_Left[0] B^nabc_Left[1] C^nabc_Left[2], and same for _Right variant.
// //    cx_assert_rt(nabc_Left + nabc_Right == nabc);
// //    cx_assert_rt((nabc_Left + nabc_Right - nabc).abs().max() == 0);
//    for (int ixyz = 0; ixyz < 3; ++ ixyz)
//       cx_assert_rt(nabc_Left[ixyz] + nabc_Right[ixyz] == nabc[ixyz]);
//    FPolynomialN const
//       &Left = *ComputePowMN(nabc_Left, Int),
//       &Right = *ComputePowMN(nabc_Right, Int);
//    // and multiply them together
//    FPolynomialNCptr
//       res = FPolynomialNCptr(new FPolynomialN(Left*Right));
//    // remember the intermediate result
//    Int[nabc] = res;
//    return res;
// }


typedef Eigen::Array<int, 3, 1> FPow3Index;

// less predicate ordering the monomials so that we can use them in a map.
struct FPow3IndexCompareLess {
   bool operator() (FPow3Index const &a, FPow3Index const &b) const
   {
      for (unsigned i = a.size()-1; i < a.size(); --i)
         if (a[i] < b[i])
            return true;
         else if (a[i] > b[i])
            return false;
      return false; // never reached, just to please the compiler.
   }
};


typedef std::map<FPow3Index, FPolynomialNCptr, FPow3IndexCompareLess>
   FPolynomialNPowerM_Map;

FPolynomialNCptr ComputePowMN(FPow3Index nabc, FPolynomialNPowerM_Map &Int)
{
   typedef FPolynomialNPowerM_Map
      FMap;
   // is this a result we already have?
   FMap::const_iterator
      itInt = Int.find(nabc);
   if (itInt != Int.end())
      return itInt->second;
   // nope. So make it now.
   FPow3Index
      // power indices of the two intermediate terms we will multiply together
      nabc_Left,
      nabc_Right;
   
//    if (nabc.maxCoeff() == 1) {
//       // all exponents are 1 or 0. reduce one which is not yet zero.
//       nabc_Left = nabc;
//       nabc_Right = FPow3Index(0,0,0);
//       for (int ixyz = 0; ixyz < 3; ++ ixyz) {
//          if (nabc[ixyz] != 0) {
//             nabc_Left[ixyz] -= 1;
//             nabc_Right[ixyz] += 1;
//             break;
//          }
//       }
//    } else {
//       if (0) {
//          nabc_Left = nabc / 2;
//       } else {
//          // these ones are for balancing even/odd numbers on average
//          nabc_Left = nabc;
//          nabc_Left[0] += 1;
//          if (nabc.sum() % 2 == 0)
//             nabc_Left[2] += 1;
//          nabc_Left /= 2;
//       }
//       nabc_Right = nabc - nabc_Left;
//    }
   // this variant: first split off either the largest or the smallest
   // of the indices in ^ijk (e.g., (i,j,k) into (i,j,0) and (0,0,k)),
   // then reduce what is left with the incremental bisection variant
   // above.
   if (nabc[0] == 0 || nabc[1] == 0 || nabc[2] == 0) {
      if (nabc.maxCoeff() == 1) {
         // all exponents are 1 or 0. reduce one which is not yet zero.
         nabc_Left = nabc;
         nabc_Right = FPow3Index(0,0,0);
         for (int ixyz = 0; ixyz < 3; ++ ixyz) {
            if (nabc[ixyz] != 0) {
               nabc_Left[ixyz] -= 1;
               nabc_Right[ixyz] += 1;
               break;
            }
         }
      } else {
         if (0) {
            nabc_Left = nabc / 2;
         } else {
            // these ones are for balancing even/odd numbers on average
            nabc_Left = nabc;
            nabc_Left[0] += 1;
            if (nabc.sum() % 2 == 0)
               nabc_Left[2] += 1;
            nabc_Left /= 2;
         }
         nabc_Right = nabc - nabc_Left;
      }
   } else {
      nabc_Left = nabc;
      nabc_Right = FPow3Index(0,0,0);
      // note: i tried both <= and >=, and splitting off the smaller one into an
      // isolated set worked better (not by very much, though. None of these
      // variants made a very large difference; +/- 30% compared to the straight
      // pow method, maybe. Which is a bit odd, because it *did* make a difference
      // when making all the intermediates... and more a factor 10-kind of
      // difference.
      if (nabc[0] <= nabc[1] && nabc[0] <= nabc[2])
         std::swap(nabc_Left[0], nabc_Right[0]);
      else if (nabc[1] <= nabc[2])
         std::swap(nabc_Left[1], nabc_Right[1]);
      else
         std::swap(nabc_Left[2], nabc_Right[2]);
   }
   
   // compute A^nabc_Left[0] B^nabc_Left[1] C^nabc_Left[2], and same for _Right variant.
//    cx_assert_rt(nabc_Left + nabc_Right == nabc);
   for (int ixyz = 0; ixyz < 3; ++ ixyz)
      cx_assert_rt(nabc_Left[ixyz] + nabc_Right[ixyz] == nabc[ixyz]);
   FPolynomialN const
      &Left = *ComputePowMN(nabc_Left, Int),
      &Right = *ComputePowMN(nabc_Right, Int);
   // and multiply them together
   FPolynomialNCptr
      res = FPolynomialNCptr(new FPolynomialN(Left*Right));
   // remember the intermediate result
   Int[nabc] = res;
   return res;
}


// return A^na * B^nb * C^nc, forming suitable intermediate results
FPolynomialN MultiplyPow(FPolynomialN const &A, size_t na, FPolynomialN const &B, size_t nb, FPolynomialN const &C, size_t nc)
{
   if (na == 0 && nb == 0 && nc == 0)
      // that's 1.
      return FPolynomialN(1, FPolynomialN::FMonomial(0,0,0));
   if (na == 1 && nb == 0 && nc == 0)
      return A;
   if (na == 0 && nb == 1 && nc == 0)
      return B;
   if (na == 0 && nb == 0 && nc == 1)
      return C;

   FPolynomialNPowerM_Map
      Int; // <- cache for intermediate results
   // store A, B, and C themselves, under their respective powers
   Int[FPow3Index(1,0,0)] = FPolynomialNCptr(new FPolynomialN(A));
   Int[FPow3Index(0,1,0)] = FPolynomialNCptr(new FPolynomialN(B));
   Int[FPow3Index(0,0,1)] = FPolynomialNCptr(new FPolynomialN(C));
   return *ComputePowMN(FPow3Index(na,nb,nc), Int); 
}




// return A*B*C, but checks in which order to multiply.
FPolynomialN Multiply(FPolynomialN const &A, FPolynomialN const &B, FPolynomialN const &C)
{
   // make an estimate of multiplication costs for the different orders of
   // multiplication (this does not account for some terms potentially
   // cancelling or combining to common terms, but should provide a reasonable
   // worst-case estimate)
   size_t
      nA = A.size(),
      nB = B.size(),
      nC = C.size();
   size_t
      ThrLarge = 100;
   (void)ThrLarge;
//    if (nA < ThrLarge || nB < ThrLarge || nC < ThrLarge) {
   if (1) {
      size_t
         // compute number of ops for first forming A*B, then, assuming no terms cancel,
         // forming from this the product (A*B)*C
         nOps_AB_C = nA * nB + (nA * nB) * nC,
         nOps_BC_A = nB * nC + (nB * nC) * nA,
         nOps_CA_B = nC * nA + (nC * nA) * nB;
      
      if (nOps_AB_C <= nOps_BC_A && nOps_AB_C <= nOps_CA_B)
         return (A*B)*C;
      else if (nOps_BC_A <= nOps_CA_B)
         return (B*C)*A;
      else
         return (C*A)*B;
   } else {
      // we will apply this with some very large polynomial expansions in the
      // symmetry unique target function computation. In this case it may make
      // sense to calculate the intermediate three intermediate products
      // A*B, B*C, C*A exactly, including the amount of combined/cancelled terms,
      // to check in which order to execute the final (and likely most costly)
      // multiplication.
      //
      // This may seem like overkill at first glance, but for polynomial
      // expansions with 10k terms each in A,B,C, almost the entire calculation
      // cost is in the last multiplication combining all three terms, and doing
      // the two additional two-term multiplications to find the right
      // multiplication order, is not a significant overhead.
      FPolynomialN
         AB = A*B,
         BC = B*C,
         CA = C*A;
      size_t
         // compute cost for scanning the product terms, only. This time
         // we already have all two-term products.
         nOps_AB_C = AB.size() * nC,
         nOps_BC_A = BC.size() * nA,
         nOps_CA_B = CA.size() * nB;
      
      if (nOps_AB_C <= nOps_BC_A && nOps_AB_C <= nOps_CA_B)
         return AB*C;
      else if (nOps_BC_A <= nOps_CA_B)
         return BC*A;
      else
         return CA*B;
   }
}




typedef std::tuple<int, int, int> FMonomialSortKey_Basic;
inline FMonomialSortKey_Basic SortKey_Basic(FMonomialN const &m) {
   return FMonomialSortKey_Basic(m[0], m[1], m[2]);
};


template<unsigned N>
bool HasNegativeParity(unsigned const *ii) {
   unsigned itransp = 0; // number of transpositions in permutation
   for (unsigned i = 1; i < N; ++i)
      for (unsigned j = 0; j < i; ++j)
         itransp += int(ii[i] < ii[j]);
   return bool(unsigned(itransp) & 1);
}

// hm... something might be a bit shady here. This one:
// make && aigg_64_d '{point-group:I; degree:[20]; initial-points: hex-grid{5;1;01c}; print:0; target-space:exact; max-iter:10; export:last.dat}'
// currently goes out of orbit if I reactive the complex SortKey function.
// But there are too many other open Qs at the moment. Will look again later.
// FIXME: investigate if the better SortKey() function is working correctly.

// this monomial sorting key is used in the export of polynomial data. It is meant
// to (1) group like terms together and (2) provide a reasonable order for a schmidt
// orthogonalization
typedef std::tuple<int, ptrdiff_t, int, ptrdiff_t> FMonomialSortKey_Eo;
FMonomialSortKey_Eo SortKey_ExportOrder(FMonomialN const &m) {
   unsigned const g = 3;
   unsigned xyz[g] = {m[0], m[1], m[2]};
   unsigned ii[g] = {0,1,2}; // permutation to arrive at xyz[0] >= xyz[1] >= xyz[2] >= ...
   unsigned mx = 0;
   for (unsigned i = 0; i < g; ++ i)
      mx += unsigned(xyz[i]); // need total degree+1 for accumulator, not max individual entry.
//       mx = std::max(unsigned(mx), unsigned(xyz[i]));
   for (unsigned i = 0; i < g - 1; ++ i)
      for (unsigned j = i+1; j < g; ++ j)
         if (xyz[i] < xyz[j]) { std::swap(ii[i], ii[j]); std::swap(xyz[i], xyz[j]); }
   if (1) {
      mx = 511; // FIXME....
      for (unsigned i = 0; i < g; ++ i)
         cx_assert_rt(xyz[i] <= mx);
   }
   unsigned iperm = 0;
   unsigned deg = 0; // monomial degree
   unsigned AccPerm = 1;
   ptrdiff_t ixyz = 0, AccXyz = 1, ixyz_safety = 0;
   for (unsigned i = g-1; i < g; -- i) {
      deg += xyz[i];
      iperm += ii[i] * AccPerm;
      AccPerm *= g;
      ixyz += xyz[i] * AccXyz;
      ixyz_safety += m[i] * AccXyz; // FIXME: ...
      AccXyz *= ptrdiff_t(mx+1);
   }
   if (HasNegativeParity<g>(ii))
      iperm += AccPerm; // put negative-sign permutations behind positive sign permutations
   // sort first by degree (large first), then by sorted(xyz) (largest powers first; e.g., (10,0,0) before (8,2,0)),
   // then by permutation sign, then by permutation
   return FMonomialSortKey_Eo(-int(deg), -ixyz, iperm, ixyz_safety);
};

#ifdef INCLUDE_ABANDONED
// typedef std::tuple<int, int, int, int, int> FMonomialSortKey_Eo;
// inline FMonomialSortKey_Eo SortKey_ExportOrder(FMonomialN const &m) {
//    unsigned ix = m[0], iy = m[1], iz = m[2];
//    unsigned iix = 0, iiy = 1, iiz = 2;
//    if (iy < iz) { std::swap(iy, iz); std::swap(iiy, iiz); } // now: iy >= iz
//    if (ix < iy) { std::swap(ix, iy); std::swap(iix, iiy); } // now: ix >= max(iy,iz), but possibly iy < iz
//    if (iy < iz) { std::swap(iy, iz); std::swap(iiy, iiz); }
//    int iperm = int(9*iix + 3*iiy + iiz);
//    unsigned itransp = int(iix < 0) + int(iiy < 1) + int(iiz < 2); // number of transpositions in permutation
//    if (bool(unsigned(itransp) & 1))
//       iperm += 3*3*3; // put negative-sign permutations behind positive sign permutations
//    return FMonomialSortKey_Eo(int(ix+iy+iz), -int(ix), -int(iy), -int(iz), iperm);
// };
#endif // INCLUDE_ABANDONED

FScalar CalcMonomialUnitSphereIntegralRaw(size_t ix, size_t iy, size_t iz);

// typedef std::tuple<double, int, ptrdiff_t, ptrdiff_t> FMonomialSortKey_Sif;
typedef std::tuple<ptrdiff_t, ptrdiff_t> FMonomialSortKey_Sif;
FMonomialSortKey_Sif SortKey_SmallIntegralsFirstOrder(FMonomialN const &m) {
   unsigned const g = 3;
   unsigned xyz_orig[g] = {m[0], m[1], m[2]};
   unsigned xyz_sorted[g] = {xyz_orig[0], xyz_orig[1], xyz_orig[2]};
   unsigned ii[g] = {0,1,2}; // permutation to arrive at xyz_sorted[0] <= xyz_sorted[1] <= xyz_sorted[2] <= ...
   unsigned deg = 0;
   if (0) {
      for (unsigned i = 0; i < g; ++ i)
         deg += unsigned(xyz_orig[i]);
   } else {
      deg = 511;
      // fix the "degree".. so that we can compare monomials of different
      // degrees directly by their packed index (need to have compatible
      // accumulators for that).
      for (unsigned i = 0; i < g; ++ i)
         cx_assert_rt(xyz_orig[i] <= deg);
   }
   for (unsigned i = 0; i < g - 1; ++ i)
      for (unsigned j = i+1; j < g; ++ j)
         if (xyz_sorted[i] > xyz_sorted[j]) { std::swap(ii[i], ii[j]); std::swap(xyz_sorted[i], xyz_sorted[j]); }
   ptrdiff_t
      ixyz_sorted = 0, ixyz_orig = 0, iAccXyz = 1;
   for (unsigned i = g-1; i < g; -- i) {
      ixyz_sorted += xyz_sorted[i] * iAccXyz;
      ixyz_orig += xyz_orig[i] * iAccXyz;
      iAccXyz *= ptrdiff_t(deg+1);
   }
   if (1) {
      // other properties being equal, put negative-sign permutations
      // behind positive-sign permutations
      ixyz_sorted *= 2;
      if (!HasNegativeParity<g>(ii))
         ixyz_sorted += 1;
   }
#if 0
   double
      SphInt(0);
   if (0) {
      // if on, explicitly calculate the sphere integral (and truncate it to doubles)
      SphInt = double(CalcMonomialUnitSphereIntegralRaw(xyz_orig[0], xyz_orig[1], xyz_orig[2]));
   }
   // otherwise: try to get the monomials with the largest *minimum* exponent
   // to the front. things like x^20 y^20 z^20 make small contributions because
   // they max out at normal direction (1,1,1), and the actual points on the sphere
   // have cartesian components < 1 at those points (because x^2+y^2+z^2 should be 1,
   // meaning that all x,y,z = 1/sqrt(3) in the worst case)
   return FMonomialSortKey_Sif(SphInt, -int(deg), -ixyz_sorted, ixyz_orig);
#else
   return FMonomialSortKey_Sif(-ixyz_sorted, ixyz_orig);
#endif
   // Note: with this ordering this comes out pretty well also without actually computing
   // the integrals--even with groups requiring mixed cartesian degrees. E.g.:
   //
   //  Target function list ():
   //      1  M(xyz)[ 6, 6, 8]  1/nfac = 4.09e-06  degn = 1
   //      2  M(xyz)[ 5, 7, 7]  1/nfac = 7.35e-06  degn = 1
   //      3  M(xyz)[ 5, 5, 9]  1/nfac = 9.81e-06  degn = 1
   //      4  M(xyz)[ 4, 8, 8]  1/nfac = 5.74e-06  degn = 1
   //      5  M(xyz)[ 4, 6,10]  1/nfac = 7.38e-06  degn = 1
   //      6  M(xyz)[ 4,10, 6]  1/nfac = 7.38e-06  degn = 1
   //      7  M(xyz)[ 4, 4,12]  1/nfac = 1.63e-05  degn = 1
   //      8  M(xyz)[ 3, 7, 9]  1/nfac = 1.48e-05  degn = 1
   //      9  M(xyz)[ 3, 9, 7]  1/nfac = 1.48e-05  degn = 1
   //     10  M(xyz)[ 3, 5,11]  1/nfac = 2.47e-05  degn = 1
   //     11  M(xyz)[ 3,11, 5]  1/nfac = 2.47e-05  degn = 1
   //     12  M(xyz)[ 3, 3,13]  1/nfac = 7.46e-05  degn = 1
   //     13  M(xyz)[ 2, 8,10]  1/nfac = 1.74e-05  degn = 1
   //     14  M(xyz)[ 2,10, 8]  1/nfac = 1.74e-05  degn = 1
   //     15  M(xyz)[ 2, 6,12]  1/nfac = 2.74e-05  degn = 1
   //     16  M(xyz)[ 2,12, 6]  1/nfac = 2.74e-05  degn = 1
   //     17  M(xyz)[ 2, 4,14]  1/nfac = 7.16e-05  degn = 1
   //     18  M(xyz)[ 2,14, 4]  1/nfac = 7.16e-05  degn = 1
   //     19  M(xyz)[ 2, 2,16]  1/nfac = 3.63e-04  degn = 1
   //     20  M(xyz)[ 1, 9, 9]  1/nfac = 6.09e-05  degn = 1
   //     21  M(xyz)[ 1, 7,11]  1/nfac = 7.62e-05  degn = 1
   //     22  M(xyz)[ 1,11, 7]  1/nfac = 7.62e-05  degn = 1
   //     23  M(xyz)[ 1, 5,13]  1/nfac = 1.53e-04  degn = 1
   //     24  M(xyz)[ 1,13, 5]  1/nfac = 1.53e-04  degn = 1
   //     25  M(xyz)[ 1, 3,15]  1/nfac = 5.39e-04  degn = 1
   //     26  M(xyz)[ 1,15, 3]  1/nfac = 5.39e-04  degn = 1
   //     27  M(xyz)[ 1, 1,17]  1/nfac = 4.45e-03  degn = 1
   //     28  M(xyz)[ 0,10,10]  1/nfac = 1.81e-04  degn = 1
   //     29  M(xyz)[ 0, 8,12]  1/nfac = 2.21e-04  degn = 1
   //     30  M(xyz)[ 0,12, 8]  1/nfac = 2.21e-04  degn = 1
   //     31  M(xyz)[ 0, 6,14]  1/nfac = 4.11e-04  degn = 1
   //     32  M(xyz)[ 0,14, 6]  1/nfac = 4.11e-04  degn = 1
   //     33  M(xyz)[ 0, 4,16]  1/nfac = 1.24e-03  degn = 1
   //     34  M(xyz)[ 0,16, 4]  1/nfac = 1.24e-03  degn = 1
   //     35  M(xyz)[ 0, 2,18]  1/nfac = 7.12e-03  degn = 1
   //     36  M(xyz)[ 0,18, 2]  1/nfac = 7.12e-03  degn = 1
   //     37  M(xyz)[ 0, 0,20]  1/nfac = 1.56e-01  degn = 1
   // ░ G = I   l = 20   #monomials[lmin: 19  lmax: 20  all: 441   req: 37    nonz-P_G = 36   ] #linearly-independent = 8
   //
   // I think I'll turn off the sphere integral and degree components for that reason.
};


bool FMonomialN_ExportOrderPred::operator() (FMonomialN const &A, FMonomialN const &B) const {
#if 1
   return SortKey_ExportOrder(A) < SortKey_ExportOrder(B);
#else
   return SortKey_Basic(A) < SortKey_Basic(B);
#endif // #if 0
}

bool FMonomialN_SmallIntegralsFirstOrder::operator() (FMonomialN const &A, FMonomialN const &B) const {
   return SortKey_SmallIntegralsFirstOrder(A) < SortKey_SmallIntegralsFirstOrder(B);
}

bool FMonomialN_BasicOrder::operator() (FMonomialN const &A, FMonomialN const &B) const {
   return SortKey_Basic(A) < SortKey_Basic(B);
}
