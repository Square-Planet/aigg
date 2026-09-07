#ifndef CX_MONOMIAL_MATH_H
#define CX_MONOMIAL_MATH_H

#include <type_traits>
#include <array>
#include <cmath>
#include <stdexcept>
#include <ostream>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "format.h"

namespace ct {
using std::size_t;

/// return number of non-negative integer tuples (i_0,i_1,...,i_{TupleLength-1}) with
/// \sum i_\alpha = L. E.g., (L+1)*(L+2)/2 for TupleLength = 3 (same as nCartY in IR).
///
/// Note: this one has an overflow check (enabled at runtime)---it is not meant
/// to be fast.
template<class uint_type = size_t>
uint_type NumberOfTuples_DegreeEqual(uint_type L, uint_type TupleLength) {
//    static_assert(sizeof(L) >= sizeof(TupleLength), "unusual combination of integer types. intended?");
   uint_type
      Num = 1, Denom = 1;
   for (uint_type i = 1; i < TupleLength; ++ i) { // <-- note loop starting at 1, not 0.
      uint_type LastNum = Num;
      Num *= L + i;
      if (!(Num >= LastNum))
         throw std::overflow_error(fmt::format("NumberOfTuples_DegreeEqual<{}>({},{}): {} integer overflow.",
               sizeof(uint_type), L, TupleLength, (std::is_signed<uint_type>::value? "signed" : "unsigned")));
         // ^-- note: signed integer overflow is undefined behavior. This check
         //     may not be triggered in this case---compiler decides. (According
         //     to the C++ standard, it would also be legal if the compiler
         //     builds in code into the executable which, if this happens,
         //     organizes an angry protesting mob in front of your house
         //     chanting "NO MORE OVERFLOW!!!" "YOU TEMPTED FATE! SIGNED INTEGER
         //     OVERFLOWER!!" and similar things).
      Denom *= i;
   }
   return Num / Denom;
}

/// return number of non-negative integer tuples (i_0,i_1,...,i_Base) with
/// \sum i_\alpha <= L. E.g., (L+1)*(L+2)*(L+3)/6 for Base = 3 (same as nCartX in IR)
template<class uint_type = size_t>
inline uint_type NumberOfTuples_DegreeUpTo(uint_type L, uint_type Base) {
   return NumberOfTuples_DegreeEqual(L, Base+1);
}


/// For a sequence of unsigned integer exponents e[0], e[1], e[2], ..., compute
/// and return the multinomial coefficient
///
///    (e[0] + e[1] + ... + e[B])!
///    ──────────────────────────
///     e[0]!  e[1]!  ...   e[B]!
///
/// Arguments:
/// - itExp: input iterator to start of exponent sequence. *itExp must
///   be an unsigned integer type.
/// - itEnd: iterator to end of exponent sequence (exclusive)
///
/// Comments:
/// - Formally, the multinomial coefficient is always an integer.
/// - Less formally, it is not at all difficult to come up with applications
///   in which the coefficient would exceed the limits of standard integral
///   types and cause overflow errors in the computation.
/// - So better use with floating point types, unless you are
///   absolutely sure about the range covered in applications...
template<class FScalar, class FInputIterator>
FScalar EvalMultinomialCoefficient(FInputIterator itExp, FInputIterator itEnd)
{
   if (itExp == itEnd)
      throw std::runtime_error("EvalMultinomialCoefficient(): there must be at least one exponent");

   // For the base functions used, see: https://www.boost.org/doc/libs/1_76_0/libs/math/doc/html/special.html
   //
   // We mean to compute the multinomial coefficient
   //
   //    (e[0] + e[1] + ... + e[B])!
   //    ──────────────────────────
   //     e[0]!  e[1]!  ...   e[B]!
   //
   // We do this incrementally: We compute the coefficient for exponent subsets
   // e[0:1], e[0:2], e[0:3] ... each time fixing up the intermediate
   // coefficient to go from e[0:b-1] to e[0:b].
   //
   // We start with the first term, which is just:
   //
   //    (e[0])!
   //    ─────── = 1
   //     e[0]!
   //
   FScalar
      Coeff(1);
   unsigned
      // lnfa = "last numerator's factorial argument".
      // This accumulates the sum e[0] + e[1] + ..., the factorial of which is
      // (implicitly) the current numerator of this->coeffs
      lnfa = *itExp;
   ++itExp;
   for (; itExp != itEnd; ++itExp) {
      // current coeff is: (b = itExp - itBeg)
      //
      //     (e[0] + e[1] + ... + e[b-1])!
      //     ────────────────────────────
      //     (e[0]!  e[1]! ...  e[b-1]!)
      //
      // So to incorporate the next base exponent, we need to multiply it by
      //
      //     (e[0] + e[1] + ... + e[b-1] + e[b])!         1
      //     ────────────────────────────────────   *   ──────
      //     (e[0] + e[1] + ... + e[b-1])!               e[b]!
      //
      //     ^-- for new numerator: cancel old          ^-- for new denominator:
      //         numerator and replace by new               it just gets one
      //         factorial (e[0] + ... + e[b])!             more e[b]! factor
      //
      // The sum e[0] + ... + e[b-1] in the denominator is what we currently
      // have in `lnfa`. So the entire factor is just
      //
      //    (lnfa + e[b])!
      //    ──────────────
      //     lnfa!  e[b]!
      //
      // ...which is another binomial coefficient.
      // Interesting. Didn't expect that. But argument looks fine and
      // I cross checked:
      //
      //    FunctionExpand[FullSimplify[Binomial[i+j,i]*Binomial[i+j+k,k]*Binomial[i+j+k+l,l]*Binomial[i+j+k+l+m,m]]]
      //
      // indeed evaluates to Γ[1+i+j+k+l+m]/(Γ[1+i] Γ[1+j] Γ[1+k] Γ[1+l] Γ[1+m]).
      // So... easier than expected. No need for rising factorials, I guess.
      unsigned eb = *itExp; // e[b]
      Coeff *= boost::math::binomial_coefficient<FScalar>(lnfa + eb, eb);
      lnfa += eb;
   }
   return Coeff;
}



template<unsigned BaseSizeT, class FScalar>
struct TMultinomialTerm {
   enum { BaseSize = BaseSizeT };
   typedef unsigned short
      exponent_type;
   typedef std::array<exponent_type, BaseSizeT>
      FExpArray;
   FExpArray
      // exponents for base terms:  base[0]^exp[0] * base[1]^exp[1] * ...
      exp;
   FScalar
      // trinomial expansion coefficient: (ix + iy + iz)!/(ix! iy! iz!)
      coeff;

   // evaluate this->coeff from this->exp.
   void _Finalize();
   void Print(std::ostream &xout) const;
public:
   exponent_type const &operator [] (size_t i) const { return exp[i]; }
   exponent_type &operator [] (size_t i) { return exp[i]; }
   size_t size() const { return exp.size(); }
   typename FExpArray::const_iterator begin() const { return exp.begin(); }
   typename FExpArray::const_iterator end() const { return exp.end(); }
};

// represents the expansion of (base[0] + base[1] + ... + base[BaseSizeT-1])^L into
// \sum_{i} t[i].coeff \prod_{a=0}^{BaseSizeT} base[a]^{t[i].exp[a]}
template<unsigned BaseSizeT, class FScalar>
struct TMultinomialExpansion1 : public std::vector<TMultinomialTerm<BaseSizeT, FScalar> > {
   enum { BaseSize = BaseSizeT };
   typedef std::vector<TMultinomialTerm<BaseSizeT, FScalar> >
      FTermList;
   typedef typename FTermList::value_type
      FTerm;
   explicit TMultinomialExpansion1(unsigned Power);

   void Print(std::ostream &xout, std::string const &ind) const;
   // return the sum of exponents for each term
   unsigned Degree() const { return m_Power; };
protected:
   unsigned m_Power;

   void GenTermsR(unsigned iBase, FTerm &PartialTerm, size_t nPowerRemaining);
};


// Stores all multinomial expansions for degrees 0,1,...,L.
template<unsigned BaseSizeT, class FScalar>
struct TMultinomialExpansionN : public std::vector<TMultinomialExpansion1<BaseSizeT, FScalar> > {
   typedef std::vector<TMultinomialExpansion1<BaseSizeT, FScalar> >
      FExpansion1List;
   typedef typename FExpansion1List::value_type
      FExpansion1;
   enum {
      BaseSize = BaseSizeT,
      InvalidPower = unsigned(-1)
   };
   // if MaxPower != InvalidPower, construct and cache multinomial expansions
   // for orders L = 0,1,...,MaxPower (inclusive),
   explicit TMultinomialExpansionN(unsigned MaxPower = InvalidPower);

   void Init(unsigned MaxPower);

   bool HaveQ(unsigned p) const { return p < m_LastPower;  }
   void Print(std::ostream &xout, std::string const &ind) const;
protected:
   unsigned
      // first power we do not have tabulated anymore; i.e., valid iPower fulfill 0
      // <= iPower < m_LastPower. Since 0 constitutes a valid power, this is also
      // the total number of expansions we have tabulated.
      m_LastPower;
};


template<unsigned BaseSizeT, class FScalar>
TMultinomialExpansionN<BaseSizeT,FScalar>::TMultinomialExpansionN(unsigned MaxPower)
   : m_LastPower(0)
{
   Init(MaxPower);
}


template<unsigned BaseSizeT, class FScalar>
void TMultinomialExpansionN<BaseSizeT,FScalar>::Init(unsigned MaxPower)
{
   if (!this->empty()) {
      // delete all elements -- clear memory
      FExpansion1List dummy;
      static_cast<FExpansion1List*>(this)->swap(dummy);
   }
   assert(this->empty());
   if (MaxPower == InvalidPower) {
      // place holder for "don't do anything if default constructed without arguments".
      m_LastPower = 0;
      return;
   }
   m_LastPower = MaxPower + 1;
   // make new elements
   this->reserve(m_LastPower);
   for (unsigned Power = 0; Power < m_LastPower; ++ Power) {
      this->emplace_back(Power);
   }
   assert(this->size() == m_LastPower);
}



template<unsigned BaseSizeT, class FScalar>
void TMultinomialTerm<BaseSizeT, FScalar>::_Finalize()
{
   static_assert(BaseSizeT >= 1, "Tried to instanciate a multinomial expansion (b[0] + b[1] + .. + b[K-1])^N with K = 0 base terms.");
   this->coeff = EvalMultinomialCoefficient<FScalar>(this->begin(), this->end());
}



template<unsigned BaseSizeT, class FScalar>
TMultinomialExpansion1<BaseSizeT, FScalar>::TMultinomialExpansion1(unsigned Power)
   : m_Power(Power)
{
   this->reserve(NumberOfTuples_DegreeEqual<size_t>(m_Power, BaseSize));

   {
      FTerm PartialTerm; // <-- contains random data. Will be filled by GenTermsR.
      GenTermsR(0, PartialTerm, Power);
   }

   // sort expansion terms into order of increasing coefficient size
   // (better to start accumulating with small ones in linear sum)
   std::stable_sort(this->begin(), this->end(),
      [&](FTerm const& A, FTerm const &B) {
         return A.coeff < B.coeff;
      });
}

template<unsigned BaseSizeT, class FScalar>
void TMultinomialExpansion1<BaseSizeT, FScalar>::GenTermsR(unsigned iBase, FTerm &PartialTerm, size_t nPowerRemaining)
{
   if (nPowerRemaining < 0)
      // already distributed too many exponents --- no valid terms can be assembled from this.
      return;
   if (iBase == BaseSizeT) {
      // exponents have been assigned to all base elements.
      // Is the sum of what was assigned equal to m_Power?
      if (nPowerRemaining == 0) {
         // that's a valid combination of exponents. Compute actual
         // coefficient and add it to the list.
         PartialTerm._Finalize();
         this->push_back(PartialTerm);
      }
   } else if (iBase + 1 == BaseSizeT) {
      // last base: whateve was left undistributed in terms of exponents
      // needs to go to this one to arrive at \sum_i pow[i] = m_Power.
      // So no loop needed.
      assert(nPowerRemaining <= m_Power);
      PartialTerm[iBase] = typename FTerm::exponent_type(nPowerRemaining);
      GenTermsR(iBase+1, PartialTerm, 0);
   } else if (iBase < BaseSizeT) {
      for (size_t i = 0; i <= size_t(nPowerRemaining); ++ i) {
         // assign `i` as exponent to current base, and recurse down to distribute
         // the remaining (nPowerRemaining - i) exponents to the other bases.
         PartialTerm[iBase] = typename FTerm::exponent_type(i);
         GenTermsR(iBase+1, PartialTerm, nPowerRemaining - i);
      }
   } else {
      // iBase > BaseSizeT not supposed to happen.
      assert(0);
   }
}


template<unsigned BaseSizeT, class FScalar>
void TMultinomialTerm<BaseSizeT, FScalar>::Print(std::ostream &xout) const
{
   xout << "(";
   for (size_t i = 0; i < BaseSizeT; ++i) {
      if (i != 0) xout << ", ";
      xout << this->exp[i];
   }
   xout << "): " << this->coeff;
}


template<unsigned BaseSizeT, class FScalar> std::ostream &operator << (std::ostream &xout, TMultinomialTerm<BaseSizeT, FScalar> const &mt)
{
   mt.Print(xout);
   return xout;
}


template<unsigned BaseSizeT, class FScalar>
void TMultinomialExpansion1<BaseSizeT, FScalar>::Print(std::ostream &xout, std::string const &ind) const
{
   xout << fmt::format("{}MultinomialTerms<{}>(L = {}) = {{\n", ind, BaseSizeT, m_Power);
   for (FTerm const &mt : *this) {
      xout << ind << "   " << mt << ",\n";
   }
   xout << ind << "}\n";
}


template<unsigned BaseSizeT, class FScalar>
void TMultinomialExpansionN<BaseSizeT, FScalar>::Print(std::ostream &xout, std::string const &ind) const
{
   for (FExpansion1 const &e : *this) {
      e.Print(xout, ind);
   }
}

} // namespace ct

#endif // CX_MONOMIAL_MATH_H
