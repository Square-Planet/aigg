#ifdef INCLUDE_ABANDONED
// make && time OMP_NUM_THREADS=1 aigg '{ point-group: icosahedral; degree: [24]; max-step: 1e-2; initial-points: hex-grid{4;1;01c} }'
// make && time OMP_NUM_THREADS=12 aigg '{ point-group: icosahedral; degree: [10,step:+2,-1,-2]; max-step: 1e-2; initial-points: hex-grid{7;1;01c}; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes} }'
// test: make && aigg '{ point-group: icosahedral; degree: [24]; max-step: 1e-2; initial-points: hex-grid{4;1;01c} }'
// (that one *should* integrate l=24 exactly)
// for full incremental search:
// make && time OMP_NUM_THREADS=12 aigg '{ point-group: icosahedral; degree: [4,step:+2,-1,-2]; max-step: 1e-2; initial-points: hex-grid{4;1;01c}; print: 0; max-it: 128 }' > /tmp/aigg.log

// ...and maybe something completely unrelated is also still wrong. This one:
// make && OMP_NUM_THREADS=12 aigg '{point-group: icosahedral; max-it: 128; screen-target-fn: no; degree: [5,step:+2,-1,-2]; print: 1; max-step: {1e8*residual; max:1e-3; min:1e-8}; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{9;4;01c}}'
// used to integrate up to L=62 (I stole the settings from a hex grid rule which did get L=62 right, it's in runs somewhere):
//       ~/dev/sphere_grids> cat ./runs/icosahedral_S9_T4_v2.dat
//       # icosahedral integration rule: npts = 1332, lmax = 62, wtsp = 1.811, res = 7.49e-29
//       # generated 2020-04-17 04:54 (UTC) by kgaigg v@Pr0gVeR (compiled: Apr 16 2020)
//       # options: '{point-group: icosahedral; max-it: 128; degree: [5,step:+2,-1,-2]; max-step: {1e8*residual; max:1e-3; min:1e-8}; export: {file: runs/icosahedral_S9_T4_v2.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{9;4;01c}}'
// - That's even a 2020-04-17 version, so certainly not very old.
// - But nowadays I can't get it beyond L=55, regardless of settings here:
//   + It's not the Scd or symmetrization for this.
//   + it's also not the prescreening of monomials, or the inclusion or exclusion of extra monomials in the simple monomial list function
//   + It could be the nFreeParameters business and screening of sigma values in MakeStep;
//     I ran through some examples with print:3, and saw some rather fishy looking
//     ThrSigB thresholds. Didn't yet check in detail, though. Major refactoring first.
//   + oops... several tests suddenly started to work again after changing maxstep from 1e-2
//     to 1e-4 8). Oh my... I did redefine the steps in terms of orthogonalized residuals
//     at some point, and this *did* change the step definitions...
//
//     E.g., even without symmetrization and simple monomial enum, this one:
//
//       make && time OMP_NUM_THREADS=12 aigg '{ point-group: icosahedral; degree: [6,step:+2,-1,-2]; thr-res:1e-27; max-step: 1e-4; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{9;1;01c}; print: 0; max-it: 128 }'
//
//     goes through to npts = 912 (lmax = 49, wtsp = 1.709, res = 1.14e-29, mxe = 4.85e-30)
//     (best: ./runs/icosahedral_S9_T1_v2.dat:# icosahedral integration rule: npts = 912, lmax = 51, wtsp = 1.230, res = 3.10e-27)
//     and the 7/1 one to l=39 (there were some 40s/41s for that, but from old time without
//     residual cross-verification, so I'm not sure about those)
//
//     And in both cases the only real trick to change was the maximum step
//     size, which with the aggressive group averaging would anyway not be that
//     big of a deal...
//
//     Anyway: some more program cleanups now, then re-check the convergence stuff again.
#endif // INCLUDE_ABANDONED


#include "Aigg.h"
#include "TargetFunctions.h"

#include "CxIo.h" // for FFileLocator
#include <fstream> // for writing group-symmetrized target function lists
#include <algorithm> // for stable_sort
#include <array>

#include <stdexcept>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "PointCloud.h"
#include "SymmetryGroups.h"
#include "SearchOptions.h"
#include "MathOps.h"

#include "CxTiming.h"
#include "CxIterTools.h"
#include "CxMonomialMath.h"


using cx::iter::enumerate;
using ct::TMultinomialExpansionN;
using ct::TMultinomialExpansion1;
using ct::TMultinomialTerm;

// Note: should probably move these guys (and the prescreen thing) into
// a new FTargetFnOptions object in SearchOptions.h etc.

const bool g_UseGroupSymmetrizedFn = false;
// bool g_UseGroupSymmetrizedFn = true;
// bool g_TrackEquivMonomials = true;
// bool g_UseGroupSymmetrizedFn = true;
const bool g_TrackEquivMonomials = false;
// bool g_UseGroupSymmetrizedFn = false;
// bool g_UseSimpleMonomialEnum = true;
// bool g_UseSimpleMonomialEnum = false;
// const bool g_UseSimpleMonomialEnum = true; // <- much faster for large rules, and with the optimal prescreening we now have, and the PtLDLtP decomp for KtSmK, it doesn't hurt at all. That's what I used for the ico-19-1-refine which is still running, but already made proved the 3812pt rule to be of order lmax=100 min.
enum FTargetFunctionNormalization {
   TARGETFN_NormalizeSymmetrized,
   TARGETFN_NormalizeRawMonomials,
   TARGETFN_DoNotNormalize
};
FTargetFunctionNormalization
//    s_TargetFnNormalization = TARGETFN_DoNotNormalize;
   s_TargetFnNormalization = TARGETFN_NormalizeRawMonomials;
//    s_TargetFnNormalization = TARGETFN_NormalizeSymmetrized;

static int
//    s_TimeGroupAveraging = 2;
   s_TimeGroupAveraging = 0;

// i think there *should* be some weighting of residuals by this, but I couldn't
// quite get it right:
// - Currently, for example, the column-scaling for the SVD
//   means that both gradients and J^t Scd J matrices get effectively scaled by
//   degeneracy^2, so that this actually has no effect.
// - Residuals are multiplied by S^{-1} to form the non-redundant function
//   subspace when making the updates. I am currently unsure if that is correct
//   or not (I think, technically, the residuals and rhs already are effectively
//   P_G-projected at this point, so I am unsure if they can scaled there by
//   a degeneracy factor)
// - Anyway... needs some more intense thinking, I would guess.
//   It is also quite possible that we really only want to scale the final
//   residuals and possibly steps by the degeneracy factor, and otherwise
//   leave the process alone.
// - Needs understanding of what *would* be happening if we did *not* filter out
//   symmetry-equivalent monomials, but left them in (still group-averaged, and
//   zeros removed, though. That can probably be checked explicitly.)
//   This setup seemed manageable:
//
//   make && time OMP_NUM_THREADS=1 aigg '{ point-group: icosahedral; degree: [14]; max-step: 1e-2; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{3;1;01c}; print: 3; max-it: 128 }' > /tmp/brk_deg.txt
//
//   has relatively small but non-trivial target fn space (with symmetry),
//   but still a non-zero residual.
//
// - in any case... atm turning this on sometimes yields incorrect results
//   (e.g., this one stops at lmax = 17 instead of 23:
//
//   make && time OMP_NUM_THREADS=12 aigg '{ point-group: icosahedral; degree: [6,step:+2,-1,-2]; max-step: 1e-2; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{5;1;01c}; print: 0; max-it: 128 }'
//
// - FIXME, IMPORTANT: it is possible that I am still missing equivalent
//   functions because my monomials are not normalized. As a result, one could
//   possibly obtain multiple linearly dependent target functions (after symmetrization!)
//   and therefore a singular Scd matrix! That *would* be a possible explanation for the
//   order dependency below...


//   FIXME: Find out what is going on and fix it.
// static bool s_MultiplyByDegeneracy = true;
static bool s_MultiplyByDegeneracy = false;

static FScalar degn_or_one(FScalar const &degeneracy) {
   if (s_MultiplyByDegeneracy)
      return degeneracy;
   else
      return FScalar(1);
}


FScalar ComputeOverlapRaw(FMonomial const &mi, FPolynomialN const &mj);
FScalar ComputeOverlapRaw(FPolynomialN const &poly_i, FPolynomialN const &poly_j);


// static std::string FmtTiming(fmt::BasicStringRef<char> Name, double fTimeInSeconds)
// {
//    return fmt::format(p1TimingFmt, fmt::format(" Time for {}:", Name), fTimeInSeconds);
// }


FTransformedPointInfo::FTransformedPointInfo(const FPointArray &UntransformedPoints_, size_t lmax_)
   : m_OrigPoints(UntransformedPoints_), lmax(lmax_), nPoints(UntransformedPoints_.cols())
{
   TransformedPoints.resize(3, nPoints);
   for (size_t ixyz = 0; ixyz != 3; ++ ixyz)
      XyzPow[ixyz].resize(1+lmax, nPoints);
}


void FTransformedPointInfo::TransformPoints(FRotationMatrix const &SymOp)
{
   assert(TransformedPoints.cols() == XyzPow[0].cols() && XyzPow[0].cols() == XyzPow[1].cols() && XyzPow[1].cols() == XyzPow[2].cols());
   assert((1+lmax) == size_t(XyzPow[0].rows()) && (1+lmax) == size_t(XyzPow[1].rows()) && (1+lmax) == size_t(XyzPow[2].rows()));

   TransformedPoints = SymOp * m_OrigPoints;
   // format: (lmax+1, nPoints)
   for (size_t ixyz = 0; ixyz != 3; ++ ixyz) {
      XyzPow[ixyz].row(0).setOnes();  // power 0 -> 1.0
      // higher powers:  pow[l] = xyz[ixyz] * pow[l-1]
      FDenseVector
         XyzForSeeds = TransformedPoints.row(ixyz);
      for (size_t l = 1; l <= lmax; ++ l)
         XyzPow[ixyz].row(l) = (XyzForSeeds.transpose().array() * XyzPow[ixyz].row(l-1).array()).matrix();
   }
}


FMonomialEvalInfo FTransformedPointInfo::TransformedPointInfo(size_t iPoint) const
{
   FMonomialEvalInfo
      r = {&XyzPow[0](0, iPoint), &XyzPow[1](0, iPoint), &XyzPow[2](0, iPoint)};
   return r;
}



// size_t DoubleFactR(ptrdiff_t n) {
//    if (n <= 1) return 1;
//    return n * DoubleFactR(n - 2);
// }
// ^- may overflow for higher lmax (starting at ~35)

// FScalar DoubleFactR(FScalar n) {
//    if (n <= FScalar(1e-8))
//       return FScalar(1);
//    else
//       return n * DoubleFactR(n - FScalar(2));
// }

// FScalar DoubleFactR(int n) {
//    assert(n >= -1); // -2 is some sort of infinity, and -3 is -1...
//    if (n == -1 || n == 0 || n == 1)
//       return 1;
//    else
//       return boost::math::double_factorial<FScalar>(unsigned(n));
// }



// return (n + d)!!/n!!. d be even.
FScalar DoubleFactRatio(int n, int d) {
   if ((!((d % 2) == 0)) || d < 0)
      throw std::invalid_argument(fmt::format("DoubleFactRatio(n,d) requires d to be even d >= 0 (arguments: n={}, d={}", n, d));
   FScalar
      out = 1;
   for (int i = n + 2; i <= n + d; i += 2)
      out *= FScalar(i);
   return out;
}

FScalar DoubleFactR(int n) {
   if (n == -1 || n == 0) return FScalar(1);
   int n0 = (n % 2 == 1)? -1 : 0;
   return DoubleFactRatio(n0, n - n0);
}



static size_t nCartX(size_t lmax_) {
   ptrdiff_t lmax = ptrdiff_t(lmax_);
   return size_t((lmax+1)*(lmax+2)*(lmax+3)/6);
}


template<class T, class Pred = std::less<T>>
void SortArguments(T &i, T &j, T &k, Pred const &pred = Pred()) {
   using std::swap;
   if (!pred(i,j)) swap(i,j); // now i <= j
   if (!pred(j,k)) swap(j,k); // now i <= k && j <= k
   if (!pred(i,j)) swap(i,j); // now i <= j && i <= k && j <= k
}


FScalar CalcMonomialUnitSphereNormSq(unsigned ix, unsigned iy, unsigned iz) {
   return CalcMonomialUnitSphereIntegralRaw(2*ix, 2*iy, 2*iz);
}

FScalar CalcMonomialUnitSphereNorm(unsigned ix, unsigned iy, unsigned iz) {
   using std::sqrt;
   return sqrt(CalcMonomialUnitSphereNormSq(ix, iy, iz));
}

FScalar CalcMonomialUnitSphereNorm(FMonomialN const &m) {
   return CalcMonomialUnitSphereNorm(m[0], m[1], m[2]);
}

FScalar CalcMonomialUnitSphereNormSq(FMonomialN const &m) {
   return CalcMonomialUnitSphereNormSq(m[0], m[1], m[2]);
}

FScalar CalcMonomialUnitSphereIntegralRaw(unsigned ix, unsigned iy, unsigned iz)
{
   if ((ix % 2) == 1 || (iy % 2) == 1 || (iz % 2) == 1)
      return FScalar(0);
   FScalar
      out;
   if (1) {
      int
         f = int(ix) - 1, // note: f,g,h all odd.
         g = int(iy) - 1,
         h = int(iz) - 1;
      SortArguments(f, g, h);
      FScalar
         IntRatio = 1;
      IntRatio *= DoubleFactRatio(f,g+1)/DoubleFactR(g); // now (f+g+1)!!/(f!! g!!)
      IntRatio *= DoubleFactRatio(f+g+1,h+1)/DoubleFactR(h); // now (f+g+h+2)!!/(f!! g!! h!!)
      IntRatio *= DoubleFactRatio(f+g+h+2,2); // now (f+g+h+4)!!/(f!! g!! h!!)
      // ^-- this is here because I got very confused and thought "IntRatio" should evaluate
      //     to an actual integer in all cases (as it would with normal factorials!).
      //     In this case some do not:
      //     aigg_64_d '{ point-group: icosahedral; degree: [4,step:+1,-1]; max-step: 1e4*residual; initial-points: hex-grid{5;1;01c}; print: {target-fn:all} }'
      //     however, I got this wrong. For some combinations (e.g., [ix iy, iz] = [14,14,0]), one *really*
      //     does get fractional ratios. Because then there are mixed even-series and odd-series
      //     double factorials, which do not cancel each other. For test:
      //     IntS2[i_,j_,k_]:=Factorial2[i+j+k+1]/(Factorial2[i-1]*Factorial2[j-1]*Factorial2[k-1])
      out = 1/IntRatio;
#ifdef INCLUDE_ABANDONED
//       for (int i = -1; i <= 10; ++i)
//          io.Write("{:2}!! = {}", i, DoubleFactR(i));
//       for (int j = 4; j <= 5; j += 1) {
//          for (int i = 0; i <= 6; i += 2) {
//             io.Write("({}+{})!!/{}!! = {}", j,i,j,DoubleFactRatio(j,i));
//          }
//       }
//       io.Flush();
//       throw std::runtime_error("absicht.");
#endif // INCLUDE_ABANDONED
   } else if (0) {
      // note: 4pi prefactor omitted (our grids are supposed to integrate to 1, not to 4pi)
      out = FScalar(
            DoubleFactR(int(ix) - int(1)) *
            DoubleFactR(int(iy) - int(1)) *
            DoubleFactR(int(iz) - int(1))
         ) / FScalar(DoubleFactR(int(ix + iy + iz + 1)));
      // ^-- for very high orders, this can overflow, even with the large floats we use.
      //    Also, boost's double-factorial says[1]:
      //
      //       """The double factorial is implemented in terms of the factorial and
      //        gamma functions using the relations: ..."""
      //
      //    So I guess we can just as well reformulate this with Gamma functions
      //    directly (just below) which can be evaluated as ratios directly.
      //
      //    [1] https://www.boost.org/doc/libs/1_76_0/libs/math/doc/html/math_toolkit/factorials/sf_double_factorial.html
   } else {
      // for full sphere integral evaluating to 1 (instead of 4π). That is our default mode.
      static FScalar const Prefactor = 1/(2*X_PI);
      // sort monomial powers such that ix >= iy >= iz. tgamma_delta_ratio(x,a) is meant
      // to be used with a << x (but at least a < x, I assume).
      //
      // See text at end of function for comments on what this does (and why)
      SortArguments(iz, iy, ix);
      FScalar
         f = (ix + 1)/FScalar(2),
         g = (iy + 1)/FScalar(2),
         h = (iz + 1)/FScalar(2);
      out = Prefactor *
         boost::math::tgamma_delta_ratio(f,g) *
         boost::math::tgamma(g) *
         boost::math::tgamma_delta_ratio(f+g,h) *
         boost::math::tgamma(h);
   }
   // Note: until now, "out" is supposed to be the inverse of a (possibly very large) integer.
   // See some ratios: Table[{i,j,Factorial2[i+j]/(Factorial2[i]*Factorial2[j])},{i,0,20,2},{j,i,20,2}] // Flatten[#,1]& // Column



#ifdef INCLUDE_4PI_SPHERE_INTEGRAL_PREFACTOR // note: OFF by default!
   static FScalar const FourPi = 4*X_PI;
   out *= FourPi;
#endif
   return out;


   // Note: For even i,j,k (otherwise we'd not be here), it is possible to
   // reformulate this into[1]:
   //
   //   (i-1)!! (j-1)!! (k-1)!! / (i + j + k + 1)!!
   //   = (1/(2 \[Pi])) (Gamma[1/2 + i/2] Gamma[1/2 + j/2] Gamma[1/2 + l/2])/(2 \[Pi] Gamma[ (i/2 + 1/2) + (j/2 + 1/2) + (k/2 + 1/2)])
   //
   // We had a very similar one for the multinomial coefficients, which we could
   // evaluate in a rather stable way using binomial coefficients (see
   // CxMonomialMath.h);
   // Unfortunately, there are no "binomial coefficients for half-integers" (at
   // least not implemented atm). But we can still do this with
   //
   //   tgamma_delta_ratio(x,a) = Γ(x)/Γ[x + a]
   //
   // (see https://www.boost.org/doc/libs/1_76_0/libs/math/doc/html/math_toolkit/sf_gamma/gamma_ratios.html).
   // So we can write the ijk part of this as (with f=i/2+1/2, g=j/2+1/2, h=k/2+1/2):
   //
   //    (Γ(f) Γ(g) Γ(h)) / Γ(f + g + h)
   //     = (Γ(f) Γ(g)) / Γ(f + g) * (Γ(f + g) Γ(h)) / Γ((f + g) + h)
   //     = Γ(g) * [Γ(f) / Γ(f + g)] * [Γ(f + g) / Γ((f + g) + h)] * Γ(h)
   //     = Γ_{ratio}(f,g) * Γ(g) * Γ_{ratio}(f+g,h) * Γ(h)
   //
   // So it doesn't even use more function calls (no idea how the implementation
   // of these relate, though).
}


// compute estimate for squared L2(S_2) 2-norm ⟨f, f⟩ = (1/4π) \int_{S_2} |f(x,y,z)|^2 \d\Omega
// of monomial functions f(x,y,z) = m_{ijk}(x,y,z) = x^i y^j z^k,
// with integration domain the 3D unit sphere S_2 (= {(x,y,z) | x^2 + y^2 + z^2 == 1}).
template<class ScalarT>
ScalarT EstimateMonomialUnitSphereNormSquaredT(unsigned ix, unsigned iy, unsigned iz)
{
   // some comments:
   // - Because |m_{ijk}|^2 = m_{2i,2j,2k}, when computing exactly, the result
   //   of this function is supposed to be identical to the result of
   //   CalcMonomialUnitSphereIntegralRaw(2*ix, 2*iy, 2*iz)
   //
   // - However, the norm (all exponents even) can be estimated more efficiently
   //   (and we need to if using this inside polynomial term screening)
   //
   // + First, the "are all the exponents even?" check will fall away.
   // + Second, note that (2n)!! == 2^n * n!:
   //
   //    In[3]:= Assuming[Element[n,Integers], FullSimplify[FunctionExpand[Factorial2[2*n]]]]
   //    Out[3]= 2^n Gamma[1+n]
   //
   //    In[4]:= Table[{Factorial[n], Gamma[1+n]}, {n,0,5}]
   //    Out[4]= {{1,1}, {1,1}, {2,2}, {6,6}, {24,24}, {120,120}}
   //
   //    In[6]:= Table[{Factorial2[2*n], 2^n*Factorial[n]}, {n,0,5}]
   //    Out[6]= {{1,1}, {2,2}, {8,8}, {48,48}, {384,384}, {3840,3840}}
   //
   // + The regular monomial-over-S_2 integral formula is
   //
   //    I_{S_2}[m_{ijk}] = ∫_{S_2} m_{ijk} dΩ = 4π (i-1)!! (j-1)!! (k-1)!! / (1+i+j+k)!! * θ(i even && j even && k even)
   //
   //   Inserting the formula above for doubled i,j,k (because it's the integral
   //   over m_{ijk}^2 for the norm), and omitting the 4π (our integral formula
   //   doens't have it either):
   //
   //    $Assumptions:=Element[{i,j,k},NonNegativeIntegers]
   //    NormSqRef[i_,j_,k_] := Factorial2[2*i-1] * Factorial2[2*j-1] * Factorial2[2*k-1] / Factorial2[1+2*i+2*j+2*k]
   //    FullSimplify[NormSqRef[i,j,l]]
   //
   //    -->    (Gamma[1/2 + i] Gamma[1/2 + j] Gamma[1/2 + l])/(2 \[Pi] Gamma[ 3/2 + i + j + l])
   //
   // + We could take this literally, but we can play a few additional tricks.
   //   Namely:
   //
   //    (*Consider the following cgk-special for a closed form approximate expression for Gamma[n] for n >= 1/2:
   //    - This expression is exact at n = 1/2, n = 3/2, and n -> Infinity, and has a maximum relative error of about 3.5e-4 on n\in[1/2,Infinity] (at n ≈ 0.633)
   //    - This version has terms arranged so that if evaluated (i+1/2) for integer i, the pow(n/E, ...) will have integer powers as-is.
   //      See GammaApproxNi for the evaluate-at-integer version (both variants are supplied with the actual n and both work with either argument---this is a continuous global fit)
   //    *)
   //    GammaApproxNh[n_] := With[
   //       { a0 = (3/2)*(E^2 + E^3 - 27)/(4 E^4 \[Pi]),
   //         a1 = ((81 - E^2 - 4 E^3)/(4 E^4 \[Pi])),
   //         a2 = (1/(2E \[Pi]))
   //       },
   //       (n/E)^(1/2 + n) / Sqrt[a0 + a1*n + a2*n^2]
   //    ];
   //    (* test: *)
   //    Block[{RelErr},
   //       RelErr[n_] := GammaApproxNh[n]/Gamma[n]-1;
   //       {
   //          Plot[With[{n=(-Log[1-x];x)}, RelErr[n+1/2]], {x,0,1}], (* relative errors on full n\in[1/2, infty] (via x\in[0,1] -> n\in[0,\infty] variable trafo)*)
   //          Plot[With[{n=x}, RelErr[n]], {x,0,1}], (* highlight of n\in[0,1] region, without variable trafo *)
   //          NMaximize[{Abs[RelErr[1/2-Log[1-x]]],x>=0&&x<=1}, x],
   //          NMaximize[{Abs[RelErr[n]],n>=1/2}, n]
   //       }
   //    ]
   //
   //    Note (not the variant used here): The mathematically equivalent full-integer variant is:
   //
   //    GammaApproxNi[n_] := With[
   //       {b0 = (3/2)*(E^2 + E^3 - 27)/(4 E^3 \[Pi]), (* b_i = E * a_i. That's all. *)
   //        b1 = (81 - E^2 - 4 E^3)/(4 E^3 \[Pi]),
   //        b2 = 1/(2\[Pi])
   //       },
   //       (n/E)^n*Sqrt[n/(b0 + b1*n + b2*n^2)]
   //    ]
   //
   //
   // + Based on the Gamma estimates, we can derive the following accurate approximation
   //   for the squared norm of a monomial x^i y^j z^k integrated over the unit sphere
   //   (max relative error < 1e-4, see below -- definitely good enough for screening):
   //   (...also, we drop the 4π--we normalize sphere integrals to 1):
   //
   //    NormSqRef[i_,j_,k_] := Factorial2[2*i-1] * Factorial2[2*j-1] * Factorial2[2*k-1] / Factorial2[1+2*i+2*j+2*k]
   //
   //    {va0,va1,va2} = {(3/2)*(E^2 + E^3 - 27)/(4 E^4 \[Pi]), ((81 - E^2 - 4 E^3)/(4 E^4 \[Pi])),  a2 = (1/(2E \[Pi]))}
   //    NormSqApprox5[i_,j_,k_] := (
   //       Block[{ti,tj,tk,t5,tijk,invtijk},
   //          ti=1+2*i;
   //          tj=1+2*j;
   //          tk=1+2*k;
   //          tijk=(3 + 2*i + 2*j + 2*k);
   //          invtijk=1/tijk;
   //          ((2 * tijk * (ti*invtijk)^(1+i) * (tj*invtijk)^(1+j) * (tk*invtijk)^(1+k) )
   //          *
   //          \[Sqrt](
   //                   (tijk^2 + 4 E \[Pi] (2*a0 + a1*tijk))
   //                   /
   //                   ((ti^2 + 4 E \[Pi] (2*a0 + a1*ti))*
   //                    (tj^2 + 4 E \[Pi] (2*a0 + a1*tj))*
   //                    (tk^2 + 4 E \[Pi] (2*a0 + a1*tk))
   //                   )
   //                )
   //          )
   //      ]
   //    );
   //
   //    Max[Flatten[Table[With[{ref=NormSqRef[i,j,k],ta=NormSqApprox4[i,j,k]/.{a0->va0,a1->va1,a2->va2}},Abs[N[(ta-ref)/ref]]],{i,0,20},{j,0,i},{k,0,j}],2]]
   //    0.000080144
   //
   // + This one needs two divisions, one sqrt, and three pow(x,n) with (integer
   //   exponents n) to get a norm estimate.
   using std::sqrt;
   using std::pow;

//    static ScalarT constexpr
//       e = boost::math::constants<ScalarT>::e(),
//       pi = boost::math::constants<ScalarT>::pi(),
//    ^--- beware! Not necessarily accurate enough. They only guarantee 34 decimal digits of precision (at least in older libraries)
   static ScalarT const
      e = ScalarT(X_EULER),
      pi = ScalarT(X_PI),
      f4epi = 4 * e * pi,
      ee = e*e;
   static ScalarT const
      // \{ a0 = (3/2)*(E^2 + E^3 - 27)/(4 E^4 \[Pi]),
      // \  a1 = ((81 - E^2 - 4 E^3)/(4 E^4 \[Pi])),
      // \  a2 = (1/(2E \[Pi]))
      // \}
      two_a0 = 3*(ee + ee*e - 27)/(4*ee*ee*pi),
      a1 = ((81 - ee - 4*ee*e)/(4*ee*ee*pi));
//       a2 = 1/(2*e*pi);  // <-- apparently not needed to compute the norm ratio
   unsigned const
      i = ix, j = iy, k = iz;
   // next up: NormSqApprox5 implementation
   ScalarT const
      ti = 1 + 2*i,
      tj = 1 + 2*j,
      tk = 1 + 2*k,
      tijk = 3 + 2*(i + j + k),
      inv_tijk = 1/ScalarT(tijk);
   ScalarT fNormSqEstimate =
      //  (2 * tijk * (ti*invtijk)^(1+i) * (tj*invtijk)^(1+j) * (tk*invtijk)^(1+k) )
      (2*tijk * pow(ti*inv_tijk, 1+i) * pow(tj*inv_tijk, 1+j) * pow(tk*inv_tijk, 1+k)) *
      //  \[Sqrt](
      sqrt(
            // (tijk^2 + 4 E \[Pi] (2*a0 + a1*tijk))
            (tijk*tijk + f4epi*(two_a0 + a1*tijk))
            / (
               // (ti^2 + 4 E \[Pi] (2*a0 + a1*ti))  + j,k
               (ti*ti + f4epi*(two_a0 + a1*ti)) *
               (tj*tj + f4epi*(two_a0 + a1*tj)) *
               (tk*tk + f4epi*(two_a0 + a1*tk))
            )
      );
#ifdef _DEBUG
   if (1) {
      // add a consistency check of the norm estimate vs the regular target integral routine
      using std::abs;
      ScalarT fNormSqRef = static_cast<ScalarT>(CalcMonomialUnitSphereIntegralRaw(2*i,2*j,2*k));
      assert(abs(fNormSqEstimate/fNormSqRef - 1) < ScalarT(1e-3));
      if (0) {
         // ...worked first try!
         io.Write("           M({:2},{:2},{:2})   fNrm2(est) = {:12.6g}   fNrm2(ref) = {:12.6g}   fRelErr = {:8.2e}", i,j,k, double(fNormSqEstimate), double(fNormSqRef), double(fNormSqEstimate/fNormSqRef - 1));
      }
   }
#endif // _DEBUG
   return fNormSqEstimate;
}


double EstimateMonomialNormSq_S2(FPolynomialN::FMonomial const &m) {
//    using std::sqrt;
//    return sqrt(EstimateMonomialUnitSphereNormSquaredT<double>(m[0], m[1], m[2]));
   using std::sqrt;
   return EstimateMonomialUnitSphereNormSquaredT<double>(m[0], m[1], m[2]);
}


FMonomial::FMonomial(size_t ix_, size_t iy_, size_t iz_)
   : ix(ix_), iy(iy_), iz(iz_), m_Degeneracy(1), m_Factor(1)
{
   _ResetNorm();
}

void FMonomial::_ResetNorm() {
   if (s_TargetFnNormalization == TARGETFN_NormalizeRawMonomials) {
      // compute raw monomial normalization factor (without possible
      // symmetrization/group averaging later on). Note: used mostly
      // implicitly as part of BaseFactor().
//       m_Factor = FScalar(1)/sqrt(CalcMonomialUnitSphereIntegralRaw(2*ix,2*iy,2*iz));
      m_Factor = FScalar(1)/CalcMonomialUnitSphereNorm(ix,iy,iz);
#ifdef INCLUDE_ABANDONED
//       m_Factor = FScalar(1)/CalcMonomialUnitSphereNormSq(ix,iy,iz); // <-- this is wrong --- just for checking if I can add any factor here
//       m_Factor = FScalar(1)/pow(CalcMonomialUnitSphereNormSq(ix,iy,iz), 2./3.); // <-- this is wrong --- just for checking if I can add any factor here
#endif // INCLUDE_ABANDONED
   } else {
      m_Factor = 1;
   }
   m_Degeneracy = 1;
}



// FScalar FMonomial::CalcUnitSphereIntegral()
// {
//    weight = 1.;
//    if ((ix % 2) == 1 || (iy % 2) == 1 || (iz % 2) == 1)
//       return 0.;
//    // note: 4π prefactor omitted (our grids are supposed to integrate to 1, not to 4π)
//    FScalar f = 1;
//    f = 4*X_PI;
//    return f * FScalar(DoubleFactR(FScalar(ix)-FScalar(1))*DoubleFactR(FScalar(iy)-FScalar(1))*DoubleFactR(FScalar(iz)-FScalar(1))) / FScalar(DoubleFactR(FScalar(ix+iy+iz+1)));
// }



FScalar FMonomial::CalcUnitSphereIntegralRaw() const
{
   return CalcMonomialUnitSphereIntegralRaw(ix, iy, iz);
//    if ((ix % 2) == 1 || (iy % 2) == 1 || (iz % 2) == 1)
//       return FScalar(0);
//    // note: 4π prefactor omitted (our grids are supposed to integrate to 1, not to 4π)
//    FScalar f = 1;
// //    f = 4*X_PI;
//    return f * FScalar(DoubleFactR(FScalar(ix)-FScalar(1))*DoubleFactR(FScalar(iy)-FScalar(1))*DoubleFactR(FScalar(iz)-FScalar(1))) / FScalar(DoubleFactR(FScalar(ix+iy+iz+1)));
}


// FScalar FMonomial::CalcUnitSphereIntegral()
// {
//    // absorb integral prefactor into monomial weight (currently contains degeneracy of monomial)
//    FScalar
//       Degeneracy = weight, // number symmetry equivalent monomials under the given point group
//       IntegralPrefactor = 1;
// //    Degeneracy = 1;
// //    IntegralPrefactor = 4*X_PI;
//    IntegralPrefactor = 1.;
//    IntegralPrefactor *= FScalar(DoubleFactR(FScalar(ix)-FScalar(1))*DoubleFactR(FScalar(iy)-FScalar(1))*DoubleFactR(FScalar(iz)-FScalar(1))) / FScalar(DoubleFactR(FScalar(ix+iy+iz+1)));
// //    std::cout << fmt::format("    monomial: {:4} {:4} {:4}  degeneracy: {:16.8f}\n", ix, iy, iz, double(weight));
//    weight = Degeneracy/IntegralPrefactor;
//
//    if ((ix % 2) == 1 || (iy % 2) == 1 || (iz % 2) == 1)
//       return FScalar(0);
//
//    if (1) {
//       return Degeneracy;
//    } else {
//       weight = Degeneracy;
//       return IntegralPrefactor * Degeneracy;
//    }
// }

FScalar FMonomial::BaseFactor() const
{
   return m_Factor * degn_or_one(m_Degeneracy);
}


FScalar FMonomial::CalcUnitSphereIntegral() const
{
//    // absorb integral prefactor into monomial weight (currently contains degeneracy of monomial)
//    FScalar
//       Degeneracy = weight; // number symmetry equivalent monomials under the given point group
// //       Degeneracy = 1;
//    weight = 1.;
//
//    if ((ix % 2) == 1 || (iy % 2) == 1 || (iz % 2) == 1)
//       return FScalar(0);
//
//    FScalar
//       IntegralPrefactor = FScalar(DoubleFactR(FScalar(ix)-FScalar(1))*DoubleFactR(FScalar(iy)-FScalar(1))*DoubleFactR(FScalar(iz)-FScalar(1))) / FScalar(DoubleFactR(FScalar(ix+iy+iz+1)));
//
//    return IntegralPrefactor * Degeneracy;
//    return CalcUnitSphereIntegralRaw();
   if (g_UseGroupSymmetrizedFn && !m_RgSymmetrized.empty()) {
      // return rhs value of the not the monomial itself, but its
      // group-symmetrized replacement
      //
      //   ⟨r⃗|P̂_G mⱼ⟩ = (1/|G|) ∑_ĝ mⱼ(R_g r⃗).
      //
      // This one goes with symmetrized overlap matrices.
      // Note: should there be a multiplication with |G| here? Not sure.
      FScalar s = 0;
      FPolynomialN::const_iterator
         itMj;
      for (itMj = this->m_RgSymmetrized.begin(); itMj != this->m_RgSymmetrized.end(); ++ itMj)
         s += itMj->second * CalcMonomialUnitSphereIntegralRaw(itMj->first[0], itMj->first[1], itMj->first[2]);
      return BaseFactor() * s;
   } else {
      return BaseFactor() * CalcUnitSphereIntegralRaw();
   }
}



FScalar FMonomial::EvalPoint(FPoint const &p) const
{
   return BaseFactor() * pow(p[0], ix) * pow(p[1], iy) * pow(p[2], iz);
}


FScalar FMonomial::EvalPoint(FPoint &dMdXyz, FPoint const &p) const
{
   // do it in the stupid way. fix later. maybe.
   if (ix == 0) dMdXyz[0] = 0; else dMdXyz[0] = BaseFactor() * ix * pow(p[0], ix-1) * pow(p[1], iy) * pow(p[2], iz);
   if (iy == 0) dMdXyz[1] = 0; else dMdXyz[1] = BaseFactor() * iy * pow(p[0], ix) * pow(p[1], iy-1) * pow(p[2], iz);
   if (iz == 0) dMdXyz[2] = 0; else dMdXyz[2] = BaseFactor() * iz * pow(p[0], ix) * pow(p[1], iy) * pow(p[2], iz-1);
   return BaseFactor() * pow(p[0], ix) * pow(p[1], iy) * pow(p[2], iz);
}


FScalar FMonomial::EvalPoint(FMonomialEvalInfo const &p) const
{
   return BaseFactor() * p.pPowX[ix] * p.pPowY[iy] * p.pPowZ[iz];
//    return p.pPowX[ix] * p.pPowY[iy] * p.pPowZ[iz];
   // ^- if we do the real overlap thing then we should not multiply those or the unit integrals by weight.
}


FScalar FMonomial::EvalPoint(FPoint &dMdXyz, FMonomialEvalInfo const &p) const
{
   if (ix == 0) dMdXyz[0] = 0; else dMdXyz[0] = BaseFactor() * ix * p.pPowX[ix-1] * p.pPowY[iy] * p.pPowZ[iz];
   if (iy == 0) dMdXyz[1] = 0; else dMdXyz[1] = BaseFactor() * iy * p.pPowX[ix] * p.pPowY[iy-1] * p.pPowZ[iz];
   if (iz == 0) dMdXyz[2] = 0; else dMdXyz[2] = BaseFactor() * iz * p.pPowX[ix] * p.pPowY[iy] * p.pPowZ[iz-1];
   return BaseFactor() * p.pPowX[ix] * p.pPowY[iy] * p.pPowZ[iz];
//    return p.pPowX[ix] * p.pPowY[iy] * p.pPowZ[iz];
   // ^- if we do the real overlap thing then we should not multiply those or the unit integrals by weight.
}


// make a list of all cartesian monomials in x,y,z of total lmin <= l <= lmax
FMonomialList MakeMonomialList(size_t lmin, size_t lmax)
{
   FMonomialList r;
   for (size_t l = lmin; l <= lmax; ++l) {
      for (size_t ix = 0; ix <= l; ++ ix)
         for (size_t iy = 0; iy <= l - ix; ++ iy) {
            size_t iz = l - ix - iy;
            r.push_back(FMonomial(ix,iy,iz));
         }
   }
   return r;
}

FMonomialList MakeSymmetryReducedMonomialList(size_t lmin, size_t lmax, FMonomialSymmetry const *pMonomialSymmetry, FPrintLevel iPrintLevel, bool PrescreenTargetFunctions)
{
   // we got a more or less closed list of what would end up in the target
   // monomials. Nevertheless, we collect them in a set first because this is
   // easy and it allows imposing a well-defined order for free.
   typedef FMonomialN_SmallIntegralsFirstOrder FMonomialOrder1;
//    typedef FMonomialN_ExportOrderPred FMonomialOrder1;
//    typedef std::less<FMonomialN> FMonomialOrder1;
   typedef std::set<FMonomialN, FMonomialOrder1>
      FMonomialN_Set;
   FMonomialN_Set
      RetainedMonomials;

   for (size_t l = lmin; l <= lmax; ++l) {
      for (size_t ix = 0; ix <= l; ++ ix) {
         for (size_t iy = 0; iy <= l - ix; ++ iy) {
            size_t iz = l - ix - iy;
            FMonomialN
               mijk(ix, iy, iz);
            if (PrescreenTargetFunctions && !pMonomialSymmetry->SelectAsCanonialReprQ(mijk))
               continue;
            RetainedMonomials.insert(mijk);
         }
      }
   }
   FMonomialList
      r;
   r.reserve(RetainedMonomials.size());
   for (FMonomialN const &m : RetainedMonomials)
      r.push_back(FMonomial(m[0], m[1], m[2]));

//    iPrintLevel = 2;
   PrintMonomialList(std::cout, r, iPrintLevel); // note: this doesn't do anything unless iPrintLevel >= 2;
   return r;
}


#ifdef INCLUDE_ABANDONED
// FMonomialList MakeSymmetryReducedMonomialList(size_t lmin, size_t lmax, FPointGroup const &PointGroup, FPrintLevel iPrintLevel, bool PrescreenTargetFunctions)
// {
//    // we got a more or less closed list of what would end up in the target
//    // monomials. Nevertheless, we collect them in a set first because this is
//    // easy and it allows imposing a well-defined order for free.
//    typedef FMonomialN_SmallIntegralsFirstOrder FMonomialOrder1;
// //    typedef FMonomialN_ExportOrderPred FMonomialOrder1;
// //    typedef std::less<FMonomialN> FMonomialOrder1;
//    typedef std::set<FMonomialN, FMonomialOrder1>
//       FMonomialN_Set;
//    FMonomialN_Set
//       RetainedMonomials;
//
//    for (size_t l = lmin; l <= lmax; ++l) {
//       for (size_t ix = 0; ix <= l; ++ ix) {
//          for (size_t iy = 0; iy <= l - ix; ++ iy) {
//             size_t iz = l - ix - iy;
//             if (PrescreenTargetFunctions && !PointGroup.CouldMonomialBeATargetFunction(ix, iy, iz))
//                continue;
//             RetainedMonomials.insert(FMonomialN(ix,iy,iz));
//          }
//       }
//    }
//    FMonomialList
//       r;
//    r.reserve(RetainedMonomials.size());
//    for (FMonomialN const &m : RetainedMonomials)
//       r.push_back(FMonomial(m[0], m[1], m[2]));
//
// //    iPrintLevel = 2;
//    PrintMonomialList(std::cout, r, iPrintLevel); // note: this doesn't do anything unless iPrintLevel >= 2;
//    return r;
// }
#endif // INCLUDE_ABANDONED


// returns total number of Cartesian monomials x^i y^j z^k
// which have for which i+j+k <= l
size_t nCartX(int l)
{
   assert(l >= -3);
   return size_t((l+1)*(l+2)*(l+3)/6);
}


// returns a unique index of a monomial x^i y^j z^k in such a way that
// indices are ordered by increasing i+j+k.
// (That means, in particular, that for all non-negative integer
// 3-tuples i,j,k, we have iCartX(i,j,k) < nCartX(i+j+k))
//
// If given, index iBaseAc will be subtracted from the actual
// monomial index (useful to exp
size_t iCartX(int i, int j, int k, size_t iBaseAc=0)
{
   assert(i + j + k >= 0);
   assert(i >= 0 && j >= 0 && k >= 0);
   size_t
      l = size_t(i + j + k);
   size_t
      // start with the total number of monomials (i,j,k) with i+j+k < l.
      idx = (l+0)*(l+1)*(l+2)/6;
   // there are (l+1)*(l+2)/2 monomials total with i+j+k = l (exactly l):
   //
   //    Sum[Sum[1,{j,0,l-i}],{i,0,l}] = (l+1)*(l+2)//2
   //
   // That applies because we have two degrees of freedom:
   // - i may go from 0 to l (inclusive)
   // - j may go from 0 to (l−i) (inclusive)
   // - but once i and j have assigned values, k is fixed by i+j+k = l.
   //
   // We assume we store the data in the following order, with 'j' data
   // up to index N-i in the fast dimension:
   //
   //   j=0  j=1  j=2 … j=l-2  j=l-1  j=l      // i=0 data
   //   j=0  j=1  j=2 … j=l-2  j=l-1           // i=1 data
   //   j=0  j=1  j=2 … j=l-2                  // i=2 data
   //   …
   //   j=0                                    // i=l data
   //
   // So, to get a starting index for the data of a given 'i', we compute
   // the data size for all 'j' entries *before* the current 'i', i.e., up to
   // amount of all data stored up to index i-1:
   //
   //    Sum[Sum[1,{j,0,l-i}],{i,0,im1}]/.im1->(i−1) = −(1/2)*i*(-3 + i - 2*l)
   //    = i*(2*l + 3 - i)/2
   //
   idx += size_t((i*(2*signed(l) + 3 - i))/2); // amount of 'j' data *before* current 'i'

   // to this we then just sum the index inside the 'j' data which is left
   idx += size_t(j);
   // UPDATE: tested; works -> see test_iac3.py

   assert(size_t(int(iBaseAc)) == iBaseAc);
   idx -= iBaseAc;
   assert(idx >= 0);
   return idx;
}





struct FGroupAveragingContext
{
   explicit FGroupAveragingContext(FPointGroup const &PointGroup);

   // (potentially) compute resident information used in the calculation of
   // large numbers of symmetrized monomials x^ix y^iy z^iz with ix+iy+iz between lmin and lmax
   virtual void Init(size_t lmin, size_t lmax) { m_lmin = lmin; m_lmax = lmax; };

   // Compute and return the symmetrized (=group-averaged) monomial
   //
   //    SymXyzn(i,j,k) = ∑_g (ĝ m_{ijk})(x,y,z)
   //                   = ∑_g m_{ijk}((R_g)⁻¹(x,y,z)).
   //
   // In here, ĝ is the group action, R_g is its representation in 3D space as a
   // (3,3)-shape matrix, and m_{ijk}(x,y,z) = x^i y^j z^k are the raw
   // Cartesian monomials.
   virtual FPolynomialN ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz) = 0;
   virtual std::string Desc() const { return "base"; };
protected:
   size_t
      m_lmin, m_lmax;
   FRotationMatrixList const
      &m_TrafoList;
   FPointGroup const
      &m_PointGroup;
   FPolynomialN::FCullCriteria
      // base culling criteria for intermediate monomials, assuming a norm scale of 1.
      // Object will be adjusted for concrete case if needed
      m_IntdCullCriteriaBase;

   // For the iRg'th group action g, return either of the values x', y', and z',
   // where x' = Rg * x   (with x = (1,0,0) for ixyz=0),
   //       y' = Rg * y   (with y = (0,1,0) for ixyz=1),
   //       z' = Rg * z   (with z = (0,0,1) for ixyz=2);
   // where x', y', z' is expressed as a polynomial in the original raw x,y,z.
   FPolynomialN const &GetRg_w(size_t ixyz, size_t iRg) { assert(ixyz < 3 && iRg < m_TrafoList.size()); size_t idx = 3*iRg + ixyz; assert(idx < m_Rg_w_as_poly.size()); return m_Rg_w_as_poly[idx]; }
   void InitRg_w();
private:
   std::vector<FPolynomialN>
      // cache for Rg-transformed base monomials x,y,z for all group elements;
      // it's an array of length `3*m_TrafoList.size()` (these would most likely
      // be required in whatever kind of method we use to compute the group
      // averages)
      m_Rg_w_as_poly;
};


// this variant uses no/little memory for storing intermediates during the
// computation of the transformed monomials x'^i y'^j z'^k. It assembles
// each such 3-index monomial independently of the others.
struct FGroupAveragingContext_Direct : public FGroupAveragingContext
{
//    using FGroupAveragingContext::FGroupAveragingContext;
   explicit FGroupAveragingContext_Direct(FPointGroup const &PointGroup, bool CacheOneIdxPowers);
   void Init(size_t lmin, size_t lmax); // override
   FPolynomialN ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz); // override
   virtual std::string Desc() const { return fmt::format("direct{{cache-1ix: {}}}", m_CacheOneIdxPowers? "on" : "off"); };
protected:
   bool
      // if set, precompute w'^k for lmin ≤ k ≤ max and w' ∈ {x',y',z'}
      m_CacheOneIdxPowers;
   std::vector<FPolynomialN>
      // if m_CacheOneIdxPowers: a ng*(1+lmax)*3 array containing w'^l entries.
      m_CachedOneIdxPowers;
   size_t iCachedOneIdxEntry(size_t ixyz, size_t iRg, size_t k) {
      assert(m_CacheOneIdxPowers);
      assert(ixyz < 3 && iRg < m_TrafoList.size() && 0 <= k && k <= m_lmax);
      return iRg + m_TrafoList.size() * (k + (m_lmax + 1)*ixyz);
   }
   FPolynomialN &CachedOneIdxPow(size_t ixyz, size_t iRg, size_t k) { return m_CachedOneIdxPowers[iCachedOneIdxEntry(ixyz, iRg, k)];  }
};


#ifdef INCLUDE_ABANDONED
// // this variant assembles the transformed monomials x'^i y'^j z'^k by
// // incrementally building from the lower monomials. This can require
// // significant amounts of intermediate memory.
// struct FGroupAveragingContext_CachedIncremental : public FGroupAveragingContext
// {
//    explicit FGroupAveragingContext_CachedIncremental(FPointGroup const &PointGroup);
//    void Init(size_t lmin, size_t lmax); // override
//    FPolynomialN ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz); // override
//    virtual std::string Desc() const { return "cached/incremental"; };
// protected:
//    size_t
//       m_iAc0; // = nCartX(m_lmin-1) ... base component of requested output cartesian monomials
//    std::vector<FPolynomialN>
//       m_avg_xyzn;
//    FPolynomialN &SYM_XYZN(int ix, int iy, int iz) { return m_avg_xyzn[iCartX(ix,iy,iz,m_iAc0)]; }
// };
//
//
// // this variant assembles the transformed monomials x'^i y'^j z'^k by
// // incrementally building from the lower monomials. This can require
// // significant amounts of intermediate memory.
// struct FGroupAveragingContext_CachedIncremental2 : public FGroupAveragingContext
// {
//    explicit FGroupAveragingContext_CachedIncremental2(FPointGroup const &PointGroup);
//    void Init(size_t lmin, size_t lmax); // override
//    FPolynomialN ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz); // override
//    virtual std::string Desc() const { return "cached/incremental_v2"; };
// protected:
//    size_t
//       m_iAc0; // = nCartX(m_lmin-1) ... base component of requested output cartesian monomials
//    std::vector<FPolynomialN>
//       m_Rg_xyzn;
// //    FPolynomialN &SYM_XYZN(int ix, int iy, int iz) { return m_avg_xyzn[iCartX(ix,iy,iz,m_iAc0)]; }
//    inline FPolynomialN &RG_XYZN(size_t iRg, size_t ix, size_t iy, size_t iz);
//
// };
//
//
// FGroupAveragingContext_CachedIncremental::FGroupAveragingContext_CachedIncremental(FPointGroup const &PointGroup)
//    : FGroupAveragingContext(PointGroup)
// {
// }
//
// void FGroupAveragingContext_CachedIncremental::Init(size_t lmin, size_t lmax)
// {
//    FGroupAveragingContext::Init(lmin, lmax);
//
// // #define SYM_XYZN(ix,iy,iz) avg_xyzn[iCartX(ix,iy,iz,iAc0)]
//    m_iAc0 = nCartX(int(m_lmin)-1);
//    m_avg_xyzn.resize(nCartX(m_lmax) - m_iAc0);
//
//    // initialize [(x^0 y^0 z^0)] to 1
//    if (lmin == 0)
//       SYM_XYZN(0,0,0) = FPolynomialN(1, FPolynomialN::FMonomial(0,0,0));
//
//    // for all group actions, evaluate transformed monomials [Rg,(x^i y^j
//    // z^k)] with i+j+k==l from transformed monomials with i+j+k == l - 1.
//    //
//    // note: can't parallelize over Rg... that would lead to (large) memory
//    // demand increasing linearly with number of threads.
//    size_t
//       ng = m_TrafoList.size();
//    for (int iRg = 0; iRg < int(ng); ++ iRg) {
//       std::vector<FPolynomialN>
//          // accumulator for Rg (x^i y^j z^j) being built incrementally for
//          // current group element Rg.
//          Rg_xyzn(nCartX(lmax));
// #define RG_XYZN(iRg,ix,iy,iz) Rg_xyzn[iCartX(ix,iy,iz,0)]
//       RG_XYZN(iRg,0,0,0) = FPolynomialN(1, FPolynomialN::FMonomial(0,0,0));
//
//       FPolynomialN const
//          &Rg_x_as_poly = GetRg_w(0, iRg),
//          &Rg_y_as_poly = GetRg_w(1, iRg),
//          &Rg_z_as_poly = GetRg_w(2, iRg);
//
//       for (size_t l = 1; l <= lmax; ++ l) {
//          #pragma omp parallel for schedule(dynamic)
//          for (int ixyz_ = 0; ixyz_ < int((l+1)*(l+2)/2); ++ ixyz_) {
//             size_t ixyz = size_t(ixyz_);
//             {
//                size_t
//                   ix,iy,iz;
//                // Regarding the indexing:
//                // - we run ix from 0 to l, and iy from 0 to l - ix. iz is fixed by ix and iy (as iz = l - ix - iz)
//                // - (ix,iy,iz) index offsets for any given ix:
//                //   FullSimplify[Sum[Sum[1,{iy,0,l-ix}],{ix,0,x-1}]] = -(1/2) x (-3 - 2 l + x)
//                // - given a ixyz, solve for the largest ix with offset < ixyz:
//                //   Assuming[ix>=0&&ixyz>=0&&l>=ix,FullSimplify[Solve[-(ix)*(-3-2*l+(ix))/2==ixyz,ix]]]
//                //   {{ix -> 3/2 + l - Sqrt[-2 ixyz + (3/2 + l)^2]}, {ix -> 3/2 + l + Sqrt[-2 ixyz + (3/2 + l)^2]}}
//                //   ...of this we want the smaller solution (larger one may yield ix with ix > l)
//                ix = size_t(1.5 + double(l) - std::sqrt(sqr(1.5 + l) - 2*ixyz));
//                iy = ixyz % (l+1-ix);
//                iz = l - ix - iy;
//
//                // so we now know the monomial powers (ix,iy,iz) of [x^ix * y^iy * z^iz].
//                // At this point we have all monomials with ix+iy+iz < l already constructed.
//                // Choose the reduction direction which yields the least amount of numerical
//                // work (in case you wonder: yes, it DOES make a difference).
//                size_t
//                   nOpsRedX = (ix > 0)? (Rg_x_as_poly.size() * RG_XYZN(iRg,ix-1,iy,iz).size()) : size_t(-1),
//                   nOpsRedY = (iy > 0)? (Rg_y_as_poly.size() * RG_XYZN(iRg,ix,iy-1,iz).size()) : size_t(-1),
//                   nOpsRedZ = (iz > 0)? (Rg_z_as_poly.size() * RG_XYZN(iRg,ix,iy,iz-1).size()) : size_t(-1);
//                if (nOpsRedX <= nOpsRedY && nOpsRedX <= nOpsRedZ) {
//                   cx_assert_rt(ix != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix-1,iy,iz) * Rg_x_as_poly;
//                } else if (nOpsRedY <= nOpsRedZ) {
//                   cx_assert_rt(iy != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix,iy-1,iz) * Rg_y_as_poly;
//                } else {
//                   cx_assert_rt(iz != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix,iy,iz-1) * Rg_z_as_poly;
//                }
//
//                // in case we reached one of the required outputs, increment the corresponding
//                // symmetrized quantity.
//                if (int(lmin) <= int(l) && int(l) <= int(lmax)) {
//                   // this is a component required for the output. add it to target location
//                   SYM_XYZN(ix,iy,iz) += RG_XYZN(iRg,ix,iy,iz);
//                   // ^- actual averaging (division by ng) done later
//                }
//             }
//          }
//       }
// #undef RG_XYZN
//    }
// }
//
//
//
//
// FPolynomialN FGroupAveragingContext_CachedIncremental::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
// {
// //    FPolynomialN
// //       AvgM;
// //    for (size_t iRg = 0; iRg != ng; ++ iRg)
// //       AvgM += RG_XYZN(iRg,ix,iy,iz);
// //    AvgM *= 1/FScalar(m_TrafoList.size());
// //    return AvgM;
//    FPolynomialN
//       AvgM = SYM_XYZN(ix,iy,iz);
//    AvgM *= 1/FScalar(m_TrafoList.size());
//    return AvgM;
// }
//
//
//
//
//
// FGroupAveragingContext_CachedIncremental2::FGroupAveragingContext_CachedIncremental2(FPointGroup const &PointGroup)
//    : FGroupAveragingContext(PointGroup)
// {
// }
//
// FPolynomialN &FGroupAveragingContext_CachedIncremental2::RG_XYZN(size_t iRg, size_t ix, size_t iy, size_t iz) {
// //    if (0) {
// //       size_t nl = m_lmax + 1;
// //       size_t ng = m_TrafoList.size();
// //       return m_Rg_xyzn[(iRg) + ng * ((ix) + nl * ((iy) + nl * (iz)))];
// //       // ^- todo: rebase with iCartX, and in post-processing delete entries for lower l?
// //    } else {
//       return m_Rg_xyzn[iRg + m_TrafoList.size() * iCartX(ix,iy,iz,0)];
// //    }
// }
//
//
// void FGroupAveragingContext_CachedIncremental2::Init(size_t lmin, size_t lmax)
// {
//    FGroupAveragingContext::Init(lmin, lmax);
//
//    size_t
//       ng = m_TrafoList.size(),
//       nl = lmax + 1;
// //    m_Rg_xyzn.resize(ng*nl*nl*nl); // <- not all of those needed, I think...
//    m_Rg_xyzn.resize(ng*nCartX(nl));
//
//    {
//       // initialize [Rg,(x^0 y^0 z^0)] to 1
//       for (size_t iRg = 0; iRg != ng; ++iRg)
//          RG_XYZN(iRg,0,0,0) = FPolynomialN(1, FPolynomialN::FMonomial(0,0,0));
//       // for all group actions, evaluate transformed monomials [Rg,(x^i y^j
//       // z^k)] with i+j+k==l from transformed monomials with i+j+k == l - 1.
//       for (size_t l = 1; l <= lmax; ++ l) {
//
//          #pragma omp parallel for schedule(dynamic)
//          for (int iRgXyz = 0; iRgXyz < int(ng) * int((l+1)*(l+2)/2); ++ iRgXyz) {
//             size_t iRg = size_t(iRgXyz) % ng;
//             size_t ixyz = size_t(iRgXyz) / ng;
//
//             FPolynomialN const
//                &Rg_x_as_poly = GetRg_w(0, iRg),
//                &Rg_y_as_poly = GetRg_w(1, iRg),
//                &Rg_z_as_poly = GetRg_w(2, iRg);
//
//             {
//                size_t
//                   ix,iy,iz;
//                // Regarding the indexing:
//                // - we run ix from 0 to l, and iy from 0 to l - ix. iz is fixed by ix and iy (as iz = l - ix - iz)
//                // - (ix,iy,iz) index offsets for any given ix:
//                //   FullSimplify[Sum[Sum[1,{iy,0,l-ix}],{ix,0,x-1}]] = -(1/2) x (-3 - 2 l + x)
//                // - given a ixyz, solve for the largest ix with offset < ixyz:
//                //   Assuming[ix>=0&&ixyz>=0&&l>=ix,FullSimplify[Solve[-(ix)*(-3-2*l+(ix))/2==ixyz,ix]]]
//                //   {{ix -> 3/2 + l - Sqrt[-2 ixyz + (3/2 + l)^2]}, {ix -> 3/2 + l + Sqrt[-2 ixyz + (3/2 + l)^2]}}
//                //   ...of this we want the smaller solution (larger one may yield ix with ix > l)
//                ix = size_t(1.5 + double(l) - std::sqrt(sqr(1.5 + l) - 2*ixyz));
//                iy = ixyz % (l+1-ix);
//                iz = l - ix - iy;
//
//                // so we now know the monomial powers (ix,iy,iz) of [x^ix * y^iy * z^iz].
//                // At this point we have all monomials with ix+iy+iz < l already constructed.
//                // Choose the reduction direction which yields the least amount of numerical
//                // work (in case you wonder: yes, it DOES make a difference).
//                size_t
//                   nOpsRedX = (ix > 0)? (Rg_x_as_poly.size() * RG_XYZN(iRg,ix-1,iy,iz).size()) : size_t(-1),
//                   nOpsRedY = (iy > 0)? (Rg_y_as_poly.size() * RG_XYZN(iRg,ix,iy-1,iz).size()) : size_t(-1),
//                   nOpsRedZ = (iz > 0)? (Rg_z_as_poly.size() * RG_XYZN(iRg,ix,iy,iz-1).size()) : size_t(-1);
//                if (nOpsRedX <= nOpsRedY && nOpsRedX <= nOpsRedZ) {
//                   cx_assert_rt(ix != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix-1,iy,iz) * Rg_x_as_poly;
//                } else if (nOpsRedY <= nOpsRedZ) {
//                   cx_assert_rt(iy != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix,iy-1,iz) * Rg_y_as_poly;
//                } else {
//                   cx_assert_rt(iz != 0);
//                   RG_XYZN(iRg,ix,iy,iz) = RG_XYZN(iRg,ix,iy,iz-1) * Rg_z_as_poly;
//                }
//             }
//          }
//
//          if (0) {
//             // clear data of intermediates for l < lmin which we do not need anymore
//             // now that we got shell 'l' fully assembled---to free up memory.
//             // (this helps quite a bit, but resulting memory requirement is still way too
//             // large, and this memory free up is actually *very* slow...)
//             if (l-1 < m_lmin) {
//                size_t l1 = l - 1;
//                for (size_t iz = 0; iz <= l1; ++ iz) {
//                   for (size_t iy = 0; iy <= l1-iz; ++ iy) {
//                      size_t ix = l1 - iz - iy;
//                      for (size_t iRg = 0; iRg != ng; ++iRg)
//                         RG_XYZN(iRg,ix,iy,iz).clear();
//                   }
//                }
//             }
//          }
//       }
//    }
// }
//
//
//
// FPolynomialN FGroupAveragingContext_CachedIncremental2::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
// {
//    FPolynomialN
//       AvgM;
//    for (size_t iRg = 0; iRg < m_TrafoList.size(); ++ iRg)
//       AvgM += RG_XYZN(iRg,ix,iy,iz);
//    AvgM *= 1/FScalar(m_TrafoList.size());
//    return AvgM;
// }
#endif // INCLUDE_ABANDONED



FGroupAveragingContext::FGroupAveragingContext(FPointGroup const &PointGroup)
   : m_lmin(0), m_lmax(0), m_TrafoList(PointGroup.GetOps()), m_PointGroup(PointGroup)
{
   InitRg_w(); // I guess we won't get around using these guys. So just make them always & now.
}


#ifdef INCLUDE_ABANDONED
// void FGroupAveragingContext::InitRg_w()
// {
//    assert(m_Rg_w_as_poly.empty());
//    m_Rg_w_as_poly.resize(3*m_TrafoList.size());
//
//    FPolynomialN::FCullCriteria
//       // coefficient threshold for neglect of intermediate monomials.
//       // should propagate from here to rest of FPolynomialN calculations
//       // automatically. Note: these propagate automatically to all subsequent
//       // calculations
// //       Cull_PolyIntd(g_ThrAlmostZero*g_ThrAlmostZero),
// //       Cull_PolyIntd(FLOAT_LITERAL(1e-10)*g_ThrAlmostZero), // <- FIXME: Needed?
// //       Cull_PolyIntd(FLOAT_LITERAL(1e-2)*g_ThrAlmostZero),
// //       Cull_PolyIntd(g_ThrVerySmall),
// //       Cull_PolyIntd(sqr(g_ThrVerySmall)),
// //       Cull_PolyIntd(pow(g_ThrAlmostZero, 4/FScalar(3))),
// //       Cull_PolyIntd(g_ThrVerySmall),
// //       Cull_PolyIntd(g_ThrAlmostZero),
//       Cull_PolyIntd(sqr(g_ThrAlmostZero)),
//       // note @ Cull_FirstTermOnly: that only applies to the linear terms
//       // (transformed Rg m_{1,0,0}, Rg m_{0,1,0} and Rg m_{0,0,1} for all
//       // Rg). Anything != 0 should be fine here---our operator matrices are
//       // exact to almost full floating point precision.
//       Cull_FirstTermOnly(g_ThrAlmostZero);
// //       Cull_FirstTermOnly(sqr(g_ThrAlmostZero)); // TODO: check if necessarily. something is wrong. this should ensure at least the linear terms stay
//    if (0) {
//       // add function to culling criteria which provide an estimate of the
//       // monomial norm as a function in function space. This can help making
//       // sense of whether polynomial terms with small coefficients
//       // are actually significant or just numerical noise.
//       FPolynomialN::FEstimateSquaredNormFn
//          pEstNormSq(EstimateMonomialNormSq_S2);
//       Cull_PolyIntd.SetNormEstimateFn(pEstNormSq);
// //       Cull_FirstTermOnly.SetNormEstimateFn(pEstNormSq);
//    }
//
//    // FIXME: an issue. The thresholds need to be at least low enough that
//    // the smallest non-zero matrix element of the R_g group actions does
//    // not simply get cancelled if raised to the l'th power for our target
//    // degree. E.g., for I that's a 0.309017:
//    //  g:[ 52]  R_g ∘ m_{1,0,0} = -0.809017 x + 0.5 y + 0.309017 z
//    // This means that, for example, if we make a l=80 grid, in which
//    // Rg m_{80,0,0} is computed by raising this Rg m_{1,0,0} expansion
//    // to the 80'th power, then this alone will already yield a coefficient of
//    // 0.309017**80 == 1.5797550753732336e-41 ....
//
//    m_IntdCullCriteriaBase = Cull_PolyIntd;
//    if (0)
//       // disable special treatment for linear terms (these don't have high
//       // monomial powers---very small coefficients therefore likely mean that
//       // the term is indeed averaged to zero)
//       Cull_FirstTermOnly = Cull_PolyIntd;
//    FPolynomialN
//       x_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(1,0,0), Cull_PolyIntd),
//       y_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,1,0), Cull_PolyIntd),
//       z_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,0,1), Cull_PolyIntd);
//
//
//    // for all group actions, evaluate transformed monomials w' = Rg * w,
//    // where w is either x, y, or z.
//    for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
//       FRotationMatrix const
//          &Rg = m_TrafoList[iRg];
//       // compute image of x' as image of (1,0,0) under group action
//       FPoint
//          Rg_x = Rg * FPoint(1,0,0).matrix(),
//          Rg_y = Rg * FPoint(0,1,0).matrix(),
//          Rg_z = Rg * FPoint(0,0,1).matrix();
//       FPolynomialN
//          &Rg_x_as_poly = m_Rg_w_as_poly[3*iRg + 0],
//          &Rg_y_as_poly = m_Rg_w_as_poly[3*iRg + 1],
//          &Rg_z_as_poly = m_Rg_w_as_poly[3*iRg + 2];
//       Rg_x_as_poly = (Rg_x[0] * x_as_poly + Rg_x[1] * y_as_poly + Rg_x[2] * z_as_poly).Purged(Cull_FirstTermOnly);
//       Rg_y_as_poly = (Rg_y[0] * x_as_poly + Rg_y[1] * y_as_poly + Rg_y[2] * z_as_poly).Purged(Cull_FirstTermOnly);
//       Rg_z_as_poly = (Rg_z[0] * x_as_poly + Rg_z[1] * y_as_poly + Rg_z[2] * z_as_poly).Purged(Cull_FirstTermOnly);
//       // reset screening criteria as basis for further computations
//       Rg_x_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
//       Rg_y_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
//       Rg_z_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
//       if (0) {
// #ifdef INCLUDE_ABANDONED
// //          io.Write("     [{0:4}] Rg ∘ m[1,0,0] = {1}", iRg, Rg_x_as_poly);
// //          io.Write( "            Rg ∘ m[0,1,0] = {1}", iRg, Rg_y_as_poly);
// //          io.Write( "            Rg ∘ m[0,0,1] = {1}", iRg, Rg_z_as_poly);
// //          io.Write("     [{0:4}] (Rg ∘ m_{{1,0,0}})(x,y,z) = {1}", iRg, Rg_x_as_poly);
// //          io.Write( "            (Rg ∘ m_{{0,1,0}})(x,y,z) = {1}", iRg, Rg_y_as_poly);
// //          io.Write( "            (Rg ∘ m_{{0,0,1}})(x,y,z) = {1}", iRg, Rg_z_as_poly);
// //          io.Write("      g:[{0:3}] (R_g ∘ m_{{1,0,0}}) = {1}", iRg, Rg_x_as_poly);
// //          io.Write(  "              (R_g ∘ m_{{0,1,0}}) = {1}", iRg, Rg_y_as_poly);
// //          io.Write(  "              (R_g ∘ m_{{0,0,1}}) = {1}", iRg, Rg_z_as_poly);
// #endif // INCLUDE_ABANDONED
//          io.Write("      g:[{0:3}]  R_g ∘ m_{{1,0,0}} = {1}", iRg, Rg_x_as_poly);
//          io.Write(  "               R_g ∘ m_{{0,1,0}} = {1}", iRg, Rg_y_as_poly);
//          io.Write(  "               R_g ∘ m_{{0,0,1}} = {1}", iRg, Rg_z_as_poly);
//       }
//    }
// }
#endif // INCLUDE_ABANDONED

void FGroupAveragingContext::InitRg_w()
{
   assert(m_Rg_w_as_poly.empty());
   m_Rg_w_as_poly.resize(3*m_TrafoList.size());

   FPolynomialN::FCullCriteria
      // coefficient threshold for neglect of intermediate monomials.
      // should propagate from here to rest of FPolynomialN calculations
      // automatically.
//       Cull_PolyIntd(g_ThrAlmostZero*g_ThrAlmostZero),
//       Cull_PolyIntd(FLOAT_LITERAL(1e-10)*g_ThrAlmostZero), // <- FIXME: Needed?
//       Cull_PolyIntd(FLOAT_LITERAL(1e-2)*g_ThrAlmostZero),
//       Cull_PolyIntd(g_ThrVerySmall),
      Cull_PolyIntd(sqr(g_ThrAlmostZero)),
      // note @ `Cull_FirstTermOnly`: that only applies to the linear terms
      // (transformed Rg m_{1,0,0}, Rg m_{0,1,0} and Rg m_{0,0,1} for all
      // Rg). Anything ≠ 0 should be fine here — our operator matrices are
      // exact to almost full floating point precision.
      Cull_FirstTermOnly(g_ThrAlmostZero);
   if (0) {
      // add function to culling criteria which provide an estimate of the
      // monomial norm as a function in function space. This can help making
      // sense of whether polynomial terms with small coefficients
      // are actually significant or just numerical noise.
      FPolynomialN::FEstimateSquaredNormFn
         pEstNormSq(EstimateMonomialNormSq_S2);
      Cull_PolyIntd.SetNormEstimateFn(pEstNormSq);
   }

   m_IntdCullCriteriaBase = Cull_PolyIntd;
   FPolynomialN
      x_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(1,0,0), Cull_PolyIntd),
      y_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,1,0), Cull_PolyIntd),
      z_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,0,1), Cull_PolyIntd);


   // for all group actions, evaluate transformed monomials w' = Rg * w,
   // where w is either x, y, or z.
   for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
      FRotationMatrix const
         &Rg = m_TrafoList[iRg];
      // compute image of x' as image of (1,0,0) under group action
      FPoint
         Rg_x = Rg * FPoint(1,0,0).matrix(),
         Rg_y = Rg * FPoint(0,1,0).matrix(),
         Rg_z = Rg * FPoint(0,0,1).matrix();
      FPolynomialN
         &Rg_x_as_poly = m_Rg_w_as_poly[3*iRg + 0],
         &Rg_y_as_poly = m_Rg_w_as_poly[3*iRg + 1],
         &Rg_z_as_poly = m_Rg_w_as_poly[3*iRg + 2];
      Rg_x_as_poly = (Rg_x[0] * x_as_poly + Rg_x[1] * y_as_poly + Rg_x[2] * z_as_poly).Purged(Cull_FirstTermOnly);
      Rg_y_as_poly = (Rg_y[0] * x_as_poly + Rg_y[1] * y_as_poly + Rg_y[2] * z_as_poly).Purged(Cull_FirstTermOnly);
      Rg_z_as_poly = (Rg_z[0] * x_as_poly + Rg_z[1] * y_as_poly + Rg_z[2] * z_as_poly).Purged(Cull_FirstTermOnly);
      // reset screening criteria as basis for further computations
      Rg_x_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
      Rg_y_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
      Rg_z_as_poly.SetCullCriteria(m_IntdCullCriteriaBase);
      if (0) {
         io.Write("      g:[{0:3}]  R_g ∘ m_{{1,0,0}} = {1}", iRg, Rg_x_as_poly);
         io.Write(  "               R_g ∘ m_{{0,1,0}} = {1}", iRg, Rg_y_as_poly);
         io.Write(  "               R_g ∘ m_{{0,0,1}} = {1}", iRg, Rg_z_as_poly);
      }
   }
}


FGroupAveragingContext_Direct::FGroupAveragingContext_Direct(FPointGroup const &PointGroup, bool CacheOneIdxPowers)
   : FGroupAveragingContext(PointGroup), m_CacheOneIdxPowers(CacheOneIdxPowers)
{
}

void FGroupAveragingContext_Direct::Init(size_t lmin, size_t lmax)
{
   FGroupAveragingContext::Init(lmin, lmax);
   if (!m_CacheOneIdxPowers) {
      // do nothing. This one computes all group averages independently of each other.
   } else {
      m_CachedOneIdxPowers.clear();
      m_CachedOneIdxPowers.resize(m_TrafoList.size() * (m_lmax + 1) * 3);
      for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
         for (size_t iw = 0; iw < 3; ++ iw) {
            // start with w^0 = 1
            CachedOneIdxPow(iw, iRg, 0) = FPolynomialN(1, FPolynomialN::FMonomial(0,0,0));
            // make other w^k by multiplying w^{k-1} to w.
            for (size_t k = 1; k <= m_lmax; ++k ) {
               CachedOneIdxPow(iw, iRg, k) = CachedOneIdxPow(iw, iRg, k-1) * GetRg_w(iw,iRg);
//                CachedOneIdxPow(iw, iRg, k) = (CachedOneIdxPow(iw, iRg, k-1) * GetRg_w(iw,iRg)).Purged(FPolynomialN::FCullCriteria(1e-2));
            }
         }
      }
   }
}



#ifdef INCLUDE_ABANDONED
// // make a copy of `PolyIn`, replace the copie's cull criteria by Cull, then return this copy.
// FPolynomialN CullReplaced(FPolynomialN const &PolyIn, FPolynomialN::FCullCriteria const &NewCull) {
// //    return PolyIn;
//    FPolynomialN
//       PolyOut = PolyIn;
//    PolyOut.SetCullCriteria(NewCull, true); // true: exec purge now.
//    return PolyOut;
// }
//
//
// FPolynomialN FGroupAveragingContext_Direct::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
// {
//    FPolynomialN
//       // accumulator for sum over Rg-transformed monomials. Will
//       // be divided by ng at end to form the average.
//       // Note: we start with empty monomial (no terms), which is a zero monomial.
//       AvgM;
//    FPolynomialN::FCullCriteria
//       // set up a specific culling criterion for the current monomial, taking
//       // into account that the monomial is not normalized and may therefore
//       // aquire an external norm scaling factor on output which should be
//       // considered when screening contributions. It can be large!
//       CullForThis(m_IntdCullCriteriaBase);
//    FScalar
//       fSourceObjectNorm = sqrt(CalcMonomialUnitSphereIntegralRaw(2*ix, 2*iy, 2*iz));
//    if (0) {
//       assert(CullForThis.m_pEstimateSquaredNormFn);
//    //    FScalar
//    //       fSourceObjectNorm = sqrt(CullForThis.m_pEstimateSquaredNormFn(FMonomialN(ix, iy, iz)));
//       CullForThis.m_ThrCoeffNeglectAbs *= fSourceObjectNorm;
//       CullForThis.SetNormEstimateFn(EstimateMonomialNormSq_S2);
//       if (1) {
//    //       CullForThis.m_pEstimateSquaredNormFn = 0;
//    //       CullForThis.m_ThrCoeffNeglectAbs *= 1e-10; // to check if lowering it *here* still does anything (...or if it is already too late)
//          // UPDATE: at this point it is too late already!
//          //  make && time aigg_64_d --preprocess '{group:I; order:20 to 40; skip:no; print:2}'
//          // fails at l=22, despite all the thresholds being so tiny. Could be prep of the Rg things,
//          // or the cached coordinate powers.
//          //
//          // UPDATE: it works for l=22 if I turn off m_CachedOneIdxPowers now. In this case
//          // the power will be recomputed each time... but at that point the new screening
//          // settings have already been set up, unlike before.
//          m_CacheOneIdxPowers = false; // FIXME: remove this.
//
//          //
//          // TODO: hm.. maybe it's just about the number of terms, too? too high powers of 1/sqrt(3) and
//          // you'd have lots of terms which are all significant, but all equally small and below thr...
//       }
//    }
//    // IDEA: try running two group avg computations one after another, with different thresholds;
//    // one which works and one which does not. Try to locate source of problem in this way.
//    AvgM.SetCullCriteria(CullForThis, false);
//    if (1) {
// //       io.Write("           GAvg/D   M({:2},{:2},{:2})   fNrm2(M) = {:12.6g}  //  fThrC(base) = {:10.4e} --> fThrC(P_G M) = {:10.4e}", ix,iy,iz, double(fSourceObjectNorm), double(m_IntdCullCriteriaBase.m_ThrCoeffNeglectAbs), double(CullForThis.m_ThrCoeffNeglectAbs));
//    }
//
// #define ADJUST_CULL(x) CullReplaced((x), CullForThis)
//    // for all group actions, evaluate transformed monomials [Rg,(x^i y^j
//    // z^k)] with i+j+k==l from transformed monomials with i+j+k == l - 1.
//    for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
//       // multiply the transformed monomials together and increment accumulator
//       // for group average
//       if (!m_CacheOneIdxPowers) {
//          AvgM += Multiply(pow(ADJUST_CULL(GetRg_w(0,iRg)), ix), pow(ADJUST_CULL(GetRg_w(1,iRg)), iy), pow(ADJUST_CULL(GetRg_w(2,iRg)), iz));
//       } else {
//          size_t l = ix + iy + iz; (void)l;
//          assert(m_lmin <= l && l <= m_lmax);
//          AvgM += Multiply(ADJUST_CULL(CachedOneIdxPow(0,iRg,ix)), ADJUST_CULL(CachedOneIdxPow(1,iRg,iy)), ADJUST_CULL(CachedOneIdxPow(2,iRg,iz)));
//       }
//    }
// #undef ADJUST_CULL
//    // divide by number of group elements.
//    FScalar f = FScalar(1)/FScalar(m_TrafoList.size());
//    if (0) {
// //       f *= FScalar(1)/sqrt(CalcMonomialUnitSphereIntegralRaw(2*ix, 2*iy, 2*iz));
// //       f *= 2;
//    }
//    AvgM *= f;
//    if (0) {
//       std::string PfType = "";
//       FScalar prefac = 1;
//       // Hm... for ico group (I), m_{4,0,0} and m_{2,2,0} definitely average to polynomials
//       // which are colinear. However, even if I do normalize the input monomials, they
//       // are still not quite identical. I guess this means they both have slightly different components
//       // in subspaces which P_G annihilates...
//       //
//       // Note: without scaling, the coeffs come out like this:
//       //
//       //       1/‖m‖ * P_G[m_{0,2,2}] = 0.0666667 x^4 + 0.133333 x^2y^2 + 0.0666667 y^4 + 0.133333 x^2z^2 + 0.133333 y^2z^2 + 0.0666667 z^4
//       //
//       //       1/‖m‖ * P_G[m_{0,0,4}] = 0.2 x^4 + 0.4 x^2y^2 + 0.2 y^4 + 0.4 x^2z^2 + 0.4 y^2z^2 + 0.2 z^4
//       //
//       // So there is a factor of 3 in difference for the raw monomials.
//       // However, this doesn't fit to the ratio between the original monomial norms:
//       //
//       // >>> from scipy.special import factorial2
//       // >>> fac2 = factorial2
//       // >>> def n(i,j,k): return fac2(2*i-1)*fac2(2*j-1)*fac2(2*k-1)/fac2(2*(i+j+k)+1)
//       // ...
//       // >>> n(4,0,0)
//       // 0.11111111111111113
//       // >>> n(4,0,0)/n(2,2,0)
//       // 11.666666666666664
//       // >>> (n(4,0,0)/n(2,2,0))**.5
//       // 3.415650255319866
//       //
//       // It's close, but not identical.
//       // I guess I can normalize the target polynomial instead. But would need to
//       // be careful with the "I can get S_{ij} = ⟨P_G[m_i]| P_G[m_j]⟩ = ⟨m_i| P_G[m_j]⟩
//       // because P_G is a projector and thereby an involution"-business...
//
//       if (0) {
//          if (0) {
//             PfType = "1/‖m‖ * ";
//             prefac = 1/fSourceObjectNorm;
//    //          prefac = 1;
//          } else {
//             // emit normalizatin for the actual averaged polynomial. I think in theory
//             // this should indiscrimenantly yield identical functions for the group averages.
//             //
//             // ...the issue is just: the group averages in I are typically not at all equal.
//             // Rather, there are some centralizer-functions, and the different monomials
//             // with non-zero averages have different degrees of overlap with them.
//             PfType = "1/‖P_G[m]‖ * ";
//             FScalar fTargetObjectNorm = sqrt(ComputeOverlapRaw(AvgM, AvgM));
//             prefac = 1/fTargetObjectNorm;
//          }
//       }
//       io.Write("\n      {}P_G[m_{{{},{},{}}}] = {}", PfType, ix, iy, iz, prefac * AvgM);
//    }
//    return AvgM;
// }

// m224 = 0.00416667 x^8 + 0.0295424 x^6 y^2 + 0.00416667 x^4 y^4 - 0.0170424 x^2 y^6 + 0.00416667 y^8 - 0.0170424 x^6 z^2 + 0.1125 x^4 y^2 z^2 + 0.1125 x^2 y^4 z^2 + 0.0295424 y^6 z^2 + 0.00416667 x^4 z^4 + 0.1125 x^2 y^2 z^4 + 0.00416667 y^4 z^4 + 0.0295424 x^2 z^6 - 0.0170424 y^2 z^6 + 0.00416667 z^8
//
// m044 = 0.00833333 x^8 + 0.0178825 x^6 y^2 + 0.075 x^4 y^4 + 0.0737842 x^2 y^6 + 0.00833333 y^8 + 0.0737842 x^6 z^2 + 0.025 x^4 y^2 z^2 + 0.025 x^2 y^4 z^2 + 0.0178825 y^6 z^2 + 0.075 x^4 z^4 + 0.025 x^2 y^2 z^4 + 0.075 y^4 z^4 + 0.0178825 x^2 z^6 + 0.0737842 y^2 z^6 + 0.00833333 z^8
//
// m025 = 0.00992486 x^8 - 0.0375 x^6 y^2 + 0.18446 x^4 y^4 + 0.24181 x^2 y^6 + 0.00992486 y^8 + 0.24181 x^6 z^2 - 0.255636 x^4 y^2 z^2 - 0.255636 x^2 y^4 z^2 - 0.0375 y^6 z^2 + 0.18446 x^4 z^4 - 0.255636 x^2 y^2 z^4 + 0.18446 y^4 z^4 - 0.0375 x^2 z^6 + 0.24181 y^2 z^6 + 0.00992486 z^8
//
// m062 = 0.0192418 x^8 + 0.12069 x^6 y^2 + 0.0447062 x^4 y^4 - 0.0375 x^2 y^6 + 0.0192418 y^8 - 0.0375 x^6 z^2 + 0.443136 x^4 y^2 z^2 + 0.443136 x^2 y^4 z^2 + 0.12069 y^6 z^2 + 0.0447062 x^4 z^4 + 0.443136 x^2 y^2 z^4 + 0.0447062 y^4 z^4 + 0.12069 x^2 z^6 - 0.0375 y^2 z^6 + 0.0192418 z^8
//
// m008 = 0.116667 x^8 + 0.538771 x^6 y^2 + 0.583333 x^4 y^4 + 0.277896 x^2 y^6 + 0.116667 y^8 + 0.277896 x^6 z^2 + 1.75 x^4 y^2 z^2 + 1.75 x^2 y^4 z^2 + 0.538771 y^6 z^2 + 0.583333 x^4 z^4 + 1.75 x^2 y^2 z^4 + 0.583333 y^4 z^4 + 0.538771 x^2 z^6 + 0.277896 y^2 z^6 + 0.116667 z^8

// FPolynomialN FGroupAveragingContext_Direct::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
// {
//    FPolynomialN
//       x_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(1,0,0)),
//       y_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,1,0)),
//       z_as_poly = FPolynomialN(1, FPolynomialN::FMonomial(0,0,1));
//
//    FPolynomialN
//       // accumulator for sum over Rg-transformed monomials. Will
//       // be divided by ng at end to form the average.
//       // Note: we start with empty monomial (no terms), which is a zero monomial.
//       AvgM;
//
//    // for all group actions, evaluate transformed monomials [Rg,(x^i y^j
//    // z^k)] with i+j+k==l from transformed monomials with i+j+k == l - 1.
//    for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
//       // For all group actions g, compute (x'^i y'^j z'^k) for i+j+k <= l, where x' = Rg * x
//       // (same with y' = Rg * y, z' = Rg * z);
//       FRotationMatrix const
//          &Rg = m_TrafoList[iRg];
//       // compute image of x' as image of (1,0,0) under group action
//       FPoint
//          Rg_x = Rg * FPoint(1,0,0).matrix(),
//          Rg_y = Rg * FPoint(0,1,0).matrix(),
//          Rg_z = Rg * FPoint(0,0,1).matrix();
//       FPolynomialN
//          Rg_x_as_poly = (Rg_x[0] * x_as_poly + Rg_x[1] * y_as_poly + Rg_x[2] * z_as_poly).Purged(),
//          Rg_y_as_poly = (Rg_y[0] * x_as_poly + Rg_y[1] * y_as_poly + Rg_y[2] * z_as_poly).Purged(),
//          Rg_z_as_poly = (Rg_z[0] * x_as_poly + Rg_z[1] * y_as_poly + Rg_z[2] * z_as_poly).Purged();
//
//       // multiply the transformed monomials together and increment accumulator
//       // for group average
//       AvgM += Multiply(pow(Rg_x_as_poly, ix), pow(Rg_y_as_poly, iy), pow(Rg_z_as_poly, iz));
// //          AvgM += MultiplyPow(Rg_x_as_poly, ix, Rg_y_as_poly, iy, Rg_z_as_poly, iz);
//
//    }
//    // divide by number of group elements.
//    AvgM *= (1/FScalar(m_TrafoList.size()));
//    return AvgM;
// }



// // make a copy of `PolyIn`, replace the copie's cull criteria by Cull, then return this copy.
// FPolynomialN CullReplaced(FPolynomialN const &PolyIn, FPolynomialN::FCullCriteria const &NewCull) {
// //    return PolyIn;
//    FPolynomialN
//       PolyOut = PolyIn;
//    PolyOut.SetCullCriteria(NewCull, true); // true: exec purge now.
//    return PolyOut;
// }


// FPolynomialN FGroupAveragingContext_Direct::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
// {
//    FPolynomialN
//       // accumulator for sum over Rg-transformed monomials. Will
//       // be divided by ng at end to form the average.
//       // Note: we start with empty monomial (no terms), which is a zero monomial.
//       AvgM(m_IntdCullCriteriaBase);
//
// // #define ADJUST_CULL(x) CullReplaced((x), CullForThis)
// #define ADJUST_CULL(x) (x)
//    // for all group actions, evaluate transformed monomials [Rg,(x^i y^j
//    // z^k)] with i+j+k==l from transformed monomials with i+j+k == l - 1.
//    for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
//       // multiply the transformed monomials together and increment accumulator
//       // for group average
//       if (!m_CacheOneIdxPowers) {
//          AvgM += Multiply(pow(ADJUST_CULL(GetRg_w(0,iRg)), ix), pow(ADJUST_CULL(GetRg_w(1,iRg)), iy), pow(ADJUST_CULL(GetRg_w(2,iRg)), iz));
//       } else {
//          size_t l = ix + iy + iz; (void)l;
//          assert(m_lmin <= l && l <= m_lmax);
//          AvgM += Multiply(ADJUST_CULL(CachedOneIdxPow(0,iRg,ix)), ADJUST_CULL(CachedOneIdxPow(1,iRg,iy)), ADJUST_CULL(CachedOneIdxPow(2,iRg,iz)));
//       }
//    }
// #undef ADJUST_CULL
//    // divide by number of group elements.
//    AvgM *= FScalar(1)/FScalar(m_TrafoList.size());
//    return AvgM;
// }
#endif // INCLUDE_ABANDONED

FPolynomialN FGroupAveragingContext_Direct::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
{
   FPolynomialN
      // accumulator for sum over Rg-transformed monomials. Will
      // be divided by ng at end to form the average.
      // Note: we start with empty monomial (no terms), which is a zero monomial.
      AvgM(m_IntdCullCriteriaBase);

   // for all group actions, evaluate transformed monomials [Rg, (x^i y^j
   // z^k)] with i+j+k==l from transformed monomials with i+j+k == l − 1.
   for (int iRg = 0; iRg < int(m_TrafoList.size()); ++ iRg) {
      // multiply the transformed monomials together and increment accumulator
      // for group average
      if (!m_CacheOneIdxPowers) {
         AvgM += Multiply(pow(GetRg_w(0,iRg), ix), pow(GetRg_w(1,iRg), iy), pow(GetRg_w(2,iRg), iz));
      } else {
         size_t l = ix + iy + iz; (void)l;
         assert(m_lmin <= l && l <= m_lmax);
         AvgM += Multiply(CachedOneIdxPow(0,iRg,ix), CachedOneIdxPow(1,iRg,iy), CachedOneIdxPow(2,iRg,iz));
      }
   }
   // divide by number of group elements.
   AvgM *= FScalar(1)/FScalar(m_TrafoList.size());
   return AvgM;
}



// this variant bypasses any actual calculations using `FPolynomialN` objects.
// Instead, it invokes a total expansion of
//
//     ∑_g m_{ix,iy,iz}([R_g r⃗]_x, [R_g r⃗]_y, [R_g r⃗]_z)
//
// for all (!) terms at once using high-order trinomial coefficients.
// It may cache some expansion coefficients, but generally treats all monomials
// independently, and memory demand should be quite reasonable.
struct FGroupAveragingContext_TrinomialExpand : public FGroupAveragingContext
{
//    using FGroupAveragingContext::FGroupAveragingContext;
   explicit FGroupAveragingContext_TrinomialExpand(FPointGroup const &PointGroup);
   void Init(size_t lmin, size_t lmax); // override
   FPolynomialN ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz); // override
   virtual std::string Desc() const { return "trinomial-expand"; };
protected:
   typedef TMultinomialExpansionN<3,FScalar>
      FTrinomialExpansionN;
   typedef TMultinomialExpansion1<3,FScalar>
      FTrinomialExpansion1;
   typedef TMultinomialTerm<3,FScalar>
      FTrinomialTerm;
   FTrinomialExpansionN
      m_TrinomialExpansions;
};


FGroupAveragingContext_TrinomialExpand::FGroupAveragingContext_TrinomialExpand(FPointGroup const &PointGroup)
   : FGroupAveragingContext(PointGroup)
{
}


void FGroupAveragingContext_TrinomialExpand::Init(size_t lmin, size_t lmax)
{
   FGroupAveragingContext::Init(lmin, lmax);
   if (s_TimeGroupAveraging >= 2) {
      ct::FTimer tInit;
      m_TrinomialExpansions.Init(lmax);
      size_t
         nt = 0;
      for (auto const &e : m_TrinomialExpansions)
         nt += e.size();
      io.WriteTiming("trinomial expansions", double(tInit), fmt::format("lmax = {}, terms = {}", lmax, nt));
   } else {
      m_TrinomialExpansions.Init(lmax);
   }

   if (0) {
      m_TrinomialExpansions.Print(std::cout, "  │  ");
   }
}





// Compute `(nRg,la+1)`-shape matrix of scalar powers of [R_g]_{i,j}. For fixed i,j
// but open group element index `g` (in rows) and open power (in cols):
static FDenseArray MakeRgijPowerArray(unsigned i, unsigned j, unsigned la, FRotationMatrixList const &Rgs) {
   FDenseArray
      RgijPows(Rgs.size(), la+1); // =  { pow(Rgs[iRg](i,j),iPow) }_{iRg,iPow}
   // set all [R_g]_{i,j}^0 to one. Technically there are some 0^0s there,
   // but `1` appears to be the right choice for those in the current application.
   RgijPows.col(0).setOnes();
   if (la >= 1) {
      // copy the base [R_g]_{i,j} elements from Rgs to the iPow=1 column.
      for (auto &&[iRg, Rg] : enumerate(Rgs))
         RgijPows(iRg,1) = Rg(i,j);
   }
   // build higher powers up to ≤ la (inclusive!) incrementally.
   for (unsigned ip = 2; ip <= la; ++ ip)
      RgijPows.col(ip) = RgijPows.col(ip-1) * RgijPows.col(1);
   return RgijPows;
};

// Compute `(nRg, nTrinom(i_α))`-shape matrix of intermediate Rg-products
//
//     M[iRg,ia] = ([Rg]_{α,0})^a[0] * ([Rg]_{α,1})^a[1] * ([Rg]_{α,2})^a[2]
//
template<class FTrinomialExpansion1>
static FDenseArray MakeRgProdsA(unsigned ixyz, FTrinomialExpansion1 const &TrinExpa1, FRotationMatrixList const &Rgs) {
   unsigned la = TrinExpa1.Degree(); // <- that is what a[0] + a[1] + a[2] will evaluate to for the terms in TrinExpa1
   FDenseArray
      // (nRg,la+1)-shape arrays of scalar powers of [R_g]_{i,j} for i=ixyz and different j
      RgijPows_ixj0 = MakeRgijPowerArray(ixyz,0, la, Rgs),
      RgijPows_ixj1 = MakeRgijPowerArray(ixyz,1, la, Rgs),
      RgijPows_ixj2 = MakeRgijPowerArray(ixyz,2, la, Rgs);
   FDenseArray
      RgProdsA(Rgs.size(), TrinExpa1.size());
   for (auto &&[ia,a] : enumerate(TrinExpa1))
      RgProdsA.col(ia) = RgijPows_ixj0.col(a[0]) * RgijPows_ixj1.col(a[1]) * RgijPows_ixj2.col(a[2]);
   return RgProdsA;
};


static bool _IsAnnihilatedByIco(unsigned ix, unsigned iy, unsigned iz) {
   bool
      ex = (ix % 2 == 0),
      ey = (iy % 2 == 0),
      ez = (iz % 2 == 0);
   bool
      AllEvenOrAllOdd = (ex == ey && ey == ez);
   if (!AllEvenOrAllOdd)
      return true;
   return false;
}


FPolynomialN FGroupAveragingContext_TrinomialExpand::ComputeGroupAveragedMonomial(size_t ix, size_t iy, size_t iz)
{
   cx_assert_rt(ix+iy+iz <= m_lmax); cx_assert_rt(m_TrinomialExpansions.HaveQ(ix+iy+iz));
   // FIXME: this still needs monomial screening and monomial averaging!

   FScalar
      GroupOrderFactor = FScalar(1)/FScalar(m_TrafoList.size());
   FDenseArray
      RgProdsA(MakeRgProdsA(0, m_TrinomialExpansions[ix], m_TrafoList)),
      RgProdsB(MakeRgProdsA(1, m_TrinomialExpansions[iy], m_TrafoList)),
      RgProdsC(MakeRgProdsA(2, m_TrinomialExpansions[iz], m_TrafoList));

   FPolynomialN::FCullCriteria
      AccCull(0,0); // <-- for a start, do not delete *any* small terms while we accumulate terms.
//       AccCull(sqr(g_ThrAlmostZero));
//       AccCull(g_ThrAlmostZero);
   FPolynomialN
      // accumulator for output Rg-transformed monomial. Note: we start with
      // empty monomial (no terms), which is a zero monomial.
      AccM(AccCull);
   FDenseArray1
      RgProdsAB(m_TrafoList.size()); // just a temp array.

   for (auto &&[ia, a] : enumerate(m_TrinomialExpansions[ix])) {
      for (auto &&[ib, b] : enumerate(m_TrinomialExpansions[iy])) {
         RgProdsAB = RgProdsA.col(ia) * RgProdsB.col(ib);
         for (auto &&[ic, c] : enumerate(m_TrinomialExpansions[iz])) {
            FPolynomialN::FMonomial
               Monomial(a[0]+b[0]+c[0], a[1]+b[1]+c[1], a[2]+b[2]+c[2]);
            if (dynamic_cast<FIcosahedralGroup const*>(&m_PointGroup) != 0) {
               // TODO: can I just skip over monomials with zero group averages?
               // Should be easy to try, but I need to patch up CouldMonomialBeATargetFunction
               if (_IsAnnihilatedByIco(Monomial[0], Monomial[1], Monomial[2]))
                  continue;
            }
            FScalar
               TrinomFactor = a.coeff * b.coeff * c.coeff;
            FScalar
               RgFactor = (RgProdsAB * RgProdsC.col(ic)).sum();
            AccM.Add(GroupOrderFactor * TrinomFactor * RgFactor, Monomial);
         }
      }
   }
   // try deleting extra terms at the end?
   if (0) {
      FPolynomialN::FEstimateSquaredNormFn
         pEstNormSq(EstimateMonomialNormSq_S2);
      m_IntdCullCriteriaBase.SetNormEstimateFn(pEstNormSq);
      m_IntdCullCriteriaBase.SetThrCoeffNeglect(g_ThrVerySmall);

      AccM.SetCullCriteria(m_IntdCullCriteriaBase);
   }
//    AccM.SetCullCriteria(FPolynomialN::FCullCriteria(0,0));
   return AccM;
}





struct FCompareMonomialFn
{
   bool operator() (FMonomial const &A, FMonomial const &B) const {
//       int la = A.ix + A.iy + A.iz, lb = B.ix + B.iy + B.iz;
//       if (la < lb) return tru
      // TODO: is this a good order?
      return FPolynomialN::FMonomial(A.ix, A.iy, A.iz) < FPolynomialN::FMonomial(B.ix, B.iy, B.iz);
   }
};


FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FPointGroup const &PointGroup, FPrintLevel iPrintLevel, FGridSearchOptions const *pOptions)
{
   ct::FTimer
      tTotal;
   bool UseScreening = pOptions? pOptions->PrescreenTargetFunctions : true;
//    FMonomialList r = MakeSymmetryReducedMonomialList(lmin, lmax, PointGroup, iPrintLevel, UseScreening);
   FMonomialList r = MakeSymmetryReducedMonomialList(lmin, lmax, PointGroup.pMonomialSymmetry(), iPrintLevel, UseScreening);
   // ^-- hmpf... already with this it makes a difference:
   //     make && aigg_64 '{point-group:O;degree:[4,step:+1,-1];initial-points:hex-grid{7;1;01c}; export:last.dat}'
   // the old one only finds l=24 (can't converge afterwards), the new one taking
   // all permuflections finds this:
   // 230  (lmax = 25, wtsp = 1.237, res = 1.83e-15, mxe = 4.11e-15, edof = 0)
   // I guess the core difference is the full incorporation of the C4 symmetries in
   // the new version. This elimininates many redaundant/competing monomials. Still not happy :/.
   if (iPrintLevel >= 0)
      io.WriteTiming("initial monomial list", double(tTotal));
   return r;
}



#ifdef INCLUDE_ABANDONED
// FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FPointGroup const &PointGroup, FPrintLevel iPrintLevel, FGridSearchOptions const *pOptions)
// {
//    if (g_UseSimpleMonomialEnum) {
// //       g_UseGroupSymmetrizedFn = false; // not supported here (well, one could add the monomials in post processing, but...)
//       ct::FTimer
//          tTotal;
//       FMonomialList r = MakeSymmetryUniqueMonomialListSimple(lmin, lmax, PointGroup, iPrintLevel, pOptions);
//       if (iPrintLevel >= 0)
//          io.WriteTiming("initial monomial list", double(tTotal));
//       if (g_UseGroupSymmetrizedFn) {
//          ct::FTimer
//             tTotal;
//          FGroupAveragingContext_Direct
//             GrpAvg(PointGroup, true);
//          GrpAvg.Init(lmin, lmax);
//          for (size_t i = 0; i != r.size(); ++ i) {
//             r[i].m_RgSymmetrized = GrpAvg.ComputeGroupAveragedMonomial(r[i].ix, r[i].iy, r[i].iz);
//             if (r[i].m_RgSymmetrized.IsZero()) {
//                r.erase(r.begin() + i);
//                -- i;
//                continue;
//             }
// //             if (r[i].m_RgSymmetrized.IsZero())
// //                r[i].m_RgSymmetrized = FPolynomialN(FScalar(1), FMonomialN(r[i].ix, r[i].iy, r[i].iz));
//          }
//          if (s_TargetFnNormalization == TARGETFN_NormalizeSymmetrized) {
//             cx_assert_rt(false); // will break the unscaling code in orthogonal target space construction
//             // this replaces the un-symmetrized monomial normalization factors
//             // by the symmetrized monomial normalization factors.
//             for (size_t i = 0; i != r.size(); ++ i)
//                r[i].m_Factor = FScalar(1)/sqrt(ComputeOverlapRaw(r[i], r[i].m_RgSymmetrized));
//          }
//          if (iPrintLevel >= 0)
//             io.WriteTiming("symmetrizing of monomials", double(tTotal));
//       }
//       return r;
//    }
//
// //    return MakeMonomialList(lmin, lmax);
//    // ^- FIXME: code is disabled atm.
//
// //    ct::FTimer
// //       tXyzPow;
//    ct::FTimer
//       tTotal;
//    FGroupAveragingContext_Direct
//       GrpAvg(PointGroup, true);
// //    FGroupAveragingContext_CachedIncremental
// //    FGroupAveragingContext_CachedIncremental2
// //       GrpAvg(PointGroup);
//
//    GrpAvg.Init(lmin, lmax);
// //    std::cout << FmtTiming("symmetrization (setup)", double(tXyzPow));
//
// //    ct::FTimer
// //       tMonomialScan;
//    size_t
//       nFnMax = nCartX(lmax) - nCartX(lmin-1);
//    std::vector<FPolynomialN>
//       GroupAveragedMonomials;
//    FMonomialList
//       r;
//    r.reserve(nFnMax);
//    GroupAveragedMonomials.reserve(nFnMax);
//
// //    typedef std::multimap<ptrdiff_t, FPolynomialN::FMonomial>
// //       FAverageResultMap;
// //    FAverageResultMap
// //       AverageResults;
//
//
//    for (size_t l = lmin; l <= lmax; ++l) {
//       #pragma omp parallel for schedule(dynamic)
//       for (int ixy_ = 0; ixy_ < int((l+1)*(l+1)); ++ ixy_) {
//          size_t ix = size_t(ixy_/(l+1));
//          size_t iy = size_t(ixy_%(l+1));
//          if (iy <= l - ix) {
//             size_t iz = l - ix - iy;
//
//             if (pOptions && pOptions->PrescreenTargetFunctions)
//                if (!PointGroup.CouldMonomialBeATargetFunction(ix, iy, iz))
//                   continue;
//             if (0) {
//                bool
//                   ex = ix % 2 == 0,
//                   ey = iy % 2 == 0,
//                   ez = iz % 2 == 0;
//                bool
//                   AllOdd = !ex && (ex == ey && ey == ez);
//                if (AllOdd)
//                   continue;
//             }
//
//             FPolynomialN
//                AvgM = GrpAvg.ComputeGroupAveragedMonomial(ix, iy, iz);
//             if (0)
//                // the monomials with uneven powers do not contribute to the target integral...
//                // but if we just delete them, we still get wrong results.
//                AvgM.EraseMonomialsWithOddPowers();
// //             if (AvgM.IsZero()) {
//             AvgM.SetCullCriteria(FPolynomialN::FCullCriteria(g_ThrAlmostZero), true);
//             bool
//                NeglectThisOne = AvgM.IsZero();
// //             AvgM.SetCullCriteria(FPolynomialN::FCullCriteria(g_ThrAlmostZero), false);
// //             FScalar
// //                MxCoeff = 0;
// //             for (FPolynomialN::const_iterator it = AvgM.begin(); it != AvgM.end(); ++ it) {
// //                MxCoeff = std::max(abs(it->second), MxCoeff);
// //             }
// //             if (IsAlmostZero(MxCoeff*MxCoeff))
// //                NeglectThisOne = true;
//             if (0) {
//                FMonomial mi(ix,iy,iz);
//                AvgM.SetCullCriteria(FPolynomialN::FCullCriteria(g_ThrAlmostZero*mi.m_Factor), true);
//                FScalar AvgNrm1 = ComputeOverlapRaw(mi, AvgM);
//                if (AvgNrm1 < 0) AvgNrm1 = 0;
// //                AvgNrm1 = mi.m_Factor * sqrt(AvgNrm1);
// //                AvgNrm1 = sqr(mi.m_Factor) * AvgNrm1;
//                AvgNrm1 = mi.m_Factor * AvgNrm1;
//                NeglectThisOne = IsAlmostZero(AvgNrm1);
//                std::cout << fmt::format("       [{:>6}]  {} M(xyz)[{:2},{:2},{:2}]  AvgNrm1 = {:8.2e}  Skip? {}\n", -1, "zer?", ix, iy, iz, double(AvgNrm1), NeglectThisOne);
//             }
//             if (NeglectThisOne) {
// //                std::cout << fmt::format("       [{:>6}]  {}: M(xyz)[{:2},{:2},{:2}]: avg = {}\n", -1, "zer", ix, iy, iz, AvgM_Copy);
// //                std::cout << fmt::format("       [{:>6}]  {}: M(xyz)[{:2},{:2},{:2}]\n", -1, "zer", ix, iy, iz);
//
//                // it seems like for even l and ico group, only monomials survive
//                // the group averaging in which either all exponents are even or all are odd.
//                // ...but even then not necessarily all of them.
//                if (0) {
//                   bool
//                      ex = ix % 2 == 0,
//                      ey = iy % 2 == 0,
//                      ez = iz % 2 == 0;
//                   bool
//                      AllEvenOrOdd = (ex == ey && ey == ez);
//                   if (AllEvenOrOdd)
//                      std::cout << fmt::format("       [{:>6}]  {}: M(xyz)[{:2},{:2},{:2}]\n", -1, "zer", ix, iy, iz);
//                }
//                // this monomial gets averaged to zero by the group action.
//                // it doesn't do anything in any rule, and can simply be
//                // omitted! With the target symmetry, any rule will integrate
//                // it exactly.
//                continue;
//             }
//             // this will break things, unless the norm factors are also introduced
//             // in many other counterparts. Also, it does not actually fix our issue
//             // of generating linearly dependent sets of target functions if Rg symmetrization
//             // is on. I think I will have to find a way of doing it without Scd
//             // and with actually determining the target function range...
//             //
//             // (maybe with a context object and SmhSolve?...)
// //             if (0) {
// //                FScalar
// //                   fNrmSq = ComputeOverlapRaw(FMonomial(ix,iy,iz), AvgM);
// //                AvgM *= FScalar(1)/sqrt(fNrmSq);
// //                FScalar
// //                   fNrmSq = ComputeOverlapRaw(AvgM, AvgM);
// //                AvgM *= FScalar(1)/sqrt(fNrmSq);
// //             }
//
//             #pragma omp critical
//             {
//                // check if this monomial is equal to one we already have after
//                // group averaging
//                ptrdiff_t
//                   iLast = -1,
//                   iThis; // <- where this one got stored.
//                for (size_t iPrevAvg = 0; iPrevAvg != GroupAveragedMonomials.size(); ++ iPrevAvg) {
// //                   FPolynomialN diff = AvgM - GroupAveragedMonomials[iPrevAvg];
// //                   if (diff.IsZero()) {
//                   if (GroupAveragedMonomials[iPrevAvg] == AvgM) {
//                      iLast = iPrevAvg;
//                      break;
//                   }
//                }
//                if (0) {
//                   iLast = -1; // FIXME: this turns off finding of equivalent monomials!
//                }
//                if (iLast != -1) {
//                   // another monomial which yields the same values on all trial
//                   // points was there (after symmetrization). Increase that one's
//                   // weight and skip this one.
//                   r[iLast].m_Degeneracy += 1;
//                   iThis = iLast;
//                } else {
//                   r.push_back(FMonomial(ix,iy,iz));
//                   if (g_UseGroupSymmetrizedFn)
//                      r.back().m_RgSymmetrized = AvgM;
//                   GroupAveragedMonomials.push_back(AvgM);
//                   iThis = ptrdiff_t(GroupAveragedMonomials.size()) - 1;
//    //                std::cout << fmt::format("       [{:>6}]  add: M(xyz)[{:2},{:2},{:2}]: wt = {}  avg = {}\n", 1+i, r[i].ix, r[i].iy, r[i].iz, r[i].weight, AvgM);
//                }
//                if (g_TrackEquivMonomials)
//                   // this: this also includes (ix,iy,iz) of the target monomial itself
//                   r[iThis].m_SymmetryEquivFn.insert(FPolynomialN::FMonomial(ix,iy,iz));
// //                   AverageResults.insert(FAverageResultMap::value_type(iThis, FPolynomialN::FMonomial(ix,iy,iz)));
// //                if (iPrintLevel >= 2) {
// //                   size_t i = r.size() - 1; char const *pAction = "add: ";
// //                   if (iLast != -1) {
// //                      i = iLast; pAction = "skip:";
// //                   }
// //                   std::cout << fmt::format("       [{:>6}]  {} M(xyz)[{:2},{:2},{:2}]: wt = {}\n", 1+i, pAction, r[i].ix, r[i].iy, r[i].iz, r[i].weight);
// //                }
//                if (iPrintLevel >= 3) {
//                   size_t i = r.size() - 1;
//                   std::string sAction = "add: ";
//                   if (iLast != -1) {
//                      i = iLast; sAction = fmt::format("skip: M(xyz)[{:2},{:2},{:2}] --avg->", ix, iy, iz);
//                   }
//                   std::cout << fmt::format("       [{:>6}]  {} M(xyz)[{:2},{:2},{:2}]: degn = {}\n", 1+i, sAction, r[i].ix, r[i].iy, r[i].iz, r[i].m_Degeneracy);
//                }
//             }
//          }
//       }
//    }
//    if (0) {
//       for (size_t i = 0; i != r.size(); ++ i)
//          r[i].m_Degeneracy = 1;
//    }
//    if (s_TargetFnNormalization == TARGETFN_NormalizeSymmetrized) {
//       // this replaces the un-symmetrized monomial normalization factors
//       // by the symmetrized monomial normalization factors. With this the
//       // target function overlap matrix will consist of exact ones.
//       //
//       // Not sure if we want or do not want this. Either variant appears
//       // to work fine
//       for (size_t i = 0; i != r.size(); ++ i)
//          r[i].m_Factor = FScalar(1)/sqrt(ComputeOverlapRaw(r[i], r[i].m_RgSymmetrized));
//       cx_assert_rt(false);
//    }
//    if (1) {
//       // FIXME: ...I wonder if this fixes the sometimes random changes in parallel runs? E.g., this one:
//       // make && time OMP_NUM_THREADS=12 aigg '{ point-group: icosahedral; degree: [6,step:+2,-1,-2]; max-step: 1e-2; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: hex-grid{5;1;01c}; print: 0; max-it: 128 }'
//       // The rules always come out as l=23 and integrating the target exactly, but with weight spreads between 1.306 to 2.xxx.
//       // Result still unclear.
//
//       // replace each equivalency class, if present, by its minimum representative under
//       // lexicographical order
//       if (g_TrackEquivMonomials) {
//          for (size_t i = 0; i != r.size(); ++ i) {
//             FMonomial::FSymmetryEquivSet const
//                &EquivSet = r[i].m_SymmetryEquivFn;
//             FMonomial::FSymmetryEquivSet::const_iterator
//                it = EquivSet.begin();
//             // this: this also includes (ix,iy,iz) of the target monomial itself
//             r[i].ix = (*it)[0];
//             r[i].iy = (*it)[1];
//             r[i].iz = (*it)[2];
//          }
//       }
//       // sort the entire set.
//       std::sort(r.begin(), r.end(), FCompareMonomialFn());
//
//       // now, the REALLY interesting question: WHY CAN THE FUNCTION ORDER MAKE A DIFFERENCE?
//    }
//    PrintMonomialList(std::cout, r, iPrintLevel); // note: this doesn't do anything unless iPrintLevel >= 2;
// //    if (iPrintLevel >= 1) {
// //       std::cout << " Target function list:\n";
// //       for (size_t i = 0; i != r.size(); ++ i) {
// // //          std::cout << fmt::format("{:>6}  M(xyz)[{:2},{:2},{:2}]  wt = {}\n", 1+i, r[i].ix, r[i].iy, r[i].iz, r[i].weight);
// //          std::string sAdditional;
// //          if (g_TrackEquivMonomials) {
// //             fmt::MemoryWriter w;
// //             w << "  eq-set = [";
// // //             std::pair<FAverageResultMap::const_iterator,FAverageResultMap::const_iterator>
// // //                itEq = AverageResults.equal_range(i);
// // //             for (FAverageResultMap::const_iterator it = itEq.first; it != itEq.second; ++ it) {
// // //                if (it != itEq.first)
// // //                   w << ", ";
// // //                FPolynomialN::FMonomial const &m = it->second;
// // // //                w << fmt::format("M(xyz)[{:2},{:2},{:2}]", m[0], m[1], m[2]);
// // //                w << fmt::format("({:2},{:2},{:2})", m[0], m[1], m[2]);
// // //             }
// //             FMonomial::FSymmetryEqivSet const
// //                &EqivSet = r[i].m_SymmetryEquivFn;
// //             for (FMonomial::FSymmetryEqivSet::const_iterator it = EqivSet.begin(); it != EqivSet.end(); ++ it) {
// //                if (it != EqivSet.begin())
// //                   w << ", ";
// //                FPolynomialN::FMonomial const &m = *it;
// // //                w << fmt::format("M(xyz)[{:2},{:2},{:2}]", m[0], m[1], m[2]);
// //                w << fmt::format("({:2},{:2},{:2})", m[0], m[1], m[2]);
// //             }
// //             w << "]";
// //             sAdditional = w.str();
// //          }
// //
// //          if (iPrintLevel >= 3)
// //             sAdditional += fmt::format("  Rg-avg = {}", r[i].m_RgSymmetrized);
// //
// //          std::cout << fmt::format("{:>6}  M(xyz)[{:2},{:2},{:2}]  wt = {}{}\n", 1+i, r[i].ix, r[i].iy, r[i].iz, r[i].weight, sAdditional);
// //       }
// //    //    if (lmax == 4) throw std::runtime_error("absicht");
// //    }
//    io.WriteTiming("symmetrization of monomials", double(tTotal));
// //    Rg_xyzn.clear();
//    return r;
// }
// #endif
#endif // INCLUDE_ABANDONED

bool HasSymmetrizedFnQ(FMonomialList const &Monomials) {
   for (FMonomial const &m : Monomials)
      if (!m.m_RgSymmetrized.empty())
         return true;
   // ^- check if *any* of the monomials has a non-zero symmetrization stored.
   //    (in this case, all the other monomials should have one, too). Note
   //    that, technically, an empty polynomial is also valid (it's a
   //    representation of zero), so this is why we do not check this in the
   //    loop below.
   return false;
}

std::string FormatLargeFloatAsInteger(FScalar const &fValue)
{
   using std::log10;
   std::stringstream ss;
   int nDigits10 = std::numeric_limits<FScalar>::digits10; // <-- that's only 15 for doubles; it indicates the number of decimal digits in integers which survive int -> float -> int unchanged.
   int npre = 1 + int(-FScalar(1e-10) + log10(fValue)); // estimate num digits to the left of period
   if (npre < 0) npre = 0;
   if (npre <= nDigits10) {
      ss.precision(nDigits10 - npre);
//       ss.fill(' ', nDigits10 - npre);
      for (size_t iw = 0; iw < size_t(nDigits10 - npre); ++iw)
         ss << " "; // <-- for decimal point alignment
      ss << std::fixed;
   } else {
      ss.precision(nDigits10 - npre);
      ss << std::scientific;
   }
   ss << fValue;
   return ss.str();
}

void PrintMonomialList(std::ostream &xout, FMonomialList const &r, FPrintLevel iPrintLevel, std::string const &Comment)
{
   if (iPrintLevel.MoreQ()) {
      if (Comment.empty())
         xout << " Target function list:\n";
      else
         xout << fmt::format(" Target function list ({}):\n", Comment);
      bool
         HaveSymmetrizedFn = HasSymmetrizedFnQ(r);
      for (size_t i = 0; i != r.size(); ++ i) {
         std::string sAdditional;
         if (g_TrackEquivMonomials) {
            fmt::MemoryWriter w;
            w << "  eq-set = [";
#ifdef INCLUDE_ABANDONED
//             std::pair<FAverageResultMap::const_iterator,FAverageResultMap::const_iterator>
//                itEq = AverageResults.equal_range(i);
//             for (FAverageResultMap::const_iterator it = itEq.first; it != itEq.second; ++ it) {
//                if (it != itEq.first)
//                   w << ", ";
//                FPolynomialN::FMonomial const &m = it->second;
// //                w << fmt::format("M(xyz)[{:2},{:2},{:2}]", m[0], m[1], m[2]);
//                w << fmt::format("({:2},{:2},{:2})", m[0], m[1], m[2]);
//             }
#endif // INCLUDE_ABANDONED
            FMonomial::FSymmetryEquivSet const
               &EquivSet = r[i].m_SymmetryEquivFn;
            for (FMonomial::FSymmetryEquivSet::const_iterator it = EquivSet.begin(); it != EquivSet.end(); ++ it) {
               if (it != EquivSet.begin())
                  w << ", ";
               FPolynomialN::FMonomial const &m = *it;
               w << fmt::format("({:2},{:2},{:2})", m[0], m[1], m[2]);
            }
            w << "]";
            sAdditional = w.str();
         }

         if (iPrintLevel.MostQ()) {
            FScalar fIntS2 = CalcMonomialUnitSphereIntegralRaw(r[i].ix, r[i].iy, r[i].iz);
            std::string ss;
            if (fIntS2 == 0)
               ss = "1/0 (by symmetry)";
            else {
               ss = FormatLargeFloatAsInteger(1/fIntS2);
               // ^-- I specifically do *not* want scientific
            }
//                fmt::format("{:.18g}", 1/fIntS2);
            sAdditional += fmt::format("  1/⟨1,M⟩ = {}", ss);
         }
         if (HaveSymmetrizedFn && iPrintLevel.AllQ())
            sAdditional += fmt::format("  Rg-avg = {}", r[i].m_RgSymmetrized);

         xout << fmt::format("{:>6}  M(xyz)[{:2},{:2},{:2}]  1/nfac = {:8.2e}  degn = {}{}\n",
            1+i, r[i].ix, r[i].iy, r[i].iz, 1/double(r[i].m_Factor), double(r[i].m_Degeneracy), sAdditional);
      }
   //    if (lmax == 4) throw std::runtime_error("absicht");
   }
}




FScalar ComputeOverlapRaw(FMonomial const &mi, FPolynomialN const &mj_RgSymmetrized)
{
   // note: since the group averaging operator, defined by
   //
   //     P̂_G f(r⃗) := |G|⁻¹ ∑_{ĝ ∈ G} f(ĝ⁻¹ r⃗)
   //
   // is a projector (i.e., P̂_G = (P̂_G)^2), for getting the overlap matrix
   // over symmetrized functions it is still sufficient to evaluate the
   // symmetrization over one of the monomials mᵢ, mⱼ only:
   //
   //     S_{ij} = ⟨P̂_G mᵢ| P̂_G mⱼ⟩
   //            = ⟨mᵢ|P̂_G P̂_G|mⱼ⟩
   //            = ⟨mᵢ|P̂_G|mⱼ⟩
   //            = ⟨mᵢ|(P̂_G mⱼ)⟩
   //
   // That is what we do here.
   FPolynomialN::const_iterator
      itMj;
   FScalar
      s = 0;
   for (itMj = mj_RgSymmetrized.begin(); itMj != mj_RgSymmetrized.end(); ++ itMj) {
      s += itMj->second * CalcMonomialUnitSphereIntegralRaw(mi.ix + itMj->first[0], mi.iy + itMj->first[1], mi.iz + itMj->first[2]);
   }
   return s;
}
#ifdef INCLUDE_ABANDONED
/* tex version: P_G f(\vec r) := \frac{1}{\abs{G}}\sum_{g\in G} f(\hat g^{-1} \vec r) */
#endif // INCLUDE_ABANDONED


FMonomialN AsMonomialN(FMonomial const &m) {
   return FMonomialN(m.ix, m.iy, m.iz);
}


FScalar ComputeOverlapRaw(FScalar f1byNumAxpyElem, FMonomialSymmetry::FMonomialCoeffList const &mi_AxprSum, FMonomialN const &mj)
{
   // note: since the group averaging operator, defined by
   //
   //     P̂_G f(r⃗) := |G|⁻¹ ∑_{ĝ ∈ G} f(ĝ⁻¹ r⃗)
   //
   // is a projector (i.e., P̂_G = (P̂_G)^2), for getting the overlap matrix
   // over symmetrized functions it is still sufficient to evaluate the
   // symmetrization over one of the monomials mᵢ, mⱼ only:
   //
   //     S_{ij} = ⟨P̂_G mᵢ| P̂_G mⱼ⟩
   //            = ⟨mᵢ|P̂_G P̂_G|mⱼ⟩
   //            = ⟨mᵢ|P̂_G|mⱼ⟩
   //            = ⟨mᵢ|(P̂_G mⱼ)⟩
   //
   // That is what we do here.
   FScalar
      s = 0;
   for (auto &&[mi,iFactor] : mi_AxprSum) {
      s += iFactor * CalcMonomialUnitSphereIntegralRaw(mi[0] + mj[0], mi[1] + mj[1], mi[2] + mj[2]);
   }
   return f1byNumAxpyElem*s;
}


FScalar ComputeOverlapRaw(FPolynomialN const &poly_i, FPolynomialN const &poly_j)
{
   FPolynomialN::const_iterator
      itMj, itMi;
   FScalar
      s = 0;
   for (itMi = poly_i.begin(); itMi != poly_i.end(); ++ itMi) {
      for (itMj = poly_j.begin(); itMj != poly_j.end(); ++ itMj) {
         s += (itMi->second * itMj->second) * CalcMonomialUnitSphereIntegralRaw(
            itMi->first[0] + itMj->first[0],
            itMi->first[1] + itMj->first[1],
            itMi->first[2] + itMj->first[2]);
      }
   }
   return s;
}


FDenseMatrix MakeOverlapMatrix(FMonomialList const &Monomials)
{
   // FIXME: This should use the symmetrized monomials, too!
   // should store them when making them, and then iterate over one side of this (due to the projector
   // property we only need it once).
   size_t
      N = Monomials.size();
   bool
      HaveSymmetrizedFn = HasSymmetrizedFnQ(Monomials);
   for (size_t i = 0; i != N; ++ i)
      // make sure we either have the symmetrized ones for all functions, or
      // for none of them.
      assert(!Monomials[i].m_RgSymmetrized.empty() == HaveSymmetrizedFn);

//    if (1) {
//       // hmpf. I think we may need to modify the `rhs` vectors if we use this as
//       // overlap matrix — i.e., get for `rhs` (`CalcUnitSphereIntegral`) also the
//       // result of the symmetrized function, not of the raw monomial. However,
//       // unless we use the symmetrized functions *everywhere*, I wonder if this
//       // would be compatible with my actual main optimization procedure…
//       HaveSymmetrizedFn = false;
//    }
   // ^- I think the answer should be "yes!". I can just use the group-averaged
   //    functions everywhere, and I think I need not even modify the core evaluator.
   //    That should all come in implicitly, I think.
   //    Just need to be careful about factors of |G|: The paper definition for
   //    K does not have the 1/|G| factors, but the group averages here do.
   //    Let's try and see what we get...
//    if (HaveSymmetrizedFn)
//       std::cout << "\n##! NOTE: using symmetrized overlap matrix\n" << std::endl;
   FDenseMatrix
      S(N,N);
   bool EnforceSymmetry = true;
//    EnforceSymmetry = false;
   for (size_t i = 0; i != N; ++ i) {
      for (size_t j = 0; j <= (EnforceSymmetry ? i : (N-1)); ++ j) {
         FMonomial const
            &mi = Monomials[i],
            &mj = Monomials[j];
         FScalar
            s = 0;
         if (!HaveSymmetrizedFn) {
            // using raw monomials? (no group averaging)
            FMonomial
               mij(mi.ix + mj.ix, mi.iy + mj.iy, mi.iz + mj.iz);
            s = mij.CalcUnitSphereIntegralRaw();
   //          s = mij.CalcUnitSphereIntegralRaw() * (mi.weight * mj.weight);
            // ^- hm.. that would be the product of the monomial with weight:
            // product of coefficients and sum of powers. not division by weights
            // or square roots or etc.
            // … actually, should there be *any* weights here? On second thought,
            // I think the answer should be "no weights in S", even if there
            // *are* weights in the calculation of points (`EvalPoint`) and rhs
            // integrals (`CalcUnitSphereIntegral`).
         } else {
//             // note: since the group averaging operator defined by
//             //
//             //     P_G f(\vec r) := \frac{1}{\abs{G}}\sum_{g\in G} f(\hat g^{-1} \vec r)
//             //
//             // is a projector (i.e., P_G=P_G^2), for getting the overlap matrix
//             // over symmetrized functions it is still sufficient to evaluate the
//             // symmetrization over one of the monomials m_i,m_j only:
//             //
//             //     S_{ij} = <P_G m_i| P_G m_j>
//             //            = <m_i|P_G P_G|m_j>
//             //            = <m_i|P_G|m_j>
//             //            = <m_i|P_G m_j>
//             //
//             // That is what we do here.
//             FPolynomialN::const_iterator
//                itMj;
//             for (itMj = mj.m_RgSymmetrized.begin(); itMj != mj.m_RgSymmetrized.end(); ++ itMj) {
// //                FMonomial
// //                   mij(mi.ix + itMj->first[0], mi.iy + itMj->first[1], mi.iz + itMj->first[2]);
// //                s += itMj->second * mij.CalcUnitSphereIntegralRaw();
//                s += itMj->second * CalcMonomialUnitSphereIntegralRaw(mi.ix + itMj->first[0], mi.iy + itMj->first[1], mi.iz + itMj->first[2]);
//             }
            s = ComputeOverlapRaw(mi, mj.m_RgSymmetrized);
         }
         s *= mi.m_Factor * mj.m_Factor;
         S(i,j) = s;
         if (EnforceSymmetry)
            S(j,i) = s;
      }
   }
   return S;
}


FDenseMatrix MakeOverlapMatrix_AxprAveraged(FMonomialList const &Monomials, FMonomialSymmetry const *pMonomialSymmetry)
{
   size_t
      N = Monomials.size();
   FDenseMatrix
      S(N,N);
   S.setZero();
   bool EnforceSymmetry = true;
//    EnforceSymmetry = false;
   FMonomialSymmetry::FMonomialCoeffList
      mi_AxprSum;
   FScalar
      fAxprScale = 1/FScalar(pMonomialSymmetry->EasySubgroup().size());
   // ^-- or would 1 be better? that would retain the group sum instead of
   //     average. Would not be a projector anymore, though...
   //     And, technically, the raw norm factors would also not be quite right.
   //     On the other hand, all of that should implicitly be fixed
   //     when doing a proper implicit orthogonalization for the actual
   //     overlap of the functions we have.
   for (size_t i = 0; i != N; ++ i) {
      FMonomial const
         &mi = Monomials[i];
      if (!pMonomialSymmetry->EvalEasySubgroupSum(mi_AxprSum, AsMonomialN(mi)))
         // mi gets annihilated by permuflection group average.
         continue;

      for (size_t j = 0; j <= (EnforceSymmetry ? i : (N-1)); ++ j) {
         FMonomial const
            &mj = Monomials[j];
         FScalar
            s = ComputeOverlapRaw(fAxprScale, mi_AxprSum, AsMonomialN(mj));
         s *= mi.m_Factor * mj.m_Factor;
         S(i,j) = s;
         if (EnforceSymmetry)
            S(j,i) = s;
      }
   }
   return S;
}




// FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FPointGroup const &PointGroup, int iPrintLevel)
// {
//    return
// }

int TestMemoryAndTargetFnScreening()
{
   FIcosahedralGroup
//    FTetrahedralGroup
//    FOctahedralGroup
      PointGroup;
   FGridSearchOptions
      Options("{prescreen-target-fn: no}");
//       Options("{prescreen-target-fn: yes}");
   FGridSearchOptions
      *pOptions = &Options;
//    pOptions = 0;

//    FRotationMatrixList const
//       &SymmetryOps = PointGroup.GetOps();
//    for (int i = 0; i < 1000; ++ i) {
//    for (int i = 0; i < 10; ++ i) {
   for (int i = 0; i < 1; ++ i) {
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(30, 30, PointGroup, 2);

      // so.. with the Ico group it seems we either get only non-zero
      // monomials with either ix,iy,iz all even (if ix+iy+iz is even)
      // or ix,iy,iz all odd (if ix+iy+iz is odd).
      // However, that probably depends on the exact axis alignment if
      // the icosahedron with respect to our x,y,z axes.
      // For now I'll just hardcode it.
      FMonomialList ml = MakeSymmetryUniqueMonomialList(20, 20, PointGroup, 2, pOptions);
      FMonomialList ml1 = MakeSymmetryUniqueMonomialList(21, 21, PointGroup, 2, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(20, 20, PointGroup, 2, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(20, 21, PointGroup, 2, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(20, 21, PointGroup, 1, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(30, 31, PointGroup, 1, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(40, 41, PointGroup, 1, pOptions);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(10, 10, PointGroup, 2, pOptions);
//       FMonomialList ml1 = MakeSymmetryUniqueMonomialList(11, 11, PointGroup, 2, pOptions);
//       FMonomialList ml1 = MakeSymmetryUniqueMonomialList(31, 31, PointGroup, 1, pOptions);
//       FMonomialList ml1 = MakeSymmetryUniqueMonomialList(41, 41, PointGroup, 1, pOptions);
   }
   return 1;
};


// void ComputeAndStoreGroupAverages(FMonomialList &r, FPointGroup const &PointGroup, FGridSearchOptions const &Options, bool DeleteZeroAverages)
void ComputeAndStoreGroupAverages(FMonomialList &r, FPointGroup const &PointGroup, bool DeleteZeroAverages)
{
//    std::cout << FmtTiming("initial monomial list", double(tTotal));
   size_t lmin = size_t(-1), lmax = 0;
   for (FMonomial &m : r) {
      lmin = std::min(lmin, m.l());
      lmax = std::max(lmax, m.l());
   }

   ct::FTimer
      tTotal;
#ifdef INCLUDE_ABANDONED
//    FGroupAveragingContext_Direct
//       GrpAvg(PointGroup, true); // `true`: cache one-index coordiante powers (i.e., x^i, y^i, z^i individually) for all Rg transforms
//    FGroupAveragingContext_TrinomialExpand
//       GrpAvg(PointGroup);
//    GrpAvg.Init(lmin, lmax);
//    for (size_t i = 0; i != r.size(); ++ i) {
//       r[i].m_RgSymmetrized = GrpAvg.ComputeGroupAveragedMonomial(r[i].ix, r[i].iy, r[i].iz);
//       if (DeleteZeroAverages && r[i].m_RgSymmetrized.IsZero()) {
//          r.erase(r.begin() + i);
//          -- i;
//          continue;
//       }
//    }
#endif // INCLUDE_ABANDONED
   FGroupAveragingContext_TrinomialExpand
      GrpAvg(PointGroup);
//    FGroupAveragingContext_Direct
//       GrpAvg(PointGroup, true); // `true`: cache one-index coordiante powers (i.e., x^i, y^i, z^i individually) for all Rg transforms
   GrpAvg.Init(lmin, lmax);
   FMonomialList
      Keep;
   Keep.reserve(r.size());
   for (size_t i = 0; i != r.size(); ++ i) {
      FPolynomialN gAvg = GrpAvg.ComputeGroupAveragedMonomial(r[i].ix, r[i].iy, r[i].iz);
      if (!(DeleteZeroAverages && gAvg.IsZero())) {
         Keep.push_back(r[i]);
         Keep.back().m_RgSymmetrized.swap(gAvg);
      }
   }
   r.swap(Keep);
#ifdef INCLUDE_ABANDONED
//    if (s_TargetFnNormalization == TARGETFN_NormalizeSymmetrized) {
//       // this replaces the un-symmetrized monomial normalization factors
//       // by the symmetrized monomial normalization factors.
//       for (size_t i = 0; i != r.size(); ++ i)
//          r[i].m_Factor = FScalar(1)/sqrt(ComputeOverlapRaw(r[i], r[i].m_RgSymmetrized));
//    }
//    std::cout << FmtTiming("symmetrizing of monomials", double(tTotal));
#endif // INCLUDE_ABANDONED
}


typedef std::vector<FPolynomialN>
   FPolynomialList;

enum FComputeGroupAveragesFlags {
   GROUPAVERAGE_ReportTiming = 0x0001,
   // if set, do not add polynomials which are manifestly zero to the output lists.
   GROUPAVERAGE_OmitZeros = 0x0002
};

template<class TGroupAveragingContext, class... _ArgTypes>
FPolynomialList ComputeGroupAveragesT(FMonomialList const &Monomials, FPointGroup const &PointGroup, unsigned Flags, _ArgTypes&&... ContextArgs)
{
   size_t lmin = size_t(-1), lmax = 0;
   for (FMonomial const &m : Monomials) {
      lmin = std::min(lmin, m.l());
      lmax = std::max(lmax, m.l());
   }

   ct::FTimer
      tTotal;

   FPolynomialList
      out;
   out.reserve(Monomials.size());

   TGroupAveragingContext
      AvgContext(PointGroup, std::forward<_ArgTypes>(ContextArgs)...);
   AvgContext.Init(lmin, lmax);
   for (FMonomial const &m : Monomials) {
      out.emplace_back(AvgContext.ComputeGroupAveragedMonomial(m.ix, m.iy, m.iz));
      if (bool(Flags & GROUPAVERAGE_OmitZeros) && out.back().IsZero()) {
         out.pop_back();
      }
   }

   if (bool(Flags & GROUPAVERAGE_ReportTiming))
      io.WriteTiming("group averaged monomials", double(tTotal), AvgContext.Desc());
   return out;
}


void _CompareGroupAveragingRoutines(FMonomialList const &Monomials, FPointGroup const &PointGroup)
{
   int
      Verbosity = 1;
   unsigned
      Flags = 0;
   if (s_TimeGroupAveraging >= 1)
      Flags |= GROUPAVERAGE_ReportTiming;
   FPolynomialList
      pl1 = ComputeGroupAveragesT<FGroupAveragingContext_Direct>(Monomials, PointGroup, Flags, true),
      pl2 = ComputeGroupAveragesT<FGroupAveragingContext_TrinomialExpand>(Monomials, PointGroup, Flags);
   char const *AvgDesc1 = "direct";
   char const *AvgDesc2 = "trinom";

   if (bool(Flags & GROUPAVERAGE_ReportTiming))
      io.WriteLine();

   std::stringstream ssPass, ssErr;
   size_t nPass = 0, nErr = 0;

   assert(pl1.size() == Monomials.size() && pl2.size() == Monomials.size());
   for (size_t im = 0; im != Monomials.size(); ++ im) {
      if (im != 0 && Verbosity >= 3)
         io.Write(Repeat(78, "─", "  ", "\n"));
      FMonomial const &m = Monomials[im];
      char const *pFmt = "      {:<8} P_G[m_{{{},{},{}}}] = {}\n";
      FPolynomialN
         &p1 = pl1[im], &p2 = pl2[im];
      if (1) {
         FPolynomialN::FCullCriteria Cull(g_ThrAlmostZero);
         p1.SetCullCriteria(Cull, true); // true: apply purge with new criteria now
         p2.SetCullCriteria(Cull, true); // true: apply purge with new criteria now
      }
      FPolynomialN
         p12_diff = p1 - p2;
      if (Verbosity == 1) {
         std::string sTerm = fmt::format("P_G[m_{{{},{},{}}}]", m.ix, m.iy, m.iz);
         if (p12_diff.IsZero()) {
            if (nPass != 0) ssPass << ", ";
            ssPass << sTerm;
            nPass += 1;
         } else {
//             io.Write("\n WARNING: Comparison failed: diff({}) = {}", sTerm, p12_diff);
            ssErr << fmt::format("\n WARNING: Comparison failed: diff({}) = {}", sTerm, p12_diff);
            nErr += 1;
         }
      }
      if (Verbosity >= 3) {
         io.Write(pFmt, AvgDesc1, m.ix, m.iy, m.iz, p1);
         io.Write(pFmt, AvgDesc2, m.ix, m.iy, m.iz, p2);
      }
      if (Verbosity >= 2) {
         io.Write(pFmt, "diff", m.ix, m.iy, m.iz, p12_diff);
      }
   }
   if (Verbosity == 1) {
      io.Write(" {:<15} {}", fmt::format("Passed({}):", nPass), ssPass.str());

      if (nErr != 0)
//          io.Write(" {:<15} {}", fmt::format("Failed({}):", nErr), ssErr.str());
         io.Write(" {:<15} {}", fmt::format("Failed({}):", nErr), "ERROR PRINT DISABLED");
   }

}


// std::string MakeTargetSpaceDefaultFileName(FPointGroup const &PointGroup, int l)
// {
//    ct::FFileLocator
//       FileLocator;
//    std::string
//       BaseName(fmt::format("minimal-target-space-S2-{}-{}.dat", PointGroup.Name(NAMETYPE_Symbol), l)),
//       SpaceFileName(FileLocator.JoinPath("data", BaseName));
//    return SpaceFileName;
// }
//
//
// typedef std::vector<FPolynomialNPtr>
//    FPolynomialList;
//
// void WriteTargetSpaceDataFile(FPointGroup const &PointGroup, int l, FPolynomialList const &TargetFns, FTargetSpacePreprocessOptions const &Options)
// {
//    std::string
//       FileName = MakeTargetSpaceDefaultFileName(PointGroup, int l);
//    std::ofstream
//       out(FileName.c_str(), std::ofstream::trunc | std::ofstream::out);
//
//    char const *pFunctionListType = "polyN";
//    out << fmt::format("target-space{{symm: {}; order: {}, type: {}, fp: {}}}\n", PointGroup.Name(NAMETYPE_Symbol), l, pFunctionListType, g_ScalarFloatSizeBits);
//    if (1) {
//       out << "# generated ";
//       if (Options.StoreIdentityInfo) {
//          std::time_t t = time(0);
//          struct std::tm *now = std::gmtime(&t);
//          if (now) {
//             out << fmt::format("{:04}-{:02}-{:02} {:02}:{:02} (UTC)", now->tm_year+1900, now->tm_mon+1, now->tm_mday, now->tm_hour, now->tm_min);
//          }
//       }
//       if (Options.StoreVersionInfo) {
//          if (Options.StoreIdentityInfo)
//             out << " ";
//          out << "by " << MakeAiggVersionString();
//       }
//       out << "\n";
//    }
// //    if (Options.StoreOptions)
// //       out << fmt::format("# options: '{}'\n", Options.MakeOptionString());
// }
//
//
// void TestSymmetrizedSpaceSize(int l, FPointGroup const &PointGroup, FGridSearchOptions const &Options)
// {
//    FGridSearchOptions const
//       *pOptions = &Options;
//    pOptions = 0;
//    int lmin = l-1;
//    if (lmin < 0) lmin = 0;
//    int lmax = l;
// //    g_UseGroupSymmetrizedFn = true;
//    g_UseGroupSymmetrizedFn = false;
//    FMonomialList ml = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, -1, pOptions);
//    g_UseGroupSymmetrizedFn = true;
//    // ^-- do I really want symmetry-unique ones only?
//    // UPDATE: I think I have to go with this, and then add the group average members up.
//    // It's so costly...
//    FMonomialList ml_all = MakeMonomialList(lmin, lmax);
// //    std::cout << fmt::format("G = {}  l = {}  #monomials = {}", PointGroup.Name(NAMETYPE_Symbol), l, ml.size()) << std::endl;
//    size_t nMonomialsOrig = ml.size();
//    ComputeAndStoreGroupAverages(ml, PointGroup, Options, true); // last: delete zero-averaged polynomials?
//    FDenseMatrix
//       S = MakeOverlapMatrix(ml);
// //    std::cout << fmt::format("| S.shape = ({},{})", S.rows(), S.cols()) << std::endl;
//    size_t nNonZero = 0;
//    FPolynomialList
//       // output list of orthogonalized significant (i.e., non-redundant) group-symmetrized
//       // target functions for the current symmetry group and input monomial list.
//       OrthTargetFns;
//    OrthTargetFns.reserve(ml.size());
//    if (S.rows() != 0) {
//       FEigh
//          Sed(S, FEigh::LargeEwFirst);
//       using std::abs;
//       for (size_t iEw = 0; iEw < size_t(Sed.ew.size()); ++ iEw)
//          if (abs(Sed.ew[iEw]) > g_ThrAlmostZero) {
//             nNonZero += 1;
//             FPolynomialNPtr
//                pPolyEv = new FPolynomialN();
//             for (size_t iComp = 0; iComp != size_t(Sed.ev.rows()); ++ iComp)
//                pPolyEv->Add(Sed.ev(iComp,iEw), ml[iComp].m_RgSymmetrized);
//             OrthTargetFns.push_back(pPolyEv);
//          }
//    }
//    std::cout << fmt::format("G = {:2}  l = {:<3}  #monomials[lmin: {:<3} lmax: {:<3} all: {:<5} req: {:<5} nonz-P_G = {:<5}] #linearly-independent = {}", PointGroup.Name(NAMETYPE_Symbol), l, lmin,lmax, ml_all.size(), nMonomialsOrig, ml.size(), nNonZero) << std::endl;
//    // ^- wait.. in I, apart from l=1, none of them have zero group averages?!
//
//    // ??next steps:
//    // - get coverage of monomials (\sum ev[:,i]^2)
//    // - for all non-zero ones, sort them by order (large order first)
//    // - project & orthogonalize to form smaller basis sets.
//    WriteTargetSpaceDataFile(PointGroup, l, OrthTargetFns);
// }
//
//
// int TestSymmetrizedSpaceSizes()
// {
//    FIcosahedralGroup PointGroup;
// //    FOctahedralGroup PointGroup(FOctahedralGroup::GROUP_Full_Oh);
//    FGridSearchOptions
//       Options("{}");
//    for (int l = 0; l <= 20; l += 1) {
//       TestSymmetrizedSpaceSize(l, PointGroup, Options);
//    }
//    return 1;
// }








int OverlapOverflowTest()
{
   io.Write("enter sym test.");
//    io.Write("before S(0,0,0).");
//    FScalar s= CalcMonomialUnitSphereIntegralRaw(0, 0, 0);
//    io.Write("after S(0,0,0).");
   auto make_even = [](size_t x) { return x & (~1); };
   for (unsigned l = 0; l < 90; ++ l) {
//       io.Write("next: = {}  make_even(l) = {}  make_even(l/2) = {}", l, make_even(l), make_even(l/2));
      FScalar
         Sxx = CalcMonomialUnitSphereIntegralRaw(2*l, 0, 0),
         Syy = CalcMonomialUnitSphereIntegralRaw(make_even(l), make_even(l), 0),
         Szz = CalcMonomialUnitSphereIntegralRaw(make_even(l), make_even(l/2), make_even(l/2));
//       io.Write("after compute, before write.");
//       io.Write("l = {:2}  Sxx = {}  Syy = {}  Szz = {}", l, Sxx, Syy, Szz);
         // ^-- hm... freezes. with aigg_256. no idea. don't want to look.
      io.Write("l = {:2}  Sxx = {}  Syy = {}  Szz = {}", l, double(Sxx), double(Syy), double(Szz));
   }
   return -1;
}



void _TestGroupAveraging1(int l, FPointGroup const &PointGroup, FGridSearchOptions const &Options)
{
   FGridSearchOptions const
      *pOptions = &Options;
   pOptions = 0;
   int lmin = l-1;
   if (lmin < 0) lmin = 0;
   int lmax = l;
   FMonomialList ml = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, -1, pOptions);
   FMonomialList ml_all = MakeMonomialList(lmin, lmax);

   _CompareGroupAveragingRoutines(ml, PointGroup);
}


int _TestGroupAveraging()
{
   FIcosahedralGroup PointGroup;
//    FOctahedralGroup PointGroup(FOctahedralGroup::GROUP_Full_Oh);
   FGridSearchOptions
      Options("{}");
//    int lmax = 10, lmin = 0;
   int lmax = 30, lmin = lmax;
//    int lmax = 20, lmin = lmax;
   for (int l = lmin; l <= lmax; l += 1) {
      io.Write(Repeat(80, "━", "\n"));
//       io.Write("░░ {} // monomials order {}", PointGroup.Name(NAMETYPE_Group), l);
      io.Write("░░ {:<74} ░░", fmt::format("{} // monomials order {}", PointGroup.Name(NAMETYPE_Group), l));
      io.Write(Repeat(80, "━", "", "\n"));
      _TestGroupAveraging1(l, PointGroup, Options);
   }
   return 1;
}



// FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FMonomialSymmetry const *pMonomialSymmetry, FPrintLevel iPrintLevel)
// {
//    // we got a more or less closed list of what would end up in the target
//    // monomials. Nevertheless, we collect them in a set first because this is
//    // easy and it allows imposing a well-defined order for free.
//    typedef FMonomialN_SmallIntegralsFirstOrder FMonomialOrder1;
// //    typedef FMonomialN_ExportOrderPred FMonomialOrder1;
// //    typedef std::less<FMonomialN> FMonomialOrder1;
//    typedef std::set<FMonomialN, FMonomialOrder1>
//       FMonomialN_Set;
//    FMonomialN_Set
//       RetainedMonomials;
//
//    for (size_t l = lmin; l <= lmax; ++l) {
//       for (size_t ix = 0; ix <= l; ++ ix) {
//          for (size_t iy = 0; iy <= l - ix; ++ iy) {
//             size_t iz = l - ix - iy;
//             FMonomialN
//                mijk(ix, iy, iz);
// //             if (!pOptions || (pOptions && pOptions->PrescreenTargetFunctions))
//             if (!pMonomialSymmetry->SelectAsCanonialReprQ(mijk))
//                continue;
//             RetainedMonomials.insert(mijk);
//          }
//       }
//    }
//    FMonomialList
//       r;
//    r.reserve(RetainedMonomials.size());
//    for (FMonomialN const &m : RetainedMonomials)
//       r.push_back(FMonomial(m[0], m[1], m[2]));
//
// //    iPrintLevel = 2;
//    PrintMonomialList(std::cout, r, iPrintLevel); // note: this doesn't do anything unless iPrintLevel >= 2;
//    return r;
// }
//
//
// FMonomialList MakeSymmetryUniqueMonomialListSimple(size_t lmin, size_t lmax, FPointGroup const &PointGroup, FPrintLevel iPrintLevel, FGridSearchOptions const *pOptions)
// {
//    // we got a more or less closed list of what would end up in the target
//    // monomials. Nevertheless, we collect them in a set first because this is
//    // easy and it allows imposing a well-defined order for free.
//    typedef FMonomialN_SmallIntegralsFirstOrder FMonomialOrder1;
// //    typedef FMonomialN_ExportOrderPred FMonomialOrder1;
// //    typedef std::less<FMonomialN> FMonomialOrder1;
//    typedef std::set<FMonomialN, FMonomialOrder1>
//       FMonomialN_Set;
//    FMonomialN_Set
//       RetainedMonomials;
//
//    for (size_t l = lmin; l <= lmax; ++l) {
//       for (size_t ix = 0; ix <= l; ++ ix) {
//          for (size_t iy = 0; iy <= l - ix; ++ iy) {
//             size_t iz = l - ix - iy;
//             if (!pOptions || (pOptions && pOptions->PrescreenTargetFunctions))
//                if (!PointGroup.CouldMonomialBeATargetFunction(ix, iy, iz))
//                   continue;
//             RetainedMonomials.insert(FMonomialN(ix,iy,iz));
//          }
//       }
//    }
//    FMonomialList
//       r;
//    r.reserve(RetainedMonomials.size());
//    for (FMonomialN const &m : RetainedMonomials)
//       r.push_back(FMonomial(m[0], m[1], m[2]));
//
// //    iPrintLevel = 2;
//    PrintMonomialList(std::cout, r, iPrintLevel); // note: this doesn't do anything unless iPrintLevel >= 2;
//    return r;
// }



int _TestPermFlections()
{
   FGridSearchOptions
      Options("{}"); // <-- without that MakeSymmetryUniqueMonomialListSimple won't do screening
//    FIcosahedralGroup PointGroup;
//    FOctahedralGroup PointGroup(FOctahedralGroup::GROUP_Full_Oh);
//    P
//    FOctahedralGroup PointGroup;
//    FMonomialSymmetry MonSym(PointGroup.GetOps());
   for (char const *pGs : {"T", "Th", "Td", "O","Oh","I", "Ih"}) {
//    for (std::string const &Gs : {"I"}) {
//    FOctahedralGroup PointGroup;
//    FMonomialSymmetry MonSym(PointGroup.GetOps());
      FPointGroupPtr pPointGroup = MakePointGroup(pGs);
//       FMonomialSymmetry MonomialSymmetry(pPointGroup->GetOps(), AIG_SPACE_AXES, AIG_SYMOP_MAX_PHASE_ANGLES, pPointGroup->Name(NAMETYPE_Symbol), FPrintLevel::More);
//       FMonomialSymmetryCptr
//          pMonomialSymmetry(new FMonomialSymmetry(pPointGroup->GetOps(), AIG_SPACE_AXES, AIG_SYMOP_MAX_PHASE_ANGLES, pPointGroup->Name(NAMETYPE_Symbol), FPrintLevel::More));
//       size_t lmin = 8, lmax=lmin+1;
//       FPrintLevel iPrintLevel(3);
//       io.Write("\n\n░░░ Group {} target function list ░░░", Gs);
//       io.Write("\n--- FPointGroup variant:\n", Gs);
//       FMonomialList
//          t1 = MakeSymmetryUniqueMonomialListSimple(lmin, lmax, *pPointGroup, iPrintLevel, &Options);
//          // ^-- that's what I c/p'd the FMonomialSymmetry variant from.
//       io.Write("\n--- FMonomialSymmetry variant:\n", Gs);
//       FMonomialList
//          t2 = MakeSymmetryUniqueMonomialList(lmin, lmax, &*pMonomialSymmetry, iPrintLevel);

   }


   return 18;
}


int SymFnTest()
{
//    return _TestPermFlections();
//    return _TestGroupAveraging();
//    return OverlapOverflowTest();
   return 0;
//    return TestSymmetrizedSpaceSizes();
   return TestMemoryAndTargetFnScreening();
   FIcosahedralGroup
      PointGroup;
   FGridSearchOptions
      Options("{}");
//    FRotationMatrixList const
//       &SymmetryOps = PointGroup.GetOps();
//    for (int i = 0; i < 1000; ++ i) {
//    for (int i = 0; i < 1; ++ i) {
// //       FMonomialList ml = MakeSymmetryUniqueMonomialList(20, 21, SymmetryOps, 1);
//       FMonomialList ml = MakeSymmetryUniqueMonomialList(30, 31, SymmetryOps, 1);
// //       FMonomialList ml = MakeSymmetryUniqueMonomialList(40, 41, SymmetryOps, 1);
//    }
   for (int l = 0; l <= 30; l += 5) {
      FMonomialList ml = MakeSymmetryUniqueMonomialList(l, l+1, PointGroup, 0, &Options);
   }
   return 1;
};



