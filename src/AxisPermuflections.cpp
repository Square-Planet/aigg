#include "AxisPermuflections.h"

#include <type_traits>
#include <array>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <sstream>

#include <functional> // for std::function and std::hash
#include <algorithm> // for equal_range and stable_sort

#include "CxSequenceIo.h"
#include "GenericGroups.h"


using groups::atoa;

// Compute for integer `a` and positive integer N, compute `a mod N`,
// with non-negative result in {0,1, …, (N-1)} (as in ℤ_N ring).
template<class index_t>
inline index_t _integer_modulo(index_t a, index_t N) {
   assert(N >= 1);
   index_t r = a % N;
   return (r >= 0)? r : (r + N);
   // ^-- warning: please do not change that to allow different integer types
   //     for `a` and `N`! Otherwise the C++ integer arithmetic promotion rules
   //     ultimately decide what gets calculated, and for several combinations
   //     of integer types that would not be "mod(a,N)" with these formulas.
}

// hm… ugly.
static phase_scaled_int_t _get_phase_factor_2(int iPhase) {
//    assert(index_valid_q(iPhase, 2));
   return (iPhase % 2 == 0)? +1 : -1;
}

template<class scalar_t>
static void _inplace_apply_phase_factor_2(scalar_t *inout, int iPhase) {
//    assert(index_valid_q(iPhase, 2));
   if (iPhase % 2 != 0)
      *inout = -*inout;
}

template<class scalar_t>
static scalar_t _multiply_arg_by_phase_factor_2(scalar_t in, int iPhase) {
   return (iPhase % 2 != 0)? -in : in;
}

FGetPhaseFactorFn *_PhaseFactorFn_Get(size_t nPhases) {
   assert(nPhases == 2); return &_get_phase_factor_2;
}
template<class scalar_t> TInplaceApplyPhaseFactorFn<scalar_t> *_PhaseFactorFn_InplaceApply(size_t nPhases) {
   assert(nPhases == 2); return &_inplace_apply_phase_factor_2<scalar_t>;
}
template<class scalar_t> TMultiplyArgByPhaseFactorFn<scalar_t> *_PhaseFactorFn_MultiplyArg(size_t nPhases) {
   assert(nPhases == 2); return &_multiply_arg_by_phase_factor_2<scalar_t>;
}

phase_scaled_int_t _phase_factor(phase_index_t iPhase, size_t nPhases) {
   return _PhaseFactorFn_Get(nPhases)(iPhase);
}

#ifdef INCLUDE_ABANDONED
// template<class scalar_t>
// TPhaseFactorFns<scalar_t> _GetPhaseFactorFns(size_t nPhases)
// {
//    assert(nPhases == 2);
//    return {&_get_phase_factor_2, &_inplace_apply_phase_factor_2<scalar_t>, &_multiply_arg_by_phase_factor_2<scalar_t>};
// }

// FGetPhaseFactorFn *_PhaseFactorFn_Get(size_t nPhases) {
//    assert(nPhases == 2);
//    return &_get_phase_factor_2;
// }
//
// template<class scalar_t> TInplaceApplyPhaseFactorFn<scalar_t> *_PhaseFactorFn_InplaceApply(size_t nPhases) {
//    assert(nPhases == 2);
//    return &_inplace_apply_phase_factor_2<scalar_t>;
// }
//
// template<class scalar_t> TMultiplyArgByPhaseFactorFn<scalar_t> *_FindFn_MultiplyArgByPhaseFactor(size_t nPhases) {
//    assert(nPhases == 2);
//    return &_multiply_arg_by_phase_factor_2<scalar_t>;
// }
#endif // INCLUDE_ABANDONED


#ifdef INCLUDE_ABANDONED
// phase_scaled_int_t _phase_factor(phase_index_t iPhase, size_t nPhases) {
//    assert(nPhases == 2);
//    return (iPhase % 2 == 0)? +1 : -1;
// }
// struct FPhaseFactorFns {
//    FApplyPhaseFactorFn scale_by_phase_factor;
//    FApplyPhaseFactorFn get_phase_factor;
// };
// FPhaseFactorFns _phase_mod_fns(size_t nPhases);
#endif // INCLUDE_ABANDONED


std::ostream &WriteSign(std::ostream &out, phase_scaled_int_t Sign) {
   if (Sign == 1)
      out << "+";
   else if (Sign == -1)
      out << "-";
   else
      out << "#SIGN/ERR (i=" << Sign << ")";
   return out;
}


std::ostream &WritePhase(std::ostream &out, phase_index_t iPhase, size_t nPhases) {
   assert(iPhase >= 0 && iPhase < nPhases);
   if (nPhases != 2) throw std::runtime_error("WritePhase(i,N) encountered unexpected N.");
   WriteSign(out, _phase_factor(iPhase, nPhases));
   return out;
}


std::ostream &operator << (std::ostream &out, FAxisPermuflection::FAxisEntry const &pfe)
{
//    WriteSign(out, pfe.Sign);
   WritePhase(out, pfe.iPhase, (pfe.nPhasesIfKnown != 0)? pfe.nPhasesIfKnown : 2);
   out << (pfe.iTargetAxis + 1);
   return out;
}

using ct::sequence_io::open_sep_close_t;
using ct::sequence_io::TSequenceWriter;

std::ostream &operator << (std::ostream &out, FAxisPermuflection const &pf) {
   return TSequenceWriter(open_sep_close_t{"(","/",")"}).WriteList(out, pf);
}

std::ostream &operator << (std::ostream &out, FAxisPermuflectionList const &L) {
   return TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(out, L);
}

std::ostream &operator << (std::ostream &out, FAxisPermuflectionSet const &L) {
   return TSequenceWriter(open_sep_close_t{"{",", ","}"}).WriteList(out, L);
}


// FAxisPermuflection::FPermOnly FAxisPermuflection::PermQ() const {
//    // well… if that isn't super clear and concise, I also don't know…
//    FPermOnly out;
//    std::transform(m_AxisOps.begin(), m_AxisOps.end(), out.begin(),
//       [&](FAxisEntry const &ae) { return ae.iTargetAxis; });
//    return out;
// }
//
//
// FAxisPermuflection::FSignsOnly FAxisPermuflection::SignsQ() const {
//    FSignsOnly  out;
//    std::transform(m_AxisOps.begin(), m_AxisOps.end(), out.begin(),
//       [&](FAxisEntry const &ae) { return ae.Sign; });
//    return out;
// }


std::ostream &operator << (std::ostream &out, FRotationMatrix const &m)
{
   out << "[";
   for (size_t iRow = 0; iRow < size_t(m.rows()); ++iRow) {
      if (iRow != 0) out << ", ";
      out << fmt::format("[{:2},{:2},{:2}]", m(iRow,0), m(iRow,1), m(iRow,2));
   }
   out << "]";
   return out;
}


std::optional<FAxisPermuflection> AsPermuflection(FRotationMatrix const &m, unsigned nAxes, unsigned nPhases)
{
//    io.Write("      {:<8} Try convert: {}", "", m);
   std::array<FAxisPermuflection::FAxisEntry, FAxisPermuflection::MaxAxes>
      PermfOps;
   assert(size_t(m.rows()) >= nAxes && size_t(m.cols()) >= nAxes);

   auto _phase_factor = _PhaseFactorFn_Get(nPhases);
   // the matrix `m` represents the mapping of the source Cartesian coordinates [x_0, x_1, x_2, …]
   // to the transformed new coordinates:
   //
   //    x_i' = \sum_{j} m_{ij} x_j
   //
   // For this mapping to be representable as Cartesian axis permuflection, for each
   // source index `j` there must be exactly one target index `i` for which
   // m_{ij} is non-zero, and this non-zero entry must be either +1 or -1.
   for (unsigned iSourceCoord = 0; iSourceCoord != nAxes; ++ iSourceCoord) {
      int
         iFoundPhase = -1;
      unsigned
         iFoundAxis = 0;
      for (unsigned iTargetCoord = 0; iTargetCoord != nAxes; ++ iTargetCoord) {
         FScalar
            mij = m(iTargetCoord, iSourceCoord);
         if (IsAlmostZero(mij))
            // source axis not mapped to this target axis.
            continue;
#ifdef INCLUDE_ABANDONED
//          int
//             iPhaseIndex(-1);
//          if (IsAlmostZero(mij))
//             // source axis not mapped to this target axis.
//             continue;
//          else if (IsAlmostEqual(mij, 1)) {
//             iPhaseIndex = 0;
//          } else if (IsAlmostEqual(mij, -1)) {
//             assert(FAxisPermuflection::nPhases % 2 == 0);
//             iPhaseIndex = FAxisPermuflection::nPhases / 2;
//          } else {
// //             io.Write("      {:<8} m[{},{}] = {}  -> not -1, or 1", "", iTargetCoord, iSourceCoord, double(mij));
//             // not one of {+1, 0, -1} --> not a permuflection.
//             break;
//          }
#endif // INCLUDE_ABANDONED
         // have we already found a m_{ij} which is a phase factor?
         if (iFoundPhase >= 0)
            // If yes, all the other entries in the column should vanish, otherwise
            // the column vector cannot have norm 1 in the L2 norm (at least provided
            // `m` is expressed in terms of orthonormal Cartesian axes)
            throw std::runtime_error("AsPermuflection(): supposedly-unitary matrix contains multiple |m_{ij}| = 1 entries in a column.");
         for (phase_index_t iTrialPhase = 0; iTrialPhase < nPhases; ++ iTrialPhase) {
            if (IsAlmostEqual(mij, _phase_factor(iTrialPhase))) {
               iFoundAxis = iTargetCoord;
               iFoundPhase = iTrialPhase;
               break;
            }
         }
         if (iFoundPhase < 0)
            // neither 0 (zero is `continue`'d above) nor one of the phase
            // factors we were looking for (iFoundPhase would be >= 0 now)
            // --> not a permuflection
            break;
      }
      if (iFoundPhase < 0)
         // not a permuflection. Either all entries in column vanish
         // or the column contains some components which are neither -1,0, nor 1.
         return {};
      else
         PermfOps[iSourceCoord] = FAxisPermuflection::FAxisEntry(iFoundPhase, iFoundAxis, nPhases);
   }
//    io.Write("      {:<8} Created Op: {}", "", PermfOps);
//    io.Write("      {:<8} Created Op: |{}| {} {} {}", "", PermfOps.size(), PermfOps[0], PermfOps[1], PermfOps[2]);
   return std::make_optional<FAxisPermuflection>(std::make_pair(PermfOps.begin(), PermfOps.begin() + nAxes));
}


// returns permuflection list of all symmetry operators in SymOps which can be
// represented as permuflections
FAxisPermuflectionList ExtractPermuflections(FRotationMatrixList const &SymOps, unsigned nAxes, unsigned nPhases)
{
   FAxisPermuflectionList
      out;
   out.reserve(SymOps.size());
   for (FRotationMatrix const &SymOp : SymOps) {
      auto MaybePermfOp = AsPermuflection(SymOp, nAxes, nPhases);
      if (bool(MaybePermfOp))
         out.emplace_back(std::move(MaybePermfOp.value()));
         // ^-- I know it doesn't actually have a move constructor, and for such
         // an all-data type it anyway wouldn't make sense. I'm just learning
         // and trying new C++ features \o/
   }
   return out;
}







inline size_t CombineHash1(size_t x, size_t y) {
   // pull both hash values through (different) basic PRNGs and add them up.
   // These are actually meant for 32bit integers, but I guess it still beats
   // just xoring x and (y + 1234567).
   y ^= y << 5; y ^= y >> 7; y ^= y << 22;   // 3 xor, 3 shift
   x = 314527869 * x + 1234567;              // 1 mul, 1 add
   return x + y;                             // 1 add
}
                                             // --> 9 ops
inline size_t CombineHash2(size_t x, size_t y) {
   return y ^ (2654435769 + x + (y >> 2) + (y << 6)); // 6 ops (1 xor, 3 add, 2 shift)
}


inline size_t CombineHashBasic(size_t x, size_t y) {
   return y ^ (314527869 * y + 1234567 + x); // 4 ops
}


// quick & dirty rounding conversion of floats to ints. just meant to work
// around the truncation vs round-to-nearest-int issue, not denorms, out of ranges,
// "round to even"/equal distribution stuff etc.
template<class int_t, class float_t>
int_t AsInt(float_t f) {
   if (std::is_signed<int_t>::value)
      return f < 0 ? int_t(f - .5) : int_t(f + .5);
   else
      return int_t(f + .5);
}


#ifdef INCLUDE_ABANDONED
// using FSortKey_RotationMatrix = std::array<int32_t, 1 + AIG_SPACE_AXES*AIG_SPACE_AXES>; // <-- std::array has lex compare built in
//
// template<class scalar_t>
// FSortKey_RotationMatrix MakeSortKey_RotationMatrix(TRotationMatrix<scalar_t> const &Rg)
// {
// //    using std::round; // <-- mostly about what should be exact zeros and ending up slightly above/below that.
//    FSortKey_RotationMatrix out;
//    typedef FSortKey_RotationMatrix::value_type signed_int_t;
//    typedef std::make_unsigned_t<signed_int_t> unsigned_int_t;
//    size_t nBits = (sizeof(signed_int_t) - 1)*CHAR_BIT; // <-- (…-1) is for sign.
//    // Rg should be unitary, so all elements would be smaller than |1| in magnitude.
//    // I'll still cut a few bits off. Remaining precision comes out 5e-7ish
//    nBits -= 3;
//    scalar_t fScale = (unsigned_int_t(1) << nBits);
// //    io.Write("data @ scale factor for hash:  nBits = {}  1/fScale = {:.2e}  fScale = {:.2e}", nBits, fScale, 1/fScale);
//    size_t HashAcc(0);
//    for (size_t i = 0; i < AIG_SPACE_AXES; ++ i)
//       for (size_t j = 0; j < AIG_SPACE_AXES; ++ j) {
//          scalar_t Scaled_Rg_ij = fScale * Rg(i,j);
//          signed_int_t IntKey_Rg_ij = AsInt<signed_int_t>(Scaled_Rg_ij));
//          out[1 + i + AIG_SPACE_AXES*j] = IntKey_Rg_ij;
//          size_t Hash_Rg_ij = std::hash(IntKey_Rg_ij);
//          HashAcc = CombineHashBasic(HashAcc, Hash_Rg_ij);
//       }
//    out[0] = HashAcc; // <-- put this first. should make the lex compare decisive with first element in most cases.
//    return out;
// }
#endif // INCLUDE_ABANDONED

// let's make this simpler… keep *only* the hash value, and pretend that is
// always going to uniquely identify the elements. Which, in practice, is nearly
// guaranteed. Looks ugly and unsafe. However: *ANY* sort of lookup structure
// for float-containing objects is inherently dangerous, because there is always
// a possibility of just landing *exactly* on the boundary between two bins. For
// actual discrete key-based approaches it simply cannot be done, on a
// conceptual level, I think. So… in that case there is little reason not to
// cut corners and simply pretend that the hashes themselves are enough.
using FSortKey_RotationMatrix = size_t;

template<class scalar_t>
struct FSortKeyGen_RotationMatrix
{
   typedef float
      truncating_scalar_t;
   std::hash<truncating_scalar_t>
      m_HashFn;
   size_t operator() (TRotationMatrix<scalar_t> const &Rg) const
   {
      size_t HashAcc1(0), HashAcc2(0);
      for (size_t i = 0; i < size_t(Rg.rows()); ++ i)
         for (size_t j = 0; j < size_t(Rg.cols()); ++ j) {
            scalar_t f = Rg(i,j);
            // shift float baseline by adding a constant. This works around issues
            // with very small non-meaningful number distinctions introduced by
            // rounding (e.g., positive and negative zero, 1e-17 vs 1e+17, etc),
            // assuming that the matrix elements in Rg have a order of magnitude
            // comparable with unity.
            f += 17.827498237;

            // truncate f to 32bit floats in order to cut off numerical roundoff
            // errors. The idea is that two floats a,b differing only by a few eps
            // of numerical roundoff errors will most likely still convert to the
            // same 32bit float value (of course there is no guarantee---sorting
            // floats into discrete brackets is intrinsically unsafe).
            static_assert(sizeof(truncating_scalar_t) < sizeof(scalar_t));
            truncating_scalar_t f_truncated(f);

            // now just feed it through std::hash's default implementation
            // and combine with the previous member hashes.
            size_t f_hash = m_HashFn(f_truncated);
            HashAcc1 = CombineHash1(HashAcc1, f_hash);
            HashAcc2 = CombineHash2(HashAcc2, f_hash);
         }
      return HashAcc1 ^ HashAcc2;
   }
};


// this one is not shady. It should always be exact, as long as (2d)!, with d
// the space dimension, fits into a size_t.
using FSortKey_AxisPermuflection = size_t;

bool constexpr _CheckIfSizeTSufficientForUniqueEncoding(size_t nAxes, size_t nPhases) {
   size_t iFactorAcc = 1;
   for (size_t iAxis = 0; iAxis < nAxes; ++ iAxis) {
      size_t iPrevFactorAcc = iFactorAcc;
      iFactorAcc *= nAxes * nPhases;
      if (iFactorAcc < iPrevFactorAcc)
         return false;
//             throw std::runtime_error("FSortKeyGen_AxisPermuflection::eval(): integer overflow in unique sort key generation. (2d)! does not fit into a size_t.");
   }
   return true;
}

struct FSortKeyGen_AxisPermuflection {
   size_t operator() (FAxisPermuflection const &Pm) const {
//       size_t HashAcc(0), iFactorAcc(1), nChoicesForOneAxis = 2 * Pm.size();
      size_t HashAccPerm(0), iFactorAccPerm(1);
      size_t HashAccSign(0), iFactorAccSign(1);
//       size_t constexpr nAxes = Pm.size(), nPhases = Pm::nPhases;
      size_t constexpr MaxAxes = FAxisPermuflection::MaxAxes, MaxPhases = FAxisPermuflection::MaxPhases;
      static_assert(_CheckIfSizeTSufficientForUniqueEncoding(MaxAxes, MaxPhases));

      // sort equal perms of different signs next to each other,
      // with neutral signs coming first.
      for (size_t iAxis = Pm.nAxes() - 1; iAxis < Pm.nAxes(); -- iAxis) {
         assert(groups::index_valid_q(Pm.TargetAxisQ(iAxis), Pm.nAxes()));
//          assert(Pm.SignQ(iAxis) == 1 || Pm.SignQ(iAxis) == -1);
//          size_t iSignChoice = (Pm.SignQ(iAxis) == -1)? 0 : 1;
         assert(groups::index_valid_q(Pm.PhaseIndexQ(iAxis), Pm.nPhases()));

         HashAccPerm += iFactorAccPerm * Pm.TargetAxisQ(iAxis);
         HashAccSign += iFactorAccSign * Pm.PhaseIndexQ(iAxis);

         iFactorAccPerm *= MaxAxes;
         iFactorAccSign *= MaxPhases;
      }
      return HashAccSign + iFactorAccSign * HashAccPerm;
   }
 };






#ifdef INCLUDE_ABANDONED
// template<class T>
// struct element_search_helper_base {
//    // by default use group elements themselves as key type.
//    typedef T key_t;
// };
//
//
// namespace element_search {
//    // by default use group elements themselves as key type.
//    template<class T> struct key_type { typedef T type; }
//    template<class T> using key_type_t = typename key_type<T>::type;
//
//    // rotation matrices
//    template<class scalar_t>
//    struct key_type<TRotationMatrix<scalar_t> > { typedef FSortKey_RotationMatrix type;  }
//    template<class scalar_t>
//    inline FSortKey_RotationMatrix MakeKey(TRotationMatrix<scalar_t> const &Rg) { return MakeSortKey_RotationMatrix(Rg); }
//
//
//    inline FAxisPermuflection &&MakeKey(FAxisPermuflection &&Pr) { return std::forward(Pr); }
//
//    template<class TMatrix>
//    struct FUniqueMatrixSet {
//       typedef std::map<
//    }
//
// }
#endif // INCLUDE_ABANDONED


char const *FormatSubgroupId(FSubGroupId id) {
   switch(id) {
      case SYMGROUP_AxisPerms_Cyclic_All:
         return "ax:π|cyclic";
      case SYMGROUP_AxisPerms_Cyclic_Even:
         return "ax:π|cyclic,even";
      case SYMGROUP_AxisPerms_Even:
         return "ax:π|even";
      case SYMGROUP_AxisPerms_All:
         return "ax:π|all";
      case SYMGROUP_AxesReflect_AnyTwo:
         return "ax:σ|even";
      case SYMGROUP_AxesReflect_Any:
         return "ax:σ|all";
      case SYMGROUP_AxesRotate4_AnyTwo:
         return "ax:c4|all";
      default:
         return "?[?!]";
   }
};

static long _LameFactorial(long n, long acc = 1) {
   if (n <= 0) return acc;
   if (std::numeric_limits<long>::max()/(2*n) < acc)
      throw std::overflow_error(fmt::format("_LameFactorial({}): stopped at integer near-overflow", n));
   return _LameFactorial(n-1, n*acc);
}

size_t SubgroupOrderExpected(FSubGroupId id, size_t nAxes) {
   size_t N = nAxes;
   size_t const NotThere = size_t(-1);
   switch(id) {
      case SYMGROUP_AxisPerms_Cyclic_All:
         return N;
      case SYMGROUP_AxisPerms_Cyclic_Even:
         if (N % 2 == 0)
            return N/2;
         return N; // for odd N, all cyclic permutations are even.
      case SYMGROUP_AxisPerms_Even:
         return (nAxes < 2)? 1 : (_LameFactorial(N)/2);
      case SYMGROUP_AxisPerms_All:
         return size_t(_LameFactorial(N));
      case SYMGROUP_AxesReflect_AnyTwo:
         // return 2^N/2 … this covers all reflection operators with
         // determinant 1, doesn't it?
         return (nAxes < 2)? 0 : size_t(1ul << size_t(N-1));
      case SYMGROUP_AxesReflect_Any:
         return size_t(1ul << N); // <-- 2^N
      case SYMGROUP_AxesRotate4_AnyTwo:
//          return "ax:c4|all";
         if (N < 2) return 0;
         if (N == 2) return 4; // it's only the cycle of the one c4
         if (N == 3) return 24; // in R^3, two orthogonal C4s generate the entirety of O.
         // general: no idea.
         return NotThere;
      default:
         return NotThere;
   }
}


// char const *FormatSubgroupId(FSubGroupId id) {
//    switch(id) {
//       case SYMGROUP_AxisPerms_Cyclic_All:
//          return "{x[░] ⇒ x[(░+k) mod N]; k∈ℤ}";
//       case SYMGROUP_AxisPerms_Cyclic_Even:
//          return "{x[░] ⇒ x[(░+k) mod N]; k∈ℤ, (k+N)∈2ℤ}";
//       case SYMGROUP_AxisPerms_Even:
//          return "{x[░] ⇒ x[π₊(░)]; π₊∈A(N)}"; // alternating group
// //          return "{π₊[x₀,x₁,…,xₙ]; π₊∈A(N)}"; // alternating group
//       case SYMGROUP_AxisPerms_All:
// //          return "{π[x₀,x₁,…,xₙ]; π∈S(N)}";
// //          return "π ∈ S_{N}";
//          return "{x[░] ⇒ x[π(░)]; π∈S(N)}"; // symmetric group
//       case SYMGROUP_AxesReflect_AnyTwo:
//          return "{x[░] ⇒ ξ[░]⋅x[░]; |ξₖ|=1, ξ₀⋅ξ₁⋅…⋅ξₙ=1}";
// //          return "{x[░] ⇒ ξ[░]⋅x[░]; |ξₖ|=1, ∏{ξₖ}=1}"
// //          return "⟪(x[i],x[j]) ⇒ (−x[i],−x[j]); j∈[N],i∈[j]⟫"
// //          return "{σ_K; K⊂{xₖ}, |K|∈2ℕ₀}";
// //          return "{σ[K]; K⊂[N], |K|∈2ℕ₀}";
//       case SYMGROUP_AxesReflect_Any:
//          return "{x[░] ⇒ ξ[░]⋅x[░]; |ξₖ|=1}";
// //          return "{σ[K]; K⊂[N]}";
//       case SYMGROUP_AxesRotate4_AnyTwo:
//          return "⟪x[░], x[▒] ⇒ −x[▒], x[░]⟫";
//       default:
//          return "?[?!]";
//       // [N] := {0,1,…,(N-1)}
//       // unicode subscripts as of 2020: 0123456789()+-=aeoxhklmnpst, and ₔ via E (real subscripts, not phonetic markers)
//    }
// };


bool FAxisPermuflectionSubgroup::ContainsAlsoCoveringSubgroupQ(FSubGroupId H) const
{
   if (!this->ContainsAsSubgroupQ(H))
      return false;
   enum TestMode{IncludesAny, IncludesAll};
   auto SubTest1 = [&](FSubGroupId ImpliedH, TestMode tm, std::vector<FSubGroupId> &&SuperGroupIds) -> bool {
      if (tm == IncludesAll) {
#ifndef _DEBUG
         if (H != ImpliedH) return false;
#endif // _DEBUG
         bool AllSuperContained = true;
         for (auto &&F : SuperGroupIds)
            if (!this->ContainsAsSubgroupQ(F))
               AllSuperContained = false;
         if (AllSuperContained)
            cx_assert_rt(this->ContainsAsSubgroupQ(ImpliedH));
         if (H != ImpliedH) return false;
         return AllSuperContained;
      } else {
         assert(tm == IncludesAny);
#ifdef _DEBUG
         // also check reverse direction: if one of the supergroups, is contained
         // as subgroup in *this, then also the subgroup of the supergroup must be
         // contained.
         for (auto &&F : SuperGroupIds)
            if (this->ContainsAsSubgroupQ(F))
               cx_assert_rt(this->ContainsAsSubgroupQ(ImpliedH));
#endif // _DEBUG
         if (H != ImpliedH) return false;
         for (auto &&F : SuperGroupIds)
            if (this->ContainsAsSubgroupQ(F))
               return true;
         return false;
      }
   };

   if (SubTest1(SYMGROUP_AxisPerms_Cyclic_Even, IncludesAny, {SYMGROUP_AxisPerms_Cyclic_All, SYMGROUP_AxisPerms_Even}))
      return true;
   if (SubTest1(SYMGROUP_AxisPerms_Cyclic_All, IncludesAny, {SYMGROUP_AxisPerms_All}))
      return true;
   if (this->nAxes() % 2 == 1)
      // for uneven space dimension, all cyclic permutations are even. Therefore
      // in this case A_N is a supergroup of the general cyclic perms.
      if (SubTest1(SYMGROUP_AxisPerms_Cyclic_All, IncludesAny, {SYMGROUP_AxisPerms_Even}))
         return true;
   if (SubTest1(SYMGROUP_AxisPerms_Even, IncludesAny, {SYMGROUP_AxisPerms_All}))
      return true;
   if (SubTest1(SYMGROUP_AxesReflect_AnyTwo, IncludesAny, {SYMGROUP_AxesReflect_Any}))
      return true;
   if (SubTest1(SYMGROUP_AxesRotate4_AnyTwo, IncludesAll, {SYMGROUP_AxesReflect_Any, SYMGROUP_AxisPerms_All}))
      // π_{ij}, σ_{i}, and σ_j together cover c4_{ij}. If we have all of the
      // former, we also have all c4s.
      return true;
   return false;
}


void FAxisPermuflectionSubgroup::_PrintSpecialSubgroups(std::string const &Caption) {
   io.WriteLine();
   std::vector<FSubGroupId> SubGroupIds;
   for (size_t iSubGroupId = 0; iSubGroupId < SYMGROUP_Count; ++ iSubGroupId) {
      FSubGroupId SubGroupId = static_cast<FSubGroupId>(iSubGroupId);
//       if (m_ContainsAsSubgroup[iSubGroupId])
      if (ContainsAsSubgroupQ(SubGroupId) && !ContainsAlsoCoveringSubgroupQ(SubGroupId))
//          SubGroupIds.push_back(FormatSubgroupId(SubGroupId));
         SubGroupIds.push_back(SubGroupId);
   }
   std::stringstream str;
//    str << "[";
//    str << "(";
   for (auto &&[iSubGroup, SubGroupId] : enumerate(SubGroupIds)) {
      if (iSubGroup != 0) str << "  ";
      str << FormatSubgroupId(SubGroupId);
      size_t n = SubgroupOrderExpected(SubGroupId, nAxes());
      if (n != size_t(-1))
         str << "(" << n << ")";
   }
//    str << ")";
//    str << "]";
//    using ct::sequence_io::open_sep_close_t;
//    ct::sequence_io::TSequenceWriter(open_sep_close_t{"[",", ","]"}).WriteList(str, SubGroupIds, &ct::sequence_io::DefaultWriteElement_Unquoted<std::string>);

//    io.Write(" ■ largest identified apf subgroups of {}: {}", m_ParentGroupName, str.str());
//    io.Write(" ■ non-covered apf subgroups of {}: {}", m_ParentGroupName, str.str());
   io.Write(" ■ {:<32}{}", fmt::format("{}({}) core axpr subgroups:", m_ParentGroupName, m_ParentGroupOrder), str.str());
   // or Axpp for axis permutation+phase-shift instead of reflection?
};


FAxisPermuflectionSubgroup::FAxisPermuflectionSubgroup(FGroupElementList_Parent const &ParentGroupElements,
      unsigned nAxes_, unsigned nPhases_, std::string const &ParentGroupName, FPrintLevel PrintLevel)
   : m_PrintLevel(PrintLevel), m_ParentGroupName(ParentGroupName), m_ParentGroupOrder(ParentGroupElements.size()),
     m_nAxes(nAxes_), m_nPhases(nPhases_)
{
   m_ContainsAsSubgroup.fill(false);

   // extract symmetry ops in ParentGroupElements which can be represented as permuflections
   m_Elements = ExtractPermuflections(ParentGroupElements, m_nAxes, m_nPhases);
   {
      // collect them into a set (to allow easy lookup of existence of elements)
      FAxisPermuflectionSet OpSet;
      for (FAxisPermuflection const &pf : m_Elements)
         OpSet.insert(pf);
      if (OpSet.size() != m_Elements.size())
         throw std::logic_error("FAxisPermuflectionSubgroup: processing and conversion of group elements"
            "revealed inconsistent initial group structure (permuflection elements got lost during list"
            "-> set conversion). That is not supposed to happen.");
      _SearchSpecialSubgroups(OpSet);
      if (m_PrintLevel.MoreQ())
         _PrintSpecialSubgroups("final subgroup list");

   }
// #ifdef _DEBUG
#if 0
   {
      typedef TRotationMatrix<double> FGroupElement;
      std::vector<FGroupElement> ParentGroupElements1d; // reduced from big multiprecision scalars to double scalars. should be enough.
      for (FRotationMatrix const &Rg : ParentGroupElements)
         ParentGroupElements1d.emplace_back(Rg);

      FGroupElement id(FGroupElement::Identity());
      auto prod = [&] (FGroupElement const &e, FGroupElement const &f) { return (e * f).eval(); };
      auto inverse = [&] (FGroupElement const &e) { return e.inverse(); };
//       auto element_key = [&] (auto const &e) { return e; }; // they already have < operators.
//       auto element_key = [&] (auto e) { return std::forward(e); }; // they already have < operators.
//       , element_key
      groups::_VerifyGroupStructureT(id, ParentGroupElements1d, prod, inverse, FSortKeyGen_RotationMatrix<double>());
   }
   {
      typedef FAxisPermuflection FGroupElement;
      FGroupElement id({{+1,+2,+3}});
      auto prod = [&] (FGroupElement const &e, FGroupElement const &f) { return e * f; };
      auto inverse = [&] (FGroupElement const &e) { return ~e; };
      auto element_key = [&] (auto const &e) { return e; }; // they already have < operators.
//       auto element_key = [&] (auto e) { return std::forward(e); }; // they already have < operators.
//       , element_key
      groups::_VerifyGroupStructureT(id, m_Elements, prod, inverse, element_key);
   }

#endif // _DEBUG
// ^--- wow.. super slow in debug mode. Hooray for expression templates…

   _MakeCosetReps(m_LeftCosetReps, COSET_Left, ParentGroupElements);
   _MakeCosetReps(m_RightCosetReps, COSET_Right, ParentGroupElements);

}


// // returns whether all elements in `OpsToQuery` are contained in `AllOps`
// static bool _AllPresent(FAxisPermuflectionList const &OpsToQuery, FAxisPermuflectionSet const &AllOps) {
//    for (FAxisPermuflection const &op : OpsToQuery)
//       if (AllOps.find(op) != AllOps.end())
//          return false;
//    return true;
// }


void FAxisPermuflectionSubgroup::_SearchSpecialSubgroups(FAxisPermuflectionSet const &AllOps)
{
   unsigned const
      nAxes = this->nAxes(),
      nPhases = this->nPhases();

   // returns whether all elements in `OpsToQuery` are contained in `AllOps`
   auto _AllPresent = [&AllOps](FAxisPermuflectionList const &OpsToQuery) -> bool {
      for (FAxisPermuflection const &op : OpsToQuery)
         if (AllOps.find(op) == AllOps.end())
            return false;
      return true;
   };
   // makes a simultaneous-axis-reflection symmetry op, with the combination of
   // axes to be reflected described by the bitfield amf := AxisMirrorFlags.
   // Concretely, the returned operator maps:
   //
   //     (x_0, x_1, x_2, …) ⟼ (x_0', x_1', x_2', …)
   //
   // with
   //
   //     x_i' := x_i * (bit-test(amf,i)? (−1) : (+1))
   //
   auto _AxisRefl = [&nAxes,&nPhases](unsigned long AxisMirrorFlags) -> FAxisPermuflection {
      return FAxisPermuflection::Reflections(AxisMirrorFlags, nAxes, nPhases);
   };
   // cyclic(n,K) returns a permutation operator which applies a cyclic
   // permutation by `n := iCyclicShift` elements to the first `K :=
   // nAxesToPermute` axis coordinates, and acts as identity on all later
   // coordinates `x_i` with `i ≥ K`.
   // E.g., cyclic(2,K) can be viewed as splicing out the first two elements
   // (x_0, x_1) from the original sequence, and re-inserting them into the
   // sequence after x_{K-1}:
   //
   //    cyclic(2,K): (x_0, x_1, x_2, x_3, …, x_{K-1},   x_K, x_{K+1}, …, x_N)
   //              ⟼ (x_2, x_3, …, x_{K-1}, x_0, x_1,   x_K, x_{K+1}, …, x_N)
   //                 ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄   ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
   //                 cyclic perm by 2 on first K elems  | other N-K elems left alone
   auto _cyclic = [&nAxes,&nPhases](int iCyclicShift, unsigned nAxesToPermute) -> FAxisPermuflection {
      return FAxisPermuflection::CyclicPerm(iCyclicShift, nAxesToPermute, nAxes, nPhases);
   };
//    auto _cyclic = [&nAxes,&nPhases](int iCyclicShift, unsigned nAxesToPermute) -> FAxisPermuflection {
//       auto cyclp = FAxisPermuflection::CyclicPerm(iCyclicShift, nAxesToPermute, nAxes, nPhases);
//       io.Write("     _cyclic({},{},[{},{}]) = {}", iCyclicShift, nAxesToPermute, nAxes, nPhases, cyclp);
//       return cyclp;
//    };

   assert(m_ContainsAsSubgroup.size() >= SYMGROUP_Count);
   m_ContainsAsSubgroup.fill(false);

   // let's go with the important ones first: cyclic all/even, then and general
   // axis permutations, then any-one/any-two reflect.
   if (nAxes < 2) {
      // with 1 axis, there is only the identity permutation, which we have (A
      // being a group etc). And I don't think 0 axes is actually a thing.
      assert(nAxes != 0);
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_All] = true;
   } else {
      assert(nAxes >= 2);
//       FAxisPermuflection
//          pCyclic1 = _cyclic(1, nAxes),
//          p01 = _cyclic(1, 2);
      // have full cyclic shift? (=cyclic shift of "first" `nAxes` axis
      // coords by 1 element)
      if (_AllPresent({_cyclic(1, nAxes)})) {
         m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_All] = true;
         // have swap of first two elements? (= cyclic shift of first 2 axis
         // coords by 1 element, other coords remaining as they are)
         if (_AllPresent({_cyclic(1, 2)})) {
            // one cyclic permutation + one adjacent-element-swap are sufficient
            // to generate the entire symmetric group S_N as their span, for any
            // N >= 2. So that's all we need to check for that.
            m_ContainsAsSubgroup[SYMGROUP_AxisPerms_All] = true;
         }
      }
   }
   // fill out some implied subgroups
   if (m_ContainsAsSubgroup[SYMGROUP_AxisPerms_All]) {
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even] = true;
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_All] = true;
   }
   if (m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_All])
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even] = true;

   if (nAxes < 3) {
      // for 1 or 2 axes the only even-sign permutation is identity, which we have.
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even] = true;
      m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even] = true;
   }
   // even if we do not have *all* cyclic permutations, in principle we may
   // still have all even-signed ones. However, if N is odd, then *all* cyclic
   // permutations have even sign (i.e., SYMGROUP_AxisPerms_Cyclic_Even and
   // SYMGROUP_AxisPerms_Cyclic_All are identical). So an extra check is only
   // warranted for even N.
   if (nAxes % 2 == 0 && !m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even]) {
      assert(nAxes >= 3);
      // have cyclic shift of `nAxes` axis coordinates by 2 elements?
      if (_AllPresent({_cyclic(2,nAxes)}))
         m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even] = true;
   }
   // if we have even cyclic axis coordinate permutations, have a look
   // if we have general even coordinate permutations, too.
   if ( m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even] &&
       !m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even])
   {
      assert(nAxes >= 3);
      // the even-signed permutation cyclic(1, nAxis)^k (with k=1 for odd nAxis,
      // or k=2 for even nAxis) combined with the always-even-signed cyclic(1,3)
      // cyclic permutation of the first three elements is sufficient to
      // generate the full subgroup S_{N,+} ⊂ S_N of all even-signed
      // permutations (with N!/2 elements for N >= 2)
      if (_AllPresent({_cyclic(1,3)}))
         m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even] = true;
   }

//    bool HaveCyclicPerm = m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even];
//    assert(!"^--- wait, what?");
   bool
      HavePerm_AllEven = m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Even],
      HavePerm_All = m_ContainsAsSubgroup[SYMGROUP_AxisPerms_All],
      HavePerm_CyclicEven = m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_Even],
      HavePerm_CyclicAll = m_ContainsAsSubgroup[SYMGROUP_AxisPerms_Cyclic_All];
//    assert(!"^--- wait, what?");
   // up next: the checks for 'every single-axis symmetry X_i' and/or 'every
   // axis-pair symmetry X_{i,j}' (with i ≠ j). If the axis-pair operator
   // X_{i,j} can be written as X_{i,j} = X_i X_j, then the property
   // 'has symmetry X_i for every axis `i`' implies that also all
   // axis-pair symmetries are there. But the reverse generally does not hold.
   //
   // E.g., in ℝ³:
   // - let  σ_x: (x,y,z) ⟼ (−x,y,z) etc denote the individual
   //   coordinate reflections, and σ_{xy} := σ_x σ_y = [(x,y,z) ⟼ (−x,−y,z)],
   //   the simultaneous reflection of two-axes, etc.
   // - If σ_x, σ_y, σ_z ∈ G, then clearly also σ_{xy}, σ_{xz}, σ_{yz} ∈ G,
   //   because the two-axis reflections can be written as products of
   //   one-axis reflections.
   // - However, σ_{xy}, σ_{xz}, σ_{yz} ∈ G (by itself) does not say anything about
   //   whether σ_x, σ_y, σ_z ∈ G — note that σ_{xy},σ_{xz},σ_{yz} all are operators
   //   with determinant (+1), while σ_x,σ_y,σ_z have determinant (−1), so there is
   //   no way to obtain the single-axis reflections from any combination of
   //   two-axis reflections alone.

   // first determine limits of `i` and `j` indices we may need to iterate over
   // to get all unique X_i and/or X_i X_j pairs, given what we know about the
   // axis permutation symmetries already.
   size_t
      // first index: iterate over i ∈ {0,1,…,iAxisEnd-1}
      iAxisEnd = nAxes;
   if (HavePerm_CyclicEven)
      // cyclic-even permutations allow us to move any element to either index 0
      // or 1, so checking i ∈ {0,1} is sufficient to get all unique `i`.
      iAxisEnd = 2;
   if (HavePerm_CyclicAll)
      // …and with full-cyclic, we can move every element to index 0. So
      // searching i ∈ {0} is enough.
      iAxisEnd = 1;
   // clamp iAxisEnd to valid range (end exclusive)
   iAxisEnd = std::min(size_t(nAxes), size_t(iAxisEnd));

   // now for the second index `j` in (i,j) pairs.
   size_t
      // `jAxiscount` defines the max number of
      // elements after `i` to iterate over:
      // j ∈ {i+1, i+2, …, min(nAxes-1,i+jAxisCount)}
      jAxisCount = nAxes-1;
   // The cyclic permutation symmetries by themselves only allow fixing one
   // index (i). But even or general generic permutation symmetries allow
   // restricting `j`, too (and, if needed, `k`, and `l`, …).
   if (HavePerm_AllEven)
      // with general even permutation symmetries, we can move any
      // element `j` at most two indices beyond `i`.
      jAxisCount = 2;
   if (HavePerm_All)
      jAxisCount = 1;

   // Let's get started. Check reflections.
   {  // First check individual axis reflections. If present, they imply the other
      // reflection subgroups.
      FAxisPermuflectionList ToTest;
      // make list of all permutation-unique single-axis reflections
      for (size_t iAxis = 0; iAxis < iAxisEnd; ++ iAxis)
         ToTest.emplace_back(_AxisRefl(1ul << iAxis));
      if (_AllPresent(ToTest)) {
         // all present. That means that literally all combinations of axis
         // reflections are group elements, because all can be written as
         // products of single-axis reflections. Nothing else to test
         // reflection-wise.
         m_ContainsAsSubgroup[SYMGROUP_AxesReflect_Any] = true;
         m_ContainsAsSubgroup[SYMGROUP_AxesReflect_AnyTwo] = true;
      }
   }
   if (!m_ContainsAsSubgroup[SYMGROUP_AxesReflect_AnyTwo]) {
      // if not already implied, check simultaneous reflections of two axes.
      FAxisPermuflectionList ToTest;
      // make list of all permutation-unique combinations of two simultaneous axis reflections
      for (size_t iAxis = 0; iAxis < iAxisEnd; ++ iAxis)
         for (size_t jAxis = iAxis+1; jAxis < std::min(size_t((iAxis+1) + jAxisCount), size_t(nAxes)); ++ jAxis)
            ToTest.emplace_back(_AxisRefl((1ul << iAxis) | (1ul << jAxis)));
      if (_AllPresent(ToTest)) {
         m_ContainsAsSubgroup[SYMGROUP_AxesReflect_AnyTwo] = true;
      }
   }

   if (m_ContainsAsSubgroup[SYMGROUP_AxesReflect_Any] && m_ContainsAsSubgroup[SYMGROUP_AxisPerms_All])
      // the Axis-C4 rotations do swaps like (x,y) --> (-y,x) (for any pair of
      // coordinates x,y). If a group G has both the full S_N of axis perms and
      // the full M_1 of axis reflections as subgroups, that is clearly
      // sufficient to get the same result as any Axis-C4.
      m_ContainsAsSubgroup[SYMGROUP_AxesRotate4_AnyTwo] = true;

   if (!m_ContainsAsSubgroup[SYMGROUP_AxesRotate4_AnyTwo]) {
      // finally, check C4 rotations on axis pairs. These allow extending even-
      // permutation symmetry to all-permutation symmetry for the subset of
      // monomials with all-even exponents (and its presence means that
      // all-odd monomials with two identical exponents get annihilated)
      FAxisPermuflectionList ToTest;
      // make list of all permutation-unique combinations of two simultaneous axis reflections
      for (size_t iAxis = 0; iAxis < iAxisEnd; ++ iAxis)
         for (size_t jAxis = iAxis+1; jAxis < std::min(size_t((iAxis+1) + jAxisCount), size_t(nAxes)); ++ jAxis)
            ToTest.emplace_back(FAxisPermuflection::i2QuarterTurn(iAxis, jAxis, m_nAxes, m_nPhases));
      if (_AllPresent(ToTest)) {
         m_ContainsAsSubgroup[SYMGROUP_AxesRotate4_AnyTwo] = true;
      }
   }

   // For the moment, I do not mean to use explicit canonicalization. Rather, I
   // think I'll just go with explicitly computing permuflection-group averages
   // of all non-clearly-annihilated monomials, and then deciding
   // representatives based on the result of the averaging. It looks like it
   // should be fast enough, all exact integer arithmetic and so on (especially
   // compared to LA calculations with the emulated high-precision floats)

   // However, even to start out, I think it makes sense to use some subgroup
   // symmetries to skip group-averaging easily identified target monomials
   // which would be annihilated in any case. Details:
   //
   // First, recall (see comment in header following FMonomialSymmetry) that:
   //
   // - Let H ⊆ G be a subgroup of G. We define the group-averaging operator
   //   P_H and the group summation operator ∑_H, both mapping
   //   the function space L²(R^n,C) onto itself, as
   //
   //        ∑_H[m_i] := ∑_{ĥ∈H} ĥ m_i
   //        P_H[m_i] := (1/|H|) ∑_H[m_i]
   //                 =  (1/|H|) ∑_{ĥ∈H} ĥ m_i
   //
   //   where m_i ∈ L²(R^n,C) is a target function.
   //
   // - Let H ⊆ G be a subgroup of G. Then:
   //
   //        P_H[m_i] = 0  implies  ∑_{ĝ ∈ G} m_i(R_g r⃗_s) w_s = 0,       (PH0)
   //
   //   by symmetry alone, irrespectively of any grid parameters (r⃗_s, w_s).
   //   Therefore if a subgroup H can be identified for which the group-average
   //   P_H annihilates a function m_i, the function m_i can be discarded from
   //   the optimization target space without consequences.
   //
   // Now, some consequences for concrete groups.
   //
   // - Let G ⊇ SYMGROUP_AxesReflect_AnyTwo. The two-axis reflections act on
   //   Cartesian monomials m_{ijk…}(x,y,z) = x^i y^j z^k… as follows:
   //
   //        (σ_{xy} m_{ijk…}) = (-1)^i (-1)^j m_{ijk…}
   //        (σ_{yz} m_{ijk…}) = (-1)^j (-1)^k m_{ijk…}
   //        (σ_{zx} m_{ijk…}) = (-1)^i (-1)^k m_{ijk…}
   //
   //   + If *EITHER* all exponents i,j,k,… are even *OR* all exponents
   //     i,j,k,… are odd, then all the phase factors on the rhs evaluate to
   //     unity. For such monomials, the subgroup sum over M_2 := {id, σ_{xy},
   //     σ_{xz}, …} evaluates to |M_2| m_{ijk} ≠ 0.
   //
   //   + However, if there is a mixture of even and odd exponents, then
   //     at least one phase factor on the rhs is (−1). Without loss of
   //     generality, let us assume this applies to σ_{xy}, i.e.:
   //
   //        (σ_{xy} m_{ijk…}) = (−1) ⋅ m_{ijk…}
   //
   //     Regarding that, note that (a) H := {id, σ_{xy}} forms a
   //     two-element subgroup H ⊂ G, and (b) its group averaging
   //     operator P_H annihilates m_{ijk…}:
   //
   //        P_H[m_{ijk…}]
   //        = (1/|H|) (∑_{h∈H} h) m_{ijk…}
   //        = (1/2) (id + σ_{xy}) m_{ijk…}
   //        = (1/2) (id m_{ijk…} + σ_{xy} m_{ijk…})
   //        = (1/2) (m_{ijk…} + (-1) m_{ijk…})
   //        = (1/2) (0) = 0
   //
   //     By eq (PH0) this implies `m_{ijk…}` can be discarded from the target
   //     space without consequences.
   //
   //   ⇒ SYMGROUP_AxesReflect_AnyTwo SUMMARY:
   //     For a group G with this subgroup, a Cartesian monomial `m_{ijk…}`
   //     will be symmetry-annihilated, and therefore can be discarded from the
   //     target space a-priorily, unless *EITHER* all exponents {i,j,k,…} are
   //     even *OR* all exponents {i,j,k,…} are odd.
   //
   // - Let G ⊇ SYMGROUP_AxesReflect_Any. The single-axis reflections act on
   //   Cartesian monomials m_{ijk…}(x,y,z) = x^i y^j z^k… as follows:
   //
   //        (σ_x m_{ijk…}) = (-1)^i m_{ijk…}
   //        (σ_y m_{ijk…}) = (-1)^j m_{ijk…}
   //        (σ_z m_{ijk…}) = (-1)^k m_{ijk…}
   //
   //   + Unless all exponents {i,j,k,…} are even, at least one of the phase
   //     factors on the rhs evaluates to (-1). Without loss of generality, let
   //     us assume that is the phase factor for σ_x. Then again we find:
   //     (a) H := {id, σ_x} forms a two-element subgroup H ⊂ G.
   //     (b) H's group averaging opeartor annihilates m_{ijk…}:
   //
   //           P_H m_{ijk…} = (1/2) (id + σ_x) m_{ijk…}
   //                        = (1/2) (m_{ijk…} - m_{ijk…}) = 0.
   //
   //   ⇒ SYMGROUP_AxesReflect_Any SUMMARY:
   //     For a group G with this subgroup, a Cartesian monomial `m_{ijk…}`
   //     can be discarded from the target space a-priorily, unless all
   //     exponents {i,j,k,…} are even.
   //
   // - @Generality: The results on SYMGROUP_AxesReflect_Any (M₁) and
   //   SYMGROUP_AxesReflect_AnyTwo (M₂) hold for general space dimensions
   //   D ≥ 2, and pose no conditions on G beyond M_i ⊂ G being a subgroup.
   //
   //   However, non-trivial generalizations could occur in complex spaces, as
   //   they can support more general phase symmetries than the straight
   //   σ_{…} axis reflections.
   //
   // - …that's almost it regarding symmetry annihilation. The only general
   //   annihilator we have left is a bit on the odd side:
   //
   // - Let G ⊇ H = SYMGROUP_AxesRotate4_AnyTwo. A symmetry op from H maps a
   //   coordinate pair (x,y) according to
   //
   //        c4_{xy} (x,y) := (−y,x),
   //
   //   leaving all other coordinates z,w,… invariant.
   //
   //   + For all-even monomials m_{ijk…}, the minus sign gets raised to an
   //     even power, and therefore cancelled — so the c4_{xy} then effectively
   //     adds a (x,y) ⟼ (y,z) axis swap symmetry generator (but that is not
   //     helpful in `ZeroAverageQ`).
   //
   //   + For all-odd monomials m_{ijk…}, this transformation does not
   //     generally afford useful conclusions… unless we have the special case
   //     of two or more exponents {i,j,k,…} coinciding.
   //     For example, c4_{xy} applied to m_{ij…} with i = j (both odd):
   //
   //        c4_{xy} m_{ij…}(x,y,…)
   //        = c4_{xy} x^i y^i …
   //        = (−y)^i x^i …   | (factor out (-1)^i = −1 for odd i)
   //        = −x^i y^i
   //
   //     So the subgroup average P_H over H := {id, c4_{xy}} will annihilate
   //     any monomial with odd+equal exponents i = j. Similarly, the C4s for
   //     the other axis pairs will symmetry-annihilate all other odd+equal-
   //     exponent pairs (e_i, e_j) with e_i = e_j.

}



groups::FCosetSide _ConvertCosetSide(FAxisPermuflectionSubgroup::FCosetSide iSide) {
   using aps = FAxisPermuflectionSubgroup;
   static_assert(size_t(aps::COSET_Left) == size_t(groups::COSET_Left));
   static_assert(size_t(aps::COSET_Right) == size_t(groups::COSET_Right));
//    assert(Side == aps::COSET_Left || Side == aps::COSET_Right);
   if (aps::COSET_Left)
      return groups::COSET_Left;
   else if (aps::COSET_Right)
      return groups::COSET_Right;
   else
      throw std::invalid_argument("_ConvertCosetSide():: unexpected coset type argument");
}


FAxisPermuflectionSubgroup::FGroupElementList_Parent const &FAxisPermuflectionSubgroup::CosetRepsQ(FCosetSide Side) const
{
   _ConvertCosetSide(Side);
   return (Side == COSET_Left)? m_LeftCosetReps : m_RightCosetReps;
}


char const *CosetSideAsString(FAxisPermuflectionSubgroup::FCosetSide Side) {
   using aps = FAxisPermuflectionSubgroup;
   return ((Side == aps::COSET_Left)? "left" : ((Side == aps::COSET_Right)? "right" : "??#unk??"));
}


void FAxisPermuflectionSubgroup::_MakeCosetReps(FGroupElementList_Parent &CosetReps, FCosetSide Side, FGroupElementList_Parent const &ParentGroupElements)
{
   // make ref to elements of the subgroup X --- this object's own
   FGroupElementList const
      &xElements = m_Elements;

   if (true && xElements.size() == ParentGroupElements.size()) {
      // with X ⊂ G and |X| = |G|, the subgroup X and the full group G must be
      // identical. We could still go through with the motions (the code below
      // works just fine in this case), …but we can also just assign id ∈ G
      // as representative of the only "coset" of X in G, and be done with it.
      return CosetReps.assign({FGroupElement_Parent::Identity()});
   }

   // make a list of parent group elements reduced from big multiprecision
   // scalars to double scalars. should be more than enough (even 32 bit floats
   // would be fine, as we're only dealing with orthogonal matrices)
   typedef TRotationMatrix<double> FGroupElement1d;
   std::vector<FGroupElement1d>
      gElements;
   for (FRotationMatrix const &Rg : ParentGroupElements)
      gElements.emplace_back(Rg.cast<double>());
   groups::FElementSearchTable<FGroupElement1d, FSortKey_RotationMatrix>
      gTable(gElements, FSortKeyGen_RotationMatrix<double>());
   // instanciate the identity element of the reduced version of the parent group
   FGroupElement1d
      gIdentity = FGroupElement1d::Identity();

   auto [iCosetReps, iElementCosets] = groups::_AssignInitialCosets(xElements, _ConvertCosetSide(Side), gElements, gTable);

   // We now got representatives of each coset, and formally this is all
   // we promised about them.
   //
   // However, there are some caveats:
   //
   // - Each coset has been assigned a unique index, but which index this
   //   is, is essentially random. In particular, identity need not be
   //   in coset #0
   //
   // - For each coset, which of its representatives has been chosen
   //   is essentially random. In particular, there is no guarantee that
   //   identity has been chosen as representative of its coset.
   //
   // - The randomness is more complicated than it looks: we use hash
   //   based lookups to make keys. And std::hash is not guaranteed to
   //   yield consistent outputs in subsequent program runs.
   //   So even if running the program multiple times in a row, on the same
   //   computer, each time may yield different choices of coset indices  and
   //   coset representatives
//    if (m_PrintLevel.MoreQ())
//       io.Write(" Sorted |G| = {} group elements {{g; g ∈ G}} into {} {} cosets of subgroup X (with |X| = {})",
//          gElements.size(), iCosetReps.size(), AsString(Side), xElements.size());

   // try some explicit beautification for special cases
   if (gElements.size() != xElements.size()) {
      // try for the best outcome: a finding a single group element
      // ĝ of which the ĝ^i can serve as coset reps for all the cosets.
      std::vector<size_t>
         iCosetReps_CyclicGen = groups::_TryFindCosetSpanningCyclicGenerator(
               iElementCosets, iCosetReps, gIdentity, gElements, gTable);

      if (!iCosetReps_CyclicGen.empty()) {
         // note: ĉ can't be in X, because then it'd be in the same coset as
         // identity, and would not generate anything.
         if (m_PrintLevel.MostQ())
            io.Write("\n    Found cyclic generator `c ∈ (G ∖ X)` with |⟨c⟩| = {} and all c^i (i=0,...,{}) in distinct cosets (of [G:X] = {})",
               iCosetReps_CyclicGen.size(), iCosetReps_CyclicGen.size()-1, iCosetReps.size());
         // take orbit of ĝ as new coset representatives. With coset #i
         // being represented by g^i.
         assert(iCosetReps_CyclicGen.size() == iCosetReps.size());
         iCosetReps.swap(iCosetReps_CyclicGen);
         // Update ⟨element⟩ ⟼ ⟨coset index⟩ mapping to reflect the new
         // coset order implied by the new representatives
         groups::_UpdateCosetMapping(iElementCosets, iCosetReps);
      } else {
         // …well, I guess we can just take the cosets in random order, and
         // random coset reps of each, and leave it at that.
      }
   }

   // assemble a list of actual parent group elements (not just reduced
   // elements, keys, or indices) for the chosen final coset representatives
   cx_assert_rt(CosetReps.empty());
   CosetReps.reserve(iCosetReps.size());
   for (size_t ic : iCosetReps)
      CosetReps.emplace_back(ParentGroupElements[ic]);

//    if (0) {
//       // Hmpf. This does not work reliably at all. I feared it looked too simple
//       // to be right.
//       //
//       // I guess at this moment I will simply skip the construction of actual generators,
//       // they are not *strictly* needed.
//       //
//       // At some later point in time it might be reasonably to get back to this,
//       // by patching up the CxPermGroup.* stuff with the python code updates —
//       // at least the permuflection groups can be trivially embedded in a 2d-
//       // element permutation group (with d = space dimension), and
//       // FStabilizerChain will happily construct very pretty subgroup chains and
//       // generator sets for anything you throw at it.
//       if (Side == 0) {
//          std::vector<groups::generator_info_t>
//             gGensInfo = _TryBuildCyclicGeneratorChain(gIdentity, gElements,
//                gTable, FPrintLevel::More, &io);
//
//          FAxisPermuflection
//             xIdentity = FAxisPermuflection::Identity(m_nAxes, m_nPhases);
//          groups::FElementSearchTable<FAxisPermuflectionSubgroup::FGroupElement, FSortKey_AxisPermuflection>
//             xTable(xElements, FSortKeyGen_AxisPermuflection());
//
//          std::vector<groups::generator_info_t>
//             xGensInfo = _TryBuildCyclicGeneratorChain(xIdentity, xElements,
//                xTable, FPrintLevel::More, &io);
//       }
//    }
   // ^-- triggers some weird build error with multi-precision floats atm.
}







#ifdef INCLUDE_ABANDONED
// // maps Cartesian exponents to total factor obtained after summing up P_f m_{ijk}
// // for all permuflections in the group. The final object is a polynomial, of
// // course. However, note that the final factors are not floats or anything this
// // time — the handling of permutations and reflections is exact, so we get well-
// // defined integers.
// typedef std::map<FMonomialN, phase_scaled_int_t>
//    FMonomialCoeffMap;
//
// template<class TArray>
// FMonomialN PermuteAxesT(FMonomialN const &ExpIn, TArray const &p) {
//    assert(ExpIn.size() == p.size());
//    std::array<FMonomialN::value_type, AIG_SPACE_AXES>
//       ExpOut;
//    assert(ExpIn.size() == ExpOut.size());
// //    std::array<unsigned, ExpIn.size()>
// //       ExpOut;
//    // ^-- interesting… haven't seen one of those in years:
//    //
//    //   AxisPermuflection.cpp: In function 'FMonomialN PermuteAxesT(const FMonomialN&, const TArray&)':
//    //   AxisPermuflection.cpp:190:37: internal compiler error: in coerce_template_parms, at cp/pt.c:8671
//    //     190 |    std::array<unsigned, ExpIn.size()>
//    //         |                                     ^
//    //   Please submit a full bug report,
//    //   with preprocessed source if appropriate.
//
//    for (unsigned i = 0; i != ExpIn.size(); ++i) {
//       // if the axis coordinates transform by permutation `p`, the exponent
//       // array by transforms by ~p, the inverse of permutation p (see desc in
//       // FMonomialSymmetry class)
//       assert(i < ExpIn.size() && p[i] < ExpIn.size());
//       ExpOut[i] = ExpIn[p[i]];
//    }
// #ifdef _DEBUG
//    {
//       FMonomialN mOut(ExpOut);
//       for (unsigned i = 0; i != ExpIn.size(); ++i)
//          assert(mOut[i] == ExpIn[p[i]]);
//    }
// #endif // _DEBUG
//    return FMonomialN(ExpOut);
// };
//
//
//
// // flatten {key --> coefficient} map into an array, removing all terms with zero
// // coefficients. Original order is preserved.
// FMonomialCoeffList AsCompressedCoeffList(FMonomialCoeffMap const &CoeffMap)
// {
//    // count number of monomial terms which survived the group summation.
//    size_t nTerms = 0;
//    for (auto const &mc : CoeffMap) // mc: monomial and coeff.
//       nTerms += (mc.second == 0)? 0 : 1;
//    if (nTerms == 0)
//       return {};
//
//    // now create the actual array as storage and collect them.
//    FMonomialCoeffList CoeffList;
//    CoeffList.reserve(nTerms);
//    for (auto const &mc : CoeffMap) {
//       if (mc.second != 0)
//          CoeffList.emplace_back(mc);
//          // ^-- theoretically we *could* get a move working, but for these
//          //     concrete types it would not actually do anything.
//          //     (FMonomialN is physically a single unsigned integer,
//          //     and we here use int coefficients (because the processing is exact!),
//          //     not the usual super-expensive emulated high-precision floats)
//    }
//    return CoeffList;
// }
//
//
// // template<class FCoeffMap, class T = std::pair<typename FCoeffMap::key_type, typename FCoeffMap::mapped_type> >
// // std::vector<T> AsCompressedCoeffList(FCoeffMap &&CoeffMap)
// // {
// //    // count number of monomial terms which survived the group summation.
// //    size_t nTerms = 0;
// //    for (auto const &mc : CoeffMap) // mc: monomial and coeff.
// //       nTerms += (mc.second == 0)? 0 : 1;
// //    if (nTerms == 0)
// //       return {};
// //
// //    // now create the actual array as storage and collect them.
// //    std::vector<T> CoeffList;
// //    CoeffList.reserve(nTerms);
// //    for (auto &&mc : CoeffMap)
// //       if (mc.second != 0)
// //          CoeffList.emplace_back(mc.first, std::move(mc.second));
// //          // ^-- @no move on first: we are probably not supposed to destroy a key
// //          //     of a still-alive std::map object.
// //    return CoeffList;
// // }


// // For the Cartesian monomial
// //
// //    m[i0,i1,i2](x0,x1,x2) = x0^{i0} x1^{i1} x2^{i2},
// //
// // computes the group sum
// //
// //    \sum_{a ∈ A} (a m)_{i0,i1,i2}(x0, x1, x2)
// //
// // over G's subgroup of permuflection symmetry operations A.
// FMonomialCoeffList ComputeAxprGroupSum(FMonomialN m, FAxisPermuflectionList const &GroupElements)
// {
//    FMonomialCoeffMap
//       out;
//    for (FAxisPermuflection const &pf : GroupElements) {
//       // apply the permutation part of the group element to the monomial,
//       // and back-tracke the axis coordinate permutation into a exponent
//       // permutation.
// //       auto const &Perm = pf.PermQ();
//       FMonomialN pm = PermuteAxesT(m, pf.PermQ());
//       // get a reference to the coefficient for ther monomial with the permuted
//       // axes. std::map's []-operator will create the map item if it does not
//       // yet exist, and in doing so default-construct (=set to zero) the
//       // coefficient. So that is actually exactly what we want.
//       int &iCoeff_pm = out[pm];
//
//       // Find out what happens to the transformed monomial's sign factor.
//       // (for explanation, see comment below this function)
//
//       auto const &Signs = pf.SignsQ();
//       assert(Signs.size() == m.size());
//
//       phase_t TotalPhase(1);
//       for (size_t i = 0; i != m.size(); ++ i) {
//          assert(Signs[i] == -1 || Signs[i] == +1);
//          // As shown above, the exponent e[i] combines with the
//          // sign factor s[i] of the same axis, therefore yielding
//          // a total axis factor of (s[i])^e[i].
//          //
//          // so collect the exponent to which the i`th source cartesian
//          // coordinate [x0, x1, …], and therefore its permuflection sign factor
//          // s[i], is raised.
//          unsigned exp_i = m[i];
//          if (exp_i % 2 == 0) {
//             // power is even, so we get (s[i])^e[i] = 1 regardless of whether
//             // s[i] is positive or negative. No sign factor from this axis.
//          } else {
//             // sign is raised to an odd power, which simply preserves it
//             // unchanged: s^e == s  if s∈{-1,1} and `e` is an odd integer. So
//             // the axis sign factor applies unchanged to the total monomial
//             // phase factor.
//             TotalPhase = (Signs[i] < 0)? -TotalPhase : TotalPhase;
//          }
//       }
//
//       // increment the coefficient, to combine with the result from earlier
//       // and/or later group elements.
//       iCoeff_pm += int(TotalPhase);
//    }
//
// //    return AsCompressedCoeffList<FMonomialCoeffMap,FMonomialCoeffMap::value_type>(std::move(out));
//    return AsCompressedCoeffList(std::move(out));
//    // TODO: maybe just do the permutation calculation in here directly?
//    // And the collection of the signs. May be easier and more straight
//    // forward (and certainly faster) than what we do now.
//    // (I'm seriously considering just classifying all monomials like this,
//    // instead of deriving and applying effective general symmetries)
// }

//    // …and while we're at it, also remember which of them were even
//    // (--> SignMask[i] = 0) and which were odd (--> SignMask[i] = 1)
//    FAxisPermuflection::FSignsOnly
//       SignMask;
//    for (size_t i = 0; i < AIG_SPACE_AXES; ++ i) {
//       auto mExp_i = m[i];
//       ExpIn[i] = mExp_i;
//       SignMask[i] = (mExp_i % 2 == 0) ? 0 : 1;
//    }
#endif // INCLUDE_ABANDONED


// For the Cartesian monomial
//
//    m[i₀,i₁,i₂](x₀,x₁,x₂) ≡ x₀^{i₀} x₁^{i₁} x₂^{i₂},
//
// computes the group sum
//
//    ∑_{â ∈ A} (â m)_{i₀,i₁,i₂}(x₀, x₁, x₂)
//
// over G's subgroup of permuflection symmetry operations A.
void ComputeAxprGroupSum(FMonomialSymmetry::FMonomialCoeffList &terms,
      FMonomialN m, FAxisPermuflectionList const &GroupElements,
      unsigned nAxes, unsigned nPhases)
{
   // there should be at least an identity element, and therefore at least one term.
   // (the function will work and correctly return 0 if that is not the case)
   assert(!GroupElements.empty());

   // clear out temp space
   terms.clear();

   // unpack the input monomial's exponents
   std::array<FMonomialN::value_type, AIG_SPACE_AXES>
      ExpIn;
   for (size_t i = 0; i < AIG_SPACE_AXES; ++ i)
      ExpIn[i] = m[i];
   assert(ExpIn.size() == m.size());
   assert(nAxes == AIG_SPACE_AXES);
   auto const _apply_phase_factor = _PhaseFactorFn_InplaceApply<phase_scaled_int_t>(nPhases);

   // apply all group operators to m, storing the results in `out`
   for (FAxisPermuflection const &pf : GroupElements) {
      // Note: see comment under this function for derivation of the formulas
      // for permuflections applied to monomials (axis permutation/sign parts)
      assert(pf.nPhases() == nPhases);

      // compute exponents of permflected monomial
      std::array<FMonomialN::value_type, AIG_SPACE_AXES>
         ExpOut;
      for (size_t i = 0; i < AIG_SPACE_AXES; ++i)
         ExpOut[i] = ExpIn[pf.TargetAxisQ(i)];
      // …and the total phase factor resulting from the combination
      // of axis reflections with the monomial's coordinate exponents
      phase_scaled_int_t TotalPhase(1);
#ifdef INCLUDE_ABANDONED
//       assert(pf.nPhases() == 2);
//       for (size_t i = 0; i < pf.nAxes(); ++i)
// //          if (ExpIn[i] % 2 != 0 && pf.SignQ(i) < 0)
//          if (ExpIn[i] % 2 != 0 && pf.PhaseIndexQ(i) != 0)
//             TotalPhase = -TotalPhase;
//       if (pf.nPhases() == 2) {
//          for (size_t i = 0; i < pf.nAxes(); ++i)
//             if (ExpIn[i] % 2 != 0 && pf.PhaseIndexQ(i) != 0)
//                TotalPhase = -TotalPhase;
//       } else {
//          for (size_t i = 0; i < pf.nAxes(); ++i)
//             _apply_phase_factor(TotalPhase, ExpIn[i]*pf.PhaseIndexQ(i));
//       }
#endif // INCLUDE_ABANDONED
      for (size_t i = 0; i < pf.nAxes(); ++i)
         _apply_phase_factor(&TotalPhase, ExpIn[i]*pf.PhaseIndexQ(i));

      // append {permuted-monomial : phase-factor} temporarily to intermediates.
      // We will combine terms after all are done.
      terms.emplace_back(ExpOut, TotalPhase);
   }

   // we now have one term for each group element. Sort intermediates, then
   // combine all terms with equivalent keys (=monomial exponents). Note that
   // after sorting, all combinable monomials will necessarily lie adjacent to
   // each other.
   std::sort(terms.begin(), terms.end());

   {
      size_t iout = 0, iin = 0;
      while (iin < terms.size()) {
         // take next element at iin as candidate
         FMonomialSymmetry::FMonomialCoeffList::value_type
            acc = terms[iin];
         iin += 1;
         // find all other sum terms for this monomial, and sum up their
         // coefficients
         for (; iin < terms.size() && terms[iin].first == acc.first; ++ iin)
            acc.second += terms[iin].second;
         // anything left?
         if (acc.second != 0) {
            // yes, add it to the retained list of output terms. This overwrites
            // some earlier input terms, but we are already done with them.
            terms[iout] = acc;
            iout += 1;
         }
      }
      // shallow-delete all terms except the ones resulting from
      // non-zero accumulation results
      terms.resize(iout);
   }
}


// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// DETAILS & DERIVATIONS: How does a permuflection-transform a monomial?
// ─────────────────────────────────────────────────────────────────────────────
//
// This is more subtle than it might appear initially, so let us have a look:
//
// o The exact operation pf represents, for each axis, is this:
//   - The source coordinate x[i] gets mapped to (Sign[i] * x[Perm[i]]);
//   - That is, the sign applies to the source coordinate #i, which by means of
//     the permutation, gets turned into the target coordinate x[Perm[i]].
//
// o We should now find out what applying the three axis transforms implies for
//   the phase factor of the Cartesian monomial. The source monomial with
//   exponents [e₀, e₁, e₂] is:
//
//       x₀^{e₀} x₁^{e₁} x₂^{e₂},
//
//   where (xᵢ)ᵢ are the Cartesian coordinates. As just described, upon
//   application of the group element Pf, each of the coordinates undergoes the
//   change:
//
//       Pf[x_i] ⟼ (s[i] * x_{p[i]}).
//
//   However, Pf, formally, acts only an Cartesian coordinates, not on the
//   exponents or the formal structure of the monomial (…although that can be
//   rearranged, it is not Pf which is doing it — we are!).
//
// o Taking these points together, we find that the symmetry operation Pf acts
//   as follows on the original monomial (notation; αᵢ≡s[i], πᵢ≡p[i]):
//
//       Pf[m]
//       = Pf[x₀^{e₀} x₁^{e₁} x₂^{e₂}]
//       = (α₀ x_{π₀})^{e₀} ⋅ (α₁ x_{π₁})^{e₁} ⋅ (α₂ x_{π₂})^{e₂}
//
//   In this, note that although the sign factors αᵢ (≡`s[i]`) now stand with different
//   Cartesian axes, the exponents which apply to them are the ones which
//   the source coordinates **originally had**: E.g., α₀ ≡ s[0] gets raised to the
//   e₀'th power, even if the former `x` coordinate, now turned into `x_{p[0]}`,
//   was replaced by a `y` or `z` coordinate.
//
// o We can now collect, combine, and reorder terms in order to find a mapping
//   of the coordinate-permuted monomial into an exponent-permuted monomial in
//   terms of the original Cartesian axes:
//
//       Pf[m] = ⟨phase-factor⟩ ⋅
//               x₀^{e_{(~π)₀}} ⋅ x₁^{e_{(~π)₁}} ⋅ x₂^{e_{(~π)₂}}
//
//   where, as before, `π` (with elements `πᵢ≡p[i]`) denotes the Cartesian axis
//   coordinate permutation (`FAxisPerm`), and `~π` its inverse (which is what
//   gets applied to the exponents!).
//
// o From the term before, we find that the total phase factor becomes:
//
//       ⟨phase-factor⟩ = (α₀)^{e₀} ⋅ (α₁)^{e₁} ⋅ (α₂)^{e₂}
//                      = (s[0])^{e_{0}} ⋅ (s[1])^{e_{1}} ⋅ (s[2])^{e_{2}}.
//
//   So while the permutation `p` controls *which* new monomial is reached, it
//   does not actually affect the phase factor of the symmetry element
//   between the source and the target monomial.




FMonomialSymmetry::FMonomialSymmetry(FRotationMatrixList const &GroupElements, unsigned nAxes, unsigned nPhases, std::string const &GroupName, FPrintLevel PrintLevel)
   : m_pAxisPermuflections(new FAxisPermuflectionSubgroup(GroupElements, nAxes, nPhases, GroupName, PrintLevel)),
     m_GroupName(GroupName), m_PrintLevel(PrintLevel)
{
#ifdef INCLUDE_ABANDONED
//    // use the permuflections on the Cartesian coordinates to derive the actual
//    // symmetry operations we can perform on monomials using the Cartesian
//    // coordinates as arguments. These symmetry operations differ based on
//    // monomial type (combinations of exponents based on even-ness)
//    m_TypeSymmetries.reserve(nMonomialTypes);
//    for (unsigned mt = 0; mt != nMonomialTypes; ++ mt) {
//       m_TypeSymmetries.emplace_back(m_AxisPermuflections, FMonomialSymmetryOpList::FMonomialType(mt), m_GroupName);
//    }
#endif // INCLUDE_ABANDONED
};


FMonomialSymmetry::FMonomialSymmetry(FAxisPermuflectionSubgroupPtr pAxisPermuflections, FPrintLevel PrintLevel)
   : m_pAxisPermuflections(pAxisPermuflections),
     m_GroupName(pAxisPermuflections->ParentGroupNameQ()),
     m_PrintLevel(PrintLevel)
{
}



bool FMonomialSymmetry::ZeroAverageQ(FMonomialN const &m) const
{
   // If a Cartesian monomial has any odd exponent, then it gets annihilated by
   // a group average over the full SO(3) rotation group. But for the finite
   // subgroups of SO(3)/O(3) this is not necessarily the case.
   //
   // Nevertheless, the check if all exponents in `m` are even OR all exponents
   // are odd is almost sufficient to cast the screening decision for the 3D
   // polyhedral groups, so start with that
   bool
      AllExpEven = true,
      AllExpOdd = true;
   for (size_t i = 0; i < m.size(); ++ i) {
      if (m[i] % 2 == 0) {
         AllExpOdd = false;
      } else {
         assert(m[i] % 2 == 1);
         AllExpEven = false;
      }
   }

   // What is this? Documentation cross-refs:
   // - The subgroups and where they occur in 3D ⟹ see comments in `FSubGroupId` (in the header)
   // - Concrete subgroup tests performed here ⟹ comments in `FAxisPermuflectionSubgroup::_SearchSpecialSubgroups`
   // - General concept of monomial screening ⟹ text after `FMonomialSymmetry` (in the header)
   if (m_pAxisPermuflections->ContainsAsSubgroupQ(SYMGROUP_AxesReflect_Any)) {
      if (!AllExpEven)
         return true; // `m` gets symmetry-annihilated ⟹ ZeroAverageQ(m) == true
   } else if (m_pAxisPermuflections->ContainsAsSubgroupQ(SYMGROUP_AxesReflect_AnyTwo)) {
      if (!(AllExpEven || AllExpOdd))
         return true; // `m` gets symmetry-annihilated ⟹ ZeroAverageQ(m) == true
   } else {
      // We technically could check for some other point groups with reflective
      // subgroups. E.g., just σ_x ∈ G is sufficient to annihilate monomials
      // m_{ijk…} with odd exponents `i`. However…
      // - …so far we have implemented only 3d polyhedral groups. And all of
      //   them have `SYMGROUP_AxesReflect_AnyTwo` as subgroup, so there was no
      //   practical need (it could easily be added to `_SearchSpecialSubgroups`
      //   if really needed, or one could just keep the `std::set` of
      //   `AxisPermuflections` and do the testing here).
      // - …for the smaller groups also the permuflection group sums are
      //   (even) more efficient. These handle such symmetries explicitly — this
      //   here is just pre-screening. So it is probably not really worth it (at
      //   least not in 3D).
   }

   if (AllExpOdd && m_pAxisPermuflections->ContainsAsSubgroupQ(SYMGROUP_AxesRotate4_AnyTwo)) {
      // This one symmetry-annihilates all monomials m_{ijk…} in which:
      // - all exponents in {i,j,k,…} =: {e_0, e_1, e_2, …} are odd,
      // - AND at least one exponent coincides with another exponent
      //   (i.e., there is some (i,j) with i ≠ j and e_i = e_j)
      // Unless both conditions apply, this symmetry affords no statements
      // regarding zero averages.
      for (size_t i = 0; i + 1 < m.size(); ++ i)
         for (size_t j = i + 1; j < m.size(); ++ j)
            if (m[i] == m[j])
               return true; // found coinciding exponents -> `m` gets zero-averaged -> ZeroAverageQ(m) == true
   }

   // Did npt find any reason why `m` should get symmetry-annihilated by symmetry.
   // Note: for Cartesian monomials and 3d polyhedral groups, the tests above
   // near-exhausive. For many groups they cover all zero average cases, and for
   // others almost all of them.
   return false;
}


bool FMonomialSymmetry::SelectAsCanonialReprQ(FMonomialN const &m) const
{
   // let's do this the easy way for a start…
   FMonomialCoeffList
      MonomialCoeffs;
   if (!EvalEasySubgroupSum(MonomialCoeffs, m))
      // ⬑ note: this internally calls `ZeroAverageQ` before actually computing
      //   anything. So we don't have to do that, too.
      return false;
   cx_assert_rt(!MonomialCoeffs.empty());
   // well… the monomial coefficients *are* sorted… by… something.
   // Which is good! Sorted by “something” is exactly what we want!
   //
   // (…it is likely the packed index of `FMonomialN`, which would yield some
   // type lexicographical sort. Seems fine, considering any canonical
   // representative selection is in principle arbitrary.)
   return m == MonomialCoeffs.front().first;
}




bool FMonomialSymmetry::EvalEasySubgroupSum(FMonomialCoeffList &MonomialCoeffs, FMonomialN m) const
{
   if (this->ZeroAverageQ(m)) {
      MonomialCoeffs.clear();
      return false;
   } else {
      ComputeAxprGroupSum(MonomialCoeffs, m, m_pAxisPermuflections->ElementsQ(), m_pAxisPermuflections->nAxes(), m_pAxisPermuflections->nPhases());
      return !MonomialCoeffs.empty();
   }
}


FPolynomialN FMonomialSymmetry::EvalEasySubgroupSum(FPolynomialN const PolyIn) const
{
   FPolynomialN
      out(PolyIn.GetCullCriteria());
   FMonomialCoeffList
      temp;
   for (auto &&[mIn,cIn] : PolyIn) {
      if (EvalEasySubgroupSum(temp, mIn)) {
         for (auto &&[mSub,cSub] : temp) {
            out.Add(cSub*cIn, mSub);
         }
      }
   }
   return out;
}



// mostly meant to get printing working quickly.
template<class FMonomialCoeffSeq>
FPolynomialN AsPolynomial(FMonomialCoeffSeq const &Coeffs) {
   FPolynomialN
      out;
   for (auto const &mc : Coeffs) {
      out.Add(FScalar(mc.second), mc.first);
   }

   // we already cancelled the redundant terms. There should be nothing
   // left which needs deleting.
   cx_assert_rt(out.size() == Coeffs.size());
   return out;
};


void pEmitRule(int iRuleConfig) {
   size_t LineWidth = 80;
   char const *pBarChar = (iRuleConfig == 0)? "─" : "━";
   io.Write(" {}", Repeat(LineWidth-1, pBarChar));
};


void AnalyzeMonomials(std::vector<FMonomialN> Monomials, FAxisPermuflectionList const &GroupElements, unsigned nAxes, unsigned nPhases, std::string const &GroupName, FPrintLevel PrintLevel)
{
   int Verbosity = int(PrintLevel);
   Verbosity = 3;

   io.Write("");
   io.Write("      {:<8}  {}", fmt::format("|A| = {}", GroupElements.size()), "Group elements in permuflection subgroup:");
   {
      size_t nEntriesPerLine = 8;
      std::vector<FAxisPermuflection> subset;
      auto Flush = [&]() {
         if (!subset.empty()) io.Write("      {:<8} {}", "// ", subset);
         subset.clear();
      };
      for (auto &&ge : GroupElements) {
         subset.push_back(ge);
         if (subset.size() % nEntriesPerLine == 0)
            Flush();
      }
      Flush();
   }


   io.Write("");
   FMonomialSymmetry::FMonomialCoeffList CoeffList;
   for (auto &&[im, m] : enumerate(Monomials)) {
#ifdef INCLUDE_ABANDONED
//       if (im != 0 && Verbosity >= 3)
//          io.Write(Repeat(78, "─", "  ", "\n"));
//       FMonomialCoeffList
//          CoeffList = ComputeAxprGroupSum(m, GroupElements);
#endif // INCLUDE_ABANDONED
      ComputeAxprGroupSum(CoeffList, m, GroupElements, nAxes, nPhases);
      char const *pFmt = "      {:<8} P_{{A}}[m_{{{},{},{}}}] = {}\n"; // axis aligned permutation subgroup?
      FPolynomialN
         p1 = AsPolynomial(CoeffList);
      if (1) {
         FPolynomialN::FCullCriteria Cull(g_ThrAlmostZero);
         p1.SetCullCriteria(Cull, true); // true: apply purge with new criteria now
      }
      if (Verbosity >= 3) {
         io.Write(pFmt, fmt::format("#m = {}",p1.size()), m[0], m[1], m[2], p1);
      }
   }
}


// std::string AsString(FMonomialSymmetry::FMonomialType mt) {
//    if (mt == FMonomialSymmetryOpList::AllExpEven)
//       return "m[2i,2j,2k] -- All exponents EVEN";
//    if (mt == FMonomialSymmetryOpList::AllExpOdd)
//       return "m[(2i+1),(2j+1),(2k+1)] -- All exponents ODD";
//    if (mt == FMonomialSymmetryOpList::NeitherAllEvenNorAllOdd)
//       return "m[i,j,k] -- Mixed exponents, including even and odd";
//    cx_assert_rt(0);
//    throw std::runtime_error("unexpected monomial type in AsString(FMonomialSymmetry::FMonomialType mt)");
// };


// // Derives the set of symmetry operations applicable to monomials
// // of the given type from special symmetry elements of the full point group
// FMonomialSymmetryOpList::FMonomialSymmetryOpList(FAxisPermuflectionSubgroup const &AxisPermuflections, FMonomialType MonomialType, std::string const &GroupName)
//    : m_MonomialType(MonomialType), m_GroupName(GroupName)
// {
// //       AllExpEven, // all monomial exponents are even. E.g., x^4 y^2 z^2 counts.
// //       AllExpOdd, // all monomial exponents are odd. E.g., x^7 y^3 z^9 counts.
// //       NeitherAllEvenNorAllOdd, // contains mixed even-ness exponents.
//    pEmitRule(1);
//    io.WriteLine();
//    io.WriteLine();
//    if (MonomialType == 0) {
//       io.WriteLine();
//       io.WriteLine();
//    }
//
//    pEmitRule(1);
//    io.Write(" ░░ Permuflection group sum analysis: ░░ {} ░░ {}", m_GroupName, AsString(MonomialType));
//    pEmitRule(0);
//
//    // hm… I actually think, currently, that we can really just take a
//    // random sample monomial of each of the forms [e0,e1,e2], [e0,e0,e2], [e0,e0,e0]
//    // (where the exponents e0,e1,e2 are mutually distinct, of course), and
//    // just collect the transforms for that. I think they would apply unchanged
//    // to any other monomial of said form.
//    //
//    // Zero exponents might need a bit of extra thinking, but since zero
//    // is also a even number which, as exponent, anyway cancel the sign factors,
//    // it might not really do anything after all. Could check.
//    //
//    //
//
//    FPrintLevel PrintLevel(4);
//    typedef std::vector<FMonomialN> FExampleList;
//    FExampleList Monomials;
//
//    if (m_MonomialType  == AllExpEven) {
//       Monomials = FExampleList{{0,2,4},{6,12,18},{4,4,14},{0,0,0},{2,2,2},{6,6,6},{16,16,16}};
//    } else if (m_MonomialType == AllExpOdd) {
//       Monomials = FExampleList{{1,3,5},{3,9,17},{3,3,17},{1,1,1},{17,17,17}};
//    } else {
//       Monomials = FExampleList{{0,3,5},{4,10,17},{2,3,17},{1,6,18},{16,18,15}};
//    }
//
// //    if (m_MonomialType  == AllExpEven)
//    AnalyzeMonomials(Monomials, AxisPermuflections.ElementsQ(), AxisPermuflections.nAxes(), AxisPermuflections.nPhases(), m_GroupName, PrintLevel);
// }


enum _InitFromSequenceMode {
   INIT_Copy,
   INIT_Convert = INIT_Copy,
   INIT_Move
};

// Set elements in `ArrayOut` from an data taken from an iterator range. Since
// `std::arrays` are fixed size, there are some choices to make. Here:
//
// - If the input range contains more elements than the array has, the
//   remaining entries in the input range will be ignored.
//
// - If the input range contains fewer elements than the array, then the
//   remaining elements will be explicitly default-initialized (for ints,
//   floats, etc this means that they will be set to zero)
//
// - Both cases currently trigger a failed assertion in debug mode. Not sure if
//   that will stay, or if the present behavior is useful enough in practice.
//
// The function assumes that whatever data types are copied are simple and small
// — it *does* copy them! (it will try to move if you explicitly tell it to, but
// in many situations that cannot be done or does not mean anything)
template<class T, size_t N, _InitFromSequenceMode Mode, class FInputIt>
void _InitArrayFromSequence(std::array<T,N> &ArrayOut, FInputIt first, FInputIt last) {
//    assert(std::distance(first,last) == ArrayOut.size());
   assert(size_t(last - first) == ArrayOut.size());
   size_t iArg = 0;
   for (FInputIt it = first; it != last; ++it, iArg += 1) {
      // too many init list elements for the array?
      if (iArg >= ArrayOut.size())
         // just ignore the extra values.
         return;
      if constexpr (Mode == INIT_Move) {
         // try to claim the input sequence items for our the new array,
         // even if that destroys them in the process.
         ArrayOut[iArg] = std::move(*it);
      } else {
         // make explicit construction/conversion. should allow T and *it to be different,
         // as long as T can be constructed from whatever *it contains.
         ArrayOut[iArg] = T(*it);
      }
   }
   // fill up remaining elements, if any. Note that (unlike `std::vector`, etc),
   // `std::array` objects have semantics closer to normal C arrays on auto
   // storage: if not explicitly initialized, expect random data for pods.
   for ( ; iArg < ArrayOut.size(); ++ iArg)
      ArrayOut[iArg] = T();
}


// Set elements in `ArrayOut` from brace-enclosed `initializer_list`.
template<class T, size_t N>
void _InitArrayFromInitList(std::array<T,N> &ArrayOut, std::initializer_list<T> InitList) {
   return _InitArrayFromSequence(ArrayOut, InitList.begin(), InitList.end(), INIT_Convert);
}
// ⬑ (cgk note to self: there appears to be no standard way of moving
// something out of an `std::initializer_list<>`. Its iterators only provide const
// access, so neither r-value refs nor `std::move` will do anything. Apparently
// the core idea of those things is that the compiler can keep the init block
// being references in static const memory which is embedded as-is in the
// executable. Certainly not what I would have expected this type to be for.)

#ifdef INCLUDE_ABANDONED
// FAxisPermuflection::FAxisPermuflection(std::initializer_list<FAxisEntry> vlArgs) {
//    _InitArrayFromInitList(m_AxisOps, vlArgs);
// }
//
// template<class InputIt>
// FAxisPermuflection::FAxisPermuflection(InputIt first, InputIt last) {
//    _InitArrayFromSequence(m_AxisOps, first, last);
// }


// FAxisPermuflection::FAxisPermuflection(std::initializer_list<FAxisEntry> vlArgs) {
//    FAxisOps AxisOps;
//    _InitArrayFromInitList(AxisOps, vlArgs);
//    _InitFromAxisEntries(AxisOps.begin(), AxisOps.end());
// }
//
// template<class InputIt>
// FAxisPermuflection::FAxisPermuflection(InputIt first, InputIt last) {
//    FAxisOps AxisOps;
//    _InitArrayFromSequence(m_AxisOps, first, last);
//    _InitFromAxisEntries(AxisOps.begin(), AxisOps.end());
// }
#endif // INCLUDE_ABANDONED

template<class InputIt>
void FAxisPermuflection::_InitFromAxisEntries(InputIt first, InputIt last)
{
//    assert(m_TargetAxes.size() == m_Signs.size() && m_Signs.size() == nAxes);
   assert(m_TargetAxes.size() == m_PhaseIndices.size() && m_PhaseIndices.size() == MaxAxes);
   size_t i = 0;
   for (InputIt it = first; it != last; ++ it) {
      if (i >= MaxAxes)
         throw std::runtime_error("FAxisPermuflection::c'tor(): too many axis entries in input.");
//       m_TargetAxes[i] = it->iTargetAxis;
//       m_Signs[i] = it->Sign;
      FAxisEntry ai(*it); // <-- this makes a copy or creates the object via c'tor from other types.
      m_TargetAxes[i] = ai.iTargetAxis;
//       m_Signs[i] = ai.Sign;
      m_PhaseIndices[i] = ai.iPhase;
      i += 1;
   }
   if (i != MaxAxes)
      throw std::runtime_error("FAxisPermuflection::c'tor(): not enough axis entries in input.");
}


FAxisPermuflection::FAxisPermuflection(FAxisOps const &AxisOps) {
   _InitFromAxisEntries(AxisOps.begin(), AxisOps.end());
}

FAxisPermuflection::FAxisPermuflection(std::initializer_list<FAxisEntry> vlArgs) {
   _InitFromAxisEntries(vlArgs.begin(), vlArgs.end());
}

FAxisPermuflection::FAxisPermuflection(std::initializer_list<int> vlArgs) {
   _InitFromAxisEntries(vlArgs.begin(), vlArgs.end());
}

template<class InputIt>
FAxisPermuflection::FAxisPermuflection(std::pair<InputIt, InputIt> range) {
   _InitFromAxisEntries(range.first, range.second);
}


FAxisPermuflection FAxisPermuflection::Identity(unsigned nAxes_, unsigned nPhases_)
{
//    assert(nAxes_ <= MaxAxes);
//    std::array<int, MaxAxes> id;
//    for (size_t i = 0; i < MaxAxes; ++i) id[i] = i + 1;
//    return FAxisPermuflection(id.begin(), id.begin() + nAxes_);
   cx_assert_rt(nAxes_ <= MaxAxes && nPhases_ <= MaxPhases);
   FAxisPermuflection out(0, nAxes_, nPhases_);
   out.m_PhaseIndices.fill(0); // all phases neutral
   for (size_t i = 0; i < MaxAxes; ++ i)
      out.m_TargetAxes[i] = i; // all axis coordinates mapped to themselves
   return out;
}


FAxisPermuflection FAxisPermuflection::CyclicPerm(int nShift_, unsigned nAxesToPermute_, unsigned nAxes_, unsigned nPhases_)
{
   cx_assert_rt(nAxes_ <= MaxAxes && nPhases_ <= MaxPhases);
   cx_assert_rt(nAxesToPermute_ <= nAxes_);
   FAxisPermuflection out(0, nAxes_, nPhases_);
   out.m_PhaseIndices.fill(0); // all phases neutral
   for (size_t i = 0; i < nAxesToPermute_; ++ i)
      // axis coord #i gets mapped to #(i + shift) [modulo nAxes]
      out.m_TargetAxes[i] = _integer_modulo<int>(i + nShift_, int(nAxesToPermute_));
   for (size_t i = nAxesToPermute_; i < MaxAxes; ++ i)
      // if any axes remain, assign identity permutations for them.
      out.m_TargetAxes[i] = i;
   return out;
}


FAxisPermuflection FAxisPermuflection::Reflections(unsigned long AxisMirrorFlags, unsigned nAxes_, unsigned nPhases_)
{
   cx_assert_rt(nAxes_ <= MaxAxes && nPhases_ <= MaxPhases);
   FAxisPermuflection out(0, nAxes_, nPhases_);

   for (size_t i = 0; i < MaxAxes; ++ i)
      out.m_TargetAxes[i] = i; // all axis coordinates mapped to themselves * phase factors
   out.m_PhaseIndices.fill(0); // start with all phases neutral

   if (nPhases_ % 2 != 0) {
      // we have an odd number of discrete phase angles. Actual reflections
      // cannot be expressed with those.
      if (AxisMirrorFlags == 0)
         // Not actually asked to reflect anything --> odd nPhases okay in this case.
         return out; // <-- should be an identity trafo at this point
      throw std::runtime_error("FAxisPermuflection::Reflections():"
         " cannot encode axis reflections with uneven number of discretized phase angles.");
   }
   phase_index_t iMirrorPhase = nPhases_ / 2;
   assert(_phase_factor(iMirrorPhase, nPhases_) == -1);

   // set up the mirror flags. One axis reflection for each bit set.
   // (note: m_PhaseIndices is already initialized to all-neutral above).
   for (size_t i = 0; i < nAxes_; ++ i)
      if (bool(AxisMirrorFlags & (1ul << i)))
         out.m_PhaseIndices[i] = iMirrorPhase;

   assert(nAxes_ <= CHAR_BIT * sizeof(AxisMirrorFlags));
   // ⬑ can "only" encode arbitrary axis flip combinations up to (bit
   // width of unsigned long combinations) number of axes. Yes, not likely
   // to trigger anytime soon.
   return out;
}


FAxisPermuflection FAxisPermuflection::i2QuarterTurn(unsigned iAxis, unsigned jAxis, unsigned nAxes_, unsigned nPhases_)
{
   cx_assert_rt(nAxes_ <= MaxAxes && nPhases_ <= MaxPhases);
   cx_assert_rt(iAxis != jAxis && iAxis < nAxes_ && jAxis < nAxes_);
   FAxisPermuflection out(0, nAxes_, nPhases_);

   // start with all coords x[k] mapped to themselves with neutral phases
   for (size_t k = 0; k < MaxAxes; ++ k)
      out.m_TargetAxes[k] = k; // x[k] ⟼ x[k] * α_k
   out.m_PhaseIndices.fill(0); // α_k = k

   if (nPhases_ % 2 != 0)
      // the C4 rotation requires a sign inversion. This requires an even number
      // of phase indices.
      throw std::runtime_error("FAxisPermuflection::i2QuarterTurn():"
         " cannot encode C4 rotation with uneven number of discretized phase angles.");
   phase_index_t iMirrorPhase = nPhases_ / 2;
   assert(_phase_factor(iMirrorPhase, nPhases_) == -1);

   // Update mapping of axes `i` and `j` for φ = 360°/4 = π/2 ccw rotation:
   //
   //   x[i] ⟼ (cos(φ) x[i] − sin(φ) x[j])
   //        ⟼ −x[j]
   //
   out.m_TargetAxes[iAxis] = jAxis;
   out.m_PhaseIndices[iAxis] = iMirrorPhase;
   //
   //   x[j] ⟼ (sin(φ) x[i] + cos(φ) x[j])
   //        ⟼ +x[i]
   //
   out.m_TargetAxes[jAxis] = iAxis;
   out.m_PhaseIndices[jAxis] = 0;
   // … and that's it already.
   return out;
#ifdef INCLUDE_ABANDONED
   // >>> import numpy as np
   // >>> def rg(i,N): f=2*np.pi*i/N; return np.array([np.cos(f), -np.sin(f), np.sin(f), np.cos(f)]).round(6)
   // >>> rg(1,4)
   // array([ 0., -1.,  1.,  0.])
#endif // INCLUDE_ABANDONED
}


#ifdef INCLUDE_ABANDONED
// bool FAxisPermuflection::operator < (FAxisPermuflection const &other) const {
//    return std::tie(this->m_TargetAxes, other.m_Signs) < std::tie(other.m_TargetAxes, this->m_Signs);
//    // ^-- reversed order of signs: we want operations containing sign flips
//    // behind non-reflecting ones, and this is comparing the signs
//    // lexicographically as actual numerical values.
// }
//
// bool FAxisPermuflection::operator == (FAxisPermuflection const &other) const {
//    return std::tie(this->m_TargetAxes,this->m_Signs) == std::tie(other.m_TargetAxes,other.m_Signs);
// }
#endif // INCLUDE_ABANDONED

bool FAxisPermuflection::operator < (FAxisPermuflection const &other) const {
   return std::tie(this->m_TargetAxes, this->m_PhaseIndices) < std::tie(other.m_TargetAxes, other.m_PhaseIndices);
   // ⬑ I think this should yield operations with the same permutation part as
   // adjacent to each other, with operations containing sign flips behind non-
   // reflecting ones
}

bool FAxisPermuflection::operator == (FAxisPermuflection const &other) const {
   return std::tie(this->m_TargetAxes, this->m_PhaseIndices) == std::tie(other.m_TargetAxes, other.m_PhaseIndices);
}



// FAxisPermuflection::FAxisEntry FAxisPermuflection::operator [] (unsigned iSourceAxis) const {
//    assert(iSourceAxis < nAxes());
//    return FAxisEntry{m_Signs[iSourceAxis], m_TargetAxes[iSourceAxis]};
// }
FAxisPermuflection::FAxisEntry FAxisPermuflection::operator [] (unsigned iSourceAxis) const {
   assert(iSourceAxis < nAxes());
   return FAxisEntry{m_PhaseIndices[iSourceAxis], m_TargetAxes[iSourceAxis], this->nPhases()};
}


// FAxisPermuflection FAxisPermuflection::multiply(FAxisPermuflection const &A, FAxisPermuflection const &B){
//    FPermOnly TargetAxes;
//    FSignsOnly Signs;
//    for (size_t i = 0; i < nAxes; ++i) {
//       // (A*B) x[i]
//       // = A (B x[i])
//       // = A(SignB[i] * x[PermB[i]])
//       // …set j := PermB[i]
//       // = SignB[i] * (A x[j])
//       // = SignB[i] SignA[j] * x[PermA[j]]
//       auto j = B.m_TargetAxes[i];
//       TargetAxes[i] = A.m_TargetAxes[j];
//       Signs[i] = A.m_Signs[j] * B.m_Signs[i];
//    }
//    return FAxisPermuflection(TargetAxes, Signs);
// }
//
//
// FAxisPermuflection FAxisPermuflection::inverse(FAxisPermuflection const &A)
// {
//    // (A*B) x[i] = SignB[i] SignA[PermB[i]] * x[PermA[PermB[i]]]
//    // --> for B to be the inverse of A:
//    //     + PermA[PermB[i]] == i must hold, so PermB must be set to the inverse
//    //       permutation of PermA
//    //     + SignB[i] SignA[PermB[i]] == 1 must hold, so we need:
//    //
//    //       SignB[i] == SignA[PermB[i]]   # use PermB = PermA^{-1}
//    //       SignB[i] == SignA[PermA^{-1}[i]]
//    //       SignB[PermA[i]] == SignA[i]
//    FPermOnly TargetAxes;
//    FSignsOnly Signs;
//    for (size_t i = 0; i < nAxes; ++i) {
//       auto j = A.m_TargetAxes[i];
//       TargetAxes[j] = i;
//       Signs[j] = A.m_Signs[i];
//    }
//    return FAxisPermuflection(TargetAxes, Signs);
// }

FAxisPermuflection FAxisPermuflection::multiply(FAxisPermuflection const &A, FAxisPermuflection const &B){
   FPermOnly TargetAxes;
   FPhaseIndicesOnly PhaseIndices;
   assert(A.nPhases() == B.nPhases() && A.nAxes() == B.nAxes());
   for (size_t i = 0; i < A.nAxes(); ++i) {
      // (A*B) x[i]
      // = A (B x[i])
      // = A(SignB[i] * x[PermB[i]])
      // … now set j := PermB[i]
      // = SignB[i] * (A x[j])
      // = SignB[i] SignA[j] * x[PermA[j]]
      auto j = B.m_TargetAxes[i];
      TargetAxes[i] = A.m_TargetAxes[j];
      PhaseIndices[i] = _integer_modulo<int>(A.m_PhaseIndices[j] + B.m_PhaseIndices[i], A.nPhases());
      assert(groups::index_valid_q(PhaseIndices[i], A.nPhases()));
   }
   return FAxisPermuflection(A.nAxes(), A.nPhases(), TargetAxes, PhaseIndices);
}


FAxisPermuflection FAxisPermuflection::inverse(FAxisPermuflection const &A)
{
   // (A*B) x[i] = SignB[i] SignA[PermB[i]] * x[PermA[PermB[i]]]
   // ⟹ for B to be the inverse of A:
   //     + PermA[PermB[i]] == i must hold, so PermB must be set to the inverse
   //       permutation of PermA
   //     + SignB[i] SignA[PermB[i]] == 1 must hold, so we need:
   //
   //       SignB[i] == SignA[PermB[i]]   # use PermB = PermA^{-1}
   //       SignB[i] == SignA[PermA^{-1}[i]]
   //       SignB[PermA[i]] == SignA[i]
   FPermOnly TargetAxes;
   FPhaseIndicesOnly PhaseIndices;
   for (size_t i = 0; i < A.nAxes(); ++i) {
      auto j = A.m_TargetAxes[i];
      TargetAxes[j] = i;
      PhaseIndices[j] = _integer_modulo(-int(A.m_PhaseIndices[i]), int(A.nPhases()));
   }
   return FAxisPermuflection(A.nAxes(), A.nPhases(), TargetAxes, PhaseIndices);
}


#if USE_EIGEN_PERM_DIAG_MATRICES_FOR_PERMFLECTIONS

template<class T>
FAxisPermuflection::TPermDiagMatrixPair<T> FAxisPermuflection::_MakePermDiagMatrixPairT() const
{
   typedef Eigen::Matrix<FPermOnly::value_type,nAxes,1> TIndexVector;
   typedef Eigen::Matrix<FSignsOnly::value_type,nAxes,1> TSignVector;
   auto MappedIndices = Eigen::Map<TIndexVector, Eigen::Unaligned>(&m_TargetAxes[0]);
   auto MappedSigns = Eigen::Map<TSignVector, Eigen::Unaligned>(&m_Signs[0]);

   return {FAxisPermMatrix(MappedIndices), TAxisSignMatrix<T>(MappedSigns)};
   // Note: These here use `PermutationMatrix` and `DiagonalMatrix`, which copy the
   //     respective vector data. But there's also `DiagonalWrapper` and
   //     `PermutationWrappper`. However, for those the docs seem a bit fishy at
   //     the moment, and I anyway wonder whether it is a fantastic idea to make
   //     matrix referencing into char arrays…
}

// template<class T>
// TRotationMatrix<T> FAxisPermuflection::AsMatrixT() const
// {
//    auto em = _MakePermDiagMatrixPairT<T>();
//    return em.P * em.D;
// }
//  ⬑ interesting: apparently one can multiply diagonal and perm matrices
//     to general matrices — … but not to each other.


// forms the operator product of R * X, both as interpreted as actions on real-space vectors
// expressed in Cartesian coordinate column-vectors.
template<class T>
FRotationMatrix operator * (TRotationMatrix<T> const &R, FAxisPermuflection const &X)
{
   auto em = X._MakePermDiagMatrixPairT<int>();
   return (R * em.P) * em.D;
}

template<class T>
FRotationMatrix operator * (FAxisPermuflection const &X, TRotationMatrix<T> const &R)
{
   auto em = X._MakePermDiagMatrixPairT<int>();
   return em.P * (em.D * R);
}

#else

template<class T>
TRotationMatrix<T> FAxisPermuflection::AsMatrixT() const
{
   FAxisPermuflection const &X = *this;
   size_t const N = X.nAxes();
   auto const _phase_factor = _PhaseFactorFn_Get(X.nPhases());
   TRotationMatrix<T> out = TRotationMatrix<T>::Zero();
   assert(size_t(out.rows()) == N && size_t(out.cols()) == N);
   for (size_t i = 0; i < N; ++ i)
//       out(m_TargetAxes[i],i) = m_Signs[i];
//       out(m_TargetAxes[i],i) = T(_phase_factor(m_PhaseIndices[i]));
      out(X.TargetAxisQ(i), i) = T(_phase_factor(X.PhaseIndexQ(i)));
   return out;
}


FRotationMatrix FAxisPermuflection::AsMatrix() const {
   return this->AsMatrixT<FRotationMatrix::Scalar>();
}


// forms the operator product of R * X, both as interpreted as actions on real-space vectors
// expressed in Cartesian coordinate column-vectors.
template<class T>
TRotationMatrix<T> operator * (TRotationMatrix<T> const &R, FAxisPermuflection const &X)
{
//    size_t constexpr N = FAxisPermuflection::nAxes;
   size_t const N = X.nAxes();
   auto const _phase_factor = _PhaseFactorFn_Get(X.nPhases());
//    auto const _multiply_by_phase_factor = _PhaseFactorFn_MultiplyArg<T>(X.nPhases());
   assert(size_t(R.cols()) == N);
   // (R X) v = R (X v); ⟹ Must make R X act on vectors v in such a way
   // that coincides with the sequence:
   //
   //    (apply X to input vector) ⟹ (use as input to regular R).
   //
   //    temp[p[i]] = s[i]*inp[i]  ⟹ out[a] = \sum_j R[a,i] temp[i]
   //
   // So we need to apply X to the input-linked axis (columns) of R.
   TRotationMatrix<T> out;
   for (size_t i = 0; i < N; ++i)
      for (size_t a = 0; a < size_t(R.rows()); ++a)
//          out(a, X.TargetAxisQ(i)) = R(a,i) * X.SignQ(i);
         out(a, X.TargetAxisQ(i)) = R(a,i) * _phase_factor(X.PhaseIndexQ(i));
//          out(a, X.TargetAxisQ(i)) = _multiply_by_phase_factor(R(a,i), X.PhaseIndexQ(i));
   return out;
}


template<class T>
TRotationMatrix<T> operator * (FAxisPermuflection const &X, TRotationMatrix<T> const &R)
{
//    size_t constexpr N = FAxisPermuflection::nAxes;
   size_t const N = X.nAxes();
   auto const _phase_factor = _PhaseFactorFn_Get(X.nPhases());
//    auto const _multiply_by_phase_factor = _PhaseFactorFn_MultiplyArg<T>(X.nPhases());
   assert(size_t(R.rows()) == N);
   // (X R) v = X (R v) ⟹ apply X to output-linked axis (rows) of R
   TRotationMatrix<T> out;
   for (size_t a = 0; a < N; ++a)
      for (size_t i = 0; i < size_t(R.cols()); ++i)
//          out(X.TargetAxisQ(a), i) = R(a,i) * X.SignQ(a);
         out(X.TargetAxisQ(a), i) = _phase_factor(X.PhaseIndexQ(a)) * R(a,i);
//          out(a, X.TargetAxisQ(i)) = _multiply_by_phase_factor(R(a,i), X.PhaseIndexQ(a));
   return out;
}

#endif // USE_EIGEN_PERM_DIAG_MATRICES_FOR_PERMFLECTIONS

