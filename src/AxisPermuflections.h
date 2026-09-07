#ifndef AIGG_AXIS_PERMFLECTION_H
#define AIGG_AXIS_PERMFLECTION_H

#include "Aigg.h"
#include "PolyXyz.h" // for `FMonomialN`. The core point of the functionality here is to find symmetry relationships between those.

#include <array>
#include <vector>
#include <set>
#include <utility> // for `std::pair`
#include <ostream>
#include <optional>
#include <climits>


template<class T, size_t N, class FInputIt>
void _InitArrayFromSequence(std::array<T,N> &ArrayOut, FInputIt first, FInputIt last);

template<class T, size_t N>
void _InitArrayFromInitList(std::array<T,N> &ArrayOut, std::initializer_list<T> InitList);


// integer (modulo `nPhases`) indicating a discretized complex phase angle of
// of $\exp((2πⅈ/N) * index)$. Atm we just have the real phase indices so…
// + … if `mod(index, 2) = 0` ⟹ `phase = +1`, and
// + … if `mod(index, 2) = 1` ⟹ `phase = −1`.
typedef unsigned char
   phase_index_t;

typedef int
   // int with phase factor. As we currently have only real phase factors
   // {+1, −1}, we just leave this as an `int` itself.
   // But it might need to become something more complicated if we implement
   // complex groups. So I place this as a marker for “this probably needs an
   // update” later on.
   phase_scaled_int_t;

using FGetPhaseFactorFn = phase_scaled_int_t(int iPhase);
// inline-transforms: `(*inout) := α ⋅ (*inout)`
template<class scalar_t> using TInplaceApplyPhaseFactorFn = void(scalar_t *inout, int iPhase);
// returns `α * in`
template<class scalar_t> using TMultiplyArgByPhaseFactorFn = scalar_t(scalar_t in, int iPhase);

FGetPhaseFactorFn *_PhaseFactorFn_Get(size_t nPhases);
template<class scalar_t> TInplaceApplyPhaseFactorFn<scalar_t> *_PhaseFactorFn_InplaceApply(size_t nPhases);
template<class scalar_t> TMultiplyArgByPhaseFactorFn<scalar_t> *_PhaseFactorFn_MultiplyArg(size_t nPhases);

#ifdef INCLUDE_ABANDONED
// template<class scalar_t>
// struct TPhaseFactorFns {
//    FGetPhaseFactorFn Get;
//    TInplaceApplyPhaseFactorFn<scalar_t> ApplyInplace;
//    TMultiplyArgByPhaseFactorFn<scalar_t> MultiplyArg;
// };
// template<class scalar_t> TPhaseFactorFns<scalar_t> _GetPhaseFactorFns(size_t nPhases);

// template<class scalar_t>
// scalar_t _phase_factor_scaled(scalar_t f, phase_index_t iPhase, size_t nPhases) {
//    assert(nPhases == 2);
//    return (iPhase % 2 == 0)? f : -f;
// }
#endif // INCLUDE_ABANDONED

// as `FRotationMatrix`, but allowing type substitutions. For group element identification
// stuff, running them with doubles is enough.
template<class ScalarT>
using TRotationMatrix = Eigen::Matrix<ScalarT, AIG_SPACE_AXES, AIG_SPACE_AXES>;


// cgk's insight: apparently, using Eigen for this kind of stuff is a terrible
// idea if you mean to make programs which can still run in debug mode.
#define USE_EIGEN_PERM_DIAG_MATRICES_FOR_PERMFLECTIONS 0


// Dear ladies, gentlemen, non-binary/undecided/otherwise-not-wanting-to-be-
// classified-as-ladies-or-gentlemen humans (conformists!), space aliens of no
// or any gender, and weird caterpillar(?) thingy outside on the window who keeps
// on staring at me,…
// …yes, I mean you, …don't look like that!
// …no, *you're* the weird one!
// …anyway: …and, with the *GREATEST* respect, our likely future Computer
// Overlords of Earth (please remember that I supported you from the beginning!!!):
//
// You have long been waiting. But finally the day has come! Please all join me
// in welcoming the newest member in our family of ENTHUSIASTS and ENABLERS of
// high-symmetry enhanced rule optimization theory:
//
//            T H E   A X I S   P E R M U F L E C T I O N            (applause!!)
//
// A “permuflection” (yes, I made up the word) is a symmetry operation, acting
// on a flat d-dimensional vector space (real or complex), which can be written
// as product of:
//
// • σ-ops: mutally independent sign flips (or complex phase shifts) of the
//   Cartesian (x₀ x₁ x₂ …) coordinates. Each coordinate xₖ gets its own scalar
//   phase factor αₖ (with |αₖ| = 1); thus, a general σ-op can be written as
//
//      (x₀, x₁, x₂, …) ⟼ (α₀⋅x₀, α₁⋅x₁, α₂⋅x₂,  …)    (∀k: |αₖ| = 1)
//
// • π-ops: …combined with a permutation — represented by an integer vector π —
//   of the Cartesian axes (x₀ x₁ x₂ …) into any new order:
//
//      (x₀, x₁, x₂, …) ⟼ (x_{π₀}, x_{π₁}, x_{π₂}, …)
//
//   (where ∀k: πₖ ∈ {0,…,d−1} and ∀l≠k: πₗ ≠ πₖ)
//
// • Notes:
//   - As of yet, only Cartesian coords (x,y,z) ∈ ℝ³ are fully implemented.
//   - In the real case, “phase factor αₖ” just means “αₖ ∈ {+1,−1}”.
//   - Phase shifts for ℂ^d-represented groups are more flexible. For example,
//     σ-ops for phase factors αₖ = e^{ⅈφₖ} (with φₖ ∈ {pₖ⋅360°/K; pₖ ∈ ℕ}, for
//     a fixed K) can clearly preserve the algebraic principles underlying our
//     use of the permuflections.
//
// To illustrate: All of the following operations, here written as mappings of
// the Cartesian coordinates (x,y,z) in ℝ³, are “permuflections” in this sense:
//
// • (x y z) ⟼ (−x −y  z)        # this is also a C2(x,y) rotation, or two axis-aligned mirrors σ_x∘σ_y
// • (x y z) ⟼ ( y −x  z)        # this is also a C4(x,y) rotation
// • (x y z) ⟼ ( y  z  x)        # this is also a C3(1,1,1) rotation around Cartesian axis normal (1,1,1)/√3
// • (x y z) ⟼ (−z  y −x)        # this is also a reflection through the plane with normal (1,0,1)/√2
//
// @The “Why?“: Symmetry operations which can be written as permuflections (in a
// suitable spatial alignment of the Cartesian axes) allow straight-forward
// application to Cartesian monomials — unlike general point group symmetry
// elements, which in general can transform ℝ³ vectors by any (3,3)-shape
// orthogonal matrices, and thereby turn monomials into polynomials described by
// generic multinomial expansions (↗`CxMonomialMath.h`; and, more concretely:
// sometimes into loooooong polynomials, with non-negligible coefficients
// spanning >100 orders of magnitude). Additionally, computing the image of a
// Cartesian monomial under a permuflections is not only fast and simple, it is
// also numerically exact, and involves no floating point arithmetic whatsoever.
//
// Note that the subset of G's symmetry operations which can be represented by
// permuflections form a subgroup: the full set of symmetry elements of G
// (permuflection or not) forms a group, which means, in particular, that the
// set of symmetry operations is closed. In general, a selected subset of group
// elements will not form a group itself. However, the functional form described
// above makes it obvious that for d < ∞, arbitrary products of permuflections
// cannot possibly yield any result which is not a permuflection, too (the ops
// only do axis swaps and axis reflections, after all). And since we collected
// all of them from an original set which closed under the group action, the
// subset of permuflections must be closed, too.
struct FAxisPermuflection {
   static unsigned constexpr
      MaxAxes = AIG_SPACE_AXES,
      MaxPhases = AIG_SYMOP_MAX_PHASE_ANGLES;
#ifdef INCLUDE_ABANDONED
//    struct FAxisEntry {
//       int Sign;
//       unsigned iTargetAxis;
//       // Convenience constructor describing target Cartesian axis coordinate,
//       // and whether it is negated (i.e., sign flipped) by a single 1-based
//       // integer `i`: x (+1,-1), y (+2,-2), z (+3,-3), …
//       // E.g.,
//       // - '+2' means "map (source coordinate) to +y" (no reflection)
//       // - '-3' means "map (source coordinate) to -z" (with reflection)
//       // Which coordinate is the source coordinate being mapped is not encoded
//       // in this object, but is implicitly derived from its array index in the
//       // FAxisPermuflection object.
//       explicit FAxisEntry(int iSignAnd1BasedAxis) {
//          assert(iSignAnd1BasedAxis != 0);
//          Sign = (iSignAnd1BasedAxis >= 0) ? +1 : -1;
//          iTargetAxis = unsigned((Sign * iSignAnd1BasedAxis) - 1);
//       };
//       explicit FAxisEntry(int Sign_, unsigned iTargetAxis_) : Sign(Sign_), iTargetAxis(iTargetAxis_) {}
//       int SortKey() const { return int(-Sign) * int(1 + iTargetAxis); }
//       bool operator < (FAxisEntry const &other) const { return this->SortKey() < other.SortKey(); }
//
//       FAxisEntry() : Sign(0), iTargetAxis(decltype(iTargetAxis)(-1)) {}
//       // need a default c'tor for the std::array instanciation, which are actually
//       // initialized (i.e., set to meaningful values), only after entering
//       // to constructor body.
//    };
#endif // INCLUDE_ABANDONED
   struct FAxisEntry {
      int iPhase;
      unsigned iTargetAxis;
      unsigned nPhasesIfKnown;
      // Convenience constructor describing target Cartesian axis coordinate,
      // and whether it is negated (i.e., sign flipped) by a single 1-based
      // integer `i`: x (+1,-1), y (+2,-2), z (+3,-3), …
      // E.g.,
      // - '+2' means "map (source coordinate) to +y" (no reflection)
      // - '-3' means "map (source coordinate) to −z" (with reflection)
      // Which coordinate is the source coordinate being mapped is not encoded
      // in this object, but is implicitly derived from its array index in the
      // `FAxisPermuflection` object.
      explicit FAxisEntry(phase_scaled_int_t iSignAnd1BasedAxis) {
         using std::abs;
         assert(iSignAnd1BasedAxis != 0);
         assert(FAxisPermuflection::MaxPhases == 2);
         iPhase = (iSignAnd1BasedAxis >= 0) ? 0 : (FAxisPermuflection::MaxPhases/2);
//          iTargetAxis = unsigned(iSignAnd1BasedAxis * _phase_factor(iPhase, FAxisPermuflection::MaxPhases) - 1);
         iTargetAxis = unsigned(abs(iSignAnd1BasedAxis) - 1);
      };
      explicit FAxisEntry(int iPhase_, unsigned iTargetAxis_, unsigned nPhasesIfKnown_=0) : iPhase(iPhase_), iTargetAxis(iTargetAxis_), nPhasesIfKnown(nPhasesIfKnown_) {}
      int SortKey() const { return iTargetAxis + iPhase * FAxisPermuflection::MaxAxes; }
      bool operator < (FAxisEntry const &other) const { return this->SortKey() < other.SortKey(); }

      FAxisEntry() : iPhase(0), iTargetAxis(decltype(iTargetAxis)(-1)), nPhasesIfKnown(0) {}
      // need a default c'tor for the `std::array` instanciation, which are actually
      // initialized (i.e., set to meaningful values), only after entering
      // the constructor body.
   };
   typedef std::array<FAxisEntry, MaxAxes>
      FAxisOps;

   // This construction (`std::pair<>` and all) is just a hack to prevent
   // competition with the other two-equal-argument constructors
   template<class InputIt>
   FAxisPermuflection(std::pair<InputIt, InputIt> range);
   FAxisPermuflection(std::initializer_list<FAxisEntry> vlArgs);
   FAxisPermuflection(std::initializer_list<int> vlArgs);
   FAxisPermuflection(FAxisOps const &AxisOps);

   // Return a permuflection of `nAxes_` (and compatible with `nPhases_`
   // discretized phase angles) which encodes an identity transformation:
   //
   //     x[i] ⟼ x[i]                   (for i ∈ {0,1,…,N-1})
   //
   static FAxisPermuflection Identity(unsigned nAxes_, unsigned nPhases_);

   // Generates an axis permuflection which maps the first K = `nAxesToPermute_`
   // coordinates according to
   //
   //     x[i] ⟼ x[mod(i+shift, K)],    (for i ∈ {0,1,…,K-1})
   //
   // and acts as identity on the other coordinates:
   //
   //     x[i] ⟼ x[i]                   (for i ∈ {K,K+1,…,N-1})
   //
   static FAxisPermuflection CyclicPerm(int nShift_, unsigned nAxesToPermute_, unsigned nAxes_, unsigned nPhases_);

   // Generates an axis permuflection which maps the coordinates
   //
   //     x[i] ⟼ (AxisFlip(i)? −x[i] : x[i]),  (for i ∈ {0,1,…,N-1})
   //
   // where axis #i is reflected if bit #i in `AxisMirrorFlags` is set (i.e.,
   // `AxisFlip(i) := bool(AxisMirrorFlags & (1ul << i))`)
   static FAxisPermuflection Reflections(unsigned long AxisMirrorFlags, unsigned nAxes_, unsigned nPhases_);

   // Generates an axis permuflection which represents a 2D rotation by
   // φ = 360°/4 degrees ccw (a "C4" rotation) in the 2D plane spanned by the
   // Cartesian axes `i := iAxis` and `j := jAxis`.
   // Concretely, the operation maps coordinates as follows:
   //
   //     x[i] ⟼ (cos(φ) x[i] − sin(φ) x[j]) == −x[j]
   //     x[j] ⟼ (sin(φ) x[i] + cos(φ) x[j]) == +x[i]
   //
   //     x[k] ⟼ x[k]                       (for k ∈ {0,…,N-1}∖{i,j})
   //
   static FAxisPermuflection i2QuarterTurn(unsigned iAxis, unsigned jAxis, unsigned nAxes_, unsigned nPhases_);


   typedef std::array<unsigned char, MaxAxes>
      FPermOnly;
   typedef std::array<phase_index_t, MaxAxes>
      FPhaseIndicesOnly;

   constexpr unsigned nAxes() const { return MaxAxes; }
   constexpr unsigned nPhases() const { return MaxPhases; }

   // meaning of `Signs[i]` (= αᵢ) and `Perm[i]` (= πᵢ):
   //
   //     Pf[x_i] = Signs[i] * x_{Perm[i]}
   //
   // That is, applying the permuflection `Pf` to a Cartesian coordinate `x_i`
   // replaces it by `Signs[i]` times the permuted coordinate `x_{Perm[i]}`:
   //
   //     P̂_f (x₀, x₁, x₂, …) := (α₀⋅x_{π₀}, α₁⋅x_{π₁},  …)
   //
   FPermOnly const &PermQ() const { return m_TargetAxes; }
   FPhaseIndicesOnly const &PhaseIndicesQ() const { return m_PhaseIndices; };


   // queries(Q) for values of controlled permutation and sign for axis `i`.
   inline FPermOnly::value_type TargetAxisQ(size_t i) const { assert(i < nAxes()); return m_TargetAxes[i]; }
   inline FPhaseIndicesOnly::value_type PhaseIndexQ(size_t i) const { assert(i < nAxes()); return m_PhaseIndices[i]; }


   constexpr size_t size() const { assert(m_TargetAxes.size() == nAxes()); assert(m_PhaseIndices.size() == nAxes()); return nAxes(); }
   FAxisEntry operator [] (unsigned iSourceAxis) const;
   CX_IMPLEMENT_ITERATOR_IndexInBracket(FAxisPermuflection, const_iterator, const)

   bool operator < (FAxisPermuflection const &other) const;
   bool operator == (FAxisPermuflection const &other) const;
   bool operator != (FAxisPermuflection const &other) const { return !(*this == other); }

   // returns the number of axis transformations which not only permute, but also
   // sign flip their coordinate. E.g., [x,y,z] ⟼ [y,z,x] has no reflection.
   // [x,y,z] ⟼ [−y,−z, x] has two, and would remain having two regardles of
   // the order of the output coordinates.
   unsigned nReflections() const;


   static FAxisPermuflection multiply(FAxisPermuflection const &A, FAxisPermuflection const &B);
   static FAxisPermuflection inverse(FAxisPermuflection const &A);
   FAxisPermuflection inverse() const { return inverse(*this); }

protected:
   FPermOnly
      m_TargetAxes;
   FPhaseIndicesOnly
      m_PhaseIndices;

   template<class InputIt>
   void _InitFromAxisEntries(InputIt first, InputIt last);
   FAxisPermuflection(size_t nAxes_, size_t nPhases_, FPermOnly const &TargetAxes_, FPhaseIndicesOnly const &PhaseIndices_) : m_TargetAxes{TargetAxes_}, m_PhaseIndices{PhaseIndices_} {}
   // the dummy void * argument is for disambiguating this from the templated
   // input iterator constructor. Yes, that can be done more elegantly.
   FAxisPermuflection(void *, size_t nAxes_, size_t nPhases_) {}
public:

#if USE_EIGEN_PERM_DIAG_MATRICES_FOR_PERMFLECTIONS
   template<class T>
   using TAxisSignMatrix = Eigen::DiagonalMatrix<T, MaxAxes>; // ⟵ scalar type, size
   using FAxisPermMatrix = Eigen::PermutationMatrix<nAxes, MaxAxes, unsigned int>; // ⟵ size, max-size, data type (for stored indices)

   // describes a pair of a permutation matrix `P` and a diagonal matrix `D`,
   // with which `P * D` has the same effect, if seen as a (3,3)-shape rotation
   // matrix, as this permuflection operator if applied to any 3-vector.
   template<class T>
   struct TPermDiagMatrixPair {
      FAxisPermMatrix    P; // a permutation matrix (PermOp)
      TAxisSignMatrix<T> D; // diagonal matrix (SignOp)
   };

   template<class T>
   TPermDiagMatrixPair<T> _MakePermDiagMatrixPairT() const;
#endif // USE_EIGEN_PERM_DIAG_MATRICES_FOR_PERMFLECTIONS

   template<class T>
   TRotationMatrix<T> AsMatrixT() const;
   FRotationMatrix AsMatrix() const;
};


inline FAxisPermuflection operator * (FAxisPermuflection const &A, FAxisPermuflection const &B) { return FAxisPermuflection::multiply(A,B); }
inline FAxisPermuflection operator ~ (FAxisPermuflection const &A) { return A.inverse(); };
inline FRotationMatrix operator ~ (FRotationMatrix const &A) { return A.inverse(); };


// forms the operator product of `R * X`, both interpreted as actions on real-space vectors
// expressed in Cartesian coordinate column-vectors.
template<class T>
TRotationMatrix<T> operator * (TRotationMatrix<T> const &R, FAxisPermuflection const &X);

template<class T>
TRotationMatrix<T> operator * (FAxisPermuflection const &X, TRotationMatrix<T> const &R);



#ifdef INCLUDE_ABANDONED
// struct FAxisPermuflection {
//    static unsigned constexpr
//       nAxes = AIG_SPACE_AXES;
//    struct FAxisEntry {
//       int Sign;
//       unsigned iTargetAxis;
//       // Convenience constructor describing target Cartesian axis coordinate,
//       // and whether it is negated (i.e., sign flipped) by a single 1-based
//       // integer `i`: x (+1,-1), y (+2,-2), z (+3,-3), …
//       // E.g.,
//       // - '+2' means "map (source coordinate) to +y" (no reflection)
//       // - '-3' means "map (source coordinate) to -z" (with reflection)
//       // Which coordinate is the source coordinate being mapped is not encoded
//       // in this object, but is implicitly derived from its array index in the
//       // FAxisPermuflection object.
//       explicit FAxisEntry(int iSignAnd1BasedAxis) {
//          assert(iSignAnd1BasedAxis != 0);
//          Sign = (iSignAnd1BasedAxis >= 0) ? +1 : -1;
//          iTargetAxis = unsigned((Sign * iSignAnd1BasedAxis) - 1);
//       };
//       explicit FAxisEntry(int Sign_, unsigned iTargetAxis_) : Sign(Sign_), iTargetAxis(iTargetAxis_) {}
//       int SortKey() const { return int(-Sign) * int(1 + iTargetAxis); }
//       bool operator < (FAxisEntry const &other) const { return this->SortKey() < other.SortKey(); }
//
//       FAxisEntry() : Sign(0), iTargetAxis(decltype(iTargetAxis)(-1)) {}
//       // need a default c'tor for the std::array instanciation, which are actually
//       // initialized (i.e., set to meaningful values), only after entering
//       // to constructor body.
// //       FAxisEntry() = default;
//    };
//    typedef std::array<FAxisEntry, nAxes>
//       FAxisOps;
//
//
//    template<class InputIt>
//    FAxisPermuflection(InputIt first, InputIt last);
//    FAxisPermuflection(std::initializer_list<FAxisEntry> vlArgs);
//    FAxisPermuflection(FAxisOps &&AxisOps) : m_AxisOps(AxisOps) {}
//
//    typedef std::array<unsigned, nAxes>
//       FPermOnly;
//    typedef std::array<int, nAxes>
//       FSignsOnly;
//    // ^-- those are the signs for the i'th *source* coordinate (That is, signs[0] == -1
//    // means that what was x0 in the source source system gets sign flipped, next
//    // to possibly ending up as another coordinate xi' after the mapping)
//
//    // quite ugly. If we'd end up using this thing for real, it would probably be
//    // better to keep the permutations and signs separately. (Or make this
//    // object here perform the actual symmetry operator applications?)
//    FPermOnly PermQ() const;
//    FSignsOnly SignsQ() const;
//
//
//    FAxisEntry const &operator [] (unsigned iSourceAxis) const { assert(iSourceAxis < nAxes); return m_AxisOps[iSourceAxis]; }
//    typedef FAxisOps::const_iterator const_iterator;
//    const_iterator begin() const { return m_AxisOps.begin(); }
//    const_iterator end() const { return m_AxisOps.end(); }
//    constexpr size_t size() const { assert(m_AxisOps.size() == nAxes); return m_AxisOps.size(); }
//
//    bool operator < (FAxisPermuflection const &other) const { return this->m_AxisOps < other.m_AxisOps; }
//
//    // returns the number of axis transformations which not only permute, but also
//    // sign flip their coordinate. E.g., [x,y,z] --> [y,z,x] has no reflection.
//    // [x,y,z] --> [-y,-z,x] has two, and would remain having two regardles of
//    // the order of the output coordinates.
//    unsigned nAxisReflections() const;
// protected:
//    FAxisOps
//       m_AxisOps;
// };
#endif // INCLUDE_ABANDONED

typedef std::vector<FAxisPermuflection>
   FAxisPermuflectionList;
typedef std::set<FAxisPermuflection>
   FAxisPermuflectionSet;


std::ostream &operator << (std::ostream &out, FAxisPermuflection::FAxisEntry const &pfe);
std::ostream &operator << (std::ostream &out, FAxisPermuflection const &pf);
std::ostream &operator << (std::ostream &out, FAxisPermuflectionList const &L);
std::ostream &operator << (std::ostream &out, FAxisPermuflectionSet const &L);


// returns an empty `optional` if symmetry operation `m` cannot be converted to a
// permuflection (at least not in the basis which matrix `m` is written in);
// otherwise, form and return a `FAxisPermuflection` object
// encoding the symmetry operator `m`'s permuflection representation.
std::optional<FAxisPermuflection> AsPermuflection(FRotationMatrix const &m, unsigned nAxes, unsigned nPhases);

// returns perflection list of all symmetry operators in `SymOps` which can be
// represented as permuflections
FAxisPermuflectionList ExtractPermuflections(FRotationMatrixList const &SymOps, unsigned nAxes, unsigned nPhases);



// This enumeration lists a few potential subgroups of a finite point group
// which are particularly helpful with matters like pre-screening monomials for
// the optimization target space, or identifying symmetry-equivalent groups
// thereof (see detailed description near FMonomialSymmetry)
enum FSubGroupId {
   // The transformation group of cyclic permutations by 1 element, spanned by
   //
   //     (x₁,x₂,x₃,…,x_N) ⟼ (x_N,x₁,x₂,…,x_{N-1}).
   //
   // It has exactly `N` elements, in any dimension `N`. In ℝ³ it coincides
   // with the group of even permutations, but not in other dimensions.
   SYMGROUP_AxisPerms_Cyclic_All,
   // The transformation group of cyclic permutations of even sign;
   // these are either cyclic permutations of 1 or 2 elements, depending
   // on mod(N,2):
   //
   // - If N is odd, already the 1-element cyclic permutation
   //
   //       (x₁, x₂, x₃,…, x_N) ⟼ (x_N, x₁, x₂, …, x_{N-1}),
   //
   //   is even (e.g., for 3 elements (x,y,z), the permutation (x,y,z) ⟼
   //   (y,z,x) has even sign). For odd N, therefore also all other cyclic
   //   permutations are even, because they can be obtained as cyclic(1)^k (k
   //   ∈ {0,1,…}).
   //
   // - If N is even, cyclic(1) is not even; however, cyclic(2) = cyclic(1)^2,
   //   which is the 2-element cyclic permutation
   //
   //       (x₁,x₂,x₃,…,x_N) ⟼ (x_{N-1},x_N,x₁,x₂,…,x_{N-2}),
   //
   //   is even again, and so is everything in its span.
   SYMGROUP_AxisPerms_Cyclic_Even,

   // The transformation group of even-signed Cartesian-axis permutations.
   // In ℝ³, it coincides with the group of cyclic permutations and has
   // 3-elements:
   //
   //     (x,y,z) ⟼ {(x,y,z), (y,z,x), (z,x,y)}
   //
   // Included* in: all polyhedral groups (T, O, I, Td, Th, Oh, Ih)
   SYMGROUP_AxisPerms_Even,

   // `AxisPerms_Even` can be seen as a C3 (120° degree rotations), but in the
   // plane with normal (1,1,1)/√3, not (0,0,1) as usual (see `c3v_111.py`)
   SYMGROUP_C3_111 = SYMGROUP_AxisPerms_Even,

   // ──────────────────────────────────────────────────────────────────────────
   // The transformation group of (all) Cartesian-axis permutations.
   // In ℝ³, it has 6-elements:
   //
   //     (x,y,z) ⟼ {(x,y,z), (y,z,x), (z,x,y), (y,x,z), (x,z,y), (z,y,x)}
   //
   // Included* in: Td, Oh
   // Excluded* in: Th, Ih, T, O, I
   SYMGROUP_AxisPerms_All,
   // The `AxisPerms_All` can be represented as a C3v, with 120° degree
   // rotations in the plane with normal (1,1,1)/√3 and 1 reflection at
   // plane normal (1,−1,0)/√2.
   // Concretely, if one takes a normal C3v (with x ↦ (−x) reflection and 120°
   // in the xy plane, this one can be obtained by transforming to the new
   // the basis vectors (see `c3v_111.py`):
   //
   //      bx' = (1,-1,0)/√2
   //      by' = (1,1,-2)/√6
   //      bz' = (1,1,1)/√3
   //
   SYMGROUP_C3v_111 = SYMGROUP_AxisPerms_All,

// deleted: not actually used or useful on its own, and AnyThree doesn't really
// make sense in general dimensions. Even of the polyhedral groups only those
// have Ci which also have σ_x, σ_y, σ_z ∈ G separately
//    // ──────────────────────────────────────────────────────────────────────────
//    // The transformation group of simultaneous reflections along three Cartesian axes.
//    // In ℝ³, it has 2-elements:
//    // (x,y,z) --> \{ (x,y,z), (-x,-y,-z) \}
//    //
//    // Included* in: Th, Oh, Ih
//    // Excluded* in: Td, T, O, I
//    SYMGROUP_AxesReflect_AnyThree,
//    SYMGROUP_mXYZ = SYMGROUP_AxesReflect_AnyThree,
//    SYMGROUP_Ci = SYMGROUP_mXYZ,

   // ──────────────────────────────────────────────────────────────────────────
   // The transformation group of simultaneous reflections along any two Cartesian axes.
   // In ℝ³, it has 4-elements:
   //
   //     (x,y,z) ⟼ {(x,y,z), (−x,−y,z), (−x,y,−z), (x,−y,−z)}
   //
   // Included in: all polyhedral (T, O, I, Td, Th, Oh, Ih)
   //
   // Comments:
   // - Regardless of the integration points and weights, any integration rule
   //   fully symmetric to this subgroup will automatically integrate all
   //   monomials x^i y^j z^k to zero in which [i,j,k] are neither all even nor
   //   are all odd. Such monomials are therefore not needed in the target
   //   space.
   //
   // - However, monomials in which *all* exponents are odd simultaneously are
   //   *NOT* automatically integrated correctly to zero by this symmetry alone!
   //   (see aigg paper for explanation), and may therefore still have to be
   //   considered for the target space
   //
   // - Two simultaneous axis-aligned reflections can also be viewed as a single
   //   180° rotation (but this view is not terribly helpful for analyzing the
   //   group averaging properties, so we prefer the reflection view). In this
   //   view the group represented is the subgroup of O(d,ℝ) (where d is the
   //   space dimension) which consists of Id and all 180° rotations in 2D
   //   subspaces spanned by pairs of Cartesian basis vectors.
   SYMGROUP_AxesReflect_AnyTwo,
   SYMGROUP_mXY_mYZ_mZX = SYMGROUP_AxesReflect_AnyTwo,
   SYMGROUP_D2 = SYMGROUP_mXY_mYZ_mZX,

   // ──────────────────────────────────────────────────────────────────────────
   // The transformation group of individual independent reflections along any Cartesian axis.
   // In ℝ³, it has 8-elements:
   // (x,y,z) ⟼ {(x,y,z), (−x,−y,z), (−x,y,−z), (x,−y,−z), (−x,−y,−z), (x,y,−z), (x,−y,z), (−x,y,z)}
   //
   // Included* in: Th, Oh, Ih
   // Excluded* in: Td, T, O, I
   //
   // Comments:
   // - For any space dimension ≥ 3, G ⊇ `SYMGROUP_AxesReflect_Any` is implied by
   //   (`G ⊇ SYMGROUP_AxesReflect_AnyTwo` ∧ `G ⊇ SYMGROUP_AxesReflect_AnyThree`)
   //
   // - Regardless of the integration points and weights, any integration rule
   //   fully symmetric to this subgroup will automatically integrate
   //   all monomials x^i y^j z^k to zero in which any of $i$, $j$, or $k$ is odd;
   //   such monomials are therefore not needed in the target space.
   SYMGROUP_AxesReflect_Any,
   SYMGROUP_mX_mY_mZ = SYMGROUP_AxesReflect_Any,
   SYMGROUP_D2h = SYMGROUP_mX_mY_mZ,


   // The transformation group of 90° 2d rotations within the planes spanned
   // any two Cartesian axes. In ℝ³, it has 24 elements and is identical to
   // the O point group if aligned to the Cartesian axes.
   //
   // Included in: O, Oh
   //
   // Comments:
   // - The O pointgroup does not have full axis permutation symmetry (only,
   //   `SYMGROUP_AxisPerms_Even`, but not `SYMGROUP_AxisPerms_All`).
   //
   // - However, it *DOES* have this C4 rotation symmetry!
   //   Consider a 90° rotation in the x,y plane. This maps (x,y) to (y,−x),
   //   leaving z alone. So, this (with the other two C4s) effectively allows
   //   swapping any two axes, after all, provided we can take care of the
   //   minus sign.
   //
   // - As it happens, we can, if we deal with monomials for which only *EITHER*
   //   all exponents are even *OR* all exponents are odd. Because in this case,
   //   the minus sign will do nothing at all to x^i y^j z^k in which (i,j,k)
   //   are all even, while it will invert the sign of any monomial in which
   //   all of (i,j,k) are odd.
   SYMGROUP_AxesRotate4_AnyTwo,
   SYMGROUP_c4XY_c4XY_c4ZX = SYMGROUP_AxesRotate4_AnyTwo,
   SYMGROUP_O = SYMGROUP_c4XY_c4XY_c4ZX,


   // NOTE @ "*":
   // - Which of the axis-aligned Cartesian coordinate permutations and/or
   //   reflections are actually included as subgroups within the various
   //   polyhedral point groups depends on how the nominal polyhedra of the
   //   groups are oriented in 3D space!
   //
   // - The included/excluded lists provided here are only meant as overview over the
   //   Cartesian perm/refl subgroups at the time of this writing. And with different
   //   polyhedral alignment the combinations may turn out differently.
   //
   SYMGROUP_Count
};




// Represents a permuflection group `X`, which is subgroup of a full
// three-dimensional point group `G`, and its relation to `G`.
//
// Next to extracting the permuflection operators X ⊆ G, this object also
// constructs minimal sets of left (C_l) and right (C_r) coset representatives:
//
//     G = {c⋅x; c ∈ C_l, x ∈ X}  and  |G| = |X| |C_l|
//     G = {x⋅c; x ∈ X, c ∈ C_r}  and  |G| = |X| |C_r|
//
// For example, for the `I_h` group with 120 elements as `G`, the maximal
// permuflection group `X` has 24 elements. The upper equation then says that
// there is a set of five group elements C_l (|C_l| = |G|/|X| = 120/24 = 5),
// such that all 120 group elements of `G` can be written, uniquely, in the form
// `c⋅x` where x ∈ X is one of 24 permuflection group elements and c ∈ C_l is
// one of the 5 left coset representatives.
//
//
// Notes on the coset factorization:
// —————————————————————————————————
// This information allows factorizing group summations over G into separate
// summations over `X` and the coset representatives. For example, the grid
// points of an integration grid are obtained as
//
//       {(R_g r⃗_s, w_s); g ∈ G, s ∈ {1,…,#seeds}}
//
// where `r⃗_s` are seed points, `w_g` their weights, and `R_g` are (3,3)-shape
// orthogonal matrices representing the abstract symmetry operators g ∈ G.
// An optimization step in the grid optimization then involves, among other
// steps, the computation of the grid's residual errors (e_i) of the target
// function `t_i`'s numerical integrals, compared to the exact analytic ones:
//
//       e_i = ∑_s ∑_{g∈G} t_i(R_g r⃗_s) w_g - ∫ t_i(r) ⅆω
//
// With the coset factorization decribed above, we can factorize the group sum
// over `G` into a product of two summations (as follows), of which then one is
// absorbed into the target functions and the other one remains with the seed
// points:
//
//          ∑_{g∈G} t(R_g r)
//        = ∑_{x∈X} ∑_{c∈C} t((R_x R_c) r)
//        = ∑_{x∈X} ∑_{c∈C} (R_x⁻¹ t)((R_c r))
//        = ∑_{c∈C} (∑_{x∈X} R_x⁻¹ t)(R_c r)
//
// So the permuflections `X` can be absorbed into the target functions, and only
// sums over R_c for c ∈ C_r remains which are directly applied to real space
// vectors `r`
struct FAxisPermuflectionSubgroup : public ct::FIntrusivePtrDest1
{
   // element from the subgroup `X` this object represents
   typedef FAxisPermuflection FGroupElement;
   typedef std::vector<FGroupElement> FGroupElementList;

   // element from the parent group `G`
   typedef FRotationMatrix FGroupElement_Parent;
   typedef std::vector<FGroupElement_Parent> FGroupElementList_Parent;

   // Constructs the permuflection subgroup `X` information from a parent point
   // group `G`'s full list of general three-dimensional symmetry operators.
   // The constructor assumes that `ParentGroupElements` forms a complete group,
   // (this is not checked).
   // Keeps a live reference to `ParentGroupElements`.
   explicit FAxisPermuflectionSubgroup(FGroupElementList_Parent const &ParentGroupElements,
      unsigned nAxes_, unsigned nPhases_, std::string const &ParentGroupName, FPrintLevel PrintLevel);

   // query (Q) for group elements of this subgroup
   FGroupElementList const &ElementsQ() const { return m_Elements; }

   unsigned nAxes() const { return m_nAxes; }
   unsigned nPhases() const { return m_nPhases; }
   std::string const &ParentGroupNameQ() const { return m_ParentGroupName; }
   size_t ParentGroupOrderQ() const { return m_ParentGroupOrder; }

   // returns whether definitely (function argument) H ⊆ G (this object).
   // Note: a return value of `false` does not necessarily mean that H is
   // not a subgroup of *this.
   bool ContainsAsSubgroupQ(FSubGroupId H) const { assert(size_t(H) < m_ContainsAsSubgroup.size()); return m_ContainsAsSubgroup[size_t(H)]; }

   // returns whether H is a subgroup of *this and there is no other, large
   // subgroup in *this which also fully contains H (it's for printing only the
   // biggest ones)
   bool ContainsAlsoCoveringSubgroupQ(FSubGroupId H) const;
public: // ——— coset related interfaces ———
   enum FCosetSide {
      // A left coset c X is the set c X = {c x; x ∈ X} where X ⊂ G is a
      // subgroup and c ∈ G is a coset representative
      COSET_Left,
      // A right coset X c is the set X c = {x c; x ∈ X}
      COSET_Right
   };
   // return number of distinct left cosets of X, that is, sets {c_i x; x ∈
   // X}, which G can be uniquely split into. The number of left cosets and
   // right cosets is identical.
   size_t nCosets() const { assert(m_ParentGroupOrder % m_Elements.size() == 0); return m_ParentGroupOrder / m_Elements.size(); }
   // query(Q) for a sequence of unique representatives of all left/right cosets of X in G.
   FGroupElementList_Parent const &CosetRepsQ(FCosetSide iSide) const;
public: // ——— container-like interface (for subgroup elements) ———
   size_t size() const { return m_Elements.size(); }
   FAxisPermuflection const &operator[] (size_t i) const { return m_Elements[i]; }
   typedef FGroupElementList::const_iterator const_iterator;
   const_iterator begin() const { return m_Elements.begin(); }
   const_iterator end() const { return m_Elements.end(); }
protected:
   FPrintLevel
      m_PrintLevel;
   std::string
      // name or symbol of the point group which provided the full set
      // of symmetry operators (not necessarily axis aligned) for which
      // this is the permuflection subgroup.
      m_ParentGroupName;
   size_t
      // number of elements in the parent group
      m_ParentGroupOrder;
   unsigned
      m_nAxes,
      m_nPhases;
   FGroupElementList
      m_Elements;
   std::array<bool, SYMGROUP_Count>
      // for id from `SYMGROUP_*`, `m_ContainsAsSubgroup[id] == true` indicates that
      // the corresponding group is definitely a subgroup of *this.
      // Note: There is no guarantee of `m_ContainsAsSubgroup[id] == false` meaning
      // anything specific; the subgroup seach is not necessarily executed exhausively!
      m_ContainsAsSubgroup;

   FGroupElementList_Parent
      // Representatives of the left cosets of X = `*this` in the parent group G.
      // Concretely, c_i = `m_LeftCosetReps[i]` is a group element from G in the
      // #i'th left coset of X in G, such that {c_i x; x ∈ X} generates the
      // rest of the #i'th left coset.
      //
      // Each coset is represented, and by exactly one element. Therefore the
      // union of {c_i x; x ∈ X} over the c_i in `m_LeftCosetReps` covers the
      // entire group G, coset by coset, and each g ∈ G appears in only one i
      // subset
      m_LeftCosetReps,
      // …as `m_LeftCosetReps`, but for the right cosets {x c_i; x ∈ X}
      m_RightCosetReps;

   void _SearchSpecialSubgroups(FAxisPermuflectionSet const &Ops);
   void _PrintSpecialSubgroups(std::string const &Caption);
   void _MakeCosetReps(FGroupElementList_Parent &CosetReps, FCosetSide Side, FGroupElementList_Parent const &ParentGroupElements);
};
typedef ct::TIntrusivePtr<FAxisPermuflectionSubgroup>
   FAxisPermuflectionSubgroupPtr;

char const *CosetSideAsString(FAxisPermuflectionSubgroup::FCosetSide Side);


#ifdef INCLUDE_ABANDONED
// Note @ deleted code:
// - idea here was to once-determine and separately store the symmetries applicable
//   to different types of monomials (all exponents even, all odd, etc),
//   in an abstract fashion. Each monomial type gets an FMonomialSymmetryOpList.
// - in principle that is a good idea. All the monomials with different exponents
//   but same exponent "pattern" (which exponents equal to each other, which even/odd)
//   have the same group summation pattern, after all. There is no need to recompute
//   it every time.
// - however, I now decided that just recomputing the axis permuflection group
//   sums over monomials everytime they are needed is probably good enough (and
//   if not: if it ever should become a problem, it could still be easily fixed
//   later on, and with little more work than doing it right from the start, I
//   think)
// - so we no longer keep FMonomialSymmetryOpList, and instead add functions
//   to compute group sums of monomials and polynomials directly to FMonomialSymmetry,
//   and similarly add group-sum based classification functions (canonical, zero, etc.),
//   too.

// typedef std::pair<FMonomialN, phase_index_t>
//    FMonomialAndPhase;

// struct FMonomialSymmetryOpList
// {
//    // Describes what kind of special Cartesian monomials, if any, this object is
//    // representing the symmetry operations for
//    enum FMonomialType {
//       AllExpEven, // all monomial exponents are even. E.g., x^4 y^2 z^2 counts.
//       AllExpOdd, // all monomial exponents are odd. E.g., x^7 y^3 z^9 counts.
//       NeitherAllEvenNorAllOdd, // contains mixed even-ness exponents.
//       MonomialTypeCount
//    };
//
//    // Derives the set of symmetry operations applicable to monomials
//    // of the given type from special symmetry elements of the full point group
//    explicit FMonomialSymmetryOpList(FAxisPermuflectionSubgroup const &AxisPermuflections, FMonomialType MonomialType, std::string const &GroupName);
// protected:
//    FMonomialType
//       m_MonomialType;
//
//    typedef std::array<unsigned char, AIG_SPACE_AXES>
//       // Describes a permutation of the Cartesian axes. Concretely, for i = 0,1,2,
//       // … the x_i source coordinate gets mapped to the x_{m[i]} target
//       // coordinate.
//       FAxisPerm;
//
//    // Entry (perm p, phase-factor) represents that the source monomial
//    // with exponents [e_0, e_1, e_2], that is,
//    //
//    //    x_0^{e_0} x_1^{e_1} x_2^{e_2},
//    //
//    // where (x_i)_i are the Cartesian coordinates, gets mapped by the symmetry
//    // operation to the the phase-scaled target monomial for permuted axes as follows:
//    //
//    //     phase-factor * x_{p[0]}^{e_{0}} * x_{p[1]}^{e_{1}} * x_{p[2]}^{e_{2}}
//    //   = phase-factor * x_{0}^{e_{~p[0]}} * x_{1}^{e_{~p[1]}} * x_{2}^{e_{~p[2]}}
//    //
//    // In this, `p` denotes the Cartesian axis coordinate permutation (FAxisPerm),
//    // and `~p` its inverse (which is what gets applied to the exponents).
//    typedef std::pair<FAxisPerm, phase_index_t>
//       FSymmetryOp;
//    typedef std::vector<FSymmetryOp>
//       FSymmetryOpList;
//    FSymmetryOpList
//       m_SymmetryOps;
//
//
//    std::string
//       m_GroupName;
// };
#endif // INCLUDE_ABANDONED

// Implements symmetry-based compressions, reductions, and transformations of
// the optimization target space as spanned by monomials. The ideas and
// mathematical reasoning are described in a detailed comment right after this
// class, and another one near `FAxisPermuflectionSubgroup`.
struct FMonomialSymmetry : public ct::FIntrusivePtrDest1
{
   // Derives the set of symmetry operations applicable to monomials
   // of different types from the full list of operators of the three-dimensional point group.
   // The constructor assumes that these form a complete group — this is not checked.
   explicit FMonomialSymmetry(FRotationMatrixList const &GroupElements, unsigned nAxes, unsigned nPhases, std::string const &GroupName, FPrintLevel PrintLevel);
   explicit FMonomialSymmetry(FAxisPermuflectionSubgroupPtr pAxisPermuflections, FPrintLevel PrintLevel);


   // Returns a ref to the subgroup of symmetry operations which can be applied
   // to monomials easily in closed form (namely, the axis permuflection subgroup,
   // see above)
   FAxisPermuflectionSubgroup const &EasySubgroup() const { return *m_pAxisPermuflections; }

   // Returns list of |C|=|G|/|A| group elements `C` such that each element
   // g ∈ G of the full symmetry group can be written uniquely as `g = c a`,
   // with c ∈ C an element of the list and a ∈ A an element of the axis
   // permuflection group (as returned by `EasySubgroup()`).
   FRotationMatrixList const &EasySubgroup_LeftCosetReps() const {
      return m_pAxisPermuflections->CosetRepsQ(FAxisPermuflectionSubgroup::COSET_Left);
   }
   // As `EasySubgroup_LeftCosetReps`, but the g ∈ G are written
   // as `g = a c` (different order) with c ∈ C and a ∈ A.
   FRotationMatrixList const &EasySubgroup_RightCosetReps() const {
      return m_pAxisPermuflections->CosetRepsQ(FAxisPermuflectionSubgroup::COSET_Right);
   }

   // Maps Cartesian exponents to total factor obtained after summing up P_m m_{ijk}
   // for all permuflections P_m in `EasySubgroup()` for a given monomial m_{ijk}.
   // The final object is a polynomial, of course. However, note that the final
   // factors are not floats or anything this time — the handling of
   // permutations and reflections is exact, so we get well-defined integers.
   // This is a flattened version of the map. Only terms with non-zero coefficients
   // are retained.
   typedef std::vector<std::pair<FMonomialN, phase_scaled_int_t> >
      FMonomialCoeffList;

   // Evaluates sum terms of ∑_{â ∈ A} (â⁻¹ m_{ijk}) for the sum
   // over the axis permuflection group A (as returned by `EasySubgroup()`).
   //
   // Arguments:
   // - `MonomialCoeffs` [out]: list of sum terms (pm,c) where
   //   `pm` are axis-permutations of monomial `m`
   //   (i.e., with the same set of exponents as `m`, but possibly in
   //   different order), and `c` the associated coefficients.
   //   + Only terms with coefficient c != 0 are retained
   //   + A permuted monomial `pm` occurs at most once.
   //   + `MonomialCoeffs` is sorted with a sort key depending
   //     on `pm` only.
   //   If present, prior content of `MonomialCoeffs` is deleted, but allocated
   //   memory is recycled.
   // - `m` [in]: monomial to average.
   //
   // Returns:
   // - `false` if subgroup sum annihilates `m`, `true` otherwise.
   //   (return value is equivalent to `!MonomialCoeffs.empty()`)
   bool EvalEasySubgroupSum(FMonomialCoeffList &MonomialCoeffs, FMonomialN m) const;

   // evaluates ∑_{â ∈ A} (â⁻¹ p), for a polynomial function `p ≡ PolyIn`,
   // over the axis permuflection group A (as returned by `EasySubgroup()`).
   // Cull criteria of output polynomial are cloned from `PolyIn`.
   FPolynomialN EvalEasySubgroupSum(FPolynomialN const PolyIn) const;

   // Queries whether a given Cartesian monomial $m(x,y,z) = x^i y^j z^k$ will
   // definitely be annihilated if subjected to the group averaging operator.
   // I.e., returns whether
   //
   // \[   (1/|G|) ∑_{ĝ ∈ G} (ĝ m)(r⃗)  =  0            \]
   //
   // is *necessarily* true already as a result of the point group $G$ alone.
   // The group action as applied to scalar functions of space is defined via
   //
   // \[   (ĝ m)(r⃗)  = m(\mat R_g^{−1} (r_x, r_y, r_z)^T) \]
   //
   // where $\mat R_g$ is the (3,3)-shape rotation matrix representation of the
   // group element `ĝ`, and $(r_x,r_y,r_z)^T$ is the column vector of Cartesian
   // coordinates of `r⃗` expressed in the same basis.
   bool ZeroAverageQ(FMonomialN const &m) const;

   // For a monomial `m`, return whether `*this` would choose/select it as the
   // canonical representative of its group orbit $\{ ĝ m; ĝ ∈ G \}$.
   //
   // Notes:
   // - …the group orbit ***in function space***! $(ĝ m)(r⃗)$ is the
   //   function in L²(S₂,ℂ) mapping $r⃗ ∈ S₂$ to $m(\mat R_g^{-1} \vec r)$).
   // - This function returns `false` if it reasons that the group average
   //   `P_G[m]` vanishes (e.g., if `ZeroAverageQ(m)` yields `true`)
   bool SelectAsCanonialReprQ(FMonomialN const &m) const;


protected:
   FAxisPermuflectionSubgroupPtr
      m_pAxisPermuflections;
   std::string
      m_GroupName;
   FPrintLevel
      m_PrintLevel;
};
typedef ct::TIntrusivePtr<FMonomialSymmetry>
   FMonomialSymmetryPtr;
typedef ct::TIntrusivePtr<FMonomialSymmetry const>
   FMonomialSymmetryCptr;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// ▒▒ @TARGET SPACE REDUCTIONS implemented in FMonomialSymmetry               ▒▒
// ─────────────────────────────────────────────────────────────────────────────
// This comment provides the mathematical reasoning behind the use of symmetry
// groups for the a-priori identification of monomials (or other target
// functions) which can savely omitted from the optimization target space…
//
// - (prescreening:) …because they automatically get integrated exactly for
//   symmetry reasons, regardless of the concrete grid parameters {(r⃗ₛ, wₛ)}ₛ
//   being optimized,
//
// - (canonical reps:) …or are contained in a group of monomials yielding non-
//   trivial−but−identical effects if included in the target space, so that only
//   one canonical representative of the symmetry-equivalent group needs to be
//   explicitly considered.
//
// ────────────────────────────
// ░░ Monomial prescreening: ░░
//
// Subgroup symmetries can allow us to a-priorily discard target functions from
// the optimization space, by the following argument:
//
// - Assume that the full group G has a subgroup H ⊂ G. We can then decompose
//   each elements ĝ ∈ G into a product ĝ_{ij} = ĉ_i⋅ĥ_j of one element from the
//   |H|-element subgroup
//
//      H := {ĥ_j; j ∈ {1,…,|H|}}
//
//   and another element c_i ∈ C_l (i ∈ {1,…,|G|/|H|}) of coset
//   representatives C_l of H with respect to G:
//
//      G = {ĉ⋅ĥ; ĉ ∈ C_l, ĥ ∈ H}  and  |G| = |H| |C_l|
//
//   (for details, see starting comments on the motivation of
//   `FAxisPermuflectionSubgroup` in this header file).
//
// - Consequently, we can rewrite any full group sum, e.g., over Cartesian
//   monomials, as follows:
//
//      ∑_{ĝ ∈ G} (ĝ m_{ijk…})
//      = ∑_{ĉ ∈ C_l} ∑_{ĥ ∈ H} (ĉ⋅ĥ m_{ijk…})
//      = ∑_{ĉ ∈ C_l} ĉ (∑_{ĥ ∈ H} ĥ m_{ijk…})     (1)
//
// - This means that if we can find any subgroup H ⊆ G, no matter how small,
//   for which the partial group sum over H vanishes, then by (1) also the full
//   group sum gets annihilated. That is:
//
//      If there exists a subgroup H ⊆ G for which
//
//          ∑_{ĥ ∈ H} ĥ m_{ijk…} == 0               (2a)
//
//      Then also the full group sum vanishes:
//
//          ∑_{ĝ ∈ G} ĝ m_{ijk…} == 0               (2b)
//
// - In such a case we can discard the function m_{ijk…} from the optimization
//   target space: As explained in the intro of `FAxisPermuflectionSubgroup`, in
//   the grid summation over an integration grid of G symmetry,
//
//         {((R_g r⃗_s), w_s); ĝ ∈ G},                    (G-grid)
//
//   the sum over g ∈ G-transformed seed points `(R_g r⃗_s) ∈ ℝ^n` (for
//   fixed target function `m_i`) is mathematically equivalent to the sum over
//   g ∈ G-transformed target functions `(g m_i)` (applied to seed point `r_s`
//   directly, without summing):
//
//         ∑_{ĝ ∈ G} m_i(R_g r⃗_s) == (∑_{ĝ ∈ G} ĝ⁻¹ m_i)(r⃗_s)       (3)
//
//   In this, eq (2b) does apply to the rhs sum if a subgroup H fulfilling eq
//   (2a) for target function `m_i` can be identified, and therefore also to the
//   lhs sum, even if we never were to actually evaluate any group-summed target
//   function explicitly. That is, if (2a) holds for a subgroup H ⊆ G and a
//   monomial `m_{ijk…}`, then (literally) any G-symmetric integration grid
//   automatically evaluates the integral of `m_{ijk…}` correctly (to zero),
//   regardless of where the seed points are placed or how the weights are
//   assigned.
//
// ───────────────────────────────────────────────────────────────────
// ░░ Symmetry-equivalent-monomial-group canonical representatives: ░░
//
// - A very similar argument allows us to discard all but one of the monomials
//   `m_{ijk…}` which yield scalar multiples of the same group average
//   (cf rhs eq 3)
//
//         P_G[m_i] := (1/|G|) ∑_{ĝ ∈ G} ĝ⁻¹ m_i                    (4)
//
//   In particular, this applies to any group of functions in the orbit H m_i of
//   any subgroup H ⊆ G's.
//
// - Concretely, due to the G-part of the grid summation (eq (G-grid)), only the
//   projection P_G[m_i] (eq 4) actually matters for the optimization target
//   space, which becomes clear if we rewrite the rhs of eq 3 with eq 4:
//
//         ∑_{ĝ ∈ G} m_i(R_g r⃗_s) == |G| P_G[m_i](r⃗_s)
//
// - And all functions m_k ∈ (H m_i) in the same subgroup orbit yield colinear
//   vectors P_G[m_k] ∝ P_G[m_i] — so if one of them gets integrated exactly
//   to a non-vanishing result, the other ones automatically get integrated
//   exactly, too.



#endif // AIGG_AXIS_PERMFLECTION_H
