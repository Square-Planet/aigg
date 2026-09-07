#ifndef AIGG_TARGET_SPACE_H
#define AIGG_TARGET_SPACE_H

#include "Aigg.h"
#include <string>
#include <list>
#include "CxIntrusivePtr.h"
#include "SymmetryGroups.h"
#include "TargetFunctions.h"
#include "MathOps.h"
#include "SearchOptions.h" // only really for print options, atm


typedef std::vector<ptrdiff_t>
   FDegreeSchedule;
typedef std::vector<FPoint>
   FPointList;
struct FGridSearchOptions;

namespace ct {
   struct string_slice;
}


struct FTargetSpacePreprocessOptions : public ct::FIntrusivePtrDest1
{
   // construct property object from a textual grid property list.
   explicit FTargetSpacePreprocessOptions(std::string const &Desc_);

   typedef std::vector<FPointGroupPtr>
      FPointGroupList;
   FPointGroupList
      PointGroups;
   FDegreeSchedule
      // process this explicit list of angular momenta
      lList;
   int
      iPrintLevel;
   FPrintOptions
      PrintOptions;
   FAccuracyPreference
      AccuracyPreference;
   bool
      StoreVersionInfo, // store info on program version (text format only)
      StoreIdentityInfo; // store info on host name, date, and user into output file (text format only)
   enum FSkipExistingOptions {
      // do not run preprocessor if a target space file for this group/order
      // delcaration is already there (regardless of that data file's meta data).
      SKIP_All,
      // do not run preprocessor if a target space file for this group/order
      // delcaration is already there and that data file has been processed with
      // a higher floating point accuracy that the current version of aigg.
      SKIP_Better,
      // as SKIP_Better, but also skip existing files generated with the *same*
      // floating point accuracy
      SKIP_BetterOrEqual,
      // do not skip anything — all target space files, if existing, will be overwritten.
      SKIP_None
   };
   FSkipExistingOptions
      // declares whether/how to skip generating data files for point groups which
      // are already there.
      SkipExisting;
   bool
      Export;
   void SetAccuracyPreference(FAccuracyPreference const &AccuracyPreference_);
};

typedef ct::TIntrusivePtr<FTargetSpacePreprocessOptions>
   FTargetSpacePreprocessOptionsPtr;


// main driver routine for running aigg in preprocess mode (i.e., this one is called from main() more-or-less directly).
int RunTargetSpacePreprocess(FTargetSpacePreprocessOptions const &Options);


struct FIoTargetFn {
   FScalar
      lambda;
   FPolynomialNPtr
      poly;
   FIoTargetFn() : lambda(0) {}
   FIoTargetFn(FScalar const &ew_, FPolynomialNPtr poly_) : lambda(ew_), poly(poly_) {}
};

typedef std::vector<FIoTargetFn>
   FIoPolynomialList;



enum FTargetFnSpaceFlags {
   FNSPACE_Orthogonal = 0x0001,
   FNSPACE_GroupAveraged = 0x0002
};

// An auxiliary data structure for input/output of basis function data. These ones
// may be loaded from data files created in `--preprocess` mode using the
// `TryReadTargetFnSpaceDataFile` function below.
struct FIoTargetFnSpace : public ct::FIntrusivePtrDest {
   std::string
      sPointGroup, // point group symbol
      sDomain; // integration domain (S2 = unit sphere in 3d)
   std::string
      SourceFile;
   int
      // this function set represents space of "monomials of degree <= this"
      // restricted to the target domain, for integration rules of the given
      // symmetry group
      Degree,
      // number after 'fp' as read from the file. (e.g., 64 = binary64 floats).
      // If this is smaller than g_ScalarFloatSizeBits, we probably should not use
      // the data from here. Other direction should be fine, though.
      FloatSizeBits;
   typedef std::vector<FScalar>
      FScalarList;
   FScalarList
      // `lambda[iFn]` is function `iFn`s symmetrized-overlap eigenvalue as read from the file.
      // (note: these are stored with only a few digits of precision — meant for info purposes,
      // not actual calculations)
      pLambda;
   typedef std::vector<FPolynomialNPtr>
      FIoFnList;
   FIoFnList
      pFns;
   unsigned
      Flags;
   explicit FIoTargetFnSpace() : Degree(-1), FloatSizeBits(-1), Flags(0) {}
};
typedef ct::TIntrusivePtr<FIoTargetFnSpace>
   FIoTargetFnSpacePtr;
typedef ct::TIntrusivePtr<FIoTargetFnSpace const>
   FIoTargetFnSpaceCptr;

// check if a data file of sufficient accuracy is present which contains a
// subspace basis for optimizing integration rules of `PointGroup` symmetry with
// order `l`.
//
// Returns 0 if either no data file for this group/order is present, there is a
// data file, but it is of lower floating point accuracy than this version of
// aigg is running, or the file is there but cannot be read/interpreted for some
// reason. If `pFailureReason` is provided and an error occurs, `*pFailureReason`
// the routine may set `*pFailureReason` to provide indications of what went wrong.
FIoTargetFnSpacePtr TryReadTargetFnSpaceDataFile(FPointGroup const &PointGroup, int l, std::string *pFailureReason=0);


struct FTargetFnListOptions {
   FPrintLevel
      PrintLevel;
   FAccuracyPreference
      AccuracyPreference;
   bool
      AllowNeglectResidualDeriv;
   explicit FTargetFnListOptions(FPrintLevel const &PrintLevel_ = FPrintLevel::Basic, FAccuracyPreference const &AccuracyPreference_ = ACCURACY_Default, bool AllowNeglectResidualDeriv_ = false);
};


// represets a list target functions for the grid optimization (i.e., functions
// which we ideally want to integrate exactly at convergence). This structure is
// as primary abstraction mechanism for supporting different types of functions
// (in particular, raw monomials and precomputed orthonormalized group-averaged
// polynomials) in the grid optimization.
struct FTargetFnList : public ct::FIntrusivePtrDest1
{
   // returns number of logical functions (--> dimension of fn dimension)
   virtual size_t nFn() const = 0;
   // order (`=lmax`): heighest `l` such that integrating all functions of this
   // list implies that all $Y^l_m$ are integrated exactly.
   virtual unsigned Order() const = 0;
   // e.g., "raw monomials"
   virtual std::string FnTypeName() const;
   // e.g., "symmetry unique"
   virtual std::string FnTypeAnnotation() const;
   // Optionally write some extra info about the target space — will be written
   // after the number of target functions. Default implementation does nothing.
   virtual void WriteInfo(ct::FLog &io) const;
//    // returns `(nFn,nFn)`-shape overlap matrix of integrals \f$ S_{ij} = ∫ fᵢ(r⃗)^* fⱼ(r⃗) \d Ω \f$
//    // FIXME: should this really be here? comments below suggest that can be handled entirely inside this class
//    // (the Scd stuff is just a linear transform of K and t, isn't it?)
//    virtual FDenseMatrix EvalOverlapMatrix() const = 0;
   // returns `(nFn)`-shape vector of analytic integrals:
   // \f$ tᵢ = ∫_{S_2} fᵢ(r⃗) ⅆΩ \f$
   virtual FDenseVector EvalTargetIntegrals() const = 0;

   // default implementation does nothing
   virtual void Print(ct::FLog &io) const;

   // Overwrite data in `(nTargetFn, nSeedPt)`-shape matrix `K(iFn,iPt)` by \f$ K(i,q) = ∑_g f_i(R_g^{-1} r⃗_q) \f$.
   // K already must have the correct shape on input — it is not allocated or resized here.
   virtual void EvalK(FDenseMatrix &K, FPointArray const &SeedPoints) = 0;
   // Evaluate contracted kernel matrix derivatives:
   //
   //     dKw[i,p] = ∑_{q} dK[i,q]/dp w[q]
   //
   // — kernel derivative contracted to fixed weights; computes `dKw`: `(m_nTargetFns, nVars)`-shape matrix
   //
   //     dKe[q,p] = ∑_{i} dK[i,q]/dp e[i]
   //
   // — kernel derivative contracted to fixed `m_rhs e[i]`; computes `dKe`: `(m_nSeeds, nVars)`-shape matrix
   virtual void EvalDerivK(FDenseMatrix &dKw, FDenseVector const &w, FDenseMatrix &dKe, FDenseVector const &e, FPointArray const &SeedPoints, FOrbitTypeRawPtrList const &OrbitTypes, FRotationGenerator const *pAi) = 0;

   size_t size() const { return this->nFn(); }

   virtual std::string MakeFnDesc(size_t iTargetFn) const;
protected:
   explicit FTargetFnList(FPointGroup const *pPointGroup_, FTargetFnListOptions const &GeneralOptions);
   FTargetFnListOptions
      m_GeneralOptions;
   FPointGroup const
      *m_pPointGroup;
   FRotationMatrixList const
      *m_pSymmetryOps;
   FPrintLevel
      m_PrintLevel;
};
typedef ct::TIntrusivePtr<FTargetFnList>
   FTargetFnListPtr;
typedef ct::TIntrusivePtr<FTargetFnList const>
   FTargetFnListCptr;


// …or should I support inline transforms? The thing is: all that Scd stuff
// *probably* does not have to go out of here. It's not really anything
// different than linearly transforming `K` and `t`, is it?
struct FTargetFnList_RawMonomials : FTargetFnList
{
   FTargetFnList_RawMonomials(size_t lmin, size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FGridSearchOptions const *pGridSearchOptions=0);

   size_t nFn() const; // override
   unsigned Order() const; // override
   std::string FnTypeName() const; // override
   std::string FnTypeAnnotation() const; // override

   // returns `(nFn,nFn)`-shape overlap matrix of integrals \f$ S_{ij} = ∫ f_i(r⃗)^* f_j(r⃗) \dΩ \f$
   FDenseMatrix EvalOverlapMatrix() const;
   FDenseVector EvalTargetIntegrals() const; // override

   void EvalK(FDenseMatrix &K, FPointArray const &SeedPoints); // override
   void EvalDerivK(FDenseMatrix &dKw, FDenseVector const &w, FDenseMatrix &dKe, FDenseVector const &e, FPointArray const &SeedPoints, FOrbitTypeRawPtrList const &OrbitTypes, FRotationGenerator const *pAi); // override
   std::string MakeFnDesc(size_t iTargetFn) const; // override
   void Print(ct::FLog &io) const; // override;
protected:
   size_t
      m_lmin, m_lmax;
   FMonomialList
      // this set of functions defines the function space which we hope to
      // integrate exactly with the numerical integration rule being optimized
      // (after it is optimized)
      m_TargetFns;
   FLinearSolverPtr
      // A stored matrix decomposition (`L*L.T` or `eigh`) which can be used to
      // solve equation systems involving the `(nTargetFn, nTargetFn)`-shape target
      // function overlap matrix `S`, with elements either $S_{ij} = ⟨f_i|f_j⟩$ or
      // $S_{ij} = ⟨f_i|P̂_G|f_j⟩$ (where P̂_G=P̂_G^2 is the group symmetrizer; see
      // comments in `TargetFunctions.cpp`). Mostly used for re-conditioning
      // equation systems by isolating length scales associated with the
      // condition number of said `S`.
      m_pScd;
};


struct ETabulatedFnListError : public std::runtime_error
{
   explicit ETabulatedFnListError(std::string const &Task, FPointGroup const *pPointGroup, int Order, std::string const &FileName, std::string const &ProblemDesc);
   explicit ETabulatedFnListError(std::string const &ProblemDesc);
};


// Tabulated group-averaged (symmetrized) orthonormalized polynomials. Note that
// these functions have the group summation "built in", so to say, so there is
// no need for group element iteration when evaluating `t`/`e`/`K`. Also, the
// functions are already orthonormalized when we store them, so there is no need
// for overlap matrices or explicit or implicit orthogonalization.
//
// Note:
// - This object can represent an irreducible orthonormal basis of group-
//   averaged functions which span exactly the minimal subspace of function we
//   require to integrate exactly to get a order-l integration rule for
//   PointGroup. This should be the best-possible target space for the grid
//   optimization.
//
// - However, instanciating an object of this class does not compute the target
//   space basis — rather, this will attempt to load the function list from a
//   pre-computed file on disk. If unsuccessful, an exception of class
//   ETabulatedFnListError will be raised.
//
// - The function list files are made by running aigg in "--preprocess" mode for
//   the target group/order. Computing the optimal basis is expensive and needs
//   to be done with high floating point accuracy. But at least in principle, it
//   would be sufficient if it were done at most once in human history for any
//   point-group/order /max-float-precision combination.
struct FTargetFnList_PreAveragedPolys : FTargetFnList
{
   FTargetFnList_PreAveragedPolys(size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FGridSearchOptions const *pGridSearchOptions=0);
   FTargetFnList_PreAveragedPolys(size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FIoPolynomialList const &pIoPolys); // debug fn!

   size_t nFn() const; // override
   unsigned Order() const; // override
   std::string FnTypeName() const; // override
   std::string FnTypeAnnotation() const; // override

   FDenseVector EvalTargetIntegrals() const; // override

   void EvalK(FDenseMatrix &K, FPointArray const &SeedPoints); // override
   void EvalDerivK(FDenseMatrix &dKw, FDenseVector const &w, FDenseMatrix &dKe, FDenseVector const &e, FPointArray const &SeedPoints, FOrbitTypeRawPtrList const &OrbitTypes, FRotationGenerator const *pAi); // override
   std::string MakeFnDesc(size_t iTargetFn) const; // override
   void WriteInfo(ct::FLog &io) const; // override
protected:
   size_t
      m_lmax;
   unsigned
      // FNSPACE_* bit field
      m_TargetFnFlags;
   FMonomialList
      // set of (raw) monomials over which the polynomials are expanded
      m_BasisFns;
//    FDenseMatrix
//       // $S_{ij} = ⟨b_i, b_j⟩$ over basis functions in m_BasisFns
//       m_BasisOverlap;
// ^-- no longer needed.
   FScalar
      // if this is a set of group-averaged functions, then there is
      // no point in actually executing the sum over group-elements in
      //
      //     K(i,q) = ∑_g f_i(R_g⁻¹ s⃗_q)
      //
      // and its derivatives (`R_g` is the group action as a rotation matrix,
      // and `s⃗_q` the `q`'th seed point).
      // If the target functions already are group-averaged, then evaluating
      // them on the `R_g⁻¹ s⃗_q` will just yield the same value for every
      // group element `ĝ`.
      //
      // However, in the `K`-matrix definition, we formally have a *sum*
      // over all group elements — and unlike group-averaging, group-summing
      // is *not* an involution. In short: even if all the sum terms are the
      // same, they still need to be there. This factor here contains the
      // number of group elements `|G|`, so that we may replace the sum in `K` by:
      //
      //     K(i,q) = |G| f_i(s⃗_q)
      //
      // Note that a different alternative to absorbing the `|G|` factor into `K`
      // is to absorb a `1/|G|` factor into the rhs vector `t`. This option has
      // as advantage that `rhs` is only built once, and it is a much smaller
      // object than `K` and it derivatives.
      m_fGroupImageFactor;
   typedef FIoTargetFnSpace::FIoFnList
      FCoFnList;
   FCoFnList
      // list of contracted functions.
      m_TargetFns;
   typedef FIoTargetFnSpace::FScalarList
      FScalarList;
   FScalarList
      m_Lambdas;
   FDenseMatrix
      // `(m_BasisFns.size(), nFn())`-shape matrix `C` of expansion coefficients of
      // the group-averaged polynomials in terms of `m_BasisFns`. Like an orbital
      // matrix, with elements $[C]_{μ,i} = C^μ_i$; i.e. `C` carries one
      // contravariant index `μ` for the basis functions (upper, first index/row
      // index) and one orthonormal index `i` for the target functions
      // (second/column index).
      m_Coeffs,
      // `= S * C`. Needed for transforming from target-fn dimension (orth)
      // to covariant basis indices (as opposed to contravariant ones, like
      // m_Coeffs does). Used in computation of gradient for contraction
      // with intermediate `dK/dx`.
      m_CoeffsCov;
   FRotationMatrixList
      // if we have pre-group-averaged functions, we do not need to do any
      // explicit symmetrization. In this case we'd just use a single identity
      // operator here. However, for the moment I'd like to retain the ability
      // to easily switch this on or off.
      m_SymmetryOpsForK;

   void InitBasisAndCoeffs(FIoTargetFnSpace::FIoFnList const &pIoFns);
};



#endif // AIGG_TARGET_SPACE_H
