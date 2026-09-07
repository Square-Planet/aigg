#ifndef TARGET_FUNCTIONS_H
#define TARGET_FUNCTIONS_H

#include "Aigg.h"
#include <set>
#include "PolyXyz.h"


struct FPointGroup;
struct FGridSearchOptions;
struct FMonomialSymmetry;

struct FMonomialEvalInfo
{
   FScalar const
      *pPowX, *pPowY, *pPowZ;
};

// FIXME: Remove these in cleanup.
const extern bool
   // If set, track which monomials are equivalent to which others
   // after averaging them over the group
   g_TrackEquivMonomials;
const extern bool
   // If set, calculate overlap matrices and rhs/error vectors using
   // explicitly group-averaged monomials, rather than raw monomials.
   //
   // Notes:
   // - This will introduce hard singularities in the overlap matrix. After
   //   group-averaging, some format target functions and their linear
   //   combinations will be strictly unnecessary.
   // - Currently quite ugly. Probably should make a different class than FMonomial
   //   for this... currently this working relies on some projector magic
   //   which makes results with Rg-summed raw monomials and explicitly
   //   group-symmetrized monomials (which are polynomials) mathematically
   //   equivalent in the residual calculation, despite the fact that they look
   //   quite different.
   g_UseGroupSymmetrizedFn;

struct FMonomial {
   size_t
      // powers for x, y, z
      ix, iy, iz;
   FScalar
      // number of symmetry-equivalent monomials — residuals should be weighted by this.
      m_Degeneracy,
      // could (optionally) be used to re-scale monomials for normalizing them.
      m_Factor;
   // FIXME:
   // - get rid of `m_RgSymmetrized` and `m_SymmetryEquivFn`. And probably also
   //   of most of the implicit symmetrization support in the `TargetFunctions.cpp`.
   //   Do not think we really need this anymore now that we do have the code
   //   to precompute and store/load orthonormal symmetrized subspace bases.
   //
   // - Figure out what is happening with `m_Degeneracy`, `m_Factor`, etc.
   //   Note that `m_Factor` is actually used — at first glance it seems as
   //   if it only affects the overlap matrices. However, it does enter
   //   `BaseFactor()`, and via that one also all other evaluation routines
   //   (in particular, `CalcUnitSphereIntegral` and `EvalPoint`)
   //
   // - Are these base normalization factors really fully and correctly
   //   accounted for in the group average calculations? And in the import/
   //   export of polynomial subspace data?
   FPolynomialN
      // = (1/|G|) ∑_g m((R_g e_x)^ix, (R_g e_y)^iy, (R_g e_z)^iz)
      m_RgSymmetrized;
   typedef std::set<FPolynomialN::FMonomial>
      FSymmetryEquivSet;
   FSymmetryEquivSet
      // if given: set of other Cartesian monomials which did average to
      // the same value of `m_RgSymmetrized` after group averaging was applied.
      // (note: not stored unless `g_TrackEquivMonomials` is set)
      m_SymmetryEquivFn;
   size_t l() const { return ix + iy + iz; }
   FMonomial(size_t ix_, size_t iy_, size_t iz_);

   // if this is not zero the sphere integral over the monomial vanishes.
   size_t nUnevenPow() const { return (ix % 2) + (iy % 2) + (iz % 2);  }

//    FScalar CalcUnitSphereIntegral() const;
   FScalar CalcUnitSphereIntegralRaw() const;
   FScalar CalcUnitSphereIntegral() const;
   FScalar EvalPoint(FPoint const &p) const;
   FScalar EvalPoint(FPoint &dMdXyz, FPoint const &p) const;
   FScalar EvalPoint(FMonomialEvalInfo const &p) const;
   FScalar EvalPoint(FPoint &dMdXyz, FMonomialEvalInfo const &p) const;

   void _ResetNorm();

   // … what was I thinking when making this ix/iy/iz thing?!
   inline size_t operator[] (size_t i) const {
      assert(i < 3);
      return (i == 0)? ix : ((i == 1)? iy : iz);
   }
protected:
   FScalar BaseFactor() const;
};

typedef std::vector<FMonomial>
   FMonomialList;


// Estimate the squared L₂(S₂) norm (based on unit sphere integral scalar product).
//
// Note: it is a `double` specialization (even if compiling in high precision
// mode), because we are anyway computing only estimates for screening!
// These are meant for use in `FCullCriteria` objects of `FPolynomialN`.
//
// With the formula used at the time of this writing, these should be accurate
// to about 1e-4 relative error for any combination of i,j,k.
double EstimateMonomialNormSq_S2(FPolynomialN::FMonomial const &m);


struct FTransformedPointInfo
{
   explicit FTransformedPointInfo(const FPointArray &UntransformedPoints_, size_t lmax_);

   FPointArray
      // seed points transformed by current symmetry operation
      TransformedPoints;
   FDenseArray
      // powers of monomial components
      //    [1, x, x^2, x^3, …, x^lmax],
      //    [1, y, y^2, y^3, …, y^lmax],
      //    [1, z, z^2, z^3, …, z^lmax],
      // precomputed for (x,y,z) of seed points transformed by current symmetry operation
      XyzPow[3];

   void TransformPoints(FRotationMatrix const &SymOp);
   FMonomialEvalInfo TransformedPointInfo(size_t iPoint) const;
protected:
   FPointArray const
      // seed points before the transformation
      &m_OrigPoints;
   size_t
      lmax, nPoints;
};


FMonomialList MakeMonomialList(size_t lmin, size_t lmax);
// FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FRotationMatrixList const &TrafoList, int iPrintLevel);
FMonomialList MakeSymmetryUniqueMonomialList(size_t lmin, size_t lmax, FPointGroup const &PointGroup, FPrintLevel iPrintLevel, FGridSearchOptions const *pOptions);
FDenseMatrix MakeOverlapMatrix(FMonomialList const &Monomials);
FDenseMatrix MakeOverlapMatrix_AxprAveraged(FMonomialList const &Monomials, FMonomialSymmetry const *pMonomialSymmetry);

// note: this function does nothing unless the passed iPrintLevel >= 2.
void PrintMonomialList(std::ostream &xout, FMonomialList const &r, FPrintLevel PrintLevel, std::string const &Comment = std::string());
void ComputeAndStoreGroupAverages(FMonomialList &r, FPointGroup const &PointGroup, bool DeleteZeroAverages);


FScalar CalcMonomialUnitSphereNormSq(unsigned ix, unsigned iy, unsigned iz);
FScalar CalcMonomialUnitSphereNorm(unsigned ix, unsigned iy, unsigned iz);
FScalar CalcMonomialUnitSphereNormSq(FMonomialN const &m);
FScalar CalcMonomialUnitSphereNorm(FMonomialN const &m);
FScalar CalcMonomialUnitSphereIntegralRaw(unsigned ix, unsigned iy, unsigned iz);


#endif // TARGET_FUNCTIONS_H
