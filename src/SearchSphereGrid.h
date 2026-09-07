#ifndef SEARCH_SPHERE_GRID_H
#define SEARCH_SPHERE_GRID_H

#include "Aigg.h"
#include "PointCloud.h"
#include "SymmetryGroups.h"
#include "TargetFunctions.h"
#include "TargetSpace.h"
#include "SearchOptions.h"
#include "MathOps.h"
#include "CxTiming.h"


struct FExportedGridStats;

static void PrepareInitialSeeds(FPointArray &SeedPoints, FPointGroup const &PointGroup);
static void CheckSeedsForCollapse(FPointArray const &SeedPoints, FPointGroup const &PointGroup);

// Test if a given spherical grid, given by SeedPoints and weights, is a valid
// integration rule of order lmax. To this end, this function computes the
// values of all monomials of order <= lmax both analytically and numerically,
// and compares the results. This function will return 'false' if the maximum
// integration error (absolute) of any monomial exceeds ThreshAbs. If given, the
// maximal absolute and relative deviations will be stored to *pMaxErrAbs and
// *pMaxErrRel
bool TestSphereGrid(FScalar *pMaxErrAbs, FScalar *pMaxErrRel, size_t lmax, FScalar ThreshAbs, FPointArray &SeedPoints, FDenseVector const &weights, FPointGroup const &PointGroup, FGridSearchOptions::FVerifyRuleOptions VerifyRuleOptions);

// Writes a data file in export format to save the final grids for
// collection/visualization/analysis
void ExportGrid(FPointArray const &SeedPoints_, FDenseVector const &LastWeights_, FScalar FinalResidual,
   int lmax, FExportOptions const &ExportOptions, FGridSearchOptions const &Options);

// Just a wrapper for constructing and running a FSphereGridOptContext object.
FScalar OptimizeSphereGrid(FPointArray &SeedPoints, size_t lmin, size_t lmax, FExportedGridStats *pGridStats, FGridSearchOptions const &Options);

// Controls the outer main loop as described in 'Options'; in most cases this
// will involve multiple calls to OptimizeSphereGrid() to incrementally optimize
// a given grid for increasingly large amounts of target functions.
void RunGridSearchSchedule(FGridSearchOptions const &Options);


struct FExportedGridStats
{
   FScalar
      wtsp,
      res;
   ptrdiff_t
      lmax;

   FExportedGridStats();
   // returns 'yes' if the grid described by *this is clearly worse than
   // the one described by 'other'.
   bool IsWorseThan(FExportedGridStats const &other) const;
   bool IsEquivalent(FExportedGridStats const &other) const;
};


// Contains problem descriptionand intermediate data used during the
// optimization of a single sphere grid for one set of target lmin/lmax. This
// object's purpose is to simplify subdividing the optimization function into
// related parts.
struct FSphereGridOptContext
{
   // note: SeedPoints is an *output* parameter!
   explicit FSphereGridOptContext(FPointArray &SeedPoints, size_t lmin, size_t lmax, FGridSearchOptions const &Options);
   
   // set grid stats to compare to before exporting a grid: if previous grid is
   // clearly better than the current one (e.g., both have same lmax, but this
   // one has negative weights, and the previous one does not) then do not
   // export the current grid (to prevent overwriting the earlier)
   void SetGridExportReference(FExportedGridStats const &RefStats);
   // returns stats of last exported grid, or 0 if there aren't any.
   FExportedGridStats const *pGetLastExportedGridStats() const;
protected:
   FGridSearchOptions const
      &m_Options;
   inline FPrintLevel PrintQ(unsigned iPrintOption) const { return m_Options.PrintOptions[iPrintOption]; }
   
   ct::FTimer
      // reset at start of Run()
      m_tThisGridTotal;

   FPointArray
      m_SeedPoints;
   size_t
      m_lmin, m_lmax;

   FPointGroup const
      &m_PointGroup;
   FRotationMatrixList const
      &m_SymmetryOps;

   
   size_t
      // number of seed points which's point group images will make up the
      // positions of the integration grid.
      m_nSeeds,
      // total number of free geometric parameters of the seed points which can
      // be optimized. Each seed point typically brings 0,1,2, or 3 parameters,
      // depending on the type of its orbit (e.g., a vertex orbit typically has
      // no free geometric parameters, while a general orbit typically has two
      // (three cartesian coordinates x,y,z, which are, however, restricted by
      // x^2 + y^2 + z^2 = 1)). This characterizes our optimization space---
      // finding these geometric parameters is the goal of this routine.
      m_nFreeParameters,
      // total number of points in the integration rule (i.e., union of all
      // orbits of the seed points under the point group).
      m_nPointsTotal;
   FOrbitTypeRawPtrList
      // array of types of orbits each seed point belongs to
      m_OrbitTypes;

   void InitSeedsAndOrbits();

#ifdef INCLUDE_ABANDONED
//    FMonomialList
//       // this set of functions defines the function space which we hope to
//       // integrate exactly with the numerical integration rule being optimized
//       // (after it is optimized)
//       m_TargetFns;
//    size_t
//       // number of target functions in the cost function
//       m_nTargetFns;
//    FLinearSolverPtr
//       // A stored matrix decomposition (L*L.T or eigh) which can be used to
//       // solve equation systems involving the (nTargetFn,nTargetFn)-shape target
//       // function overlap matrix S, with elements either S_{ij} = <f_i|f_j> or
//       // S_{ij} = <f_i|P_G|f_j> (where P_G=P_G^2 is the group symmetrizer; see
//       // comments in TargetFunctions.cpp). Mostly used for re-conditioning
//       // equation systems by isolating length scales associated with the
//       // condition number of said S.
//       m_pScd;
//    FDenseVector
//       // vector of exact sphere integrals of the target functions.
//       m_rhs;
#endif // INCLUDE_ABANDONED
   FTargetFnListPtr
      // this set of functions defines the function space which we hope to
      // integrate exactly with the numerical integration rule being optimized
      // (after it is optimized)
      m_pTargetFns;
   size_t
      // number of target functions in the cost function
      m_nTargetFns;
   FDenseVector
      // vector of exact sphere integrals of the target functions.
      m_rhs;

   void InitTargetFns();

   
   bool
      m_Converged;
   FScalar
      m_FinalResidual;
   FDenseVector
      m_LastWeights;
   ptrdiff_t
      m_nExtraDofForWeights;
   FRotationGenerator
      // generators of rotation around the x,y, and z axis
      m_A[3];
   FLinearSolverType
      m_KtkdSolverType;
   FScalar
      // minimum residual encountered during the iterations.
      // (TODO: check if still used; was previously used in calculation of step
      // length restrictions)
      m_fMinRes;

   // lmax/wtsp of the last grid which was exported after a (successful)
   // optimization; meant for complex schedules in order to not override a
   // clearly better earlier optimizied grid in a later downward search.
   FExportedGridStats
      m_LastExportedGridStats;
   bool
      m_DidExportSomething;


   void RunIterations();
   void ComputeWeightsResidualAndJacobian(FDenseVector &weights, FDenseVector &residual, FDenseMatrix &dK_dp_alpha);
   void ComputeStep(FScalar &CurMaxStep, FDenseVector &dx, FDenseVector const &residual, FDenseMatrix const &Jacobian, size_t iIt, FScalar fMaxStep_DynamicFactor);
   void UpdateSeedPositions(FDenseVector const &dx, FScalar const &fResidual, size_t iIt);
   
   // returns as final residual the maximum numerical integration error re-
   // *calculated for all* monomials in range l=0...lmax.
   FScalar FinishAndExportData();
   
   // project a 3*nSeeds vector of rotation parameters in such a way that the
   // corresponding rotations are compatible with the orbit types.
   void ProjectRotationsToOrbits(FDenseVector &dx);
   void CheckOrbitTypes();

public:   
   // returns final residual of the optimization
   FScalar Run();
};


#endif // SEARCH_SPHERE_GRID_H
