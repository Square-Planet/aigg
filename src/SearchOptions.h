#ifndef SEARCH_OPTIONS_H
#define SEARCH_OPTIONS_H

#include "Aigg.h"
#include <string>
#include <list>
#include <array>
#include "CxIntrusivePtr.h"
#include "CxPrintLevel.h"
#include "SymmetryGroups.h"

typedef std::vector<ptrdiff_t>
   FDegreeSchedule;
typedef std::vector<FPoint>
   FPointList;

namespace ct {
   struct string_slice;
}


enum FPrintFlags {
   PRINT_General,
   PRINT_TargetFn,
   PRINT_Step,
   PRINT_Timing,
   PRINT_Seeds,
   PRINT_NumPrintSettings
};

using ct::FPrintLevel;

struct FPrintOptions : public ct::TPrintOptions<PRINT_NumPrintSettings>
{
   typedef ct::TPrintOptions<PRINT_NumPrintSettings> FBaseClass;
   FPrintOptions();
   static unsigned ReadOptionId(ct::string_slice const &slName);
   // clears all print options; then sets more common non-specialty/debug options to StandardLevel
   void SetDefaults(FPrintLevel StandardLevel = FPrintLevel::Basic);
   void SetArgs(ct::string_slice const &Args);
};
   

// Data type for specifying how the program should behave if deviating from
// the default behavior may have significant runtime performance impacts
// when given a choice to either maximize accuracy or taking computational
// shortcuts
enum FAccuracyPreference {
   // Allow performing the operations with approximations or slightly decreased
   // accuracy, typically-given/expected-to-be-given conditions as assumptions,
   // or reduced diagnostic power, if doing so may yield significant reductions
   // in computational cost.
   //
   // For example, when verifying unit sphere integration rules with this focus,
   // the program will assume that all the grid points indeed lie on the sphere
   // (so only monomials of degree lmax and lmax-1 must be tested), that the grid
   // points have full point group symmetry, and that the FPointGroup object has
   // a correct implementation of the prescreening function for monomials. And
   // using these assumptions, it will prune the list of monomials which are
   // explicitly tested in the verification to the minimal subset which is
   // sufficient to verify the rule IF all of those conditions hold.
   ACCURACY_PreferSpeed = 0,
   // Do not deviate from the (rather conservative) defaults in choices betwen
   // accuracy and speed.
   //
   // For example, when verifying unit sphere integration rules with this focus,
   // the program will still assume that the integration grid as full point group
   // symmetry and FPointGroup's prescreening is correctly implemented, but it
   // will test all Cartesian monomials from order 0 up to lmax.
   ACCURACY_Default = 1,
   // Perform the operations with maximum accuracy, minimum assumptions, or
   // maximal diagnostic power, even if doing so comes at a significant
   // computational cost and the accuracy/coverage increase is not likely large
   // compared to ACCURACY_Default.
   //
   // For example, when verifying unit sphere integration rules with this focus,
   // the program will test every single Cartesian monomial up to degree lmax.
   ACCURACY_Maximize = 2
};
#ifdef INCLUDE_ABANDONED
// enum FComputePriority {
//    COMPUTE_PreferSpeed,
//    COMPUTE_NoPriority,
//    COMPUTE_PreferAccuracy
// };
#endif // INCLUDE_ABANDONED
namespace option_names { bool AccuracyPreferenceQ(ct::string_slice const &Name); }
// ^-- namespace just because the function names get sooo long without it
FAccuracyPreference ParseAccuracyPreference(ct::string_slice const &Key, ct::string_slice const &Value);


struct FExportOptions : public ct::FIntrusivePtrDest1
{
   std::string
      FileName;
   std::string
      // e.g., text, python, fortran, c++
      Format;
   bool
      // if not append: overwrite.
      AppendToFile;
   int
      // approx number of decimal digits to maximally export.
      // If larger than what FScalar can represent, it will be truncated.
      // -1 means "whatever "
      Precision;
   enum FPointExportDesc {
      // export only non-redundant information
      EXPORT_SeedsAndWeightsOnly,
      // export all points of the orbits
      EXPORT_FullOrbits
   };
   FPointExportDesc
      PointsToExport;
   bool
      StoreOptions, // echo options on command line used to create this grid into output file (text format only)
      StoreVersionInfo, // store info on program version (text format only)
      StoreIdentityInfo; // store info on host name, date, and user into output file (text format only)
   typedef std::list<std::string>
      FCommentList;
   FCommentList
      // a list of extra comment lines to append to the description of the output grid.
      // (in addition to the base grid desc and the comments described by Store*).
      //
      // Extra comments are normally added in one of two ways:
      // - With the { ...; export{ ...; comment: "this is an extra comment"}}
      // - When using aigg to refine an existing rule with initial points given
      //   via aigg grid file, i.e.,
      //   { ...; initial-points: read-grid-file{...}}
      //   as quotes of the original grid comments.
      ExtraComments;

   explicit FExportOptions(std::string const &Desc_);
};

typedef ct::TIntrusivePtr<FExportOptions>
   FExportOptionsPtr;


// Control details of the optimization procedure in the itereative grid
// optimization.
struct FOptimizeControlOptions : public ct::FIntrusivePtrDest1 {
   size_t
      MaxIt;
   FScalar
      fThreshResidual;
   FScalar
      // actual maximum step length:
      //   max(fMaxStep_Min, min(fMaxStep_Max, fMaxStep_ResidualFactor * fResidualRms))
      // if you want steps NOT scaling with the residual, set
      // fMaxStep_Max == fMaxStep_Min == (whatever step you want)
      fMaxStep_Max,
      fMaxStep_Min,
      fMaxStep_ResidualFactor,
      // each iteration, multiply fMaxStep with this. (typically used
      // to *enlarge* the maximum step with slowly). Setting this to 1.0
      // will disable the mechanism.
      fMaxStep_IterFactor,
      // after a bad step, revert the step, multiply current maximum step by this,
      // and then rety.
      fMaxStep_DynamicControl_FactorAfterBadStep,
      // number of iterations during which the maximum step width will gradually
      // reach its target value again following a bad step.
      fMaxStep_DynamicControl_NumIterToReset,
      Damping;
   bool
      // if set, absolute thresholds (fMaxStep_Max, fMaxStep_Min) will be scaled
      // by 1/sqrt(4pi/N) before being applied. 4pi/N is the average ``area'' represented
      // by a grid point, and its square root can be seen as a measure of distance.
      MaxStep_AdjustByNpts,
      // if set, revert bad steps, and resume with decreased maximum step with.
      MaxStep_DynamicControl;
   ptrdiff_t
      // if > 0, check for progress in optimization during the last nIterProgressCheck
      // iterations, and if there is not any, first turn off dynamic step length control.
      // And if that also does not help, abort the optimization.
      nIterProgressCheck,
      // if set, rather than immediately terminating the optimization after a
      // "no-progress" note, switch dynamic step length control between on and
      // off up to this number of times. This is dirty, but it can release
      // convergence locks on saddle points.
      ProgressCheck_DynControl_nMaxSwitch;

   explicit FOptimizeControlOptions();
   
   // returns 'true' if the option was recognized as belonging to convergence control,
   // false otherwise. FullOptionsDesc_ should be the full set of control options;
   // may be used for error reporting on input errors.
   bool TryParseAndSetOption(ct::string_slice Name, ct::string_slice Value, std::string const &FullOptionsDesc_);
};
typedef ct::TIntrusivePtr<FOptimizeControlOptions>
   FOptimizeControlOptionsPtr;

namespace option_names {
   enum FOptionNameFlags {
      AllowPlural = 0x0001
   };
   // tests if `Name` corresponds to an option name we support for defining the
   // active point group/list of point groups (e.g., "point-group", "symmetry",
   // "sym").
   // Recognizes flag option_names::AllowPlural; if set, it will return `true`
   // for names like 'point-groups:', 'symmetries:', etc.
   bool PointGroupQ(ct::string_slice const &Name, unsigned Flags=0);
   bool PrintQ(ct::string_slice const &Name);
   // tests if `Name` corresponds to an option we recognize for defining target
   // orders of integration rules (e.g., "order", "degree", "lmax") and/or
   // incremental optimization schedules
   // thereof (e.g., "schedule"). Recognizes flag option_names::AllowPlural
   bool RuleOrderQ(ct::string_slice const &Name, unsigned Flags=0);
}


struct FGridSearchOptions : public ct::FIntrusivePtrDest1
{
   // construct property object from a textual grid property list.
   explicit FGridSearchOptions(std::string const &Desc_);

   FPointGroupPtr
      pPointGroup;
   FOptimizeControlOptionsPtr
      pStepOptions;

   FAccuracyPreference
      // sets general guidelines on whether the program should focus on
      // obtaining the maximum obtainable accuracy in the current run
      // or do reasonable (and typically near-unnoticeable) approximations
      // where reasonable if it thinks this may significantly improve
      // runtime performance.
      //
      // Setting this switch via 'mode: fast' or 'mode: accurate' on input
      // or via SetAccuracyPreference() triggers a number of subsequent option
      // changes (see SetAccuracyPreference function) which are meant to
      // provide reasonable combinations
      AccuracyPreference;
   void SetAccuracyPreference(FAccuracyPreference const &AccuracyPreference_);

   bool
      // if set, use some explicit information about each point group to
      // select the set of target functions to optimize for at each l.
      // If not, all monomials in the target range will be considered.
      PrescreenTargetFunctions,
      // try to bypass the use of singular value decompositions as linear
      // algebra core routines. These can be very accurate and elegant,
      // but in high precision arithmetic they also can be slow.
      AvoidSvds,
      // if set, allow neglecting the dKe term in the jacobian computation.
      // The trust-region Newton updates typically work fine anyway, and it allows
      // bypassing a potentially cumbersome covariant projector calculation if
      // using tabulated polynomials to span the target space.
      NeglectResidualDerivs;
   
   ptrdiff_t
      // l at which to start optimizing grid points. If -1, ignore (start at l=1 for positive
      // optimization, then do schedule).
      lStart,
      // after each successful (positive lStep) or unsuccessful (negative lStep) optimization,
      // increase/decrease current l by this amount.
      // If 0, do not increase l, but leave it as it is.
      lStep;
   FDegreeSchedule
      // if lStep = 0: process this explicit list of angular momenta
      // otherwise: do incremental optimization first, and then process
      // this list *relative* to the last optimization point done (either
      // the first unsuccessful one (lStep > 0), or the first successful
      // one (lStep < 0).
      lList;
   std::string
      // if given, export the data of the initial point set in terms of
      // barycentric coordinates to this file in text format. This includes:
      // - the 2D barycentric coordinates themselves,
      // - the 3D Cartesian coordinates they have been translated to
      // - the Cartesian coordinates of points v0, v1, vc of the fundamental
      //   region used in the translation
      // - the formal order of vertices (01c/1c0/c01 etc.)
      FileName_BarycentricGuess;
   
   FTriangle::FBarycentricCoordType
      // decides whether to use planar or spherical geometry for
      // barycentric coordinate conversions.
      BarycentricGeomType;
   FPointList
      InitialPoints;
   FExportOptionsPtr
      pExportOptions;
   FScalar
      // grids will be considered as converged and (possibly export-worthy) if
      // the maximum error over the full set of target monomials (mxe) for the
      // current Lmax is smaller than this. Note that this quantity is
      // technically distinct from the residual threshold, which determines how
      // hard to converge the data set (...and therefore should be smaller than
      // this).
      //
      // The program may accept grids even if the optimization is not fully
      // converged, provided that the maximum error over monomials is smaller
      // than this.
      fThreshConsiderGridAsExact;
   int
      iPrintLevel;
   FPrintOptions
      PrintOptions;
   enum FVerifyRuleOptions {
      // test literally all monomials from 0 .. lmax to determine mxe.
      // (note: can get slow for very large rules)
      VERIFY_AllMonomials,
      // test only symmetry unique monomials which are not known to vanish due
      // to group symmetry (should give identical result to VERIFY_AllMonomials,
      // unless the FPointGroup::CouldMonomialBeATargetFunction() screening
      // function of the used point group is broken.).
      VERIFY_SymmetryMonomialsOnly,
      // as VERIFY_SymmetryMonomialsOnly, but additionally test only from
      // lmax-1 to lmax, instead of from 0 to lmax.
      VERIFY_MinimalMonomialSetOnly
   };
   FVerifyRuleOptions
      // Controls the cross-validation of converged integration used to
      // determine the 'mxe' value of the exported grids.
      VerifyRuleOptions;
   std::string MakeOptionString() const;

   enum FTargetFnType {
      TARGETFN_RawMonomials,
      TARGETFN_QuasiOrthMonomials,
      TARGETFN_TabulatedPolynomials
   };
   FTargetFnType
      TargetFnType;
   
   FScalar fThreshResidual() const { return pStepOptions->fThreshResidual; }
   FScalar MaxIt() const { return pStepOptions->MaxIt; }
protected:
   void ProcessStartingPoints(std::string const &StartingPointsDecl_);
   void MakeInitialBarycentricGrid(std::string const &Type, ptrdiff_t s, ptrdiff_t t, FScalar OffsetM, FScalar OffsetN, FScalar AreaDistortion, FScalar AreaDistortionSnap, std::string const &Order, std::string const &FullDeclForExport_);
   std::string
      m_OriginalOptionString;
   // meant to read point group and lmax from given file, where the file contains
   // a grid previously exported by aigg.
   void ReadMetaInfoFromAiggGridFileHeader(std::string FileName_);

   FExportOptions::FCommentList
      // set of comments on the grid generation which will be spliced into
      // pExportOptions->ExtraOptions at end of input processing.
      m_ExtraComments;
};

typedef ct::TIntrusivePtr<FGridSearchOptions>
   FGridSearchOptionsPtr;


#endif // SEARCH_OPTIONS_H
