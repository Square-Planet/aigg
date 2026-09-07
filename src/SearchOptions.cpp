#include "Aigg.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string/replace.hpp>
#include "SearchOptions.h"
#include "CxParse1.h"
#include "PointCloud.h"


static FScalar eps = 1e-8;


FPrintOptions::FPrintOptions()
{
   SetDefaults();
}


void FPrintOptions::SetDefaults(FPrintLevel StandardLevel)
{
   SetAll(FPrintLevel::Off);
   Set(PRINT_General, StandardLevel);
   Set(PRINT_TargetFn, StandardLevel);
   Set(PRINT_Step, StandardLevel);
   Set(PRINT_Timing, StandardLevel);
}


unsigned FPrintOptions::ReadOptionId(ct::string_slice const &slName) {
   if (slName.startswith("target-fn") || slName == "fn" || slName.startswith("function") || slName == "monomials" || slName == "target-space")
      return PRINT_TargetFn;
   else if (slName.startswith("step") || slName == "update")
      return PRINT_Step;
   else if (slName.startswith("timing") || slName == "time")
      return PRINT_Timing;
   else if (slName == "base" || slName == "general" || slName == "common")
      return PRINT_General;
   else if (slName == "seeds")
      return PRINT_Seeds;
   else
      return FPrintOptions::InvalidOption;
}

void FPrintOptions::SetArgs(ct::string_slice const &Args) { 
   FBaseClass::SetArgs(Args, ReadOptionId);
}


char const *FmtYesNo(bool b) {
   if (b)
      return "yes";
   else
      return "no";
}



FExportOptions::FExportOptions(std::string const &Desc_)
   : FileName(""), Format("text"), AppendToFile(false), Precision(-1),
     PointsToExport(EXPORT_SeedsAndWeightsOnly),
     StoreOptions(true), StoreVersionInfo(true), StoreIdentityInfo(true)
{
   if (!Desc_.empty() && Desc_[0] != '{') {
      // not a property list. interpret as filename and leave other options alone.
      FileName = Desc_;
   } else {
      ct::FPropertyListStr
         Props(Desc_);
      for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
         if (it->Name == "file" || it->Name == "filename" || it->Name == "file-name") {
            FileName = it->Content.to_str();
         } else if (it->Name == "format") {
            Format = it->Content.to_str();
            if (!(Format == "text"))
               throw std::runtime_error(fmt::format("FExportOptions: export format '{}' not (yet?) supported. Know: 'text'.", it->Content.to_str()));
         } else if (it->Name == "prec" || it->Name == "digits" || it->Name == "precision") {
            if (it->Content == "full" || it->Content == "max")
               // we're asked to export all digits which are representable in FScalars.
               // just set number of digits to something very large. FmtExportFloat
               // will truncate it to a sensible level.
               Precision =  0x7fffffff;
            else
               Precision = it->Content.to_int();
         } else if (it->Name == "append") {
            AppendToFile = it->Content.to_bool();
            if (AppendToFile)
               throw std::runtime_error("FExportOptions: append option not yet supported. FIXME.");
         } else if (it->Name == "points") {
            if (it->Content == "all") {
               PointsToExport = EXPORT_FullOrbits;
            } else if (it->Content == "seeds") {
               PointsToExport = EXPORT_SeedsAndWeightsOnly;
            } else {
               throw std::runtime_error(fmt::format("FExportOptions: points option '{}' not recognized. Know: 'all' and 'seeds'.", it->Content.to_str()));
            }
         } else if (it->Name == "options") {
            StoreOptions = it->Content.to_bool();
         } else if (it->Name == "version") {
            StoreVersionInfo = it->Content.to_bool();
         } else if (it->Name == "identity") {
            StoreIdentityInfo = it->Content.to_bool();
         } else if (it->Name == "comment") {
            ExtraComments.push_back(it->Content.to_str(ct::string_slice::TOSTR_EvalStringLiteral));
         } else {
            throw std::runtime_error(fmt::format("FExportOptions: export option '{}' not recognized. Know 'file', 'format', 'append', 'points'.", it->Content.to_str()));
         }
      }
   }
   if (!FileName.empty()) {
      std::cout << fmt::format("! Will export final grid to file '{}' "
         "{{format: {}; append: {}; points: {}; options: {}; version: {}; identity: {}}}", FileName, Format,
         FmtYesNo(AppendToFile), (PointsToExport==EXPORT_FullOrbits)? "all" : "seeds",
         FmtYesNo(StoreOptions), FmtYesNo(StoreVersionInfo), FmtYesNo(StoreIdentityInfo))
         << std::endl;
   }
}


static void ParseMaxStep(FScalar &ScaleOrStep, bool &RelativeToResidual, ct::string_slice s)
{
   ct::split_result sr;
   s.split(sr, '*');
   if (sr.empty())
      throw std::runtime_error("max-step, if given, cannot be empty. Expected <float> or <float>*residual.");
   if (sr.size() > 2)
      throw std::runtime_error(fmt::format("max-step format must be <float> or '<float>*residual', but got '{}'. Too many '*'s.", s.to_str()));
   if (sr.size() == 2) {
      if (sr[1] == "res" || sr[1] == "residual") {
         RelativeToResidual = true;
         ScaleOrStep = FScalar(sr[0].to_float());
      } else {
         throw std::runtime_error(fmt::format("failed to understand max-step declaration '{}'. Must be <float> or '<float>*residual'.", s.to_str()));
      }
   } else {
      RelativeToResidual = false;
      ScaleOrStep = s.to_float();
   }
}


FOptimizeControlOptions::FOptimizeControlOptions()
   : MaxIt(2000), 
     fMaxStep_Max(5e-4), fMaxStep_Min(1e-8), fMaxStep_ResidualFactor(3e4),
     fMaxStep_IterFactor(1.01), fMaxStep_DynamicControl_FactorAfterBadStep(0.1), fMaxStep_DynamicControl_NumIterToReset(8.0),
     Damping(0.), MaxStep_AdjustByNpts(true), MaxStep_DynamicControl(false),
     nIterProgressCheck(100), ProgressCheck_DynControl_nMaxSwitch(3)
{
#ifdef INCLUDE_ABANDONED
// #ifdef SCALAR_IS_FLOAT_128
//    fThreshResidual = FScalar(1e-20Q);
// #else
//    fThreshResidual = 1e-11;
// #endif
// #ifdef SCALAR_IS_FLOAT_128
//    fThreshResidual = FScalar(1e-25Q);
// #else
//    fThreshResidual = 1e-12;
// #endif
#endif // INCLUDE_ABANDONED
   fThreshResidual = FScalar(1e1) * g_ThrVerySmall;
   // ^- this means 'set something reasonable unless the user has provided an explicit value'
}


bool FOptimizeControlOptions::TryParseAndSetOption(ct::string_slice Name, ct::string_slice Value, std::string const &FullOptionsDesc_)
{
   // hack around leaving the "it-> ..." from when this was moved from a larger option parsing function.
   ct::FPropertyStr
      PropStr;
   PropStr.Name = Name;
   PropStr.Content = Value;
   ct::FPropertyStr const
      *it = &PropStr;
      
   if (it->Name == "thrres" || it->Name == "thr-res" || it->Name == "thr" || it->Name == "thresh") {
      // threshold for convergence; by default applies to rmsd of error of
      // orthogonalized target function integrals (res).
      // From this the threshold for mxe (maximum absolute error of integration rule
      // over *all* monomials up to lmax) is also derived.
      fThreshResidual = FScalar(it->Content.to_float());
   } else if (it->Name == "maxit" || it->Name == "max-it" || it->Name == "max-iter" || it->Name == "maxiter") {
      MaxIt = it->Content.to_int();
   } else if (it->Name == "progress-check" || it->Name == "progress-check-it" || it->Name == "it-progress-check") {
      nIterProgressCheck = it->Content.to_int();
   } else if (it->Name == "progress-check-dyncontrol-switch" || it->Name == "progress-check-switch") {
      ProgressCheck_DynControl_nMaxSwitch = it->Content.to_int();
   } else if (it->Name == "damping" || it->Name == "damp") {
      Damping = it->Content.to_float();
   } else if (it->Name == "max-step-adjust-by-npts" || it->Name == "maxstep-adjust-by-npts" || it->Name == "msabn" || it->Name == "mabn") {
      MaxStep_AdjustByNpts = it->Content.to_bool();
   } else if (it->Name == "maxstep" || it->Name == "max-step") {
      // options:
      //   - a single float -> use this as fixed max-step
      //   - <float>*residual
      //   - <float>*residual;min:<float>;max:<float>
      if (it->Content.empty())
         throw std::runtime_error(fmt::format("error in grid declaration '{}'. {}, if given, cannot be empty.", FullOptionsDesc_, it->Name.to_str()));
      FScalar
         ScaleOrStep;
      bool
         RelativeToResidual;
      if (it->Content[0] != '{') {
         // single float or float*residual case.
         ParseMaxStep(ScaleOrStep, RelativeToResidual, it->Content);
      } else {
         ct::FPropertyListStr
            StepProps(it->Content, 0, ct::PROPLIST_AllowEntriesWithoutName);
         for (ct::FPropertyListStr::const_iterator itStepProp = StepProps.begin(); itStepProp != StepProps.end(); ++ itStepProp) {
            if (itStepProp->Name == "") {
               ParseMaxStep(ScaleOrStep, RelativeToResidual, itStepProp->Content);
            } else if (itStepProp->Name == "min") {
               fMaxStep_Min = FScalar(itStepProp->Content.to_float());
            } else if (itStepProp->Name == "max") {
               fMaxStep_Max = FScalar(itStepProp->Content.to_float());
            } else {
               throw std::runtime_error(fmt::format("max-step: entry '{}' of property list '{}' not recognized.", itStepProp->Name.to_str(), it->Content.to_str()));
            }
         }
      }
      if (RelativeToResidual) {
         fMaxStep_ResidualFactor = ScaleOrStep;
      } else {
         fMaxStep_ResidualFactor = 0.;
         fMaxStep_Min = ScaleOrStep;
         fMaxStep_Max = ScaleOrStep;
      }
   } else if (it->Name == "max-step-it-factor" || it->Name == "max-step-iter-factor" || it->Name == "msif") {
      fMaxStep_IterFactor = it->Content.to_float();
   } else if (it->Name == "dynamic" || it->Name == "dynamic-step-control") {
      MaxStep_DynamicControl = it->Content.to_bool();
   } else if (it->Name == "dynamic-bad-step-factor" || it->Name == "step-length-factor-after-bad-step" || it->Name == "step-adjust-after-bad-step") {
      fMaxStep_DynamicControl_FactorAfterBadStep = it->Content.to_float();
   } else if (it->Name == "dynamic-relax" || it->Name == "dynamic-relax-it" || it->Name == "iter-to-relax-after-bad-step" || it->Name == "relax-iter-after-bad-step") {
      fMaxStep_DynamicControl_NumIterToReset = it->Content.to_float();
   } else {
      // option not recognized.
      return false;
   }
   // option was recognized.
   return true;
}



void FGridSearchOptions::SetAccuracyPreference(FAccuracyPreference const &AccuracyPreference_)
{
   this->AccuracyPreference = AccuracyPreference_;
   size_t iAccu = size_t(this->AccuracyPreference);
   cx_assert_rt(iAccu >= 0 && iAccu <= 2);
   assert(ACCURACY_PreferSpeed == 0 && ACCURACY_Default == 1 && ACCURACY_Maximize == 2);
   {
      FVerifyRuleOptions const vro[] = { VERIFY_MinimalMonomialSetOnly, VERIFY_SymmetryMonomialsOnly, VERIFY_AllMonomials };
      this->VerifyRuleOptions = vro[iAccu];
   }
   {
      FPrintLevel const tfp[] = { FPrintLevel::Off, FPrintLevel::Basic, FPrintLevel::More };
      this->PrintOptions[PRINT_TargetFn] = tfp[iAccu];
   }
   {
      bool const nrd[] = { true, false, false };
      this->NeglectResidualDerivs = nrd[iAccu];
   }
}


namespace option_names {
   bool AccuracyPreferenceQ(ct::string_slice const &Name) {
      return Name == "prefer" || Name == "compute-mode" || Name == "accuracy";
   }
}

FAccuracyPreference ParseAccuracyPreference(ct::string_slice const &Key, ct::string_slice const &Value)
{
   if (Key == "prefer") {
      if (Value == "accuracy")
         return ACCURACY_Maximize;
      else if (Value == "neither")
         return ACCURACY_Default;
      else if (Value == "speed" || Value == "performance")
         return ACCURACY_PreferSpeed;
      else
         throw std::runtime_error(fmt::format("input error in accuracy-vs-speed preference: value in '{}': '{}' not recognized; should be one 'accuracy', 'speed', or 'neither'", Key, Value));
   } else if (Key == "compute-mode") {
      if (Value == "accurate")
         return ACCURACY_Maximize;
      else if (Value == "default" || Value == "normal")
         return ACCURACY_Default;
      else if (Value == "fast")
         return ACCURACY_PreferSpeed;
      else
         throw std::runtime_error(fmt::format("input error in accuracy-vs-speed preference: value in '{}': '{}' not recognized; should be one 'accurate', 'fast', or 'default'", Key, Value));
   } else if (Key == "accuracy") {
      if (Value == "max" || Value == "best")
         return ACCURACY_Maximize;
      else if (Value == "default" || Value == "normal" || Value == "good")
         return ACCURACY_Default;
      else if (Value == "fair")
         return ACCURACY_PreferSpeed;
      else
         throw std::runtime_error(fmt::format("input error in accuracy-vs-speed preference: value in '{}': '{}' not recognized; should be one 'best', 'good', or 'fair'", Key, Value));
   } else {
      throw std::runtime_error(fmt::format("programming error: option '{}' (set to '{}') was not recognized as a accuracy-vs-speed preference option; option name should be 'prefer:', 'compute-mode:', or 'accuracy:'", Key, Value));
   }
}


namespace option_names {
   bool PointGroupQ(ct::string_slice const &Name, unsigned Flags) {
      if (Name == "point-group" || Name == "pointgroup" || Name == "group" || Name == "sym" || Name == "symmetry")
         return true;
      if (bool(Flags & AllowPlural) && (Name == "point-groups" || Name == "groups" || Name == "symmetries"))
         return true;
      return false;
   }

   bool PrintQ(ct::string_slice const &Name) {
      return Name == "print" || Name == "print-level";
   }

   bool RuleOrderQ(ct::string_slice const &Name, unsigned Flags) {
      if (Name == "order" || Name == "l" || Name == "lmax" || Name == "schedule" || Name == "degree")
         return true;
      if (bool(Flags & AllowPlural) && (Name == "orders" || Name == "ls" || Name == "degrees"))
         return true;
      return false;
   }
}


FGridSearchOptions::FGridSearchOptions(std::string const &Desc_)
   : pPointGroup(new FIcosahedralGroup()), 
     pStepOptions(new FOptimizeControlOptions()),
     AccuracyPreference(ACCURACY_Default),
     PrescreenTargetFunctions(true),
     AvoidSvds(false),
     NeglectResidualDerivs(false),
     lStart(-1), lStep(+1),
     iPrintLevel(0),
     VerifyRuleOptions(VERIFY_AllMonomials)
{
   BarycentricGeomType = FTriangle::GEOMTYPE_Spherical;
   TargetFnType = TARGETFN_QuasiOrthMonomials;
//    BarycentricGeomType = FTriangle::GEOMTYPE_Planar; // FIXME: remove this
   m_OriginalOptionString = Desc_;
   fThreshConsiderGridAsExact = FScalar(-1);

   ct::FPropertyListStr
      Props(Desc_);
   for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
      if (option_names::PointGroupQ(it->Name)) {
         pPointGroup = MakePointGroup(it->Content.to_str());
#ifdef INCLUDE_ABANDONED
//          if (it->Content == "icosahedral" || it->Content == "ico")
//             pPointGroup = new FIcosahedralGroup();
//          else if (it->Content == "tetrahedral" || it->Content == "tetra")
//             pPointGroup = new FTetrahedralGroup();
//          else if (it->Content == "octahedral" || it->Content == "octa")
//             pPointGroup = new FOctahedralGroup();
//          else
//             throw std::runtime_error(fmt::format("point group '{}' not recognized. Know 'icosahedral', 'tetrahedral', 'octahedral'.", it->Content.to_str()));
#endif // INCLUDE_ABANDONED
      } else if (it->Name == "barycentric" || it->Name == "baryc" || it->Name == "barycentric-coords" || it->Name == "baryc-coords") {
         if (it->Content == "planar" || it->Content == "plane")
            BarycentricGeomType = FTriangle::GEOMTYPE_Planar;
         else if (it->Content == "spherical" || it->Content == "sph")
            BarycentricGeomType = FTriangle::GEOMTYPE_Spherical;
         else
            throw std::runtime_error(fmt::format("error in specification of barycentric coordinate type '{}': must be either 'planar' or 'spherical', but was set to '{}'", it->Name.to_str(), it->Content.to_str()));
      } else if (option_names::PrintQ(it->Name)) {
         ptrdiff_t iSinglePrintLevel;
         if (it->Content.try_convert_to_int(&iSinglePrintLevel)) {
            iPrintLevel = iSinglePrintLevel;
            PrintOptions.SetDefaults(iPrintLevel);
         } else {
            PrintOptions.SetArgs(it->Content);
            iPrintLevel = int(PrintOptions[PRINT_General]);
         }
      } else if (it->Name == "target-fns" || it->Name == "target-fn" || it->Name == "fns" || it->Name == "target" || it->Name == "target-space") {
         if (it->Content == "simple" || it->Content == "raw-monomials")
            TargetFnType = TARGETFN_RawMonomials;
         else if (it->Content == "approx" || it->Content == "quasi-orth")
            TargetFnType = TARGETFN_QuasiOrthMonomials;
         else if (it->Content == "exact" || it->Content == "poly" || it->Content == "polynomials" || it->Content == "tabulated")
            TargetFnType = TARGETFN_TabulatedPolynomials;
         else
            throw std::runtime_error(fmt::format("error in specification of target function space '{}': options '{}' not recognized", it->Name.to_str(), it->Content.to_str()));
      } else if (it->Name == "verify") {
         if (it->Content == "all" || it->Content == "full")
            VerifyRuleOptions = VERIFY_AllMonomials;
         else if (it->Content == "sym" || it->Content == "symmetry-required" || it->Content == "medium" || it->Content == "intermediate")
            VerifyRuleOptions = VERIFY_SymmetryMonomialsOnly;
         else if (it->Content == "minimal")
            VerifyRuleOptions = VERIFY_MinimalMonomialSetOnly;
         else
            throw std::runtime_error(fmt::format("error in specification of integration rule verification option '{}': must be either 'full' or 'symmetry-required', or 'minimal', but was set to '{}'", it->Name.to_str(), it->Content.to_str()));
      } else if (option_names::AccuracyPreferenceQ(it->Name)) {
         SetAccuracyPreference(ParseAccuracyPreference(it->Name, it->Content));
      } else if (it->Name == "thr-exact" || it->Name == "thr-export") {
         fThreshConsiderGridAsExact = FScalar(it->Content.to_float());
      } else if (it->Name == "screen-target-fn" || it->Name == "prescreen-target-fn") {
         PrescreenTargetFunctions = it->Content.to_bool();
      } else if (it->Name == "avoid-svds") {
         AvoidSvds = it->Content.to_bool();
      } else if (option_names::RuleOrderQ(it->Name, option_names::AllowPlural)) {
         // Definition of a grid optimization target integration order and/or optimization schedule.
         // There are several supported versions of this:
         // - a single integer (in this case we will process only this particular lmax)
         // - a list of either integers, or [<lstart>,step:+1,<l1>,<l2>,...] specifying
         //   sequences of target integration rule orders (lmax).
         // In the latter case, this defines an incremental optimization schedule, in which in each
         // step the result from the last optimization is ued as starting guess for the next.
         if (it->Content.empty())
            throw std::runtime_error(fmt::format("error in grid declaration '{}'. {}, if given, cannot be empty.", Desc_, it->Name.to_str()));
         if (it->Content[0] != '[') {
            // lmax is not a list. Assume single integer case.
            lStep = 0;
            if (it->Content != "cur")
               lStart = it->Content.to_int();
            lList.clear();
         } else {
            // starts with '[' -> list of stuff case. See what we got.
            assert(it->Content[0] == '[');
            ct::long_split_result
               sr;
            it->Content.split_list(sr);
            ptrdiff_t lStartOriginal = lStart; // might have come from earlier def or read from grid file
            lStep = 0;
            lStart = -1;
            lList.clear();
            for (size_t il = 0; il != sr.size(); ++ il) {
               // current entry is a incremental schedule declaration?
               if (sr[il].startswith("step")) {
                  ct::split_result
                     StepSplit;
                  sr[il].split(StepSplit,':');
                  if (StepSplit.size() != 2)
                     throw std::runtime_error(fmt::format("expected: 'step:<integer>', but got '{}'", sr[il].to_str()));
                  lStep = StepSplit[1].to_int();
                  if (!(lList.empty() || lList.size() == 1))
                     throw std::runtime_error(fmt::format("error in grid order schedule '{}': must be "
                        "either [int,int,int,int,...] or [int,step:int,int,int,...]. But got already {} "
                        "grid order integers before occurrence of 'step'.", it->Content.to_str(), lList.size()));
                  assert(lList.empty() || lList.size() == 1);
                  if (lList.size() == 1) {
                     // first integer in incremental schedule declaration is lmin... where to
                     // start the incremental search before the rest of the schedule.
                     lStart = lList[0];
                     lList.clear();
                  }
               } else {
                  // schedule entry does not start with 'step'. Assume integer and append to list.
                  if (sr[il] == "cur") {
                     lList.push_back(lStartOriginal);
//                      std::cout << fmt::format("\n\n ---- lStartOriginal = {}\n\n", lStartOriginal);
                     // ^- note: depending on argument order, this either may or
                     //    may not have been set by this point if input data is
                     //    read from a grid file. In the latter case, we leave
                     //    lStart at -1 here (that's its default value) and
                     //    later replace it with the grid file data.
                  } else
                     lList.push_back(sr[il].to_int());
               }
            }
            if (lStart == -1 && lStep == 0 && lList.size() == 1) {
               // first integer in incremental schedule declaration is lmin... where to
               // start the incremental search before the rest of the schedule.
               lStart = lList[0];
               lList.clear();
            }
         }
      } else if (it->Name == "starting-points" || it->Name == "initial-points" || it->Name == "start" || it->Name == "points" || it->Name == "guess") {
         // this is a bit messy... put it into an extra function.
         ProcessStartingPoints(it->Content.to_str());
      } else if (it->Name == "export") {
         pExportOptions = new FExportOptions(it->Content.to_str());
      } else if (pStepOptions->TryParseAndSetOption(it->Name, it->Content, Desc_)) {
         // was handled as step control options. done here.
      } else {
         throw std::runtime_error(fmt::format("grid search property '{}' not recognized.", it->Name.to_str()));
      }
   }
   
   if (fThreshConsiderGridAsExact == -1) {
      // decide on a "good enough" threshold for considering a grid as
      // exact of order L.
      
      // By default, something close to integrating double precision data
      // exactly appears acceptable (note that mxe is exported; this grid may
      // still be used as initial guess for a refine operation to get it to the
      // full double precision if needed)
      fThreshConsiderGridAsExact = FScalar(1e-14);
      
      // however, this number *should* be larger than the residual threshold.
      // So if a loose residual threshold is set, increase this one, too.
      FScalar fExtra(1e3);
      if (fExtra*fThreshResidual() > fThreshConsiderGridAsExact)
         fThreshConsiderGridAsExact = fExtra*fThreshResidual();
   }

   if (pExportOptions.get()) {
      // splice extra comments which we possibly could have encountered into
      // the export options.
      pExportOptions->ExtraComments.splice(pExportOptions->ExtraComments.end(), m_ExtraComments);
   }
}


static FPointList ReadPointListFromFile(std::string const &FileName)
{
   std::ifstream
      inp(FileName.c_str());
   std::string
      line;
   std::getline(inp, line); // skip first line.
   if (inp.bad() || inp.fail() || !inp.good())
      throw std::runtime_error(fmt::format("failed to read from input file '{}'", FileName));
   FPointList
      PointList;
   while (std::getline(inp, line)) {
      ct::long_split_result
         sr;
      ct::string_slice
         line_as_slice(line);
      line_as_slice.trim();
      if (line_as_slice.startswith("#") || line_as_slice.startswith("//"))
         continue;
      line_as_slice.split(sr, ' ');
//       std::cout << fmt::format(": {} [{}]", line, sr.size()) << std::endl;
      if (sr.size() < 3)
         break;
      try {
         FScalar
            x = ParseFloat(sr[0].to_str()),
            y = ParseFloat(sr[1].to_str()),
            z = ParseFloat(sr[2].to_str());
         PointList.push_back(FPoint(x,y,z));
      } catch (std::runtime_error &e) {
         if (std::string(e.what()).find("floating point number") != std::string::npos)
            // read something which was not a float. end of list.
            break;
         else
            throw e;
      }
   }
   if (PointList.empty())
      throw std::runtime_error(fmt::format("Could not read any points from  file '{}'. Expected format: (first line ignored), then <float> <float> <float> <irest ignored> per line.", FileName));
   std::cout << fmt::format("! Read {} initial points from file '{}'", PointList.size(), FileName) << std::endl;
   return PointList;
}



std::string FGridSearchOptions::MakeOptionString() const
{
   std::string s = m_OriginalOptionString;
   // if data came from file, it might have contained newline characters and/or tabs.
   // replace these to get a one-line string for export purposes.
   boost::replace_all(s, "\n", " ");
   boost::replace_all(s, "\t", " ");
   boost::replace_all(s, "  ", " ");
   return s;
}


         
std::string FmtVec3_DoublePrec(FVector3 const &v, int n = -1)
{
   if (n==6)
      return fmt::format("[{:+.6f}, {:+.6f}, {:+.6f}]", double(v[0]), double(v[1]), double(v[2]));
   return fmt::format("[{:+.16e}, {:+.16e}, {:+.16e}]", double(v[0]), double(v[1]), double(v[2]));
}

std::string FmtTriangle_DoublePrec(FTriangle const &t)
{
//    return fmt::format("{{v0={},v1={},v2={}}}", FmtVec3_DoublePrec(t.v0), FmtVec3_DoublePrec(t.v1), FmtVec3_DoublePrec(t.v2));
   return fmt::format("[{}, {}, {}]", FmtVec3_DoublePrec(t.v0), FmtVec3_DoublePrec(t.v1), FmtVec3_DoublePrec(t.v2));
}



// tests whether a point given in terms of the barycentric coordinates of one
// triangle (BaseRegion) lies inside another triangle (FundamentalRegion).
struct FTestInsideFundamentalRegionPred
{
   explicit FTestInsideFundamentalRegionPred(FTriangle const &BaseRegion, FTriangle const &FundamentalRegion, FTriangle::FBarycentricCoordType GeomType);

   bool operator()(FVector3 const& vBarycPt_BaseRegion) const;
protected:
   FTriangle
      m_BaseRegion,
      m_FundamentalRegion;
   FTriangle::FBarycentricCoordType
      m_GeomType;
   bool
      m_NeedConverstion;
};

FTestInsideFundamentalRegionPred::FTestInsideFundamentalRegionPred(FTriangle const &BaseRegion, FTriangle const &FundamentalRegion, FTriangle::FBarycentricCoordType GeomType)
   : m_BaseRegion(BaseRegion), m_FundamentalRegion(FundamentalRegion), m_GeomType(GeomType)
{
   m_NeedConverstion = BaseRegion != FundamentalRegion;
}

   
bool FTestInsideFundamentalRegionPred::operator()(FVector3 const& vBarycPt_BaseRegion) const
{
   if (!m_NeedConverstion) {
      // base region and fundamental region identical. Test directly in terms of base region.
      return m_BaseRegion.IsInside_Barycentric(vBarycPt_BaseRegion, eps);
   } else {
      // base region and fundamental region differ. So first convert
      // barcycentric coordinates given in terms of the base region into
      // fundamental region, and then test there.
      return m_FundamentalRegion.IsInside_Cartesian(m_BaseRegion.MakeCartesian(vBarycPt_BaseRegion, m_GeomType), eps, m_GeomType);
   }
}


static void MakeHexGridPointList(FPointList &BarycPoints, std::string const &Type, ptrdiff_t s, ptrdiff_t t, FScalar OffsetM, FScalar OffsetN, std::string const &Order, FTestInsideFundamentalRegionPred const *pIsInsideFundamentalRegionFn)
{
   BarycPoints.clear();
   for (ptrdiff_t m_ = -100; m_ <= 100; ++ m_) {
      for (ptrdiff_t n_ = -100; n_ <= 100; ++ n_) {
         FScalar
            m = FScalar(m_) + OffsetM,
            n = FScalar(n_) + OffsetN;
         // compute barycentric coordinates with respect to P0, P1, Pc (see pg. 14)
         FScalar
            tau0,
            tau1,
            tauc;
         if (Type == "hex-grid") {
            tau0 = (s*s + (s-t)*m + s*t + t*t - (s+2*t)*n)/FScalar(s*s + s*t + t*t),
            tau1 = ((2*m + n)*s + (m-n)*t)/FScalar(s*s + s*t + t*t),
            tauc = 3*(n*t - m*s)/FScalar(s*s + s*t + t*t);
         } else {
            throw std::runtime_error(fmt::format("FGridSearchOptions::MakeInitialBarycentricGrid: grid type '{}' not recognized.", Type));
         }

         if (Order == "01c") {
         } else if (Order == "10c") {
            std::swap(tau0, tau1);
         } else if (Order == "0c1") {
            std::swap(tau1, tauc);
         } else if (Order == "c10") {
            std::swap(tau0, tauc);
         } else if (Order == "1c0") {
            std::swap(tauc, tau1);
            std::swap(tau1, tau0);
         } else if (Order == "c01") {
            std::swap(tauc, tau0);
            std::swap(tau0, tau1);
         } else {
            throw std::runtime_error(fmt::format("barycentric coordinate permutation '{}' not recognized. Should be permutation of '0', '1', and 'c'", Order));
         }

         bool inside;
         if (pIsInsideFundamentalRegionFn) {
            inside = (*pIsInsideFundamentalRegionFn)(FVector3(tau0, tau1, tauc));
            // ^- note: we can't necessarily just test for 01c region here
            // directly, because there is at least one group (Th) for which 01c
            // is *smaller* than the fundamental region.
            //
            // FIXME: check if this is really okay with spherical coordinates
            // and Th. I would almost expect that the answer is "no".
            // Should at least add a switch to prevent using spherical
            // barycentric coordinates in combination with Th in this case.
         } else {
            // For all other groups, here only testing if all tau's (in terms of
            // the 01c region) are >= 0 is sufficent.
            inside = tau0 >= -eps && tau1 >= -eps && tauc >= -eps;
         }
         if (!inside)
            continue;
         // NOTE:
         // - as of yet, this point list is not yet final!
         // - Additional tests for inclusion in the FundamentalRegion and having
         //   only points which are symmetry-unique are done in
         //   MakeInitialBarycentricGrid.
         //   (...but only on the reduced set of points lying inside the 01c
         //   region which we let out of here)
         BarycPoints.push_back(FPoint(tau0, tau1, tauc));
      }
   }
}

#ifdef INCLUDE_ABANDONED
// static FPoint Baryc_v01c_from_v012(FPoint vBarycPt_v012)
// {
//    // we need to generate points with respect to the barycentric coordinates
//    // 0,1,c; for an icosahedron, 0 and 1 are vertices and c is a face center:
//    //
//    //   m_Vertex0 = v0;
//    //   m_Vertex1 = v1;
//    //   m_VertexC = (FScalar(1)/FScalar(3))*(v0+v1+v2);
//    //
//    // but for a subdiv triangle (with area distortion), we'd like to generate
//    // with three vertices of an equilateral triangle (012), not with 01c.
//    FScalar f = 1/FScalar(3);
//    return FPoint(vBarycPt_v012[0] - f*vBarycPt_v012[2], vBarycPt_v012[1] - f*vBarycPt_v012[2], (FScalar(1) + 2*f)*vBarycPt_v012[2]);
// }
#endif // INCLUDE_ABANDONED

#ifdef INCLUDE_ABANDONED
// static void MakeTriangleSubdivGrid(FPointList &BarycPoints, std::string const &Type, ptrdiff_t s, ptrdiff_t t, FScalar OffsetM, FScalar OffsetN, std::string const &Order, FTestInsideFundamentalRegionPred const &IsInsideFundamentalRegionFn)
#endif // INCLUDE_ABANDONED
static void MakeTriangleSubdivGrid(FPointList &BarycPoints, std::string const &Type, ptrdiff_t s, ptrdiff_t t, FScalar OffsetM, FScalar OffsetN, std::string const &Order)
{
   // make && aigg '{point-group: icosahedral; max-it: 128; degree: [4,step:+2,-1,-2]; max-step: {1e8*residual; max:2e-3; min:1e-8}; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: tri-subd{2;0;01c}}'
   // ^- doesn't work yet.
   // make && aigg '{point-group: icosahedral; max-it: 128; degree: [4,step:+2,-1,-2]; max-step: {1e8*residual; max:1e-3; min:1e-8}; damping:0.0;  export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: tri-subd{3;0;01c}}'
   // ^- this one kind of does now, with the modified subdiv fix....
   //
   // Issues:
   // - wrong orbit length... that might be a orbit with a seed of the wrong type.
   //   some of the non-general orbit types require a specific set of point of the
   //   symmetry unique set as input; and we might get the wrong one.
   // - seed collapse.... that is odd. shouldn't happen?

   BarycPoints.clear();
   int
      N = int(s);
//    for (unsigned iTau01 = 0; iTau01 <= N; ++iTau01) {
   for (int iTau01 = -abs(t); iTau01 <= N + abs(t); ++iTau01) {
      int
         iTau2 = N - iTau01;
      for (int iTau0 = 0; iTau0 <= iTau01; ++iTau0) {
         int
            iTau1 = iTau01 - iTau0;
         FScalar
            f = 1/FScalar(N);
         FPoint
            vBarycPt_v012(f*iTau0, f*iTau1, f*iTau2);
         if (1) {
            FScalar
               g = (FScalar(t) + OffsetN)/FScalar(2*N);
            g *= iTau0 / FScalar(N-1);
            vBarycPt_v012 += FPoint(0,g,-g);
         }
#ifdef INCLUDE_ABANDONED
//          if (1) {
//             vBarycPt_v012 = FPoint(f*iTau0, f*((ptrdiff_t(iTau1)+t+100*(N+1))%(N+1)), f*iTau2);
//             vBarycPt_v012 /= vBarycPt_v012.sum();
//          }
//          if (!IsInsideFundamentalRegionFn(vBarycPt_v012))
//             continue;
         // ^- test moved to main MakeInitialBarycentricGrid function.
#endif // INCLUDE_ABANDONED
         std::cout << fmt::format("--- triang subdiv baryc: itau0 = {}  itau1 = {}  itau2 = {}  N = {}  -> vBaryc = {}  sum(vBaryc) = {}",
            iTau0, iTau1, iTau2, N, FmtVec3_DoublePrec(vBarycPt_v012), vBarycPt_v012.sum()) << std::endl;
         BarycPoints.push_back(vBarycPt_v012);
      }
   }
}

#ifdef INCLUDE_ABANDONED
// #ifdef INCLUDE_OPTIONALS
// static FScalar fAreaTraf1D(FScalar taux, FScalar p, int role)
// {
//    if (role == 1)
//       taux = 1 - taux;
//    
//    // shouldn't I want a derivative of wx at F[1]? Let's try that one.
//    // Solve[F[0]==0&&F[1]==1&&Derivative[1][F][0]==1&&Derivative[1][F][1]==p,{a0,a1,a2,a3}]
//    // {{a0->0,a1->1,a2->1-p,a3->-1+p}}
//    // Plot[(t+(1-p)*t^2-(1-p)*t^3)/.p->2,{t,0,1}]
//    // Plot[t*(1+(1-p)*t*(1-t))/.p->2,{t,0,1}]
//    // --> formula look okay, I guess
//    //
//    // UPDATE:
//    // - That is not really right, because I am also scaling the other points...
//    // - hmpf... for T even the v0-v1 distortion is too large to ignore
//    //   see:
//    //   make && aigg '{ point-group: tetrahedral; degree: [2]; max-step: 1e-3; export: {file: /tmp/test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: save{_initial_points.dat}; initial-points: hex-grid{8;0;01c;dist:1.0}; print: 0; max-it: 4096; dynamic:off; progress-check:0; damping:0.2 }'
//    //   and combine with visualize-rule.sh /tmp/test_run.dat
//    FScalar
//       fSign = (taux < 0)? -1 : 1;
//    FScalar
//       t = fSign * taux,
//       new_taux = fSign * (t*(1 + (1-p)*t*(1 - t)));
// //    //          new_taux = fSign * (1-pow(1-t,p));
// //    //          new_taux = fSign * (t*(1 + t*(1 - tc)));  // 
// //       // F[t_]:=a0 + a1*t + a2*t^2+a3*t^3
// //       // Solve[F[0]==0&&F[1]==1&&Derivative[1][F][0]==1&&Derivative[1][F][1]==0,{a0,a1,a2,a3}]
// //       // {{a0->0,a1->1,a2->1,a3->-1}}
// //       new_taux = fSign * (t*(1 + t*(1 - tc)));
//    
//    if (role == 1)
//       new_taux = 1 - new_taux;
//    
//    return new_taux;
// }
// #endif // INCLUDE_OPTIONALS
// 
// 
// #ifdef INCLUDE_OPTIONALS
// // UPDATE:
// // - I think two things need to happen. First, if doing aligned grids (like hex-grid(n,0) or
// //   tri-subd{n,0}, adding a twist towards the center of the face would likely be a very
// //   good idea. Could be done inside here:
// //
// //   - Multiply point positions by ExpM(TwistAngle * A)
// //     with A being A = Ax*vCof[0] + Ay*vCof[1] + Az*vCof[2]
// //   - Angle is 0 at vertex and gets up to TwistAngle near center of face.
// //     increase should probably be some fraction of either linear angle
// //
// //           f = 1 - angle(vCur,vCof) / angle(vCur,v0)
// //   
// //     or some sort of power or area transform of that
// //     (e.g., ratio of small-circle circumferences with normal vCof at vCur vs v0?)
// //
// // - The barycentric transforms... they should be non-linear, but I should not just
// //   guess like atm. I think I should look into actual spherical geometry, and try
// //   to use that. Specify ``triangles'' in terms of great circles and so on, with the
// //   (v0,v1,v2) region being enclosed by them (if looked at from the top indivudally,
// //   they'd still all look as lines, and the symmetry ops would transform them
// //   accordingly).
// //
// //   Some guesses, part of which might not be too complicated to implement:
// //   - Specify barycentric coordinates in terms of solid angle weight of the
// //     opposite ``spherical triangle'' (instead of area of real triangle, as
// //     normal barycentric coordinates do). That one /may/ be quite complicated.
// //   - Specify barycentric coordinates in terms of linear angle fractions
// //     between v0/v1/v2 and vCur:
// //
// //           alpha0 = angle(vCur,v0)    (or angle(vCur,vCoe0) with coe0 = coe12?)
// //           alpha1 = angle(vCur,v1)
// //           alpha2 = angle(vCur,v2)
// //
// //     and then the barycentric weights somehow specifying the relative ratios between
// //     these angles... (or 1-barycentric weights specifying the ratios? Would
// //     need some thinking)
// //
// //     + ...with angle(u,v) = acos(dot(u,v)/(norm(u)*norm(v))) etc
// //     + ...and optimizing the Cartesian coordinates to match the fractions
// //       would likely be a non-linear process; but could start at
// //       Cartesian barycentric projected.
// #endif //INCLUDE_OPTIONALS
// 
// #ifdef INCLUDE_OPTIONALS
// void ApplyAreaDistortion(FPoint &vSeed, FScalar AreaDistortion, FTriangle const &eq_v012, FScalar SnapDistance)
// {
//    if (AreaDistortion != 0.) {
//       // convert to barycentric coordinates in terms of 012 triangle
// //       FPoint vBaryc_012 = eq_v012.MakeBarycentric(vSeed);
//       // compute area distortion at seed point
//       FScalar wx = 1/vSeed.squaredNorm();
//       // ^- approx linear dimension rescaling for sphere area mapping if points are equidistant in 2d barycentric space
//       // hm... what do I do with it? don't want to think this through to the end atm.
//       
// //       FPoint
// //          vc = ((eq_v012.v0 + eq_v012.v1 + eq_v012.v2)/FScalar(3));
// //       FScalar wc = 1/vc.squaredNorm();
// //       
// //       // wx = f*wc + (1-f)*1
// //       // -> f = (wx-1)/(wc-1)
// //       if (wc != 1) {
// //          FScalar
// //             f = (wx-1)/(wc-1);
// //          f = 0;
// //          if (0) {
// //             vSeed = f*vc + (1-f)*vSeed;
// //          } else {
// //             // that does the same thing. I should be scaling 1-(max(b)-min(b)) instead?
// //             FPoint
// //                vBaryc_c(1/FScalar(3), 1/FScalar(3), 1/FScalar(3));
// //             vSeed = eq_v012.MakeCartesian(f*vBaryc_c + (1-f)*vBaryc_012);
// //          }
// //          
// //             
// //             
// //          std::cout << fmt::format(" --- ApplyAreaDistortion:\n     | wc = {:.8f}  wx = {:.8f}  f = {:.8f}\n",
// //             double(wc),double(wx), double(f));
// //       }
// //       SnapDistance /= (eq_v012.v0 - eq_v012.v1).norm();
//       
//       FPoint
//          vc = ((eq_v012.v0 + eq_v012.v1 + eq_v012.v2)/FScalar(3));
//       FScalar wc = 1/vc.squaredNorm();
//       if (1) {
// //          wx = sqrt(wx);
// //          wc = sqrt(wc);
// //          FScalar fp = 0.5;
//          FScalar fp = 1.;
//          wx = pow(wx,fp);
//          wc = pow(wc,fp);
//       }
//       
//       // wx = f*wc + (1-f)*1
//       // -> f = (wx-1)/(wc-1)
//       FScalar
//          f = (wx-1)/(wc-1);
//       f = 1 - pow(1 - f, AreaDistortion);
// //       vSeed = f*vc + (1-f)*vSeed;
//          FPoint vBaryc_012 = eq_v012.MakeBarycentric(vSeed);
// //          FScalar r = vBaryc_012.maxCoeff() - vBaryc_012.minCoeff(); // if this is 0, the point is on the center. if it is 1, it is on a vertex.
// //          FVector3 vBaryc_New012 = vBaryc_012 / f;
// //          vBaryc_New012 -= ((vBaryc_New012.sum() - 1)/3)*FVector3(1,1,1);
//          FPoint
//             vBaryc_c(1/FScalar(3), 1/FScalar(3), 1/FScalar(3));
//          FVector3 vBaryc_New012 = f*vBaryc_c + (1-f)*vBaryc_012;
//          
//          // hack to move points close to edges back onto the edge so that they do not get duplicated.
//          for (size_t i = 0; i < 3; ++ i) {
// //             if (abs(vBaryc_012[i]) < 1e-6 && abs(vBaryc_New012[i]) < SnapDistance) {
// //             std::cout << fmt::format(" --- ApplyAreaDistortion:     | SNAP {:.8f} --> {:.8f} --> 0.0 (snap dist: {:.8f})\n",
// //                double(vBaryc_012[i]), double(vBaryc_New012[i]), double(SnapDistance));
// //                vBaryc_New012[i] = 0;
// //                vBaryc_New012 /= vBaryc_New012.sum();
// //             }
//             // snap edges
//             if (abs(vBaryc_012[i]) < 1e-6 && abs(vBaryc_New012[i]) < SnapDistance) {
//                vBaryc_New012[i] = 0;
//             }
//             // snap diagonals
//             int i1 = (i+1) % 3, i2 = (i+2) % 3;
//             if (abs(vBaryc_012[i1] - vBaryc_012[i2]) < 1e-6 && abs(vBaryc_New012[i1] - vBaryc_New012[i2]) < SnapDistance) {
//                FScalar c = (vBaryc_New012[i1] + vBaryc_New012[i2])/2;
//                vBaryc_New012[i1] = c;
//                vBaryc_New012[i2] = c;
//             }
//             vBaryc_New012 /= vBaryc_New012.sum();
//          }
//          
//          vSeed = eq_v012.MakeCartesian(vBaryc_New012);
//             
// //          std::cout << fmt::format(" --- ApplyAreaDistortion:\n     | wc = {:.8f}  wx = {:.8f}  f = {:.8f}  sum(vBarycNew) = {:.8f}\n",
// //             double(wc),double(wx), double(f), double(vBaryc_New012.sum()));
// 
//       
// //       FPoint
// //          ve01 = ((eq_v012.v0 + eq_v012.v1)/FScalar(2)),
// //          ve02 = ((eq_v012.v0 + eq_v012.v2)/FScalar(2)),
// //          ve12 = ((eq_v012.v1 + eq_v012.v2)/FScalar(2));
// //       FScalar
// //          we01 = 1/ve01.squaredNorm(),
// //          we02 = 1/ve02.squaredNorm(),
// //          we12 = 1/ve12.squaredNorm();
// //       FScalar
// //          f01 = (wx-1)/(we01-1),
// //          f02 = (wx-1)/(we02-1),
// //          f12 = (wx-1)/(we12-1);
// //       
// //       vSeed = (f01*(1-f02)*(1-f12))*ve01 + (f02*(1-f01)*(1-f12))*ve02 + (f12*(1-f01)*(1-f02))*ve12 + (1-((f01*(1-f02)*(1-f12)) + (f02*(1-f01)*(1-f12)) + (f12*(1-f01)*(1-f02))))*vSeed;
//    }
// //    throw std::runtime_error("absicht!/1");
// }
// #endif // INCLUDE_OPTIONALS
// 
// 
// #ifdef INCLUDE_OPTIONALS
// // void ApplyAreaDistortion(FPoint &vSeed, FScalar AreaDistortion, FTriangle const &eq_v012)
// // {
// //    if (AreaDistortion != 0.) {
// //       // convert to barycentric coordinates in terms of 012 triangle
// //       FPoint vBaryc_012 = eq_v012.MakeBarycentric(vSeed);
// //       // compute area distortion at seed point
// //       FScalar wx = 1/vSeed.squaredNorm();
// //       // ^- approx linear dimension rescaling for sphere area mapping if points are equidistant in 2d barycentric space
// //       // hm... what do I do with it? don't want to think this through to the end atm.
// //       
// //       FScalar
// //          p = wx * AreaDistortion;
// //       FVector3
// //          vBaryc_New012(0,0,0);
// // 
// // //          p = pow(p,-3/FScalar(2));
// //          p = 1/(p*p);
// // //          p = 1/p;
// //          vBaryc_New012[0] =  fAreaTraf1D(vBaryc_012[0],p,0);
// //          vBaryc_New012[1] =  fAreaTraf1D(vBaryc_012[1],p,0);
// //          vBaryc_New012[2] =  fAreaTraf1D(vBaryc_012[2],p,0);
// //          vBaryc_New012 /= vBaryc_New012.sum();
// //          
// //       
// //       std::cout << fmt::format(" --- ApplyAreaDistortion:\n     | vBaryc/old = {}\n     | vBaryc/new = {}\n     | delta = {}\n     | wx = {:.8f}\n",
// //          FmtVec3_DoublePrec(vBaryc_012,6), FmtVec3_DoublePrec(vBaryc_New012,6), FmtVec3_DoublePrec(vBaryc_New012 - vBaryc_012,6),
// //          double(wx));
// //       // convert back
// //       vSeed = eq_v012.MakeCartesian(vBaryc_New012);
// //    }
// // //    throw std::runtime_error("absicht!/1");
// // }
// 
// 
// // void ApplyAreaDistortion(FPoint &vSeed, FScalar AreaDistortion, FTriangle const &eq_v012)
// // {
// //    if (AreaDistortion != 0.) {
// //       // make a new triangle 01c from the equilateral region with c = (v0+v1+v2)/2;
// //       FTriangle eq_v01c{eq_v012.v0, eq_v012.v1, (1/FScalar(3))*(eq_v012.v0 + eq_v012.v1 + eq_v012.v2)};
// //       // convert to barycentric coordinates in terms of this triangle
// //       FPoint vBaryc_01c = eq_v01c.MakeBarycentric(vSeed);
// //       // increase/decrease weight relative to center point (normally one would want to increase number of points there)
// //       FScalar wx = sqrt((eq_v01c.v0.norm() + eq_v01c.v1.norm()) / (2*eq_v01c.v2.norm()));
// //       // ^- approx linear dimension rescaling for sphere area mapping if points are equidistant in 2d barycentric space
// //       // hm... what do I do with it? don't want to think this through to the end atm.
// //       
// // //       // F[t_]:=a0 + a1*t + a2*t^2+a3*t^3
// // //       // Solve[F[0]==0&&F[1]==1&&Derivative[1][F][0]==1&&Derivative[1][F][1]==0,{a0,a1,a2,a3}]
// // //       // {{a0->0,a1->1,a2->1,a3->-1}}
// // //       // -> so F[t] = t + t^2 - t^3 has F[0]=0, F[1]=1, and F[t] \approx t near t=0 and F[t] \approx 0 near t=1.
// // //       FScalar fSign = 1;
// // //       if (vBaryc_01c[2] < 0) fSign = -1;
// // //       FScalar
// // //          tc = fSign * vBaryc_01c[2],
// // //          tc1 = fSign * (tc * (1 + tc * (1 - tc)));
// // // //       FScalar
// // // //          d = wx*(tc1 - vBaryc_01c[2]); // actual difference in tauc.
// // //       FVector3
// // //          delta = tc1*wx*AreaDistortion*(FVector3(0, 0, 1) - vBaryc_01c);
// //       
// //       // shouldn't I want a derivative of wx at F[1]? Let's try that one.
// //       // Solve[F[0]==0&&F[1]==1&&Derivative[1][F][0]==1&&Derivative[1][F][1]==p,{a0,a1,a2,a3}]
// //       // {{a0->0,a1->1,a2->1-p,a3->-1+p}}
// //       // Plot[(t+(1-p)*t^2-(1-p)*t^3)/.p->2,{t,0,1}]
// //       // Plot[t*(1+(1-p)*t*(1-t))/.p->2,{t,0,1}]
// //       // --> formula look okay, I guess
// //       //
// //       // UPDATE:
// //       // - That is not really right, because I am also scaling the other points...
// //       // - hmpf... for T even the v0-v1 distortion is too large to ignore
// //       //   see:
// //       //   make && aigg '{ point-group: tetrahedral; degree: [2]; max-step: 1e-3; export: {file: /tmp/test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes}; initial-points: save{_initial_points.dat}; initial-points: hex-grid{8;0;01c;dist:1.0}; print: 0; max-it: 4096; dynamic:off; progress-check:0; damping:0.2 }'
// //       //   and combine with visualize-rule.sh /tmp/test_run.dat
// // 
// //       FScalar
// //          p = wx * AreaDistortion;
// // //       p *= p;
// // //       p = 1/(p*p);
// // // //       p = pow(p,3);
// //       FVector3
// //          vBaryc_New01c(0,0,0);
// // 
// //       if (1) {
// // //          p = pow(p,-3/FScalar(2));
// // //          p = 1/(p*p);
// //          p = 1/p;
// //          vBaryc_New01c[0] =  fAreaTraf1D(vBaryc_01c[0],p,1); // no... anything which does not treat v0/v1/v2 equal cannot possibly be right.
// //          vBaryc_New01c[1] =  fAreaTraf1D(vBaryc_01c[1],p,1);
// //          vBaryc_New01c[2] =  fAreaTraf1D(vBaryc_01c[2],p,0);
// //          vBaryc_New01c /= vBaryc_New01c.sum();
// //       } else {
// //          p = 1/(p*p); // <- just a guess. Maybe has something to do with the f * vBaryc_01c[0]; scaling below? But shouldn't I get a (2/3) or something power like that?
// //          vBaryc_New01c = FPoint(0, 0, fAreaTraf1D(vBaryc_01c[2], p, 0));
// //          // rescale t0 and t1 to be compatible with new tc (that should be equivalent to mapping
// //          // new(t0+t1) = (1 - F(1-(t0+t1))) to find the new sum of t0+t1)
// //          FScalar
// //             old_t01 = vBaryc_01c[0] + vBaryc_01c[1],
// //             new_t01 = 1 - vBaryc_New01c[2]; // sum of new_t0 + new_t1 + new_tc should be 1.
// //          if (old_t01 != 0) {
// //             FScalar f = new_t01/old_t01;
// //             vBaryc_New01c[0] = f * vBaryc_01c[0];
// //             vBaryc_New01c[1] = f * vBaryc_01c[1];
// //          }
// //       }
// //          
// //       
// //       std::cout << fmt::format(" --- ApplyAreaDistortion:\n     | vBaryc/old = {}\n     | vBaryc/new = {}\n     | delta = {}\n     | wx = {:.8f}  tc = {:.8f}  new_tc = {:.8f}\n",
// //          FmtVec3_DoublePrec(vBaryc_01c,6), FmtVec3_DoublePrec(vBaryc_New01c,6), FmtVec3_DoublePrec(vBaryc_New01c - vBaryc_01c,6),
// //          double(wx), double(vBaryc_01c[2]), double(vBaryc_New01c[2]));
// //       vBaryc_01c = vBaryc_New01c;
// //       // convert back
// //       vSeed = eq_v01c.MakeCartesian(vBaryc_01c);
// //    }
// // //    throw std::runtime_error("absicht!/1");
// // }
// #endif // INCLUDE_OPTIONALS
// 
// 
// #ifdef INCLUDE_OPTIONALS
// // void ApplyAreaDistortion(FPoint &vSeed, FScalar AreaDistortion, FTriangle const &eq_v012)
// // {
// //    if (AreaDistortion != 0.) {
// //       int
// //          N = 7;
// //       FScalar
// //          fWeightSum = 0;
// //       FVector3
// //          vBarycPtSum(0,0,0); // will accumulated weighted sum of barycentric coordinates of vBarycSeed;
// //       for (int iTau01 = 0; iTau01 <= N; ++iTau01) {
// //          int
// //             iTau2 = N - iTau01;
// //          for (int iTau0 = 0; iTau0 <= iTau01; ++iTau0) {
// //             int
// //                iTau1 = iTau01 - iTau0;
// //             FScalar
// //                f = 1/FScalar(N);
// //             FPoint
// //                vBarycGridPt(f*iTau0, f*iTau1, f*iTau2);
// //             
// //             FVector3
// //                vCartGridPt = eq_v012.MakeCartesian(vBarycGridPt);
// //             FScalar
// //                fNormRsq = 1/vCartGridPt.squaredNorm();
// //          }
// //       }
// //       // hm, nope.
// //       
// //    }
// // }
// #endif // INCLUDE_OPTIONALS
#endif // INCLUDE_ABANDONED

void FGridSearchOptions::MakeInitialBarycentricGrid(std::string const &Type, ptrdiff_t s, ptrdiff_t t, FScalar OffsetM, FScalar OffsetN, FScalar AreaDistortion, FScalar AreaDistortionSnap, std::string const &Order, std::string const &FullDeclForExport_)
{
   std::cout << fmt::format(" Initial seed points -- barycentric coordinates (s={}, t={})",s,t) << std::endl;
   std::cout << "\n   SEED         TAU0             TAU1           TAUC" << std::endl;
   FMappedPointCloud
      PointCloud(pPointGroup->GetOps());
   FPointList
      PointList,
      BarycPoints;
   FTriangle
      // triangle in terms of which the barycentric coordinates are specified.
      BaseRegion,
      // triangle which tiles the polyhedron surface
      FundamentalRegion = pPointGroup->GetTriangle(FPointGroup::REGION_Fundamental);
   if (Type == "hex-grid") {
      // this one is set up to make coordinates given in terms of 01c barycentrics
      BaseRegion = pPointGroup->GetTriangle(FPointGroup::REGION_01c);
      FTestInsideFundamentalRegionPred
         TestInsideFn(BaseRegion, FundamentalRegion, BarycentricGeomType),
         *pTestInsideFn = &TestInsideFn;
      if (pPointGroup->IsFundamentialRegionIncludedIn01c())
         // test only needed if fundamental region does not lie inside 01c.
         // otherwise, testing with the 01c coordinates the function builds
         // internally is sufficent.
         pTestInsideFn = 0;
      MakeHexGridPointList(BarycPoints, Type, s, t, OffsetM, OffsetN, Order, pTestInsideFn);
   } else if (Type == "tri-subd") {
      // this one is set up to make coordinates given in terms of 012 barycentrics
      BaseRegion = pPointGroup->GetTriangle(FPointGroup::REGION_012);
      MakeTriangleSubdivGrid(BarycPoints, Type, s, t, OffsetM, OffsetN, Order);
#ifdef INCLUDE_ABANDONED
//       FTestInsideFundamentalRegionPred
//          TestInsideFn(BaseRegion, FundamentalRegion);
//       MakeTriangleSubdivGrid(BarycPoints, Type, s, t, OffsetM, OffsetN, Order, TestInsideFn);
#endif // INCLUDE_ABANDONED
   } else {
      throw std::runtime_error(fmt::format("FGridSearchOptions::MakeInitialBarycentricGrid: grid type '{}' not recognized.", Type));
   }
   std::stringstream
      ExportData; // gets to accumulate text description of initial guess in case FileName_BarycentricGuess was set.
   if (!FileName_BarycentricGuess.empty()) {
      ExportData << fmt::format("# setup of initial points via barycentric guess: ");
      ExportData << fmt::format("#! sym = {}\n", pPointGroup->Name(NAMETYPE_Symmetry));
      ExportData << fmt::format("#! def = {}\n", FullDeclForExport_);
      ExportData << "# Notes:\n";
      ExportData << "# - barycentric coords are defined in terms of the base-region.\n";
      ExportData << "# - fundamental-region denotes a triangle tiling the polyhedron suface. n";
      ExportData << "# - coordinates given as [v0, v1, v2] with v... = [x,y,z]\n";
      ExportData << fmt::format("#! base-region = {}\n", FmtTriangle_DoublePrec(BaseRegion));
      ExportData << fmt::format("#! fund-region = {}\n", FmtTriangle_DoublePrec(FundamentalRegion));
      ExportData << fmt::format("#! v012-region = {}\n", FmtTriangle_DoublePrec(pPointGroup->GetTriangle(FPointGroup::REGION_012)));
      ExportData << "# Possible rejection reasons:\n";
      ExportData << "# - 0x01 -> point lies outside fundamental region\n";
      ExportData << "# - 0x02 -> point is symmetry-equivalent to a previous point\n";
      ExportData << "#";
      ExportData << fmt::format("# {:^4}  {:^8}    {:^15}  {:^15}  {:^15}       {:^15}  {:^15}  {:^15}\n",
         "SEED", "REJECT?", "TAU0-BASE", "TAU1-BASE", "TAU2-BASE", "X-POS", "Y-POS", "Z-POS");
//          "SEED", "INCLUDE", "TAU0-BASE", "TAU1-BASE", "TAU1-BASE", "X-POS", "Y-POS", "Z-POS");
   }
   for (FPointList::const_iterator itBarycPt = BarycPoints.begin(); itBarycPt != BarycPoints.end(); ++itBarycPt) {
      // vSeed = tau0*P0 + tau1*P1 + tauc*Pc
      FPoint
         vBarycPt = *itBarycPt;
      // add distortion...
      // WARNING: should be done with v012 barycentric coordinates, not with v01c!
#ifdef INCLUDE_ABANDONED
//       vBarycPt /= vBarycPt[0] + vBarycPt[1] + vBarycPt[2];
//       FPoint
//          vSeed = pPointGroup->ConvertBarycentricToCartesian(vBarycPt);
//          GeomType = FTriangle::GEOMTYPE_Planar;         
#endif // INCLUDE_ABANDONED
      FPoint
         vSeed = BaseRegion.MakeCartesian(vBarycPt, BarycentricGeomType);
      bool
         // this may happen when using spherical barycentric coordinates
         // and supplying points which cannot be represented.
         CoordConversionFailed = false;
      if (BarycentricGeomType == FTriangle::GEOMTYPE_Spherical) {
         FPoint bt = BaseRegion.MakeBarycentric(vSeed, BarycentricGeomType);
         FPoint btc = BaseRegion.MakeCartesian(bt, BarycentricGeomType);
//          std::cout << fmt::format("-- CONVERT: BASE REGION\n  vBaryc = {} --> vCart = {}\n  vBaryc = {} --> vCart = {} (after back conversion)\n  diff norms:  baryc: {:8.2e}   cart: {:8.2e}", FmtVec3_DoublePrec(vBarycPt,6), FmtVec3_DoublePrec(vSeed,6), FmtVec3_DoublePrec(bt,6), FmtVec3_DoublePrec(btc,6), double((vBarycPt-bt).norm()), double((vSeed-btc).norm())) << std::endl;
         if (!IsAlmostZero((btc - vSeed).norm()))
            CoordConversionFailed = true;
      }
#ifdef INCLUDE_ABANDONED      
// #ifdef INCLUDE_OPTIONALS
//       if (AreaDistortion != 0.) {
// //          FScalar SnapDistance = 0.5*sqrt(4*X_PI/BarycPoints.size()); // <- FIXME: we do not actually know how many points there should be at the end!
//          FScalar SnapDistance = AreaDistortionSnap*sqrt(FScalar(1)/BarycPoints.size()); // <- FIXME: we do not actually know how many points there should be at the end!
// //          FScalar SnapDistance = 100.; // snap everything to type of orbit of source triangle
// //          FScalar SnapDistance = 0.;
//          ApplyAreaDistortion(vSeed, AreaDistortion, pPointGroup->GetTriangle(FPointGroup::REGION_012), SnapDistance);
//       }
// #endif // INCLUDE_OPTIONALS
#endif // INCLUDE_ABANDONED
      // check if we had the seed point before.
      ptrdiff_t
         iSeed = ptrdiff_t(PointList.size());
      bool
         // is the point inside the fundamental region?
         Outside = CoordConversionFailed || !FundamentalRegion.IsInside_Cartesian(vSeed, eps, BarycentricGeomType);
      bool
         // did we already add a symmetry equivalent of this point?
         AlreadyThere = false;
      if (!Outside) {
         FMappedPointCloud::FPointEntry const
            &e = PointCloud.Insert(vSeed, iSeed);
   //       std::cout << fmt::format("{} {:>4}  {:15.8f}  {:15.8f}  {:15.8f}", (e.iPoint == iSeed)?".":"X", PointList.size(), double(vBarycPt[0]), double(vBarycPt[1]), double(vBarycPt[2])) << std::endl;
         AlreadyThere = (e.iPoint != iSeed);
      }

      if (!FileName_BarycentricGuess.empty()) {
         unsigned
            Flags = 0x01*unsigned(Outside) + 0x02*unsigned(AlreadyThere);
         ptrdiff_t
            iSeedId = ptrdiff_t(PointList.size());
         if (AlreadyThere || Outside)
            iSeedId = -1; // this one will not be included.
#ifdef INCLUDE_ABANDONED
//          char const
//             *pIncludeDesc = 0;
//          if (Flags == 0)
//             pIncludeDesc = "yes";
//          else if (Flags == 1)
//             pIncludeDesc = "no(out)";
//          else if (Flags == 2)
//             pIncludeDesc = "no(sym)";
//          else if (Flags == 3)
//             pIncludeDesc = "no(s,o)";
#endif // INCLUDE_ABANDONED
         vBarycPt = FundamentalRegion.MakeBarycentric(vSeed, BarycentricGeomType);// FIXME: remove this (debug stuff)
         ExportData << fmt::format("  {:>4}  {:8}    {:15.8f}  {:15.8f}  {:15.8f}       {:15.8f}  {:15.8f}  {:15.8f}\n",
//             iSeedId, pIncludeDesc,
            iSeedId, Flags,
            double(vBarycPt[0]), double(vBarycPt[1]), double(vBarycPt[2]),
            double(vSeed[0]), double(vSeed[1]), double(vSeed[2]));
      }
      if (AlreadyThere || Outside)
         // either outside of fundamental region (in this case this is not meant
         // to be added) or an image of this is already in the list (in this case,
         // don't add it again)
         continue;
      PointList.push_back(vSeed);
      std::cout << fmt::format("  {:>4}  {:15.8f}  {:15.8f}  {:15.8f}", PointList.size(), double(vBarycPt[0]), double(vBarycPt[1]), double(vBarycPt[2])) << std::endl;
      
      // should also probably check other variants of the points being there before already...
   }
   if (!FileName_BarycentricGuess.empty()) {
      std::cout << fmt::format("* Exporting initial guess data to '{}'.", FileName_BarycentricGuess) << std::endl;
      std::ofstream
         ExportFile(FileName_BarycentricGuess.c_str());
      ExportFile << ExportData.str() << std::endl;
   }
//    throw std::runtime_error("absicht!/3");
#ifdef INCLUDE_ABANDONED
//    for (ptrdiff_t m_ = -100; m_ <= 100; ++ m_) {
//       for (ptrdiff_t n_ = -100; n_ <= 100; ++ n_) {
//          FScalar
//             m = FScalar(m_) + OffsetM,
//             n = FScalar(n_) + OffsetN;
//          // compute barycentric coordinates with respect to P0, P1, Pc (see pg. 14)
//          FScalar
//             eps = 1e-8,
//             tau0,
//             tau1,
//             tauc;
//          if (Type == "hex-grid") {
//             tau0 = (s*s + (s-t)*m + s*t + t*t - (s+2*t)*n)/FScalar(s*s + s*t + t*t),
//             tau1 = ((2*m + n)*s + (m-n)*t)/FScalar(s*s + s*t + t*t),
//             tauc = 3*(n*t - m*s)/FScalar(s*s + s*t + t*t);
//          } else {
//             throw std::runtime_error(fmt::format("FGridSearchOptions::MakeInitialBarycentricGrid: grid type '{}' not recognized.", Type));
//          }
// 
//          if (Order == "01c") {
//          } else if (Order == "10c") {
//             std::swap(tau0, tau1);
//          } else if (Order == "0c1") {
//             std::swap(tau1, tauc);
//          } else if (Order == "c10") {
//             std::swap(tau0, tauc);
//          } else if (Order == "1c0") {
//             std::swap(tauc, tau1);
//             std::swap(tau1, tau0);
//          } else if (Order == "c01") {
//             std::swap(tauc, tau0);
//             std::swap(tau0, tau1);
//          } else {
//             throw std::runtime_error(fmt::format("barycentric coordinate permutation '{}' not recognized. Should be permutation of '0', '1', and 'c'", Order));
//          }
// 
//          bool inside = tau0 >= -eps && tau1 >= -eps && tauc >= -eps;
//          if (!inside)
//             continue;
// //                if (std::abs(tau1-1.) < eps)
// //                   // we already have P0 -- don't need P1, too.
// //                   inside = false;
// //                if (std::abs(tau0) < eps && std::abs(tauc-1.) >= eps)
// //                   // only one of the edges P0--PC and P1--PC is unique. Skip points lying on P1--PC.
// //                   // (unless we have PC itself (tauc=1), which should stay in.
// //                   inside = false;
//          // vSeed = tau0*P0 + tau1*P1 + tauc*Pc
//          FPoint
//             vSeed = pPointGroup->ConvertBarycentricToCartesian(FPoint(tau0, tau1, tauc));
// 
//          // check if we had the seed point before.
//          ptrdiff_t
//             iSeed = ptrdiff_t(PointList.size());
//          FMappedPointCloud::FPointEntry const
//             &e = PointCloud.Insert(vSeed, iSeed);
//          if (e.iPoint != iSeed)
//             // there was a symmetry equivalent point already. Don't add it again.
//             continue;
//          assert(inside);
//          PointList.push_back(vSeed);
//          std::cout << fmt::format("  {:>4}  {:15.8f}  {:15.8f}  {:15.8f}", PointList.size(), double(tau0), double(tau1), double(tauc)) << std::endl;
//          // should also probably check other variants of the points being there before already...
//       }
//    }
#endif // INCLUDE_ABANDONED
   InitialPoints = PointList;
}


void FGridSearchOptions::ReadMetaInfoFromAiggGridFileHeader(std::string FileName_)
{
   std::ifstream
      inp(FileName_.c_str());

   // process header line to extract core data (point-group, target degree, and
   // possibly number of points)
   {
      std::string
         HeaderLine;
      std::getline(inp, HeaderLine); // read first line.
      if (inp.bad() || inp.fail() || !inp.good())
         throw std::runtime_error(fmt::format("failed to read from input file '{}'", FileName_));

      // example for header line format:
      // # icosahedral integration rule: npts = 312, lmax = 29, wtsp = 1.628, res = 7.61e-140
      ct::string_slice
         slHeaderLine(HeaderLine);
      ct::split_result
         srDefArgs;
      slHeaderLine.split(srDefArgs, ':');
      if (srDefArgs.size() != 2)
         throw std::runtime_error(fmt::format("First line of file '{}' does not conform to expected aigg grid file format. No ':' in # ... rule: k0 = v0, k1 = v1, ...'. Actual header line is:\n'{}'", FileName_, HeaderLine));
      
      // isolate point group name.
      ct::split_result
         srDef;
      srDefArgs[0].split(srDef, ' ');
      if (srDef.size() != 4)
         throw std::runtime_error(fmt::format("First line of file '{}' does not conform to expected aigg grid file format. Expected '# <point-group> integration rule: ...'. Actual header line is:\n'{}'", FileName_, HeaderLine));
      std::string
         sPointGroup = srDef[1].to_str();
      this->pPointGroup = MakePointGroup(sPointGroup);
      
      // isolate key=value pairs on args line.
      ct::split_result
         srArgs;
      srDefArgs[1].split(srArgs, ',');
      for (size_t iArg = 0; iArg != srArgs.size(); ++ iArg) {
         ct::string_slice
            slArg = srArgs[iArg];
         ct::split_result
            srKeyVal;
         slArg.split(srKeyVal, '=');
         if (srKeyVal.size() != 2)
            throw std::runtime_error(fmt::format("First line of file '{}' does not conform to expected aigg grid file format. Error in '...rule: k0 = v0, k1 = v1, ... pair enumeration. Actual header line is:\n'{}'", FileName_, HeaderLine));
         if (srKeyVal[0] == "lmax") {
            lStart = srKeyVal[1].to_int();
            // note: this leaves lStep and lList alone. This would allow the
            // schedule to be specified before the initial grid data is read, and
            // here would just overwrite the lStart and leave the rest of the
            // schedule alone.
         }
            
      }
   }

   // process further lines in order to replicate the original grid description
   // (apart from the header line) as comments to the refined grid's output desc.
   {
      std::string
         Line;
      std::string
         sCommentType = "# ",
         sGridQuote = "| ";
      FExportOptions::FCommentList
         InheritanceComments;
      for (;;) {
         std::getline(inp, Line);
         if (inp.fail())
            // failed to read the line---e.g., due to encountering an eof.
            break;
         ct::string_slice
            slLine(Line);
         if (!slLine.try_skip_literal(sCommentType.c_str()))
            // stop once we found the first line which is not part of the aigg-
            // generated grid description at the start of the file.
            break;
         
         if (slLine.try_skip_literal(sGridQuote.c_str())) {
            // input grid was already a refinement of something else. In this
            // case, do not quote the refinement description, but only the
            // (innermost) original grid description
            InheritanceComments.clear();
            sCommentType += sGridQuote;
         }
         
         InheritanceComments.push_back(sGridQuote + slLine.to_str());
      }
      if (!InheritanceComments.empty()) {
         m_ExtraComments.push_back("comments from original input grid for refinement:");
         m_ExtraComments.splice(m_ExtraComments.end(), InheritanceComments);
      }
   }   
}


void FGridSearchOptions::ProcessStartingPoints(std::string const &StartingPointsDecl_)
{
   ct::FPropertyListStr
      Props(ct::string_slice(StartingPointsDecl_), 0, ct::PROPLIST_AllowEntriesWithoutName);
   if (Props.Name == "read-barycentric" || Props.Name == "read-cartesian" || Props.Name == "read-aigg" || Props.Name == "read-aigg-grid-file" || Props.Name == "read-grid-file") {
      if (Props.size() != 1 || !(Props[0].Name == "" || Props[0].Name == "file"))
         throw std::runtime_error(fmt::format("initial points: expected {}{{<filename>}}", Props.Name.to_str()));
      InitialPoints = ReadPointListFromFile(Props[0].Content);
      if (0) {
         for (size_t iPt = 0; iPt != InitialPoints.size(); ++ iPt)
            std::cout << fmt::format("   rd: {:8} -> {:12.6f}  {:12.6f}  {:12.6f}\n", iPt+1, double(InitialPoints[iPt][0]), double(InitialPoints[iPt][1]), double(InitialPoints[iPt][2]));
      }
      if (Props.Name == "read-barycentric") {
         FTriangle v01c = pPointGroup->GetTriangle(FPointGroup::REGION_01c);
         for (size_t iPoint = 0; iPoint != InitialPoints.size(); ++ iPoint) {
//             InitialPoints[iPoint] = pPointGroup->ConvertBarycentricToCartesian(InitialPoints[iPoint]);
            InitialPoints[iPoint] = v01c.MakeCartesian(InitialPoints[iPoint], BarycentricGeomType);
         }
      }
      if (Props.Name == "read-aigg" || Props.Name == "read-aigg-grid-file" || Props.Name == "read-grid-file") {
         ReadMetaInfoFromAiggGridFileHeader(Props[0].Content);
      }
   } else if (Props.Name == "hex-grid" || Props.Name == "rect-grid" || Props.Name == "tri-subd") {
      std::string
         order = "01c";
      ptrdiff_t
         s = 1, t = 1;
      FScalar
         OffsetM = 0, OffsetN = 0, AreaDistortion = 0;
      FScalar
         AreaDistortionSnap = 0.25;
      size_t
         iUnnamedArg = 0;
      for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
         if (it->Name != "")
            // once there was a named argument, do no longer allow unnamed ones to follow
            iUnnamedArg = size_t(-1);
//          std::cout << fmt::format("-- ProcessStartingPoints: '{}' -> {}\n", it->Name.to_str(), it->Content.to_str());
         if (it->Name == "s" || iUnnamedArg == 0) {
            s = it->Content.to_int();
         } else if (it->Name == "t" || iUnnamedArg == 1) {
            t = it->Content.to_int();
         } else if (it->Name == "order" || iUnnamedArg == 2) {
            order = it->Content.to_str();
            if (!(order == "01c" || order == "0c1" || order == "10c" || order == "c01" || order == "1c0" || order == "c10"))
               throw std::runtime_error("barycentric grid generation error: coordinate order must be a permutation of of '01c' (for vertex 0, 1, C, respectively)");
         } else if (it->Name == "offset-m" || iUnnamedArg == 3) {
            OffsetM = FScalar(it->Content.to_float());
         } else if (it->Name == "offset-n" || iUnnamedArg == 4) {
            OffsetN = FScalar(it->Content.to_float());
         } else if (it->Name == "dist") {
            AreaDistortion = FScalar(it->Content.to_float());
         } else if (it->Name == "snap") {
            AreaDistortionSnap = FScalar(it->Content.to_float());
         } else if (it->Name == "save" || it->Name == "export") {
            FileName_BarycentricGuess = it->Content.to_str();
         } else {
            if (iUnnamedArg > 4)
               throw std::runtime_error(fmt::format("too many ({}) unnamed hex-grid generator properties.", iUnnamedArg));
            else
               throw std::runtime_error(fmt::format("hex-grid generator property '{}' not recognized", it->Name.to_str()));
         }
         if (it->Name == "" && iUnnamedArg != size_t(-1))
            iUnnamedArg += 1;
      }
      MakeInitialBarycentricGrid(Props.Name.to_str(), s, t, OffsetM, OffsetN, AreaDistortion, AreaDistortionSnap, order, StartingPointsDecl_);
   } else if (Props.Name == "save" || Props.Name == "export") {
      // needs extra "initial-points:save{file.dat}" before the actual 'initial-points:'
      // option which specifies the actual initial guess properties. It's a bit ugly,
      // but good enough for now.
      FileName_BarycentricGuess = Props[0].Content.to_str();
   } else {
      throw std::runtime_error(fmt::format("initial point declaration '{}' not recognized. Should be "
         "'read-barycentric', 'tri-subd', 'read-cartesian', 'hex-grid', or 'rect-grid'", Props.Name.to_str()));
   }
}


