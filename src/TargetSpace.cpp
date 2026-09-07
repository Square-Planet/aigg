// make && OMP_NUM_THREADS=12 aigg_256 --preprocess '{point-group:[T,O,I]; order:1 to 50; skip:no}'
// make && OMP_NUM_THREADS=12 aigg_256 --preprocess '{point-group:all; order:1 to 50; skip:no}'
#ifdef INCLUDE_ABANDONED
// make && OMP_NUM_THREADS=8 aigg_128 --preprocess '{point-group:I; order:1 to 40; skip:no}' > /tmp/prep-trinom.log && cat /tmp/prep-trinom.log  | grep unit
// make && aigg_64_d '{point-group:I; degree:[29]; initial-points: hex-grid{5;1;01c}; print:0; target-space:exact; export:last.dat}'
// make && aigg_256 --preprocess '{point-group:all; order:1 to 20}'
// make && aigg_64_d --preprocess '{point-group:I; order:1 to 20}'
// make && aigg_64_d --preprocess '{point-group:I; order:[10]}'
// make && aigg_256 --preprocess '{point-group:[T,Td,Th,O,Oh,I,Ih]; order:1 to 20}'
//
//
// Some ramblings:
//
// Which data to generate/store?
// - not sure if I really want to make the full non-redundant set for *all* the groups.
//   Maybe just T/O/I by themselves are fine---not that a target function set for any
//   group G will automatically be sufficient for any group including G as a subgroup.
//   (it may just contain some functions more than strictly needed)
//
// How to represents target function spaces?
// - it could be sensible to make an abstract class for the target functions,
//   and refer to them via pointer only. That could unify raw monomials, group-
//   averaged monomials, group-averaged orthonormal polynomials, and possibly
//   other functions.
//
// - The reason is that, ultimately, we do not really need much function-type
//   specific work. We only need stuff to make:
//   + Overlap matrix
//   + rhs vector 't'
//   + K matrix
//   + K-matrix derivatives inline contracted to w-type or e-type quantities.
//   But most of the functionality should be equal for different types of functions.
#endif // INCLUDE_ABANDONED
#include "Aigg.h"

#include <ctime>
#include <fstream>
#include <sstream>
#include "CxParse1.h"
#include "CxIo.h"
#include "CxTiming.h"
#include "CxOpenMpProxy.h"

#include "TargetSpace.h"
#include "TargetFunctions.h"
#include "SearchOptions.h"
#include "MathOps.h"

#include "LocalizeVectors.h"


static std::string const
   s_AllKnownPointGroups = "T Td Th O Oh I Ih";
static bool s_TestHalfMxm = false;

static bool const
   // if set, remove `m_Factor` and `m_Degeneracy` from `FMonomial`
   // instances before computing and exporting data from `--preprocess`
   // (this will turn the monomials into raw, unnormalized Cartesian
   // monomials) and similarly remove them from data loaded for function
   // expansions.
   s_CancelPrimitiveRenormFactors = false;
//    s_CancelPrimitiveRenormFactors = true; // FIXME: put this back. It's for the localization test.

FScalar ComputeOverlapRaw(FPolynomialN const &poly_i, FPolynomialN const &poly_j);


std::string FormatTabulatedFnListError(std::string const &Task, FPointGroup const *pPointGroup, int Order, std::string const &FileName, std::string const &ProblemDesc)
{
   std::stringstream ss;
   ss << fmt::format("data-error{{task: {}; symmetry: {}; order: {}; fp: {}}}",
         Task, pPointGroup->Name(NAMETYPE_Symbol), Order, g_ScalarFloatSizeBits);
#ifdef INCLUDE_ABANDONED
//    ss << ": failed to load precomputed target function list";
#endif // INCLUDE_ABANDONED
   ss << ": failed to load precomputed target space basis functions";
   if (!FileName.empty())
      ss << fmt::format(" from '{}'", FileName);
   if (!ProblemDesc.empty())
      ss << fmt::format(": {}", ProblemDesc);
   return ss.str();
}

ETabulatedFnListError::ETabulatedFnListError(std::string const &Task, FPointGroup const *pPointGroup, int Order, std::string const &FileName, std::string const &ProblemDesc)
   : std::runtime_error(FormatTabulatedFnListError(Task, pPointGroup, Order, FileName, ProblemDesc))
{}

ETabulatedFnListError::ETabulatedFnListError(std::string const &ProblemDesc)
   : std::runtime_error(ProblemDesc)
{}



FTargetFnListOptions::FTargetFnListOptions(FPrintLevel const &PrintLevel_, FAccuracyPreference const &AccuracyPreference_, bool AllowNeglectResidualDeriv_)
   : PrintLevel(PrintLevel_), AccuracyPreference(AccuracyPreference_), AllowNeglectResidualDeriv(AllowNeglectResidualDeriv_)
{}


void FTargetSpacePreprocessOptions::SetAccuracyPreference(FAccuracyPreference const &AccuracyPreference_)
{
   this->AccuracyPreference = AccuracyPreference_;
   // … if anything needs more specific setup regarding detailed approximation/check/whatever options, that would go here.
}


FTargetSpacePreprocessOptions::FTargetSpacePreprocessOptions(std::string const &Desc_)
{
//    PointGroups.push_back(new FIcosahedralGroup());
   StoreVersionInfo = true;
   StoreIdentityInfo = true;
   Export = true;
   SkipExisting = SKIP_BetterOrEqual;
   iPrintLevel = 0;
   SetAccuracyPreference(ACCURACY_Default);

   ct::FPropertyListStr
      Props(Desc_);
   for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
      if (option_names::PointGroupQ(it->Name, option_names::AllowPlural)) {
         PointGroups.clear();
         if (it->Content[0] == '[' || it->Content == "all") {
            // a list of groups given?
            ct::long_split_result sr;
            if (it->Content == "all") {
               ct::string_slice(s_AllKnownPointGroups).split(sr, ' ');
            } else {
               it->Content.split_list(sr);
            }
            for (size_t il = 0; il != sr.size(); ++ il)
               PointGroups.push_back(MakePointGroup(sr[il].to_str()));
         } else {
            // not a list? Assume a single group.
            PointGroups.push_back(MakePointGroup(it->Content.to_str()));
         }
      } else if (option_names::RuleOrderQ(it->Name, option_names::AllowPlural)) {
         // there are a bunch of versions of this:
         // - a single integer (in this case we will process only this particular `lmax`)
         // - a list of either integers, or `[⟨lstart⟩,step:+1,⟨l1⟩,⟨l2⟩,…]`
         if (it->Content.empty())
            throw std::runtime_error(fmt::format("error in preprocess declaration '{}'. {}, if given, cannot be empty.", Desc_, it->Name.to_str()));
         if (it->Content[0] == '[') {
            // starts with '[' -> list of stuff case. See what we got.
            assert(it->Content[0] == '[');
            lList.clear();
            ct::long_split_result
               sr;
            it->Content.split_list(sr);
            for (size_t il = 0; il != sr.size(); ++ il)
               lList.push_back(sr[il].to_int());
         } else if (it->Content.find(" to ") != it->Content.last) {
            ct::long_split_result
               sr;
            it->Content.split(sr, ' ');
            if (sr.size() != 3 || sr[1] != "to")
               throw std::runtime_error(fmt::format("order declaration '{}': expected '<start> to <end>' but found '{}'.", it->Name.to_str(), it->Content.to_str()));
            lList.clear();
            int lStart = sr[0].to_int(), lEnd = sr[2].to_int();
            for (int l = lStart; l <= lEnd; ++ l)
               lList.push_back(l);
         } else {
            // lmax is not a list or a "⟨a⟩ to ⟨b⟩" type declaration. Assume single integer case.
            lList.clear();
            lList.push_back(it->Content.to_int());
         }
      } else if (option_names::AccuracyPreferenceQ(it->Name)) {
         SetAccuracyPreference(ParseAccuracyPreference(it->Name, it->Content));
      } else if (it->Name == "print" || it->Name == "print-level") {
         ptrdiff_t iSinglePrintLevel;
         if (it->Content.try_convert_to_int(&iSinglePrintLevel)) {
            iPrintLevel = iSinglePrintLevel;
            PrintOptions.SetDefaults(iPrintLevel);
         } else {
            PrintOptions.SetArgs(it->Content);
            iPrintLevel = int(PrintOptions[PRINT_General]);
         }
      } else if (it->Name == "skip" || it->Name == "skip-existing") {
         if (it->Content == "all" || it->Content == "any" )
            SkipExisting = SKIP_All;
         else if (it->Content == "better" || it->Content == "gt")
            SkipExisting = SKIP_Better;
         else if (it->Content == "better-or-equal" || it->Content == "yes" || it->Content == "on" || it->Content == "ge")
            SkipExisting = SKIP_BetterOrEqual;
         else if (it->Content == "no" || it->Content == "none" || it->Content == "off")
            SkipExisting = SKIP_None;
         else
            throw std::runtime_error(fmt::format("'{}': declaration '{}' not recognized. Should be one of: 'all', 'better', 'better-or-equal' (or 'yes'), 'none' (or 'no').", it->Name.to_str()));
      } else if (it->Name == "export" || it->Name == "write") {
         Export = it->Content.to_bool();
      } else {
         throw std::runtime_error(fmt::format("target space preprocess property '{}' not recognized.", it->Name.to_str()));
      }
   }
}


enum FDefaultFileNameOptions {
   FILENAME_OnlyBaseName = 0x001,
};
std::string MakeTargetSpaceDefaultFileName(FPointGroup const &PointGroup, int l, unsigned Flags = 0)
{
   std::string
      // subspace-
      BaseName(fmt::format("basis-S2-{}-l{}.dat", PointGroup.Name(NAMETYPE_Symbol), l));
   if (bool(Flags & FILENAME_OnlyBaseName)) {
      return BaseName;
   } else {
      ct::FFileLocator
         FileLocator;
      std::string
         FullName(FileLocator.JoinPath("data", BaseName));
      return FullName;
   }
}



typedef std::map<FMonomialN, FScalar, FMonomialN_ExportOrderPred>
   FPolynomialTermsExport;

// - "poly": polynomial in x,y,z in terms of explicit x,y,z coefficients
// - "orth": functions are orthonormalized
// - "mean": these functions already are group-averaged (so there is no need for another
//           sum over the group orbits in the computation of K and its derivatives)
static char const *s_DefaultFunctionListType = "[poly,orth,mean]";
//    char const *pFunctionListType = "[poly,orth,group-averaged]";
//    char const *pFunctionListType = "[poly,orth,group-avg]";
//    char const *pFunctionListType = "[poly,orth,group-mean]";
static char const *s_DefaultTaskStr = "S2"; // S2 = 2-dimensional unit sphere (in R^3 space)

void WriteTargetFnSpaceDataFile(FPointGroup const &PointGroup, int l, FIoPolynomialList const &TargetFns, FTargetSpacePreprocessOptions const &Options)
{
   std::string
      FileName = MakeTargetSpaceDefaultFileName(PointGroup, l);
   std::ofstream
      out(FileName.c_str(), std::ofstream::trunc | std::ofstream::out);
   char const
      *pFunctionListType = s_DefaultFunctionListType;
   bool
      // controls whether the output basis will first be normalized before their
      // coefficients will be written.
      // For example, this controls whether coefficients apply to raw monomials
      // or normalized monomials: In the latter case, a term
      //
      //     t_{n} = … + c^{n}_{ijk} x^i y^j z^k + …
      //
      // in the logical expansion will be replaced by
      //
      //     t_{n} = … + (c^{n}_{ijk}/f_{norm}) f_{norm} x^i y^j z^k + …
      //
      // with `f_{norm} = 1/sqrt(NormSq{x^i t^j z^k})` the normalization factor,
      // and therefore the coefficient will be written as `(c^{n}_{ijk}/f_{norm})`.
      //
      // This construction is technically redundant, but might make the emitted
      // data files less confusing. E.g., without the normalization, basis-S2-I-l39.dat
      // appears to contain a basis function with coefficients ranging from 1e-15 to 1e+15,
      // and whether this actually means anything is not obvious at first glance.
      bNormalizeBasisFn = true;
   auto CoeffForExport = [&](FScalar const &CoeffAsStored, FMonomialN const &m) {
      if (!bNormalizeBasisFn)
         return CoeffAsStored;
      else
         // `c * m` should stay the same. so if we normalize `m` by replacing it by
         // `m_norm = m/‖m‖`, `c` gets multiplied by `‖m‖` such that
         // `c_norm * m_norm = c * m` remains intact.
         return CoeffAsStored * CalcMonomialUnitSphereNorm(m);
   };
   if (1) {
      out << "# generated ";
      if (Options.StoreIdentityInfo) {
         out << MakeCurrentTimeString();
      }
      if (Options.StoreVersionInfo) {
         if (Options.StoreIdentityInfo)
            out << " ";
         out << "by " << MakeAiggVersionString();
      }
      out << "\n";
   }
   out << fmt::format("target-fn-space{{domain: {}; point-group: {}; degree: {}; fp: {}}} = [\n",
            s_DefaultTaskStr, PointGroup.Name(NAMETYPE_Symbol), l, g_ScalarFloatSizeBits);

   {
      // emit generic header for function lists
      if (1) {
         out << "fn-list{";
         out << "len: " << TargetFns.size();
#ifdef INCLUDE_ABANDONED
//          out << "; dom: " << s_DefaultTaskStr; // domain
//          out << "; sym: " << PointGroup.Name(NAMETYPE_Symbol);
//          out << "; deg: " << l;
//          out << "; type: " << "pom"; // polynomial, orthonormal, group-averaged
#endif // INCLUDE_ABANDONED
         out << "; type: " << pFunctionListType; // polynomial, orthonormal, group-averaged
         out << "} = [\n";
      }
#ifdef INCLUDE_ABANDONED
//       std::string ind1 = " ";
//       std::string ind2 = "  ";
#endif // INCLUDE_ABANDONED
      char const *ind1 = "";
      char const *ind2 = "";
      {
         unsigned nDigits10 = std::numeric_limits<FScalar>::max_digits10;
         if (nDigits10 == 0)
            nDigits10 = unsigned(.99 + 0.30103 * g_ScalarFloatMantissaBits); // 0.30103: that's about ln(2)/ln(10).
         out.precision(nDigits10);
         out << std::scientific;
//          out << std::right;

         for (size_t iFn = 0; iFn != TargetFns.size(); ++ iFn) {
            FIoTargetFn const
               &TargetFn = TargetFns[iFn];
            FPolynomialN const
               &poly = *TargetFn.poly;
            out << ind1 << "";
            out << "poly{";
            out << "len: " << poly.size();
            out << "; gens: xyz " << FPolynomialN::nGen;
            // ^-- for some groups, e.g., ico, it could make sense to employ a basis
            //     of monomials of coordinates from a non-orthogonal basis frame.
            //     Haven't implemented or tried yet, though.
            if (1) {
               out << "; lambda: " << fmt::format("{:.6g}", TargetFn.lambda);
            }
            if (1) {
               out << "; ethr: " << fmt::format("{:.2e}", poly.GetCullCriteria().GetThrCoeffNeglect());
            }
            if (1) {
               out << "; coeffs: " << (bNormalizeBasisFn? "normed" : "raw");
            }
            out << "} = [\n";
            if (1) {
               // reorder terms by monomial types.
               FPolynomialTermsExport
                  terms;
               for (FPolynomialN::const_iterator it = poly.begin(); it != poly.end(); ++ it)
                  terms[it->first] = it->second;
               for (FPolynomialTermsExport::const_iterator it = terms.begin(); it != terms.end(); ++ it) {
                  FMonomialN const &m = it->first;
                  out << ind2;
                  for (size_t igen = 0; igen != FPolynomialN::nGen; ++ igen) {
                     out << int(m[igen]) << " ";
                  }
                  out << CoeffForExport(it->second, m);
                  out << "\n";
               }
            } else {
               for (FPolynomialN::const_iterator it = poly.begin(); it != poly.end(); ++ it) {
                  out << ind2;
                  FMonomialN const &m = it->first;
                  for (size_t igen = 0; igen != FPolynomialN::nGen; ++ igen) {
                     out << int(m[igen]) << " ";
                  }
                  out << CoeffForExport(it->second, m);
                  out << "\n";
               }
            }
   //                   out << (igen == 0? "":" ") << int(m[igen]);
            out << ind1 << "]" << ((iFn == TargetFns.size() - 1)? "" : ",") << "\n";
         }
      }
      out << "]\n"; // fn-list
   }
   out << "]\n"; // target-space
}



void PrintOverlapSpectrum(ct::FLog &io, FEigh const &Sed, std::string const &Desc, std::string const &EwSymbol, FScalar const &ThreshEw, unsigned nRetained)
{
   // hm... maybe add it to FLog? Could add a few default rule levels and add a SetupRule(0, .., .., ..) function.
   auto pEmitRule = [&](int iRuleConfig) {
//       size_t LineWidth = 120;
      size_t LineWidth = 80;
      char const *pBarChar = (iRuleConfig == 0)? "─" : "━";
      io.Write(" {}", Repeat(LineWidth-1, pBarChar));
   };
   io.WriteLine();
   size_t
      BasisSize = size_t(Sed.lambda.size());
   cx_assert_rt(nRetained <= BasisSize);
   pEmitRule(1);
   io.Write(" ░░ List of {0} eigenvalues {1} above or near |{1}| ≥ {2:8.2e}:    ╎ K = {3}, N = {4}",
      Desc, EwSymbol, double(ThreshEw), BasisSize, nRetained);
   pEmitRule(0);
   auto pEmitEw = [&](size_t iEw) {
      char const
         *pEwFmt = "          {}{:<5} = {:11.2e}{}",
         *pNote = "";
      if (iEw < nRetained)
//          pNote = "  -->  accept"; // "✔✘"
         pNote = "    [✔]"; // "✔✘"
      io.Write(pEwFmt, EwSymbol, iEw, double(Sed.lambda[iEw]), pNote);
   };
   size_t
      nLeading = 0,
      nBelowThresh = 0;
   FScalar
      // print some extra, to see if we got a clear cutoff or if things got murky.
//       ThreshEwPrint = 1e-1 * ThreshEw;
//       ThreshEwPrint = 1e-10 * ThreshEw;
      ThreshEwPrint = 0;
   for (size_t iEw = 0; iEw < BasisSize; ++ iEw) {
//       if (abs(Sed.lambda[iEw]) > ThreshEw) {
// ^-- large negative ones would also be a serious issue. And we'd
//     certainly want to know about their presence. However, I think the
//     print of trailing eigvals below should do the trick just fine.
      if (nBelowThresh >= 4 || !(Sed.lambda[iEw] >= ThreshEwPrint))
         break;
      if (!(Sed.lambda[iEw] >= ThreshEw))
         nBelowThresh += 1;
      pEmitEw(iEw);
      nLeading = iEw + 1;
   }
   if (nLeading < nRetained) {
      io.Write("\n ...{} additional retained eigenvalues:\n", nRetained - nLeading);
      for (size_t iEw = nLeading; iEw < nRetained; ++ iEw)
         pEmitEw(iEw);
      nLeading = nRetained;
   }
   {
      size_t
         iFirst = BasisSize,
         nTrailing = std::min(size_t(nRetained), size_t(4));
      if (nTrailing > iFirst)
         nTrailing = iFirst;
      iFirst -= nTrailing;
      if (iFirst < nLeading) // don't re-print ones we've already done on the other side.
         iFirst = nLeading;
      nTrailing = BasisSize - iFirst;
      if (nTrailing != 0) {
         io.Write("\n ...trailing {} eigenvalues:\n", nTrailing);
         for (size_t iEw = iFirst; iEw < BasisSize; ++ iEw)
            pEmitEw(iEw);
      }
   }
   pEmitRule(0);
//    io.WriteLine();
//    io.Write("mystery?");
   io.Write("");
}



// clear out the renormalization (m_Factor) and weight (m_Degeneracy) factors inside the monomials to unity.
void ClearBaseFactors(FMonomialList &ml) {
   for (size_t i = 0; i != ml.size(); ++ i) {
      ml[i].m_Factor = 1;
      ml[i].m_Degeneracy = 1;
   }
};


// void PreprocessTargetFunctions(FPointGroup const &PointGroup, int l, FTargetSpacePreprocessOptions const &Options)
// {
//    FGridSearchOptions
//       DummyOptions(fmt::format("{{point-group: {}}}", PointGroup.Name(NAMETYPE_Symbol)));
//    FGridSearchOptions const
//       *pGridSearchOptions(&DummyOptions);
// //       = &Options;
// //    pGridSearchOptions = 0; // FIXME: <-- This turns of prescreening of target functions in MakeSymmetryUniqueMonomialList!!
//    int lmin = l-1;
//    if (lmin < 0) lmin = 0;
//    int lmax = l;
// //    g_UseGroupSymmetrizedFn = true;
//    g_UseGroupSymmetrizedFn = false;
//    FMonomialList ml = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, -1, pGridSearchOptions);
// //    g_UseGroupSymmetrizedFn = true;
//    // ^-- do I really want symmetry-unique ones only?
//    // UPDATE: I think I have to go with this, and then add the group average
//    // members up. Otherwise too expensive.
//    FMonomialList ml_all = MakeMonomialList(lmin, lmax);
//    #pragma omp critical
//    {
//       if (s_CancelPrimitiveRenormFactors) {
//             std::cout << "!WARNING: clearing out monomial list renorm and weight factors" << std::endl;
//          ClearBaseFactors(ml);
//          ClearBaseFactors(ml_all);
//       }
//       if (0) {
//          std::cout << "!WARNING: retaining ALL monomials in at [lmax-1,lmax]" << std::endl;
//          ml = ml_all;
//       }
//    }
//
//    size_t nMonomialsOrig = ml.size();
//    ComputeAndStoreGroupAverages(ml, PointGroup, true); // last: delete zero-averaged polynomials?
//    FDenseMatrix
//       S = MakeOverlapMatrix(ml);
// //    std::cout << fmt::format("G = {}  l = {}  #monomials = {}", PointGroup.Name(NAMETYPE_Symbol), l, ml.size()) << std::endl;
// //    std::cout << fmt::format("| S.shape = ({},{})", S.rows(), S.cols()) << std::endl;
//    size_t nNonZero = 0;
//    FIoPolynomialList
//       // output list of orthogonalized significant (i.e., non-redundant) group-symmetrized
//       // target functions for the current symmetry group and input monomial list.
//       OrthTargetFns;
//    OrthTargetFns.reserve(ml.size());
//    if (S.rows() != 0) {
//       FEigh
//          Sed(S, FEigh::LargeEwFirst);
//       using std::abs;
//       for (size_t iEw = 0; iEw < size_t(Sed.lambda.size()); ++ iEw)
//          if (abs(Sed.lambda[iEw]) > g_ThrAlmostZero) {
//             nNonZero += 1;
//             FPolynomialNPtr
//                pPolyEv = new FPolynomialN();
//             pPolyEv->SetThrCoeffNeglect(g_ThrAlmostZero, true); // true: AutoPurge on.
//             for (size_t iComp = 0; iComp != size_t(Sed.V.rows()); ++ iComp) {
//                // TODO: maybe use a slightly less ugly way to orthogonalize the functions.
//                pPolyEv->Add(Sed.V(iComp,iEw)/sqrt(Sed.lambda[iEw]) * ml[iComp].m_Factor, ml[iComp].m_RgSymmetrized);
//                #pragma omp critical
//                {
//                   if (0) {
//                      std::cout << fmt::format("   G = {:2}  l = {:<3}  iew = {:<3}  ew = {:18.12f}  iterm: {}  fthr(p) = {:.2e}  fthr(comp) = {:.2e}",
//                         PointGroup.Name(NAMETYPE_Symbol), l, iEw, Sed.lambda[iEw], iComp, pPolyEv->GetThrCoeffNeglect(), ml[iComp].m_RgSymmetrized.GetThrCoeffNeglect())
//                      << std::endl;
//                   }
//                }
//             }
//             pPolyEv->SetThrCoeffNeglect(g_ThrAlmostZero, true); // true: AutoPurge on.
//             pPolyEv->Purge();
//             OrthTargetFns.push_back(FIoTargetFn(Sed.lambda[iEw], pPolyEv));
//          }
//    }
//    #pragma omp critical
//    {
//       std::cout << fmt::format("G = {:2}  l = {:<3}  #monomials[lmin: {:<3} lmax: {:<3} all: {:<5} req: {:<5} nonz-P_G = {:<5}] #linearly-independent = {}", PointGroup.Name(NAMETYPE_Symbol), l, lmin,lmax, ml_all.size(), nMonomialsOrig, ml.size(), nNonZero) << std::endl;
//    }
//
//    if (Options.Export)
//       WriteTargetFnSpaceDataFile(PointGroup, l, OrthTargetFns, Options);
// }
// // template<class FScalar, class FDenseMatrix>
// // void LocalizeVectorsT(FDenseMatrix &CVec, size_t nRows, size_t nCols, size_t const *pGroupOffsets, size_t nGroups, FLocalizeOptions const &Opt)

#define S_TIME_SECTION(a,b,c) FTimeSection1 IR_UNIQUE_NAME(a,b,c)

struct FTimeSection1 {
   explicit FTimeSection1(ct::FLog *pLog, bool on, std::string message) : m_pLog(pLog), m_Message(message) {
      if (!on) { m_pLog = 0; return; }
   }
//    explicit FTimeSection1(ct::FLog *pLog, bool on, fmt::BasicStringRef<Char> format, fmt::ArgList args) : m_pLog(pLog) {
//       if (!on) { m_pLog = 0; return; }
//       fmt::MemoryWriter w;
//       fmt::BasicFormatter<Char>(w).format(format, args);
//       m_Message = w.str();
//    }
   ~FTimeSection1() {
      if (!m_pLog) return;
      double t = double(m_tSection);
      {
         m_pLog->WriteTiming(m_Message, t);
      }
   }
protected:
   ct::FLog *m_pLog;
   std::string m_Message;
   ct::FTimer m_tSection;
};


std::tuple<FIoPolynomialList, std::string> ConstructOrthogonalTargetFns(ct::FLog &io, int l, FMonomialList &ml, FPointGroup const &PointGroup, FPrintLevel plInfo, FPrintLevel plTiming)
{
   S_TIME_SECTION(&io, plTiming.BasicQ(), "orth polynomial basis (total)");
   using std::abs;
//    size_t nMonomialsOrig = ml.size();
   {
      S_TIME_SECTION(&io, plTiming.MoreQ(), "group average of monomials");
      ComputeAndStoreGroupAverages(ml, PointGroup, true); // last: delete zero-averaged polynomials?

      if (0) {
         // with the gammas this is really slow. for the moment i put the iterations.
         // that should also be reasonably stable. I just wonder *why* this is faster...
         FPolynomialN::FCullCriteria
            CullCriteria(g_ThrVerySmall, 0, FPolynomialN::FEstimateSquaredNormFn(EstimateMonomialNormSq_S2));
         for (auto &&mi : ml)
            mi.m_RgSymmetrized.SetCullCriteria(CullCriteria, true);
      }
   }
   FDenseMatrix
      S;
   {
//       S_TIME_SECTION(&io, plTiming.MoreQ(), "<m_i|P_G|m_j> overlap matrix ({0},{0})", ml.size());
      S_TIME_SECTION(&io, plTiming.MoreQ(), "<m_i|P_G|m_j> overlap matrix");
      S = MakeOverlapMatrix(ml);
   }
//    std::cout << fmt::format("G = {}  l = {}  #monomials = {}", PointGroup.Name(NAMETYPE_Symbol), l, ml.size()) << std::endl;
//    std::cout << fmt::format("| S.shape = ({},{})", S.rows(), S.cols()) << std::endl;
   if (0) {
      FScalar fErr = Rmsd(S - S.transpose());
      io.WriteInfoExpf("S(ml) symmetry violation/unit", double(fErr));
   }
   cx_assert_rt(size_t(S.rows()) == ml.size());
   size_t
      BasisSize = size_t(S.rows()),
      nNonZero = 0;
   FIoPolynomialList
      // output list of orthogonalized significant (i.e., non-redundant) group-symmetrized
      // target functions for the current symmetry group and input monomial list.
      OrthTargetFns;
   std::string
      sLocalizationReport;
   if (BasisSize != 0) {
      FDenseVector
         OrthEw;
      FDenseMatrix
         OrthCoeff;
      {
         S_TIME_SECTION(&io, plTiming.MoreQ(), "overlap spectral decomp");
         FEigh
            Sed(S, FEigh::LargeEwFirst);
         FScalar
            ThreshOverlap = g_ThrAlmostZero;
         cx_assert_rt(size_t(Sed.lambda.size()) == BasisSize);
         for (size_t iEw = 0; iEw < size_t(BasisSize); ++ iEw) {
            if (!(Sed.lambda[iEw] >= ThreshOverlap))
               break;
            nNonZero += 1;
         }

         if (1) {
            PrintOverlapSpectrum(io, Sed, "S_{{ij}} = ⟨m_i|P_g|m_j⟩", "λ", ThreshOverlap, nNonZero);
         }

         OrthEw = Sed.lambda.head(nNonZero);
         // set `C := V * diag(1/sqrt(lambda))`
         OrthCoeff = Sed.V.leftCols(nNonZero) * OrthEw.array().rsqrt().matrix().asDiagonal();
      }
//       io.WriteInfoExpf("Sed/V unitarity violation", double(RmsdFromIdentity(Sed.V.transpose() * Sed.V)));
//       io.WriteInfoExpf("Orth/C unitarity violation", double(RmsdFromIdentity(OrthCoeff.transpose() * S * OrthCoeff)));
#ifdef INCLUDE_ABANDONED
      if (0) {
         // FIXME: no… vectors are not orthonormal in the right way. doesn't work like this :(
         std::stringstream sout;
         size_t nRows = size_t(OrthCoeff.rows());
         std::vector<size_t> pGroupOffsets(nRows+1);
         for (size_t i = 0; i < nRows+1; ++ i)
            pGroupOffsets[i] = i;
         ct_loc_t::FLocalizeOptions LocOpt;
         LocOpt.pxout = &sout;
         LocOpt.Verbosity = 3;
         ct_loc_t::LocalizeVectorsT<FScalar,FDenseMatrix>(OrthCoeff, nRows, nNonZero, &pGroupOffsets[0], nRows, LocOpt);
         sLocalizationReport = sout.str();
         // update the scalar characteristic values
         OrthEw = (OrthCoeff.transpose() * (S * OrthCoeff).eval()).diagonal();
      }
#endif // INCLUDE_ABANDONED

      cx_assert_rt(size_t(OrthCoeff.rows()) == BasisSize && size_t(OrthCoeff.cols()) == nNonZero && ml.size() == BasisSize);

      {  S_TIME_SECTION(&io, plTiming.AllQ(), "form polynomials");

         FPolynomialN::FCullCriteria
            CullCriteria(g_ThrAlmostZero); // Note: changing this to `sqr(g_ThrAlmostZero)` does not appear to have any remaining influence at this point.
   //       CullCriteria = FPolynomialN::FCullCriteria(0,0);
         if (1) {
   //          CullCriteria = FPolynomialN::FCullCriteria(g_ThrVerySmall, 0, FPolynomialN::FEstimateSquaredNormFn(EstimateMonomialNormSq_S2));
   //          CullCriteria = FPolynomialN::FCullCriteria(g_ThrAlmostZero, 0, FPolynomialN::FEstimateSquaredNormFn(EstimateMonomialNormSq_S2));
            using std::sqrt;
            CullCriteria = FPolynomialN::FCullCriteria(sqrt(std::numeric_limits<FScalar>::epsilon()), 0, FPolynomialN::FEstimateSquaredNormFn(EstimateMonomialNormSq_S2));
            // ^-- with the normalization, the coefficients mostly come out scattering around unity.
         }
         OrthTargetFns.reserve(BasisSize);
   //       io.WriteInfoExpf("S/1(ml) vs S/n(ml) unit", double(Rmsd(S - MakeOverlapMatrix(ml))));
         for (size_t iEw = 0; iEw < nNonZero; ++ iEw) {
            FPolynomialNPtr
               pPolyEv = new FPolynomialN();
            for (size_t iComp = 0; iComp != size_t(OrthCoeff.rows()); ++ iComp) {
               // TODO: maybe use a slightly less ugly way to orthogonalize the functions?
               pPolyEv->Add(OrthCoeff(iComp,iEw) * ml[iComp].m_Factor, ml[iComp].m_RgSymmetrized);
#ifdef INCLUDE_ABANDONED
//                #pragma omp critical
//                {
//                   if (0) {
//                      std::cout << fmt::format("   G = {:2}  l = {:<3}  iew = {:<3}  ew = {:18.12f}  iterm: {}  fthr(p) = {:.2e}  fthr(comp) = {:.2e}",
//                         PointGroup.Name(NAMETYPE_Symbol), l, iEw, OrthEw[iEw], iComp, pPolyEv->GetThrCoeffNeglect(), ml[iComp].m_RgSymmetrized.GetThrCoeffNeglect())
//                      << std::endl;
//                      // Sed.lambda[iEw]
//                   }
//                }
#endif // INCLUDE_ABANDONED
            }
            pPolyEv->SetCullCriteria(CullCriteria, true); // true: AutoPurge on (i.e., do it on call automatically if necessary)
            // ^--- note: must be set AFTER the accumulation. Because PolyN
            //      updates cull criteria to whatever is the stronger one, and
            //      we may have things set up in such a way that these are very
            //      hard (or even exactly zero, retaining all terms) for the
            //      input polynomials.
            OrthTargetFns.emplace_back(FIoTargetFn(OrthEw[iEw], pPolyEv));
         }
      }
   }
   return {OrthTargetFns, sLocalizationReport};
}



FDenseMatrix MakeOverlapMatrix(std::vector<FPolynomialNPtr> const &pl)
{
   size_t
      N = pl.size();
   FDenseMatrix
      S(N,N);
   for (auto const &[i,pi] : enumerate(pl)) {
      for (auto const &[j,pj] : enumerate(pl)) {
         if (j < i) continue; // it's symmetric. We make just one half.
         FScalar sij = ComputeOverlapRaw(*pi, *pj);
         S(i,j) = sij;
         S(j,i) = sij;
      }
   }
   return S;
}

void EvalOrthogonalityViolation(FIoPolynomialList const &pIoPolys, std::string const &Desc) {
   std::vector<FPolynomialNPtr>
      pl;
   for (auto &&pIoPoly : pIoPolys)
      pl.push_back(pIoPoly.poly);
   FDenseMatrix
      // these should already be completely orthonormal as they come — without any
      // additional averaging or whatever.
      S_TargetFn = MakeOverlapMatrix(pl);
   io.WriteInfoExpf("Raw io polys unitarity loss", double(RmsdFromIdentity(S_TargetFn)), Desc);
}


static bool s_ClearWarningEmitted = false;

void PreprocessTargetFunctions(ct::FLog &io, FPointGroup const &PointGroup, int l, FTargetSpacePreprocessOptions const &Options)
{
   FPrintLevel
      plTargetFn = Options.PrintOptions[PRINT_TargetFn],
      plTiming = Options.PrintOptions[PRINT_Timing];
   S_TIME_SECTION(&io, plTiming.MoreQ(), "pre-process (total)");
   io.Write("");
   io.WriteProgramIntro(fmt::format("PRE-PROCESS TARGET SPACE: for {}, integration rules order {}", PointGroup.Name(NAMETYPE_Group), l));
   FGridSearchOptions
      DummyOptions(fmt::format("{{point-group: {}}}", PointGroup.Name(NAMETYPE_Symbol)));
   FGridSearchOptions const
      *pGridSearchOptions(&DummyOptions);
//       = &Options;
//    pGridSearchOptions = 0; // FIXME: <-- This turns of prescreening of target functions in MakeSymmetryUniqueMonomialList!!
   int lmin = l-1;
   if (lmin < 0) lmin = 0;
   int lmax = l;
//    g_UseGroupSymmetrizedFn = true;
//    g_UseGroupSymmetrizedFn = false;
   cx_assert_rt(!g_UseGroupSymmetrizedFn);
   FMonomialList ml, ml_all;
   {
      S_TIME_SECTION(&io, plTiming.MoreQ(), "basis monomial list");
      ml = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, -1, pGridSearchOptions);
      for (FMonomial &m : ml) m._ResetNorm();
   }
   size_t nMonomialsOrig = ml.size();
//    g_UseGroupSymmetrizedFn = true;
   // ^-- do I really want symmetry-unique ones only?
   // UPDATE: I think I have to go with this, and then add the group average
   // members up. Otherwise too expensive.
   {
      S_TIME_SECTION(&io, plTiming.MoreQ(), "full monomial list");
      ml_all = MakeMonomialList(lmin, lmax);
   }
   {
      if (s_CancelPrimitiveRenormFactors) {
         if (!s_ClearWarningEmitted)
            io.Write("!WARNING: clearing out monomial list renorm and weight factors");
         ClearBaseFactors(ml);
         ClearBaseFactors(ml_all);
         s_ClearWarningEmitted = true;
      }
      if (0) {
         io.Write("!WARNING: retaining ALL monomials in at [lmax-1,lmax]");
         ml = ml_all;
      }
   }

   auto [OrthTargetFns, sLocalizationReport] = ConstructOrthogonalTargetFns(io, l, ml, PointGroup, plTargetFn, plTiming);
   {
      io.Write("░ G = {:2}  l = {:<3}  #monomials[lmin: {:<3} lmax: {:<3} all: {:<5} req: {:<5} nonz-P_G = {:<5}] #linearly-independent = {}", PointGroup.Name(NAMETYPE_Symbol), l, lmin,lmax, ml_all.size(), nMonomialsOrig, ml.size(), OrthTargetFns.size());
      if (!sLocalizationReport.empty())
         io.Write("  ^-- Localization:\n{}\n", sLocalizationReport);
   }

   if (0)
      // gives same result as the other two below (but much slower)
      EvalOrthogonalityViolation(OrthTargetFns, fmt::format("G = {}, l = {} (SLOWW!!)", PointGroup.Name(NAMETYPE_Symbol), l));

   FPrintLevel
//       PrintLevel_TestReloadData = FPrintLevel::All;
      PrintLevel_TestReloadData = FPrintLevel::Most;
   FTargetFnListOptions
      ReloadDataTestOptions(PrintLevel_TestReloadData, Options.AccuracyPreference);
   if (0 && PrintLevel_TestReloadData) {
      {
         // try instanciating the target function list object directly from here, bypassing the load & save in-between.
         FTargetFnListPtr pTargetFns(new FTargetFnList_PreAveragedPolys(l, &PointGroup, ReloadDataTestOptions, OrthTargetFns));
         io.WriteCount("Target space basis functions", pTargetFns->size(), pTargetFns->FnTypeName());
         pTargetFns->WriteInfo(io);
      }
   }

   if (Options.Export) {
      S_TIME_SECTION(&io, plTiming.BasicQ(), "data export");
      WriteTargetFnSpaceDataFile(PointGroup, l, OrthTargetFns, Options);
   }

   if (1 && PrintLevel_TestReloadData) {
      S_TIME_SECTION(&io, plTiming.BasicQ(), "instanciate target space");
      {
         // try re-loading the file we just exported.
         FTargetFnListPtr
            pTargetFns(new FTargetFnList_PreAveragedPolys(l, &PointGroup, ReloadDataTestOptions, &DummyOptions));
         io.WriteCount("Target space basis functions", pTargetFns->size(), pTargetFns->FnTypeName());
         pTargetFns->WriteInfo(io);
      }
   }
//    io.Write("");
}



bool GetDataFileLine(std::istream &in, std::string &Line, ct::string_slice &slLine)
{
   if (!std::getline(in, Line))
      return false;
   slLine = ct::string_slice(Line.begin(), Line.end());
   slLine.trim();
   // comment line?
   if (slLine[0] == '#')
      // yes. truncate it — the test in the loops will be for `slLine.empty()` to ignore.
      slLine.last = slLine.first;
   return true;
}


struct FDataListReader {
   explicit FDataListReader(char const *pListName_, size_t nListEntriesDefault_ = 0) : m_pListName(pListName_), m_nListEntriesDefault(nListEntriesDefault_) {}
   void Read(ct::string_slice slHeaderLine, std::istream &in);

protected:
   virtual bool ReadProperty(ct::string_slice const &Name, ct::string_slice const &Content) = 0;
   // called after all input properties are processed. If this returns 'false',
   // the reading is aborted and no further `ReadEntry()`/data split is performed.
   // Default implementation returns 'true' unconditionally.
   virtual bool CheckProperties();
   virtual void ReadEntry(ct::string_slice slHeaderLine, std::istream &in) = 0;

private:
   char const
      *m_pListName;
   size_t
      // number of list entries to read if no `{len: …}` property is specified.
      m_nListEntriesDefault;
};


bool FDataListReader::CheckProperties()
{
   return true;
}


void FDataListReader::Read(ct::string_slice slHeaderLine, std::istream &in)
{
   if (!slHeaderLine.remove_if_at_end(" = ["))
      throw ct::FInputError(fmt::format("expected list '{}{{...}} = [...]' but found '{}'", m_pListName, slHeaderLine.to_str()));
   ct::FPropertyListStr
      Props(slHeaderLine, m_pListName);
   size_t
      nEntriesToRead = m_nListEntriesDefault;
   for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
      if (it->Name == "len") {
         nEntriesToRead = it->Content.to_size();
      } else if (!ReadProperty(it->Name, it->Content)) {
         throw ct::FInputError(fmt::format("unrecognized {} option '{}' (set to '{}')", m_pListName, it->Name.to_str(), it->Content.to_str()));
      }
   }
   if (!CheckProperties())
      return;
   std::string
      Line;
   ct::string_slice
      slLine;
   size_t
      iEntriesRead = 0;
   for (; iEntriesRead < nEntriesToRead; ) {
      if (!GetDataFileLine(in, Line, slLine))
         throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read {} list entry '{}' of {}", m_pListName, iEntriesRead+1, nEntriesToRead));
      if (slLine.empty())
         continue;
      ReadEntry(slLine, in);
      iEntriesRead += 1;
   }
   if (iEntriesRead != nEntriesToRead)
      throw ct::FInputError(fmt::format("while reading list data of '{}': insufficient elements provided (declared length = {}, actual length = {})", m_pListName, nEntriesToRead, iEntriesRead));
   if (!GetDataFileLine(in, Line, slLine))
      throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read end-of-list marker ']'"));
   if (!(slLine == "]" || slLine == "],"))
      throw ct::FInputError(fmt::format("expected end-of-list-marker ']' but found '{}'", slLine.to_str()));
}



std::pair<FPolynomialN::FMonomial, FScalar> ReadMonomialN(ct::string_slice slLine, std::istream &in, bool CoeffsAreNormalized)
{
   char const *pListName = "monomial";
   ct::long_split_result sr;
   slLine.split(sr, ' ');
   if (sr.size() != 4)
      throw ct::FInputError(fmt::format("expected {} '<x-exp> <y-exp> <z-exp> <coeff>' '{}'", pListName, slLine.to_str()));
   std::pair<FPolynomialN::FMonomial, FScalar>
      res(FPolynomialN::FMonomial(sr[0].to_int(), sr[1].to_int(), sr[2].to_int()), ParseFloat(sr[3].to_str()));
   if (CoeffsAreNormalized)
      res.second /= CalcMonomialUnitSphereNorm(res.first);
   return res;
}


struct FReaderPolynomialN : public FDataListReader {
   FPolynomialNPtr
      m_Poly;
   FScalar
      m_lambda;
   bool
      m_CoeffsAreNormalized;

   explicit FReaderPolynomialN() : FDataListReader("poly"), m_Poly(new FPolynomialN), m_lambda(0), m_CoeffsAreNormalized(false) {}

   bool ReadProperty(ct::string_slice const &Name, ct::string_slice const &Content) {
      if (Name == "gens") {
         if (Content != "xyz 3")
            throw ct::FInputError(fmt::format("polynomial type '{}: {}' not supported", Name.to_str(), Content.to_str()));
      } else if (Name == "lambda") {
         m_lambda = FScalar(Content.to_float());
      } else if (Name == "ethr") {
//          if (0) {
//             m_Poly->SetThrCoeffNeglect(FScalar(Content.to_float()));
//             // ^-- Hmpf. that won't really work if cross-referencing data sets with different
//             //     float precisions… (note that the thresholds are sticky and propagate
//             //     to derived calculation results using those polys…)
//          }
      } else if (Name == "coeffs") {
         if (Content == "raw")
            m_CoeffsAreNormalized = false;
         else if (Content == "normed")
            m_CoeffsAreNormalized = true;
         else
            throw ct::FInputError(fmt::format("polynomial coefficient type '{}: {}' not supported", Name.to_str(), Content.to_str()));
      } else
         return false; // not recognized
      return true;
   }

   void ReadEntry(ct::string_slice slHeaderLine, std::istream &in) {
      std::pair<FPolynomialN::FMonomial, FScalar>
         res = ReadMonomialN(slHeaderLine, in, m_CoeffsAreNormalized);
      m_Poly->Add(res.second, res.first);
   }
};

struct FReaderFnList : public FDataListReader {
   FIoTargetFnSpacePtr
      m_pSpace;
   std::string
      m_sFnClass;

   explicit FReaderFnList(FIoTargetFnSpacePtr &pSpace) : FDataListReader("fn-list"), m_pSpace(pSpace) {}

   bool ReadProperty(ct::string_slice const &Name, ct::string_slice const &Content) {
      if (Name == "type") {
         ct::long_split_result
            sr;
         Content.split_list(sr);
         for (size_t il = 0; il != sr.size(); ++ il) {
            if (sr[il] == "orth") {
               m_pSpace->Flags |= FNSPACE_Orthogonal;
            } else if (sr[il] == "mean") {
               m_pSpace->Flags |= FNSPACE_GroupAveraged;
            } else if (sr[il] == "poly") {
               // that's the only type we support reading atm — Cartesian
               // polynomials. So we do not actually do anything here.
               m_sFnClass = sr[il];
            } else {
               throw ct::FInputError(fmt::format("unrecognized function type fragment '{}' in '{}: {}'", sr[il], Name.to_str(), Content.to_str()));
            }
         }
      } else
         return false;
      return true;
   }

   void ReadEntry(ct::string_slice slHeaderLine, std::istream &in) {
      if (m_sFnClass != "poly")
         throw ct::FInputError(fmt::format("function class '{}' empty or not supported", m_sFnClass));
      FReaderPolynomialN
         res;
      res.Read(slHeaderLine, in);
      m_pSpace->pFns.push_back(res.m_Poly);
      m_pSpace->pLambda.push_back(res.m_lambda);
   }
};


struct FReaderTargetFnSpace : public FDataListReader {
   FIoTargetFnSpacePtr
      m_pSpace;
   FPointGroup const
      *m_pPointGroup;
   int
      m_TargetDegree;
   std::string
      // if != 0, will set to a description of the problem if `m_pSpace` is set to
      // 0 to indicate incorrect/unusable data.
      *m_pErrorDesc;

   explicit FReaderTargetFnSpace(FIoTargetFnSpacePtr &pSpace, FPointGroup const *pPointGroup, int TargetDegree, std::string *pErrorDesc = 0)
   : FDataListReader("target-fn-space", 1), m_pSpace(pSpace), m_pPointGroup(pPointGroup), m_TargetDegree(TargetDegree), m_pErrorDesc(pErrorDesc)
   {}

   bool ReadProperty(ct::string_slice const &Name, ct::string_slice const &Content) {
      if (Name == "point-group") {
         m_pSpace->sPointGroup = Content.to_str();
      } else if (Name == "domain") {
         m_pSpace->sDomain = Content.to_str();
      } else if (Name == "degree") {
         m_pSpace->Degree = int(Content.to_int());
      } else if (Name == "fp") {
         m_pSpace->FloatSizeBits = int(Content.to_int());
      } else {
         // ignore properties we do not recognize? Maybe it's a file from a later version.
         return false;
      }
      return true;
   }

   bool CheckProperties() {
      bool Ok = true;
      if (ptrdiff_t(m_pSpace->FloatSizeBits) < ptrdiff_t(g_ScalarFloatSizeBits)) {
         // data may tabulated, but not with sufficient precision for the current run.
         Ok = false;
         if (m_pErrorDesc)
            *m_pErrorDesc = fmt::format("insufficient precision in tabulated data (fp = {}, but need fp >= {})", m_pSpace->FloatSizeBits, g_ScalarFloatSizeBits);
      } if (ptrdiff_t(m_pSpace->Degree) != ptrdiff_t(m_TargetDegree)) {
         // requested integration order (degree) and actually tabulated degree do not match.
         Ok = false;
         if (m_pErrorDesc)
            *m_pErrorDesc = fmt::format("data file declares unexpected function subspace (deg = {}, expected deg = {})", m_pSpace->Degree, m_TargetDegree);
      } if (m_pSpace->sPointGroup != m_pPointGroup->Name(NAMETYPE_Symbol)) {
         // hm… this one may need some adjustments. Using target spaces from subgroups for larger groups is actually fine.
         Ok = false;
         if (m_pErrorDesc)
            *m_pErrorDesc = fmt::format("data file declares unexpected symmetry (group = {}, expected group = {})", m_pSpace->sPointGroup, m_pPointGroup->Name(NAMETYPE_Symbol));
      } if (!Ok)
         m_pSpace = 0;
      return Ok;
   }

   void ReadEntry(ct::string_slice slHeaderLine, std::istream &in) {
      FReaderFnList
         res(m_pSpace);
      res.Read(slHeaderLine, in);
   }
};

#ifdef INCLUDE_ABANDONED
// std::pair<FPolynomialNPtr, FScalar> ReadPolynomialN(ct::string_slice slHeaderLine, std::istream &in)
// {
//    char const *pListName = "poly";
//    FPolynomialNPtr
//       pPoly(new FPolynomialN);
//    FScalar
//       lambda(0);
//    if (!slHeaderLine.remove_if_at_end(" = ["))
//       throw ct::FInputError(fmt::format("expected list '{}{{...}} = [...]' but found '{}'", pListName, slHeaderLine.to_str()));
//    ct::FPropertyListStr
//       Props(slHeaderLine, pListName);
//    size_t
//       nEntriesToRead = 0;
//    for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
//       if (it->Name == "len") {
//          nEntriesToRead = it->Content.to_size();
//       } else if (it->Name == "gens") {
//          if (it->Content != "xyz 3")
//             throw ct::FInputError(fmt::format("polynomial type '{}: {}' not supported", it->Name.to_str(), it->Content.to_str()));
//       } else if (it->Name == "lambda") {
//          lambda = FScalar(it->Content.to_float());
//       } else if (it->Name == "ethr") {
//          if (0) {
//             pPoly->SetThrCoeffNeglect(FScalar(it->Content.to_float()));
//             // ^-- nah.. that won't really work if cross-referencing data sets with different
//             //     float precisions... (note that the thresholds are sticky and propagate
//             //     to derived calculation results using those polys...)
//          }
//       } else {
//          // ignore properties we do not recognize? Maybe it's a file from a later version.
//          throw ct::FInputError(fmt::format("unrecognized {} option '{}' (set to '{}')", pListName, it->Name.to_str(), it->Content.to_str()));
//       }
//    }
//    std::string
//       Line;
//    ct::string_slice
//       slLine;
//    size_t
//       iEntriesRead = 0;
//    for (; iEntriesRead < nEntriesToRead; ) {
//       if (!GetDataFileLine(in, Line, slLine))
//          throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read {} list entry '{}' of {}", pListName, iEntriesRead+1, nEntriesToRead));
//       if (slLine.empty())
//          continue;
//       std::pair<FPolynomialN::FMonomial, FScalar>
//          res = ReadMonomialN(slLine, in);
//       pPoly->Add(res.second, res.first);
//       iEntriesRead += 1;
//    }
//    if (iEntriesRead != nEntriesToRead)
//       throw ct::FInputError(fmt::format("while reading list data of '{}': insufficient elements provided (declared length = {}, actual length = {})", pListName, nEntriesToRead, iEntriesRead));
//    if (!GetDataFileLine(in, Line, slLine))
//       throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read end-of-list marker ']'"));
//    if (!(slLine == "]" || slLine == "],"))
//       throw ct::FInputError(fmt::format("expected end-of-list-marker ']' but found '{}'", slLine.to_str()));
//    return std::pair<FPolynomialNPtr, FScalar>(pPoly, lambda);
// }
//
// void ReadFnList(FIoTargetFnSpacePtr pSpace, ct::string_slice slHeaderLine, std::istream &in)
// {
//    char const *pListName = "fn-list";
//    if (!slHeaderLine.remove_if_at_end(" = ["))
//       throw ct::FInputError(fmt::format("expected list '{}{{...}} = [...]' but found '{}'", pListName, slHeaderLine.to_str()));
//    ct::FPropertyListStr
//       Props(slHeaderLine, pListName);
//    size_t
//       nEntriesToRead = 0;
//    ct::string_slice
//       slFnClass;
//    for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
//       if (it->Name == "len") {
//          nEntriesToRead = it->Content.to_size();
//       } else if (it->Name == "type") {
//          ct::long_split_result
//             sr;
//          it->Content.split_list(sr);
//          for (size_t il = 0; il != sr.size(); ++ il) {
//             if (sr[il] == "orth") {
//                pSpace->Flags |= FNSPACE_Orthogonal;
//             } else if (sr[il] == "mean") {
//                pSpace->Flags |= FNSPACE_GroupAveraged;
//             } else if (sr[il] == "poly") {
//                // that's the only type we support reading atm---cartesian
//                // polynomials. so we don't actually do anything here.
//                slFnClass = sr[il];
//             } else {
//                throw ct::FInputError(fmt::format("unrecognized function type fragment '{}' in '{}: {}'", sr[il], it->Name.to_str(), it->Content.to_str()));
//             }
//          }
//       } else {
//          // ignore properties we do not recognize? Maybe it's a file from a later version.
//          throw ct::FInputError(fmt::format("unrecognized {} option '{}' (set to '{}')", pListName, it->Name.to_str(), it->Content.to_str()));
//       }
//    }
//    std::string
//       Line;
//    ct::string_slice
//       slLine;
//    size_t
//       iEntriesRead = 0;
//    pSpace->pFns.reserve(nEntriesToRead);
//    for (; iEntriesRead < nEntriesToRead; ) {
//       if (!GetDataFileLine(in, Line, slLine))
//          throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read {} list entry '{}' of {}", pListName, iEntriesRead+1, nEntriesToRead));
//       if (slLine.empty())
//          continue;
//       std::pair<FPolynomialNPtr, FScalar>
//          res = ReadPolynomialN(slLine, in);
//       pSpace->pFns.push_back(res.first);
//       pSpace->pLambda.push_back(res.second);
//       iEntriesRead += 1;
//    }
//    if (iEntriesRead != nEntriesToRead)
//       throw ct::FInputError(fmt::format("while reading list data of '{}': insufficient elements provided (declared length = {}, actual length = {})", pListName, nEntriesToRead, iEntriesRead));
//    if (!GetDataFileLine(in, Line, slLine))
//       throw ct::FInputError(fmt::format("I/O issue (eof?) while attempting to read end-of-list marker ']'"));
//    if (!(slLine == "]" || slLine == "],"))
//       throw ct::FInputError(fmt::format("expected end-of-list-marker ']' but found '{}'", slLine.to_str()));
// }


// FIoTargetFnSpacePtr TryReadTargetFnSpaceDataFile(FPointGroup const &PointGroup, int l)
// {
//    FIoTargetFnSpacePtr
//       pSpace(new FIoTargetFnSpace);
//    std::string
//       FileName = MakeTargetSpaceDefaultFileName(PointGroup, l);
//    std::ifstream
//       in(FileName.c_str());
//    if (!in.good())
//       // failed to open the file (if one exists)
//       return 0;
//    std::string
//       Line;
//    int
//       iPhase = 0;
//    try {
//       for ( ; ; ) {
//          ct::string_slice slLine;
//          if (!GetDataFileLine(in, Line, slLine))
//             break;
//          if (slLine.empty())
//             continue; // ignore empty lines (and comment lines)
//          if (slLine == "]" || slLine == "],") {
//             iPhase -= 1;
//             if (iPhase < 0)
//                throw ct::FInputError("unbalanced data lists ([...,...])");
//             continue;
//          }
//
//          if (iPhase == 0) {
//             ct::string_slice
//                slTargetSpaceDecl = slLine;
//             if (!slTargetSpaceDecl.remove_if_at_end(" = ["))
//                throw ct::FInputError(fmt::format("expected list '<id>{...} = [...]' but found '{}'", slTargetSpaceDecl.to_str()));
//             ct::FPropertyListStr
//                Props(slTargetSpaceDecl, "target-fn-space");
//             for (ct::FPropertyListStr::iterator it = Props.begin(); it != Props.end(); ++ it) {
//                if (it->Name == "point-group") {
//                   pSpace->sPointGroup = it->Content.to_str();
//                } else if (it->Name == "domain") {
//                   pSpace->sDomain = it->Content.to_str();
//                } else if (it->Name == "degree") {
//                   pSpace->Degree = int(it->Content.to_int());
//                } else if (it->Name == "fp") {
//                   pSpace->FloatSizeBits = int(it->Content.to_int());
//                } else {
//                   // ignore properties we do not recognize? Maybe it's a file from a later version.
//                   throw ct::FInputError(fmt::format("unrecognized target-fn-space option '{}' (set to '{}')", it->Name.to_str(), it->Content.to_str()));
//                }
//             }
//             if (ptrdiff_t(pSpace->FloatSizeBits) < ptrdiff_t(g_ScalarFloatSizeBits))
//                // data may tabulated, but not with sufficient precision for the current run.
//                return 0;
//             if (ptrdiff_t(pSpace->Degree) != ptrdiff_t(l))
//                // requested integration order (degree) and actually tabulated degree do not match.
//                return 0;
//             iPhase += 1;
//             continue;
//          } else if (iPhase == 1) {
//             ReadFnList(pSpace, slLine, in);
//          }
//       }
//    } catch (ct::FInputError const &e) {
//       throw ct::FInputError(fmt::format("while reading target space data file '{}': {}", FileName, e.what()));
//    }
//    return pSpace;
// }
#endif // INCLUDE_ABANDONED


FIoTargetFnSpacePtr TryReadTargetFnSpaceDataFile(FPointGroup const &PointGroup, int l, std::string *pFailureReason)
{
   FIoTargetFnSpacePtr
      pSpace(new FIoTargetFnSpace);
   std::string
      FileName = MakeTargetSpaceDefaultFileName(PointGroup, l),
      ErrorMsg("generic error");
   pSpace->SourceFile = FileName;
   std::ifstream
      in(FileName.c_str());
   if (!in.good()) {
      // failed to open the file (if one exists)
      if (pFailureReason)
         *pFailureReason = fmt::format("data file missing ('{}')",
            MakeTargetSpaceDefaultFileName(PointGroup, l, FILENAME_OnlyBaseName));
      return 0;
   }
   try {
      std::string
         Line;
      for ( ; ; ) {
         ct::string_slice slLine;
         if (!GetDataFileLine(in, Line, slLine))
            break;
         if (slLine.empty())
            continue; // ignore empty lines (and comment lines)
         FReaderTargetFnSpace
            res(pSpace, &PointGroup, l, &ErrorMsg);
         res.Read(slLine, in);
         if (res.m_pSpace == 0) {
            pSpace = 0;
            if (pFailureReason) {
               *pFailureReason = fmt::format("data file '{}': {}",
                  MakeTargetSpaceDefaultFileName(PointGroup, l, FILENAME_OnlyBaseName),
                  ErrorMsg);
            }
         }
         break;
      }
   } catch (ct::FInputError const &e) {
//       throw ETabulatedFnListError("∫(S2)", &PointGroup, l, fmt::format("while reading target space data file '{}': {}", FileName, e.what()));
      throw ETabulatedFnListError("∫(S2)", &PointGroup, l, FileName, e.what());
   }
   return pSpace;
}



// main driver routine for running aigg in preprocess mode (i.e., this one is called from `main()` more-or-less directly).
int RunTargetSpacePreprocess(FTargetSpacePreprocessOptions const &Options)
{
   if (0) {
      for (size_t iDegreeId = 0; iDegreeId != Options.lList.size(); ++ iDegreeId) {
         int l = Options.lList[iDegreeId];
         for (size_t iPointGroup = 0; iPointGroup != Options.PointGroups.size(); ++ iPointGroup) {
            FPointGroup const *pPointGroup = &*Options.PointGroups[iPointGroup];
            FIoTargetFnSpacePtr
               pFnSet = TryReadTargetFnSpaceDataFile(*pPointGroup, l);
            if (pFnSet) {
               std::cout << fmt::format("target-fns(lmax = {}, group = '{}') = {{nfn: {}; fp: {}; flags: {:03x}}}\n",
                  l, pPointGroup->Name(NAMETYPE_Symbol), pFnSet->pFns.size(), pFnSet->FloatSizeBits, pFnSet->Flags);
               size_t iFn = 0;
               for (FPolynomialNPtr const &pFn : pFnSet->pFns) {
                  std::cout << fmt::format("[{}] f(x,y,z) = ", iFn);
                  pFn->Print(std::cout);
                  std::cout << "\n";
                  iFn += 1;
               }
            }
         }
      }
      return 0;
   }

   typedef std::tuple<FPointGroup const *, int>
      FPreprocessJob;
   std::vector<FPreprocessJob>
      Jobs;

   for (size_t iDegreeId = 0; iDegreeId != Options.lList.size(); ++ iDegreeId) {
      int l = Options.lList[iDegreeId];
      for (size_t iPointGroup = 0; iPointGroup != Options.PointGroups.size(); ++ iPointGroup) {
         FPointGroup const *pPointGroup = &*Options.PointGroups[iPointGroup];
         if (Options.SkipExisting != FTargetSpacePreprocessOptions::SKIP_None) {
            // we're supposed to check what we already have in this case. Try to
            // load the existing data file, if present.
            FIoTargetFnSpacePtr
               pFnSet = TryReadTargetFnSpaceDataFile(*pPointGroup, l);
            bool
               IsAlreadyOk = false;
            if (pFnSet) {
               // if `pFnSet` returns != 0, then some sort of target space data set
               // for the space was already there. In this case, we should skip
               // re-generating the data file, depending on the options in SkipExisting:
               if (Options.SkipExisting == FTargetSpacePreprocessOptions::SKIP_All) {
                  // "any existing data set which works is good enough"
                  IsAlreadyOk = true;
               } else if (Options.SkipExisting == FTargetSpacePreprocessOptions::SKIP_Better) {
                  // don't overwrite data files earlier generated with higher precision
                  if (ptrdiff_t(pFnSet->FloatSizeBits) > ptrdiff_t(g_ScalarFloatSizeBits))
                     IsAlreadyOk = true;
               } else if (Options.SkipExisting == FTargetSpacePreprocessOptions::SKIP_BetterOrEqual) {
                  // don't overwrite data files earlier generated with higher or equal precision
                  if (ptrdiff_t(pFnSet->FloatSizeBits) >= ptrdiff_t(g_ScalarFloatSizeBits))
                     IsAlreadyOk = true;
               }
               if (IsAlreadyOk) {
                  std::cout << fmt::format("G = {:2}  l = {:<3}  preprocessing skipped (data exists: '{}' with fp = {}). Use --preprocess '{{...; skip:no}}' to force rebuild.", pPointGroup->Name(NAMETYPE_Symbol), pFnSet->Degree, pFnSet->SourceFile, pFnSet->FloatSizeBits) << std::endl;
                  continue;
               }
            }
         }
         // if we reached this place then either we were told to always
         // overwrite the data files, or either no working earlier data files
         // exist or they are not sufficient for what we were asked to do.
//          PreprocessTargetFunctions(*pPointGroup, l, Options);
         Jobs.push_back(FPreprocessJob(pPointGroup, l));
      }
   }


#pragma omp parallel for schedule(dynamic)
   for (int iJob = 0; iJob < int(Jobs.size()); ++ iJob) {
      FPreprocessJob const
         &Job = Jobs[size_t(iJob)];
      FPointGroup const
         *pPointGroup = std::get<0>(Job);
      int l = std::get<1>(Job);

      std::ostringstream io_buff;
      ct::FLogStdStream ThreadLocalLog(io_buff);
      ct::FLog *pLog = (omp_get_max_threads() == 1)? &io : &ThreadLocalLog;
      PreprocessTargetFunctions(*pLog, *pPointGroup, l, Options);
      if (pLog == &ThreadLocalLog) {
         #pragma omp critical
         io.w << io_buff.str();
         io.Flush();
      }
   }

   return 0;
}



FTargetFnList::FTargetFnList(FPointGroup const *pPointGroup_, FTargetFnListOptions const &GeneralOptions)
   : m_GeneralOptions(GeneralOptions), m_pPointGroup(pPointGroup_), m_pSymmetryOps(&pPointGroup_->GetOps()), m_PrintLevel(GeneralOptions.PrintLevel)
{}


void FTargetFnList::WriteInfo(ct::FLog &io) const {
   (void)io;
}

void FTargetFnList::Print(ct::FLog &io) const {
   (void)io;
}


std::string FTargetFnList::FnTypeName() const {
   return "functions";
}

std::string FTargetFnList::FnTypeAnnotation() const {
   return std::string();
}


FTargetFnList_RawMonomials::FTargetFnList_RawMonomials(size_t lmin, size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FGridSearchOptions const *pGridSearchOptions)
   : FTargetFnList(pPointGroup, GeneralOptions), m_lmin(lmin), m_lmax(lmax)
{
   // Compute set of target functions which our integration rule is to evaluate
   // exactly, and compute their overlap matrix and a suitable decomposition of
   // it to solve conditioning equations with.
   m_TargetFns = MakeSymmetryUniqueMonomialList(m_lmin, m_lmax, *m_pPointGroup, int(m_PrintLevel), pGridSearchOptions);
   // ^--- hmpf.. function screening turned off if no options object provided. Stupid.
   //      should fix some other day.

   if (pGridSearchOptions && pGridSearchOptions->TargetFnType != FGridSearchOptions::TARGETFN_RawMonomials) {
      // use at least a bit of an estimate of the monomial overlaps for
      // renormalization.
      FDenseMatrix
         S = this->EvalOverlapMatrix();
      if (g_UseGroupSymmetrizedFn) {
         // Group-averaging causes hard singularities in the function space,
         // because the parts of the functions which lie outside the space
         // projected to by
         //
         //    P̂_G = ∑_g ∑_j |R_g mⱼ⟩⟨mⱼ|
         //
         // are gone. Now, *mathematically* they are *not needed* for the
         // integration rule optimization (see paper).
         //
         // However, their absence also causes some problems. In particular,
         // the Choleksy solver can't deal with this — this needs a rank-
         // revealing solver which can explicitly determine the non-singular
         // space.
         // So use eigh-solver by default in this case.
         m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Eigh, true);
//          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_LDLT);
      } else {
         // if not symmetrizing, we will get a wide spectrum of S eigenvalues;
         // the overlap matrix will have a large condition number for larger L.
         // However, none of the functions will be strictly linearly dependent,
         // and so the entire space must be handled accurately. Choleksy
         // decompositions with triangular Gaussian elimination tend to do this
         // rather well — better than the Smh version certainly.
         // So use Scd-solver by default in this case.
//          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Cholesky);
         m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_LDLT, true);
      }
//          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Cholesky, true); // FIXME: put the others back.
   }
}

void FTargetFnList_RawMonomials::Print(ct::FLog &io) const {
   std::stringstream str;
   PrintMonomialList(str, m_TargetFns, m_PrintLevel, "target monomial list");
   io.w << str.str();
}



FDenseMatrix FTargetFnList_RawMonomials::EvalOverlapMatrix() const
{
   if (1)
      return ::MakeOverlapMatrix(m_TargetFns);
   // ^-- definitely makes a prettier rule with
   //      make && aigg_64 '{point-group:O;degree:[4,step:+1,-1];initial-points:hex-grid{7;1;01c}; export:last.dat}'
   //     might need to check for singularities.
   //     And maybe max step too small? This one does not make the lmax=32 rule with the default 1e-4:
   //    make && aigg_64 '{point-group:O;degree:[4,step:+1,-1];initial-points:hex-grid{9;1;01c};max-step:1e-2; export:last.dat}'
   else
      return ::MakeOverlapMatrix_AxprAveraged(m_TargetFns, m_pPointGroup->pMonomialSymmetry());
}


FDenseVector FTargetFnList_RawMonomials::EvalTargetIntegrals() const
{
   FDenseVector
      rhs(size());
   for (size_t i = 0; i != m_TargetFns.size(); ++ i)
      rhs[i] = m_TargetFns[i].CalcUnitSphereIntegral();
   if (s_TestHalfMxm) {
      // test if `HalfMxm1` is the inverse of `HalfSolve1`, as it is supposed to. so
      // transform `t` (=rhs) vector to orth basis, and then back, and then
      // compare to original rhs vector.
      FDenseVector
         rhs_Orth = m_pScd->HalfSolve1(rhs),
         rhs_NonOrth = m_pScd->HalfMxm1(rhs_Orth);
      std::cout << fmt::format(" rmsd(rhs,rhs_Orth) = {:8.2e}   rmsd(rhs,rhs_NonOrth) = {:8.2e}\n", Rmsd(rhs - rhs_Orth), Rmsd(rhs - rhs_NonOrth));
   }
   if (m_pScd)
      // transform to quasi-orthogonal basis
      return m_pScd->HalfSolve1(rhs);
   else
      return rhs;
}

size_t FTargetFnList_RawMonomials::nFn() const {
   return m_TargetFns.size();
}

unsigned FTargetFnList_RawMonomials::Order() const {
   return m_lmax;
}

std::string FTargetFnList_RawMonomials::FnTypeName() const {
   if (g_UseGroupSymmetrizedFn)
      return "symmetrized monomials";
   else
      return "monomials";
}

std::string FTargetFnList_RawMonomials::FnTypeAnnotation() const {
   return "symmetry unique";
}


void FTargetFnList_RawMonomials::EvalK(FDenseMatrix &K, FPointArray const &SeedPoints)
{
   size_t
      nTargetFns = m_TargetFns.size(),
      nSeeds = SeedPoints.cols();
   cx_assert_rt(SeedPoints.rows() == 3);
   // K[i,q] = ∑_g f_i(R_g r⃗_q): (i: target fn index, q: seed index)
   cx_assert_rt(size_t(K.rows()) == nTargetFns && size_t(K.cols()) == nSeeds);
   FTransformedPointInfo
      // contains information about seed points transformed by current symmetry operation.
      // In particular, powers of its x,y,z coordinates.
      TrafodSeeds(SeedPoints, m_lmax);
   K.setZero();
   for (size_t iSymOp = 0; iSymOp != m_pSymmetryOps->size(); ++ iSymOp) {
      TrafodSeeds.TransformPoints((*m_pSymmetryOps)[iSymOp]);
      for (size_t iTargetFn = 0; iTargetFn != nTargetFns; ++ iTargetFn) {
         for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
            K(iTargetFn, iSeed) += m_TargetFns[iTargetFn].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
         }
      }
   }
   if (s_TestHalfMxm) {
      FDenseMatrix
         K_Orth = m_pScd->HalfSolve1(K),
         K_NonOrth= m_pScd->HalfMxm1(K_Orth);
      std::cout << fmt::format(" rmsd(K,K_Orth) = {:8.2e}   rmsd(K,K_NonOrth) = {:8.2e}\n", Rmsd(K - K_Orth), Rmsd(K - K_NonOrth));
   }
   if (m_pScd)
      // transform to quasi-orthogonal basis
      K = m_pScd->HalfSolve1(K).eval();
}

void FTargetFnList_RawMonomials::EvalDerivK(FDenseMatrix &dKw, FDenseVector const &weights, FDenseMatrix &dKe, FDenseVector const &err, FPointArray const &SeedPoints, FOrbitTypeRawPtrList const &OrbitTypes, FRotationGenerator const *pAi)
{
   // TO FIX: m_A (maybe make this totally global? or supply as args to this fn?), m_OrbitTypes.
   size_t
      nTargetFns = m_TargetFns.size(),
      nSeeds = SeedPoints.cols();
   size_t
      nVars = 3*nSeeds;

   cx_assert_rt(SeedPoints.rows() == 3);
   cx_assert_rt(OrbitTypes.size() == nSeeds && size_t(weights.size()) == nSeeds && size_t(err.size()) == nTargetFns);
   // dKw[i,p] = ∑_{q} dK[i,q]/dp w[q] — kernel derivative contracted to fixed weights.
   cx_assert_rt(size_t(dKw.rows()) == nTargetFns && size_t(dKw.cols()) == nVars);
   // dKe[q,p] = ∑_{i} dK[i,q]/dp e[i] — kernel derivative contracted to fixed target-fn error vectors e[i]
   cx_assert_rt(size_t(dKe.rows()) == nSeeds && size_t(dKe.cols()) == nVars);
   dKw.setZero();
   dKe.setZero();

   FDenseVector
      err_NonOrth;
   if (m_pScd)
      // transform input error vector back orthogonalized basis to from raw
      // monomial basis, so that we can contract it to raw K directly (as
      // opposed to computing and transforming the bigger object K, and
      // contracting it to the original orthogonal-basis error vector)
      err_NonOrth = m_pScd->HalfMxm1(err);
   else
      err_NonOrth = err;

   FTransformedPointInfo
      // contains information about seed points transformed by current symmetry operation.
      // In particular, powers of its x,y,z coordinates.
      TrafodSeeds(SeedPoints, m_lmax);
   for (size_t iSymOp = 0; iSymOp != m_pSymmetryOps->size(); ++ iSymOp) {
      FRotationMatrix const
         &Rg = (*m_pSymmetryOps)[iSymOp],
         Rg_Ax = Rg * pAi[0],
         Rg_Ay = Rg * pAi[1],
         Rg_Az = Rg * pAi[2];
      FPointArray
         Rg_Ax_pts = Rg_Ax * SeedPoints,
         Rg_Ay_pts = Rg_Ay * SeedPoints,
         Rg_Az_pts = Rg_Az * SeedPoints;
      TrafodSeeds.TransformPoints(Rg);

      // project out rotations which would lead out of the given orbit type.
      m_pPointGroup->ProjectRotationToOrbitType(Rg_Ax_pts, Rg_Ay_pts, Rg_Az_pts,
         TrafodSeeds.TransformedPoints, &OrbitTypes[0]);

      // evaluate contracted derivatives of K for current group operation.
      for (size_t iTargetFn = 0; iTargetFn != nTargetFns; ++ iTargetFn) {
         for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
            FPoint
               dMdXyz,
               Rg_Ax_sq = Rg_Ax_pts.col(iSeed),
               Rg_Ay_sq = Rg_Ay_pts.col(iSeed),
               Rg_Az_sq = Rg_Az_pts.col(iSeed);
            m_TargetFns[iTargetFn].EvalPoint(dMdXyz, TrafodSeeds.TransformedPointInfo(iSeed));
            FPoint
               dKiq_dp_alpha(0,0,0);
            for (size_t beta = 0; beta != 3; ++ beta) {
               dKiq_dp_alpha[0] += dMdXyz[beta] * Rg_Ax_sq[beta];
               dKiq_dp_alpha[1] += dMdXyz[beta] * Rg_Ay_sq[beta];
               dKiq_dp_alpha[2] += dMdXyz[beta] * Rg_Az_sq[beta];
            }
            size_t iVar = 3*iSeed;
            for (size_t beta = 0; beta != 3; ++ beta) {
               dKw(iTargetFn, iVar + beta) += dKiq_dp_alpha[beta] * weights[iSeed];
               dKe(iSeed, iVar + beta) += dKiq_dp_alpha[beta] * err_NonOrth[iTargetFn];
            }
         }
      }
   }
   if (m_pScd)
      // transform to quasi-orthogonal basis
      dKw = m_pScd->HalfSolve1(dKw).eval();
}


std::string FTargetFnList::MakeFnDesc(size_t iTargetFn) const
{
   char const *pFmt = "F[{:}]";
   return fmt::format(pFmt, iTargetFn);
}

std::string FTargetFnList_RawMonomials::MakeFnDesc(size_t iTargetFn) const
{
   FMonomial
      m = m_TargetFns.at(iTargetFn);
//    char const *pFmt = "[{},{},{}]";
//    char const *pFmt = "M(xyz)[{:2},{:2},{:2}]";
   char const *pFmt = "M[{:},{:},{:}]";
   return fmt::format(pFmt, m.ix, m.iy, m.iz);
}

std::string FTargetFnList_PreAveragedPolys::MakeFnDesc(size_t iTargetFn) const
{
   char const *pFmt = "{}{}[{:.6f}]";
   return fmt::format(pFmt,
      (bool(m_TargetFnFlags & FNSPACE_GroupAveraged)? "Pg" : ""),
      (bool(m_TargetFnFlags & FNSPACE_Orthogonal)? "E" : ""),
      this->m_Lambdas.at(iTargetFn));
}


FTargetFnList_PreAveragedPolys::FTargetFnList_PreAveragedPolys(size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FIoPolynomialList const &pIoPolys)
  : FTargetFnList(pPointGroup, GeneralOptions), m_lmax(lmax), m_fGroupImageFactor(1)
{
   // WARNING: this is a debug interface function!
   m_SymmetryOpsForK = pPointGroup->GetOps();
   this->m_TargetFnFlags = FNSPACE_Orthogonal | FNSPACE_GroupAveraged;
   FIoTargetFnSpace::FIoFnList
      pIoFns;
   pIoFns.reserve(pIoPolys.size());
   for (auto &&iop : pIoPolys) {
      pIoFns.push_back(iop.poly);
      this->m_Lambdas.push_back(iop.lambda);
   }
//    this->m_Lambdas.resize(pIoFns.size(), 1);
   InitBasisAndCoeffs(pIoFns);
}


FTargetFnList_PreAveragedPolys::FTargetFnList_PreAveragedPolys(size_t lmax, FPointGroup const *pPointGroup, FTargetFnListOptions const &GeneralOptions, FGridSearchOptions const *pGridSearchOptions)
   : FTargetFnList(pPointGroup, GeneralOptions), m_lmax(lmax), m_fGroupImageFactor(1)
{
   // Compute set of target functions which our integration rule is to evaluate
   // exactly, and compute their overlap matrix and a suitable decomposition of
   // it to solve conditioning equations with.
//    m_TargetFns = MakeSymmetryUniqueMonomialList(m_lmin, m_lmax, *m_pPointGroup, int(m_PrintLevel), pGridSearchOptions);
   // ^--- hmpf.. function screening turned off if no options object provided. Stupid.
   //      should fix some other day.
   std::string
      ErrorMsg("unknown error");
   FIoTargetFnSpacePtr
      pFnSet = TryReadTargetFnSpaceDataFile(*pPointGroup, lmax, &ErrorMsg);
   if (!pFnSet) {
      throw ETabulatedFnListError(ErrorMsg);
   }

   if (!bool(pFnSet->Flags & FNSPACE_GroupAveraged) || pPointGroup->Name(NAMETYPE_Symbol) != pFnSet->sPointGroup) {
      // input functions are not pre-averaged. At least not for the current
      // group. In this case we'll have to retain the regular averaging. Make a
      // copy of the original ops.
      m_SymmetryOpsForK = pPointGroup->GetOps();
   } else {
      // using pre-averaged functions. Retain a symmetry operator list, but just
      // put in a single identity matrix as symmetry transformation.
      m_SymmetryOpsForK.clear();
      m_SymmetryOpsForK.push_back(FRotationMatrix::Identity());
      // see notes on `m_fGroupImageFactor` at its point of declaration in the class
      m_fGroupImageFactor = pPointGroup->size();
      // UPDATE: hmpf... can't completely ignore the group, apparently. Even if all the
      // contributions of the target function images under the group would give the same
      // value, that will still result in a constant factor!
      // This leaves us with the following options: (a) scale K-contributions up,
      // (b) scale rhs vector down, (c) leave the symmetry operation list alone
      // (at least for testing)
      if (0) {
         m_SymmetryOpsForK = pPointGroup->GetOps();
         m_fGroupImageFactor = 1;
         // ^-- FIXME: remove these once things work.
      }
   }
//    this->m_pSymmetryOps = 0; // just to provoke a crash in case I missed adjusting one...
   this->m_Lambdas = pFnSet->pLambda;
   this->m_TargetFnFlags = pFnSet->Flags;

   InitBasisAndCoeffs(pFnSet->pFns);
}


void FTargetFnList_PreAveragedPolys::InitBasisAndCoeffs(FIoTargetFnSpace::FIoFnList const &pIoFns)
{
   m_TargetFns = pIoFns;

   // Make a `std::map` to collect all monomials occurring in any of the imported
   // polynomial subspace basis functions. We also use the map's predicate order
   // to define the order of our basis functions.
   typedef std::map<FMonomialN, size_t, FMonomialN_ExportOrderPred>
      FIoMonomialMap;
   FIoMonomialMap
      iMonomials;
   // ^--[side note: We use the export sorting predicate here (rather than the
   //     actual `FPolynomialN` one): the export version is slower, but this does
   //     not matter here, and I think the order it defines should be better
   //     suited to support numerical stability (in the export order, monomials
   //     with small integrals come first, all structurally closely related
   //     monomials come in groups together, etc)]
   for (FPolynomialNPtr const &pIoPoly : m_TargetFns) {
      for (FPolynomialN::value_type const &mf : *pIoPoly) {
         FMonomialN const &m = mf.first;  // mf: (monomial, factor)
         iMonomials[m] = 0;
      }
   }

   // now that we got the monomials, we'll just go through the entire
   // (ordered) list again, and assign indices based on the map's
   // sorting predicate.
   {
      size_t ii = 0;
      for (FIoMonomialMap::iterator it = iMonomials.begin(); it != iMonomials.end(); ++ it) {
         it->second = ii;
         ii += 1;
      }
   }
   // `iMonomials` now maps `FMonomialN` objects to the designated indices of basis
   // functions in our subspace basis.

   // Next, convert the `FMonomialN` objects to the more cumbersome `FMonomial`
   // objects from `TargetFunctions.h` (because those have some useful support
   // functions for evaluating sphere integrals and stuff)
   cx_assert_rt(m_BasisFns.empty());
   m_BasisFns.reserve(iMonomials.size());
   for (FIoMonomialMap::value_type const &mi : iMonomials) { // mi: (monomial, index)
      FMonomialN const
         &m = mi.first;
      FMonomial mconv(m[0], m[1], m[2]);
      if (s_CancelPrimitiveRenormFactors) {
         // cancel the auto-normalization factors? alternatively we can
         // reapply them to the coefficient array.
         mconv.m_Factor = 1;
         mconv.m_Degeneracy = 1;
      }
      assert(m_BasisFns.size() == mi.second);
      m_BasisFns.push_back(mconv);
//       cx_assert_rt(m_BasisFns.back().m_Factor == 1);
   }
// #ifdef INCLUDE_ABANDONED
   if (m_PrintLevel.AllQ()) {
      // debug code: consistency check of indexing
      for (FIoMonomialMap::value_type const &mi : iMonomials) { // mi: (monomial, index)
         FMonomialN const
            &m = mi.first;
         FMonomial const
            &mconv = m_BasisFns[mi.second];
         cx_assert_rt(mconv.ix == m[0] && mconv.iy == m[1] && mconv.iz == m[2]);
      }
      PrintMonomialList(std::cout, m_BasisFns, 2, "monomial basis functions in FTargetFnList_PreAveragedPolys");
      // ^- note: 2 is `iPrintLevel`. The function does not do anything unless `iPrintLevel >= 2`
   }
// #endif // INCLUDE_ABANDONED

   // collect polynomial coefficients.
   m_Coeffs = FDenseMatrix::Zero(m_BasisFns.size(), m_TargetFns.size());

   for (size_t iTargetFn = 0; iTargetFn != m_TargetFns.size(); ++ iTargetFn) {
      FPolynomialN const *pIoPoly = &*m_TargetFns[iTargetFn];
      for (FPolynomialN::value_type const &mf : *pIoPoly) { // mf: (monomial, factor)
         FMonomialN const
            &m = mf.first;
         // find basis function for this monomial (the basis functions
         // are also types of monomials, but `FMonomial`, not `FMonomialN`)
         size_t
            iBasisFn = iMonomials.at(m);
         FMonomial const
            *pBasisFn = &m_BasisFns[iBasisFn];
         assert(m_Coeffs(iBasisFn, iTargetFn) == 0);
         m_Coeffs(iBasisFn, iTargetFn) = mf.second / pBasisFn->m_Factor;
      }
   }


   {  S_TIME_SECTION(&io, m_PrintLevel.BasicQ(), "t-space covariant projectors");
      // need that for back-transforming target-fn basis error vector
      // to covariant basis, so that I can contract it to `K`.
      // (storing just `(S * C)` instead of `S` in addition to `C` is more
      // compact — can make a difference for large optimizations with
      // large floating point types. They need *lots* of memory.)
      m_CoeffsCov = MakeOverlapMatrix(m_BasisFns) * m_Coeffs;
   }
#ifdef INCLUDE_ABANDONED
//       FDenseMatrix
//          BasisOverlap;
//       BasisOverlap = MakeOverlapMatrix(m_BasisFns);
//    if (0) {
//       FLinearSolverPtr
//          pScd = MakeLinearSolver(m_BasisOverlap, LINSOLVE_PositiveSymmetric_Eigh, false);
//       //LINSOLVE_PositiveSymmetric_Cholesky
//       m_Coeffs = pScd->HalfSolve1(m_Coeffs).eval();
//       m_BasisOverlap.setIdentity();
//    }
//
//    if (1) {
//       // Actually... since m_BasisFns is a FMonomialList object, can't we just
//       // do this? ...
//       FDenseMatrix
// //          Sbb = MakeOverlapMatrix(m_BasisFns),
//          Stt = m_Coeffs.transpose() * (m_BasisOverlap * m_Coeffs).eval();
//       std::cout << fmt::format("[--- DEBUG: contracted overlap matrix {{lmax: {}; group: {}; nbas: {}; nfn: {}; rmsd(id): {:8.2e}}}:\n",
//                    this->m_lmax, this->m_pPointGroup->Name(NAMETYPE_Symbol), m_BasisFns.size(), m_TargetFns.size(), RmsdFromIdentity(Stt))
//          << Stt << "\n]---" <<std::endl;
//       // TODO: actually... if the input functions are not orthogonal, with
//       // this one could easily *make* them orthogonal. Just compute Schmidt
//       // trafo with Cholesky solver and absorb it into m_Coeffs...
//    }
//    if (0) {
//       // these two appear to agree just fine. at least atm (all the m_Factors etc are gone)
//       size_t nTargetFn = m_TargetFns.size();
//       FDenseMatrix
//          S1(nTargetFn, nTargetFn);
//       for (size_t i = 0; i < nTargetFn; ++ i)
//          for (size_t j = 0; j < nTargetFn; ++ j)
//             S1(i,j) = ComputeOverlapRaw(*m_TargetFns[i], *m_TargetFns[j]);
//       std::cout << fmt::format("[--- DEBUG: {{lmax: {}; group: {}; nbas: {}; nfn: {}}}:\n", this->m_lmax, this->m_pPointGroup->Name(NAMETYPE_Symbol), m_BasisFns.size(), m_TargetFns.size())
// //          << rhs_TargetFn
// //          << "\n]--- origin primitive t vector ---[\n" << rhs_BasisFn
// //          << "\n]--- C matrix ---[\n" << m_Coeffs
// //          << "\n]--- S matrix (C.T S1 C) ---[\n" << (m_Coeffs.transpose() * (m_BasisOverlap * m_Coeffs).eval())
//          << fmt::format("\n]--- S matrix (poly/direct) // rmsd(id) = {:8.2e}---[\n", RmsdFromIdentity(S1)) << S1
//          << "\n]---" << std::endl;
//    }
#endif // INCLUDE_ABANDONED
}


void FTargetFnList_PreAveragedPolys::WriteInfo(ct::FLog &io) const {
   if (m_PrintLevel.BasicQ()) {
      FDenseMatrix
         // compute target function set overlap matrix
         S_TargetFn = m_CoeffsCov.transpose() * m_Coeffs;
#ifdef INCLUDE_ABANDONED
//       io.WriteInfoExpf("Target basis non-orthogonality", double(RmsdFromIdentity(S_TargetFn)), "rmsd[S,I]");
//       io.WriteInfoExpf("Target basis non-orthogonality", double(RmsdFromIdentity(S_TargetFn)));
//       io.WriteInfoExpf("Target basis unitarity metric", double(RmsdFromIdentity(S_TargetFn)));
//       io.WriteInfoExpf("Target basis non-unitarity", double(RmsdFromIdentity(S_TargetFn)));
//       io.WriteInfoExpf("Target basis lost unitarity", double(RmsdFromIdentity(S_TargetFn)));
#endif // INCLUDE_ABANDONED
      std::string SummaryDesc;
      if (1)
         // to simplify grepping for 'unitarity loss...'
//          SummaryDesc = fmt::format("sym {}, lmx {}, tspa {} {})", m_pPointGroup->Name(NAMETYPE_Symbol), m_lmax, m_TargetFns.size(), m_BasisFns.size());
         SummaryDesc = fmt::format("⌯{0}, ℓ{1}, ℙ[{3} ↦ {2}]", m_pPointGroup->Name(NAMETYPE_Symbol), m_lmax, m_TargetFns.size(), m_BasisFns.size());
      io.WriteInfoExpf("Target basis unitarity loss", double(RmsdFromIdentity(S_TargetFn)), SummaryDesc);
   }
}

size_t FTargetFnList_PreAveragedPolys::nFn() const {
   return m_TargetFns.size();
}

unsigned FTargetFnList_PreAveragedPolys::Order() const {
   return m_lmax;
}

std::string FTargetFnList_PreAveragedPolys::FnTypeName() const {
//    return "group-averaged polynomials";
//    return fmt::format("expanded in {} monomials", m_BasisFns.size());
   return fmt::format("expanded over {} monomials", m_BasisFns.size());
}

std::string FTargetFnList_PreAveragedPolys::FnTypeAnnotation() const {
   return "exact minimal basis";
}


FScalar CalcMonomialUnitSphereIntegralRaw(size_t ix, size_t iy, size_t iz);

FDenseVector FTargetFnList_PreAveragedPolys::EvalTargetIntegrals() const
{
   FDenseVector
      rhs_BasisFn(m_BasisFns.size());
   for (size_t iBasisFn = 0; iBasisFn != m_BasisFns.size(); ++ iBasisFn) {
      rhs_BasisFn[iBasisFn] = m_BasisFns[iBasisFn].CalcUnitSphereIntegral();
#ifdef INCLUDE_ABANDONED
      if (0) { // debug code to deal with `FMonomial::m_Factor` being or not being unity
         FMonomial const &m = m_BasisFns[iBasisFn];
         cx_assert_rt(rhs_BasisFn[iBasisFn] == CalcMonomialUnitSphereIntegralRaw(m.ix, m.iy, m.iz));
      }
#endif // INCLUDE_ABANDONED
   }

   FDenseVector
      rhs_TargetFn = m_Coeffs.transpose() * rhs_BasisFn;
#ifdef INCLUDE_ABANDONED
//    if (0) {
//       std::cout << fmt::format("[--- DEBUG: contracted t vector {{lmax: {}; group: {}; nbas: {}; nfn: {}}}:\n", this->m_lmax, this->m_pPointGroup->Name(NAMETYPE_Symbol), m_BasisFns.size(), m_TargetFns.size())
//          << rhs_TargetFn
//          << "\n]--- origin primitive t vector ---[\n" << rhs_BasisFn
//          << "\n]--- C matrix ---[\n" << m_Coeffs
//          << "\n]--- S matrix ---[\n" << m_BasisOverlap
//          << "\n]---" << std::endl;
//    }
#endif // INCLUDE_ABANDONED
   return (1/m_fGroupImageFactor) * rhs_TargetFn;
}


void FTargetFnList_PreAveragedPolys::EvalK(FDenseMatrix &K_TargetFn, FPointArray const &SeedPoints)
{
   size_t
      nBasisFn = m_BasisFns.size(),
      nTargetFns = m_TargetFns.size(),
      nSeeds = SeedPoints.cols();
   assert(SeedPoints.rows() == 3);
   // K[i,q] = ∑_g f_i(R_g r⃗_q): (i: target fn index, q: seed index)
   cx_assert_rt(size_t(K_TargetFn.rows()) == nTargetFns && size_t(K_TargetFn.cols()) == nSeeds);

   FDenseMatrix
      K_BasisFn(nBasisFn, nSeeds);
   K_BasisFn.setZero();
   FTransformedPointInfo
      // contains information about seed points transformed by current symmetry operation.
      // In particular, powers of its x,y,z coordinates.
      TrafodSeeds(SeedPoints, m_lmax);
   for (size_t iSymOp = 0; iSymOp != m_SymmetryOpsForK.size(); ++ iSymOp) {
      TrafodSeeds.TransformPoints(m_SymmetryOpsForK[iSymOp]);
      for (size_t iBasisFn = 0; iBasisFn != nBasisFn; ++ iBasisFn) {
         for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
            K_BasisFn(iBasisFn, iSeed) += m_BasisFns[iBasisFn].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
         }
      }
   }
#ifdef INCLUDE_ABANDONED
//    for (size_t iSymOp = 0; iSymOp != m_SymmetryOpsForK.size(); ++ iSymOp) {
//       FRotationMatrix const &Rg = m_SymmetryOpsForK[iSymOp];
//       for (size_t iBasisFn = 0; iBasisFn != nBasisFn; ++ iBasisFn) {
//          for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//             K_BasisFn(iBasisFn, iSeed) += m_BasisFns[iBasisFn].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
//          }
//       }
//    }
#endif // INCLUDE_ABANDONED
   K_TargetFn = m_Coeffs.transpose() * K_BasisFn;
}

void FTargetFnList_PreAveragedPolys::EvalDerivK(FDenseMatrix &dKw_TargetFn, FDenseVector const &weights, FDenseMatrix &dKe, FDenseVector const &e_TargetFn, FPointArray const &SeedPoints, FOrbitTypeRawPtrList const &OrbitTypes, FRotationGenerator const *pAi)
{
   size_t
      nBasisFn = m_BasisFns.size(),
      nTargetFns = m_TargetFns.size(),
      nSeeds = SeedPoints.cols();
   size_t
      nVars = 3*nSeeds;

   assert(SeedPoints.rows() == 3);
   assert(OrbitTypes.size() == nSeeds && size_t(weights.size()) == nSeeds && size_t(e_TargetFn.size()) == nTargetFns);
   // dKw[i,p] = ∑_{q} dK[i,q]/dp w[q] — kernel derivative contracted to fixed weights.
   cx_assert_rt(size_t(dKw_TargetFn.rows()) == nTargetFns && size_t(dKw_TargetFn.cols()) == nVars);
   // dKe[q,p] = ∑_{i} dK[i,q]/dp e[i] — kernel derivative contracted to fixed target-fn error vectors e[i]
   cx_assert_rt(size_t(dKe.rows()) == nSeeds && size_t(dKe.cols()) == nVars);
   FDenseMatrix
      dKw_BasisFn(nBasisFn, nVars);
   dKw_BasisFn.setZero();
   dKe.setZero();

   FDenseVector
      // transform input error vector to monomial basis so that we
      // can contract it to raw K directly. Where does that come from? ->
      // Let et = e_TargetFn, eb = e_BasisFn, and C = m_Coeffs:
      //
      //          et = C^T eb;        // multiply by C from left
      //        C et = C C^T eb       // multiply by S from left
      //      S C et = S (C C^T) eb   // use C C^T = S^{-1} (at least in the subspace spanned by the C)
      //      S C et = eb             // cancel terms.
      //
      // Note: S * C = m_CoeffsCov = (BasisOverlap * m_Coeffs).eval()
      e_BasisFn = m_CoeffsCov * e_TargetFn;

   if (0) { // for cross-check of m_BasisOverlap construction
      FSvd SvdC(m_Coeffs.transpose());
      SvdC.Compress(g_ThrAlmostZero);
      e_BasisFn = SvdC.Solve(e_TargetFn);
   }

   FTransformedPointInfo
      // contains information about seed points transformed by current symmetry
      // operation. In particular, powers of its x,y,z coordinates.
      TrafodSeeds(SeedPoints, m_lmax);
   for (size_t iSymOp = 0; iSymOp != m_SymmetryOpsForK.size(); ++ iSymOp) {
      FRotationMatrix const
         &Rg = m_SymmetryOpsForK[iSymOp],
         Rg_Ax = Rg * pAi[0],
         Rg_Ay = Rg * pAi[1],
         Rg_Az = Rg * pAi[2];
      FPointArray
         Rg_Ax_pts = Rg_Ax * SeedPoints,
         Rg_Ay_pts = Rg_Ay * SeedPoints,
         Rg_Az_pts = Rg_Az * SeedPoints;
      TrafodSeeds.TransformPoints(Rg);

      // project out rotations which would lead out of the given orbit type.
      m_pPointGroup->ProjectRotationToOrbitType(Rg_Ax_pts, Rg_Ay_pts, Rg_Az_pts, TrafodSeeds.TransformedPoints, &OrbitTypes[0]);

      // evaluate contracted derivatives of `K` for current group operation.
      for (size_t iBasisFn = 0; iBasisFn != nBasisFn; ++ iBasisFn) {
         for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
            FPoint
               dMdXyz,
               Rg_Ax_sq = Rg_Ax_pts.col(iSeed),
               Rg_Ay_sq = Rg_Ay_pts.col(iSeed),
               Rg_Az_sq = Rg_Az_pts.col(iSeed);
            m_BasisFns[iBasisFn].EvalPoint(dMdXyz, TrafodSeeds.TransformedPointInfo(iSeed));
            FPoint
               dKiq_dp_alpha(0,0,0);
            for (size_t beta = 0; beta != 3; ++ beta) {
               dKiq_dp_alpha[0] += dMdXyz[beta] * Rg_Ax_sq[beta];
               dKiq_dp_alpha[1] += dMdXyz[beta] * Rg_Ay_sq[beta];
               dKiq_dp_alpha[2] += dMdXyz[beta] * Rg_Az_sq[beta];
            }
            size_t iVar = 3*iSeed;
            for (size_t beta = 0; beta != 3; ++ beta) {
               dKw_BasisFn(iBasisFn, iVar + beta) += dKiq_dp_alpha[beta] * weights[iSeed];
               dKe(iSeed, iVar + beta) += dKiq_dp_alpha[beta] * e_BasisFn[iBasisFn];
            }
         }
      }
   }
   // hmpf… all of them aligned as TN. In high performance 64bit-float
   // arithmethic that would be rather unfortunate regarding performance on
   // modern platforms. I guess for the emulated floats neither choice matters
   // much.
   dKw_TargetFn = m_Coeffs.transpose() * dKw_BasisFn;
}

