#ifdef INCLUDE_ABANDONED
// make && time OMP_NUM_THREADS=12 aigg '{ max-it: 10; point-group: icosahedral; degree: [10,step:+2,-1,-2]; max-step: 1e-2; initial-points: hex-grid{13;1;01c}; export: {file: test_run.dat; format: text; points: seeds; options: yes; version: yes; identity: yes} }' > /tmp/aigg.log
// ^- the call which blows up residuals, maybe. Not sure.


// make && aigg '{ point-group: icosahedral; degree: [4,step:+1,-1]; max-step: 1e-3; initial-points: hex-grid{9;1;01c} }'

// c++ -O3 -DNDEBUG -Wall -I$HOME/Programs/eigen OptSphereGrid.cpp SymmetryGroups.cpp format.cpp -o opt_sphere_grid && ./opt_sphere_grid
#endif // INCLUDE_ABANDONED
#include "Aigg.h"
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdlib.h> // for atoi

#ifdef INCLUDE_ABANDONED
// // --------------------------------------------------
// // FIXME: remove these and replace by FEighOrSvd from MathOps.h
// //        once program works again
// //        (also should unify MakeStep and MakeStepSymmetric while we're at it....)
// // --------------------------------------------------
// #include <Eigen/Cholesky>
// #include <Eigen/SVD>
// // #define EIGEN_SVD_ALGO JacobiSVD // straight up jacobi SVD
// #define EIGEN_SVD_ALGO BDCSVD // bidirectional divide&conquer SVD (faster)
// // --------------------------------------------------
#endif // INCLUDE_ABANDONED

#include "SearchSphereGrid.h"


static ptrdiff_t g_nDampedSigSteps = 10;


#ifdef INCLUDE_ABANDONED
// // #define EIGEN_SVD_ALGO JacobiSVD // straight up jacobi SVD
// #define EIGEN_SVD_ALGO BDCSVD // bidirectional divide&conquer SVD (faster)
// 
// // hm... this one NaNs with BDCSVD... (NaNs in U and/or V in SVD in update):
// // make && aigg '{point-group: tetrahedral; degree: [1,step:+1,-1,+1,-1]; max-step: 0.005; initial-points: hex-grid{6;0;01c}; export: {file: runs/tetrahedral_6_0.dat; format: text; points: seeds; options: yes; version: yes; identity: yes} }'
// // but maybe the problem really lies in J matrix computation?
// 
// 
// // FScalar Rmsd(FScalar const *p, size_t n) {
// //    if (n == 0)
// //       return FScalar(0);
// //    FScalar
// //       f = 0;
// //    for (size_t i = 0; i != n; ++ i)
// //       f += p[i]*p[i];
// //    return sqrt(f/FScalar(n));
// // }
// 
// 
// FScalar Rmsd(FDenseVector const &v) {
//    if (v.size() == 0)
//       return 0.;
//    return sqrt(v.squaredNorm() / FScalar(v.size()));
// }
// 
// 
// FScalar Rmsd(FDenseVector const &v, size_t nDegreesOfFreedom) {
//    if (v.size() == 0)
//       return 0.;
//    if (nDegreesOfFreedom == 0)
//       nDegreesOfFreedom = 1;
//    return sqrt(v.squaredNorm() / FScalar(nDegreesOfFreedom));
// }
// 
// 
// // some trust-region Newton stuff c/p'd from MicroScf's geometry optimizer...
// static FScalar NtRmsd(FScalar const *p, size_t n, size_t nDegreesOfFreedom) {
//    if (n == 0)
//       return FScalar(0);
//    FScalar
//       f = 0;
//    for (size_t i = 0; i != n; ++ i)
//       f += p[i]*p[i];
//    if (nDegreesOfFreedom == 0)
//       nDegreesOfFreedom = 1;
//    return sqrt(f/FScalar(nDegreesOfFreedom));
// }
//    
// 
// struct FStepContext {
//    FStepContext(FScalar const *pEw_, FScalar const *gi_, size_t n_, FScalar fTargetStep_, size_t nDegreesOfFreedom_, FScalar fThrSig_)
//       : pEw(pEw_), gi(gi_), n(n_), fTargetStep(fTargetStep_), nDegreesOfFreedom(nDegreesOfFreedom_), fThrSig(fThrSig_)
//    {}
// 
//    void MakeStep(FScalar *gi_InvEw);
// protected:
//    FScalar const
//       *pEw, *gi;
//    size_t
//       n;
//    FScalar
//       fTargetStep;
//    size_t
//       nDegreesOfFreedom;
//    FScalar
//       fThrSig;
// 
//    // returns step rmsd
//    FScalar MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping);
//    void BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal);
// };
// 
// 
// void FStepContext::MakeStep(FScalar *gi_InvEw)
// {
//    FScalar
//       fRightPos = 0.,
//       fRightVal = MakeStepForDamping(gi_InvEw, fRightPos) - fTargetStep;
//    if (fRightVal < 0.)
//       // step with zero damping is already small enough.
//       return;
// 
//    FScalar
//       MaxEw = pEw[n-1];
//    if (MaxEw < pEw[0])
//       MaxEw = pEw[0]; // in case of singular values as input (sorted descending) instead of eigenvalues (sorted ascending)
//    FScalar
//       fLeftPos = 1e4*MaxEw*FScalar(n),
//       fLeftVal = MakeStepForDamping(gi_InvEw, fLeftPos) - fTargetStep;
//    assert(fLeftVal < 0.);
// 
//    return BisectR(gi_InvEw, fLeftPos, fLeftVal, fRightPos, fRightVal);
// }
// 
// 
// void FStepContext::BisectR(FScalar *gi_InvEw, FScalar fLeftPos, FScalar fLeftVal, FScalar fRightPos, FScalar fRightVal)
// {
//    for (size_t iIt = 0; iIt < 0x10000; ++ iIt) {
//       FScalar
//          fCenPos = 0.5 * (fLeftPos + fRightPos),
//          fCenVal = MakeStepForDamping(gi_InvEw, fCenPos) - fTargetStep;
//       if (abs(fCenVal - fTargetStep) < 1e-8 || abs(fCenVal) < 1e-10)
//          return;
//       if (fCenVal < 0.) {
//          fLeftPos = fCenPos;
//          fLeftVal = fCenVal;
//       } else {
//          fRightPos = fCenPos;
//          fRightVal = fCenVal;
//       }
//    }
//    // FIXME: ...this thing has no error reporting functionality...
// }
// 
// 
// FScalar FStepContext::MakeStepForDamping(FScalar *gi_InvEw, FScalar Damping)
// {
//    for (size_t i = 0; i < n; ++ i) {
//       if (abs(pEw[i]) > fThrSig)
//          gi_InvEw[i] = -gi[i]/(pEw[i] + Damping);
//       else
//          gi_InvEw[i] = 0.;
//    }
//    return NtRmsd(gi_InvEw, n, nDegreesOfFreedom);
// }
// 
// 
// static ptrdiff_t g_nDampedSigSteps = 10;
// 
// static void PrintVector(FDenseVector const &ew, std::string Title)
// {
//    std::cout << fmt::format("      {}", Title) << std::endl;
//    for (size_t i = 0; i != size_t(ew.size()); ++ i) {
//       std::cout << fmt::format("    {:>8d}  {:24.12f}  (={:8.2e})\n", i+1, double(ew[i]), double(ew[i]));
//    }
// }
// 
// static void PrintMatrix(FDenseMatrix const &M, std::string Title)
// {
//    std::cout << fmt::format("      {}", Title) << std::endl;
//    for (size_t iRow = 0; iRow != size_t(M.rows()); ++ iRow) {
//       std::cout << fmt::format("    {:>8d}", iRow+1);
//       for (size_t iCol = 0; iCol != size_t(M.cols()); ++ iCol)
//          std::cout << fmt::format("  {:24.12f}", double(M(iRow,iCol)));
//       io.WriteLine();
//    }
// }
#endif // INCLUDE_ABANDONED


FScalar MakeMaxStepWidth(FDenseVector const &g, size_t nFreeParameters, size_t iIt, size_t nPointsTotal, FScalar MaxStep_DynamicControlFactor, FOptimizeControlOptions const &StepOptions)
{
   FScalar
      fMaxStep_Scale = 1.;
   if (StepOptions.MaxStep_AdjustByNpts) {
      FScalar
         // square root of this should be a measure for the average distance
         // between points. I think theoretically that should be a reasonable
         // scale of reference for relative step widths.
         fAverageAreaPerPoint = 4*X_PI / FScalar(nPointsTotal);
//       fMaxStep_Scale *= 5.*5. * sqrt(fAverageAreaPerPoint);
      // ^- the 25. is for compatibility with earlier threshold values, and probably should not be there.
      fMaxStep_Scale *= 5. * sqrt(fAverageAreaPerPoint);
   }
   FScalar
      CurMaxStep;
   // note: at this point, 'g' would be the error vector, with one entry for each target function.
   // the issue is that this representation is not invariant to the number of functions (in case
   // there are redundant subspaces in the spanned function set), so we supply the number of 
   // actually optimized parameters (nFreeParameters) as N for the Rmsd.
   CurMaxStep = std::min(fMaxStep_Scale*StepOptions.fMaxStep_Max, StepOptions.fMaxStep_ResidualFactor * Rmsd(g,nFreeParameters));
   CurMaxStep = std::max(fMaxStep_Scale*StepOptions.fMaxStep_Min, CurMaxStep);
//    CurMaxStep = std::max(fMaxStep_Scale*StepOptions.fMaxStep_Min, std::min(fMaxStep_Scale*StepOptions.fMaxStep_Max, StepOptions.fMaxStep_ResidualFactor * Rmsd(g)));
   
//    CurMaxStep *= pow(0.98, iIt);
//    CurMaxStep *= pow(1.01, iIt);
//    CurMaxStep *= FScalar(1) + FScalar(1)*cos(2*X_PI*FScalar(iIt)/FScalar(5.4));
//    CurMaxStep *= pow(0.998, iIt);
//    CurMaxStep *= pow(FScalar(1.01), FScalar(iIt));
   // ^- hm.. starting with a small step and *increasing* it may not be such a horrible
   //    idea. could be worth an own option.
   // maybe even add random step length amplitudes?
   // (or finally just do the stupid DIIS?)
   CurMaxStep *= pow(StepOptions.fMaxStep_IterFactor, FScalar(iIt));
   CurMaxStep *= MaxStep_DynamicControlFactor;
   
   return CurMaxStep;
}

#ifdef INCLUDE_ABANDONED
// // - FIXME: clean up this stuff before release.
// // - FIXME: The ThrSigB stuff is *still* broken. At least with symmetrized
// //   functions and ico grids it definitely got the limit wrong at several
// //   occasions. I actually think this could be one of the reasons the code
// //   currently does not work for large order rules.
// // - We probably want program options to control this...
// static void MakeStep(FDenseVector &dx, FDenseMatrix const &J, FDenseVector const &g, FScalar &CurMaxStep, size_t iIt, FScalar &fMinRes, size_t nFreeParameters, size_t nPointsTotal, FGridSearchOptions const &Options)
// {
//    Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
//       svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    // ^- TODO: move SVD and eigh calculation code to MathOps.cpp? Very slow to
//    //    compile, and in that file there is already another SVD/Eigh.
//    //    (I mean just making an interface class containing U sig V matrices,
//    //    and keeping the actual *compute* routines in MathOps.cpp; not
//    //    actually moving this function which deals with the resulting
//    //    matrix decompositions)
//    FDenseVector
//       ew = svd.singularValues();
//    FDenseVector
//       gi = svd.matrixU().transpose() * g,
//       gi_InvEw = 0.*ew;
//    assert(gi.size() == ew.size());
//    FScalar
//       fThrSigB;
//    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
//       // take exactly as many singular values as we have free parameters in the optimization.
//       fThrSigB = ew[nFreeParameters-1]*FScalar(.99);
//    else
//       // not enough parameters to fix the number like this. take as many as we can.
//       fThrSigB = 1e-14Q;
//    
//    if (g_UseGroupSymmetrizedFn) fThrSigB = 1e-20Q; // FIXME: Test code! That should most likely not be here.
// //    FScalar
// //       fThrSigB = 1e-8Q;
// //    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
// //       // take exactly as many singular values as we have free parameters in the optimization.
// //       fThrSigB = ew[nFreeParameters-1]*FScalar(.99);
// //    // in the starting iterations (with very small monomial set sizes) we may not have enough
// //    // things to fix all the possible free parameters...
// // //    fThrSigB = std::max(fThrSigB, FScalar(1e-12Q));
// // //    fThrSigB = std::max(fThrSigB, FScalar(1e-8Q));
// //    fThrSigB = std::max(fThrSigB, FScalar(1e-12Q));
//    // ^- with 1e-12Q the grid Number of monomials               457  (symmetry unique)
//    //    for  make && ./opt_sphere_grid '{point-group: icosahedral; degree: [14,step:+2,-1]; max-step: 0.005; initial-points: hex-grid{8;1;01c}; export: {file: /tmp/grid0.dat; format: text; points: seeds}}' > /tmp/grid-run.txt
//    //    breaks. Hmpf. Maybe I should just increase all the small but non-zero ews to a certain target value, instead of completely ignoring them?
// 
//    bool PrintDetails = (Options.iPrintLevel >= 2);
//    if (PrintDetails) {
//       PrintVector(g, fmt::format("GRADIENT // ORIGINAL BASIS"));
//       PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
//       PrintVector(gi, fmt::format("GRADIENT // SINGULAR BASIS"));
// //       std::cout << fmt::format("      SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)) << std::endl;
// //       for (size_t i = 0; i != size_t(ew.size()); ++ i) {
// //          std::cout << fmt::format("    {:>8d}  {:24.12f}  (={:8.2e})\n", i+1, double(ew[i]), double(ew[i]));
// //       }
// //       FScalar fThrSigB1 = 1e-3Q;
// //       for (size_t i = 0; i != size_t(ew.size()); ++ i)
// //          if (ew[i] < fThrSigB1)
// //             ew[i] = fThrSigB1;
//          
// //       FScalar r1 = gi.norm()/sqrt(FScalar(gi.size()));
// //       for (size_t i = 0; i != size_t(ew.size()); ++ i)
// //          if (ew[i] < 1e-12Q) {
// //             ew[i] = 0;
// //          } else {
// //             FScalar
// //                f = FScalar(10*gi.size()),
// //                si = abs(gi[i]/ew[i])/f;
// //             if (si > r1)
// //                ew[i] *= si/r1;
// //                
// //          }
//    }
// //    PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
// //    std::cout << fmt::format("   sigma values: sum={:12.6f}   sum(s^2)={:12.6f}    sum(sqrt(s))={:12.6f}\n", double(ew.array().sum()), double(ew.array().pow(2).sum()), double(ew.array().sqrt().sum()));
// //    double fMaxStepScale = 
// //    CurMaxStep = std::max(Options.fMaxStep_Min, std::min(Options.fMaxStep_Max, Options.fMaxStep_ResidualFactor * Rmsd(g)));
//    CurMaxStep = MakeMaxStepWidth(g, nFreeParameters, nPointsTotal, Options);
// 
//    FStepContext(&ew[0], &gi[0], ew.size(), CurMaxStep, nFreeParameters, fThrSigB).MakeStep(&gi_InvEw[0]);
// 
//    if (PrintDetails) PrintVector(gi_InvEw, fmt::format("UPDATE // SINGULAR BASIS"));
// 
//    dx = svd.matrixV() * gi_InvEw;
//    if (PrintDetails) PrintVector(dx, fmt::format("UPDATE // ORIGINAL BASIS"));
//    
//    (void)fMinRes; // suppress unused warning
// }
// 
// 
// static void MakeStepSymmetric(FDenseVector &dx, FDenseMatrix const &J, FDenseVector const &g, FScalar &CurMaxStep, size_t iIt, FScalar &fMinRes, size_t nFreeParameters, size_t nPointsTotal, FGridSearchOptions const &Options)
// {
//    Eigen::SelfAdjointEigenSolver<FDenseMatrix>
//       sd(J); // compute spectral decomposition
//    // ^- TODO: move eigh decomposition code to MathOps.cpp? Very slow to compile.
//    FDenseVector
//       ew = sd.eigenvalues();
//    FDenseMatrix const
//       &ev = sd.eigenvectors(); 
//    FDenseVector
//       gi = ev.transpose() * g,
//       gi_InvEw = 0.*ew;
//    assert(gi.size() == ew.size());
// //    FScalar
// //       fThrSigB = 1e-14Q;
// //    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
// //       // take exactly as many singular values as we have free parameters in the optimization.
// //       // (note: no -1; if we have 1 free parameter, we need to take the last eigenvalue)
// //       fThrSigB = ew[ew.size()-nFreeParameters]*FScalar(.99);
// //    fThrSigB = std::max(fThrSigB, FScalar(1e-14Q));
//    FScalar fThrSigB = 0.;
//    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
//       // take exactly as many singular values as we have free parameters in the optimization.
//       fThrSigB = ew[nFreeParameters-1]*FScalar(.99);
// //    else
// //       // not enough parameters to fix the number like this? take as many as we can.
// //       fThrSigB = FScalar(1e-20Q);
//    fThrSigB = std::max(FScalar(1e-14Q), fThrSigB);
// //    fThrSigB = std::max(FScalar(1e-18Q), fThrSigB);
//       
// //    bool PrintDetails = false;
//    bool PrintDetails = (Options.iPrintLevel >= 2);
//    if (PrintDetails) {
//       PrintVector(g, fmt::format("GRADIENT // ORIGINAL BASIS"));
//       PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
//       PrintVector(gi, fmt::format("GRADIENT // SINGULAR BASIS"));
//    }
// //    PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
// //    std::cout << fmt::format("   nFreeParameters = {}   sigma values: sum={:12.6f}   sum(s^2)={:12.6f}    sum(sqrt(s))={:12.6f}\n", nFreeParameters, double(ew.array().sum()), double(ew.array().pow(2).sum()), double(ew.array().abs().sqrt().sum()));
// //    CurMaxStep = std::max(Options.fMaxStep_Min, std::min(Options.fMaxStep_Max, Options.fMaxStep_ResidualFactor * Rmsd(g)));
//    CurMaxStep = MakeMaxStepWidth(g, nFreeParameters, nPointsTotal, Options);
// 
//    FStepContext(&ew[0], &gi[0], ew.size(), CurMaxStep, nFreeParameters, fThrSigB).MakeStep(&gi_InvEw[0]);
// //    for (size_t i = 0; i != size_t(ew.size()); ++i) {
// //       if (abs(ew[i]) > fThrSigB)
// //          gi_InvEw[i] = -gi[i] / ew[i];
// //       else
// //          gi_InvEw[i] = 0;
// //    }
// //    CurMaxStep = 0;
// 
//    if (PrintDetails) PrintVector(gi_InvEw, fmt::format("UPDATE // SINGULAR BASIS"));
// 
//    dx = ev * gi_InvEw;
//    if (PrintDetails) PrintVector(dx, fmt::format("UPDATE // ORIGINAL BASIS"));
//    
//    (void)fMinRes; // suppress unused warning
// }
#endif // INCLUDE_ABANDONED


// - FIXME: clean up this stuff before release.
// - FIXME: The ThrSigB stuff is *still* broken. At least with symmetrized
//   functions and ico grids it definitely got the limit wrong at several
//   occasions. I actually think this could be one of the reasons the code
//   currently does not work for large order rules.
// - We probably want program options to control this...
static void MakeStep(FDenseVector &dx, FDenseMatrix const &J, FDenseVector const &g, FEighOrSvd::FDecompType SolverType, size_t nFreeParameters, FScalar &CurMaxStep, size_t iIt, size_t nPointsTotal, FScalar MaxStep_DynamicControlFactor, FGridSearchOptions const &Options)
{
   // compute spectral decomp or SVD of Jacobian
   FEighOrSvdPtr
      pJd = new FEighOrSvd(J, SolverType);
   FDenseVector
      // Collect eigenvalues/singular values (both ordered large to small).
      // Make a copy: need them in linear continuous memory for FStepContext.
      vals = pJd->vals();
   FDenseVector
      // not: system is   g \approx J * dx; so g gets multiplied by U first.
      gi = pJd->Ut() * g,
      gi_InvEw = FDenseVector::Zero(vals.size());
   assert(gi.size() == vals.size());
   FScalar fThrSigB = 0.;
   if (0) {
      if (nFreeParameters != 0 && size_t(vals.size()) > nFreeParameters)
         // take exactly as many singular values as we have free parameters in the optimization.
         fThrSigB = vals[nFreeParameters-1]*FScalar(.99);
   //    fThrSigB = std::max(FScalar(1e-14Q), fThrSigB);
   //    fThrSigB = std::max(FScalar(1e-18Q), fThrSigB);
      fThrSigB = std::max(FScalar(1e2)*g_ThrAlmostZero, fThrSigB);
   } else {
      // the set above used to be my default setup. For the moment I will try it
      // with a global threshold, hoping that the regularization style change
      // might have fixed the issues seen before in numerically complicated
      // high-order cases. TODO: reinvestigate if this indeed worked. Initial
      // trials looked like either choice did not make much of a difference.
      // (at least with the SVD renormalization style). And, for example, this:
      // make && aigg_64 '{point-group: I; order:[4,step:+1,-1]; initial-points: hex-grid{7;5;01c}; export: last.dat}'
      // comes out very nice-looking in both cases, making a:
      // 1092  (lmax = 56, wtsp = 1.670, res = 2.42e-10, mxe = 9.33e-11, edof = 0)
      // grid which look very regular in the visualization. Even the
      // make && aigg_64 '{point-group: I; order:[4,step:+1,-1]; initial-points: hex-grid{9;5;01c}; export: last.dat}'
      // grid seems to work and look fine. Makes
      // 1512  (lmax = 66, wtsp = 1.846, res = 1.12e-09, mxe = 8.57e-10, edof = 0)
      fThrSigB = g_ThrAlmostZero;
//       fThrSigB = g_ThrVerySmall;
      // UPDATE: I think when I wrote this, the core never actually went through
      // an actually asymmetric update, because it worked on J.T S^{-1} J....
   }

//    if (g_UseGroupSymmetrizedFn)
//       fThrSigB = 1e-20Q; // FIXME: Test code! That should most likely not be here.
      
   bool PrintDetails = (Options.PrintOptions[PRINT_Step].MostQ());
   if (PrintDetails) {
      // count number of singular values (and therefore update directions) which
      // FStepContext would use for informative purposes. Note: at the time of
      // this writing, this abs(vals[i]) > ThrSigB is the exact test used in
      // FStepContext to determine whether a direction should be included..
      size_t
         nEffectiveRank = 0;
      for (size_t i = 0; i < size_t(vals.size()); ++ i)
         if (abs(vals[i]) > fThrSigB)
            nEffectiveRank += 1;
      std::cout << fmt::format("      UPDATE // nFreeParameters = {}  nEffectiveRank = {}  nEw = {}  SolverType = {}\n", nFreeParameters, nEffectiveRank, vals.size(), int(SolverType));
      PrintVector(g, fmt::format("GRADIENT // ORIGINAL BASIS"));
      PrintVector(vals, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
      PrintVector(gi, fmt::format("GRADIENT // SINGULAR BASIS"));
   }
   CurMaxStep = MakeMaxStepWidth(g, nFreeParameters, iIt, nPointsTotal, MaxStep_DynamicControlFactor, *Options.pStepOptions);

#ifdef INCLUDE_ABANDONED
//    FStepContext(&vals[0], &gi[0], vals.size(), CurMaxStep, nFreeParameters, fThrSigB).MakeStep(&gi_InvEw[0]);
#endif // INCLUDE_ABANDONED
//    FEwType
// //       RegularizationStyle = EWTYPE_Eigh;
//       RegularizationStyle = EWTYPE_Svd;
//    // ^-- hm... using the same step lengths, the eigh thing (formally for the
//    //     minimum of symmetric-positive quadratic form) does appear to give
//    //     faster convergence also in the svd case for non-definite problems. But
//    //     from a theoretical perspective, it is on the shady side. I leave it
//    //     with the Tikhonov regularization (EWTYPE_Svd) for now.
   FEwType
      RegularizationStyle = EWTYPE_Eigh;
   if (SolverType == FEighOrSvd::DECOMP_Svd)
      RegularizationStyle = EWTYPE_Svd;
   ComputeTrustRegionStep<FScalar>(&gi_InvEw[0], &vals[0], RegularizationStyle, &gi[0],
      vals.size(), CurMaxStep, nFreeParameters, fThrSigB);

   if (PrintDetails) PrintVector(gi_InvEw, fmt::format("UPDATE // SINGULAR BASIS"));

   dx = pJd->V() * gi_InvEw;
   if (PrintDetails) PrintVector(dx, fmt::format("UPDATE // ORIGINAL BASIS"));
   
//    (void)fMinRes; // suppress unused warning
   (void)iIt; // suppress unused warning
}

#ifdef INCLUDE_ABANDONED
// static void MakeStep(FDenseVector &dx, FDenseMatrix const &J, FDenseVector const &g, FEighOrSvd::FDecompType SolverType, size_t nFreeParameters, FScalar &CurMaxStep, size_t iIt, size_t nPointsTotal, FGridSearchOptions const &Options)
// {
//    Eigen::SelfAdjointEigenSolver<FDenseMatrix>
//       sd(J); // compute spectral decomposition
//    // ^- TODO: move eigh decomposition code to MathOps.cpp? Very slow to compile.
//    FDenseVector
//       ew = sd.eigenvalues();
//    FDenseMatrix const
//       &ev = sd.eigenvectors(); 
//    FDenseVector
//       gi = ev.transpose() * g,
//       gi_InvEw = 0.*ew;
//    assert(gi.size() == ew.size());
// //    FScalar
// //       fThrSigB = 1e-14Q;
// //    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
// //       // take exactly as many singular values as we have free parameters in the optimization.
// //       // (note: no -1; if we have 1 free parameter, we need to take the last eigenvalue)
// //       fThrSigB = ew[ew.size()-nFreeParameters]*FScalar(.99);
// //    fThrSigB = std::max(fThrSigB, FScalar(1e-14Q));
//    FScalar fThrSigB = 0.;
//    if (nFreeParameters != 0 && size_t(ew.size()) > nFreeParameters)
//       // take exactly as many singular values as we have free parameters in the optimization.
//       fThrSigB = ew[nFreeParameters-1]*FScalar(.99);
// //    else
// //       // not enough parameters to fix the number like this? take as many as we can.
// //       fThrSigB = FScalar(1e-20Q);
//    fThrSigB = std::max(FScalar(1e-14Q), fThrSigB);
// //    fThrSigB = std::max(FScalar(1e-18Q), fThrSigB);
//       
// //    bool PrintDetails = false;
//    bool PrintDetails = (Options.iPrintLevel >= 2);
//    if (PrintDetails) {
//       PrintVector(g, fmt::format("GRADIENT // ORIGINAL BASIS"));
//       PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
//       PrintVector(gi, fmt::format("GRADIENT // SINGULAR BASIS"));
//    }
// //    PrintVector(ew, fmt::format("SIGMA VALUES  (THR={:8.2e})", double(fThrSigB)));
// //    std::cout << fmt::format("   nFreeParameters = {}   sigma values: sum={:12.6f}   sum(s^2)={:12.6f}    sum(sqrt(s))={:12.6f}\n", nFreeParameters, double(ew.array().sum()), double(ew.array().pow(2).sum()), double(ew.array().abs().sqrt().sum()));
// //    CurMaxStep = std::max(Options.fMaxStep_Min, std::min(Options.fMaxStep_Max, Options.fMaxStep_ResidualFactor * Rmsd(g)));
//    CurMaxStep = MakeMaxStepWidth(g, nFreeParameters, nPointsTotal, Options);
// 
//    FStepContext(&ew[0], &gi[0], ew.size(), CurMaxStep, nFreeParameters, fThrSigB).MakeStep(&gi_InvEw[0]);
// //    for (size_t i = 0; i != size_t(ew.size()); ++i) {
// //       if (abs(ew[i]) > fThrSigB)
// //          gi_InvEw[i] = -gi[i] / ew[i];
// //       else
// //          gi_InvEw[i] = 0;
// //    }
// //    CurMaxStep = 0;
// 
//    if (PrintDetails) PrintVector(gi_InvEw, fmt::format("UPDATE // SINGULAR BASIS"));
// 
//    dx = ev * gi_InvEw;
//    if (PrintDetails) PrintVector(dx, fmt::format("UPDATE // ORIGINAL BASIS"));
// }
#endif // INCLUDE_ABANDONED


FExportedGridStats::FExportedGridStats()
{
//    m_LastExportedWtsp = -1e99;
//    m_LastExportedResidual = 1e99;
//    m_LastExportedLmax = -1;
   wtsp = -1e99;
   res = 1e99;
   lmax = -1;
}


bool FExportedGridStats::IsWorseThan(FExportedGridStats const &other) const
{
   // note: we should only ever get here if both grids are converged (i.e., in
   // principle acceptable) or if 'other' is the default-initialized
   // FExportedGridStats object with very bad stats. For this reason there is no
   // "either grid not converged?"-residual comparison here.
   if (other.lmax < lmax)
      return false;
   if (lmax < other.lmax)
      return true;
   // if both grids exactly integrate the same lmax, and the earlier grid had
   // negative weights, then select the current, newer one (even if it, too, has
   // negative weights).
   if (other.wtsp < 0.)
      return false;
   if (wtsp >= 0. && other.wtsp >= 0.) {
      // both grids have positive weights. In this case, take the one with lower
      // weight spread. But add a tiny bit of leeway to bias towards the current,
      // newer grid.
      FScalar f = 1.001;
      if (f*other.wtsp < wtsp)
         return true;
      if (wtsp < f*other.wtsp)
         return false;      
   }
   // otherwise: select the newer grid.
   return false;
}

bool FExportedGridStats::IsEquivalent(FExportedGridStats const &other) const
{
   if (this->IsWorseThan(other))
      return false;
   if (other.IsWorseThan(*this))
      return false;
   return true;
}


FSphereGridOptContext::FSphereGridOptContext(FPointArray &SeedPoints, size_t lmin, size_t lmax, FGridSearchOptions const &Options)
   : m_Options(Options), m_SeedPoints(SeedPoints), m_lmin(lmin), m_lmax(lmax),
     m_PointGroup(*Options.pPointGroup),
     m_SymmetryOps(m_PointGroup.GetOps())
{
//    Run();
   m_fMinRes = 1e99;

//    m_LastExportedWtsp = -1e99;
//    m_LastExportedResidual = 1e99;
//    m_LastExportedLmax = -1;
   m_LastExportedGridStats = FExportedGridStats(); // inits to 'this is a clearly very bad grid'
   m_DidExportSomething = false;
}


FScalar FSphereGridOptContext::Run()
{
   m_tThisGridTotal.Reset();
   InitSeedsAndOrbits();
   InitTargetFns();
   io.WriteLine();
   io.WriteTiming("initalization", double(m_tThisGridTotal));
   RunIterations();
   return FinishAndExportData();
}


// auxiliary structure for sorting seed points of an integration rule.
// Mainly for decorative purposes (sorting the seeds makes it easier to compare if
// two rules with the same number of points are indeed equal)
struct FSeedSortEntry {
   FPoint
      p;
   FOrbitType const
      *pOrbitType;
   FSeedSortEntry() {}
   FSeedSortEntry(FPoint p_, FOrbitType const *pOrbitType_) : p(p_), pOrbitType(pOrbitType_) {}
   bool operator < (const FSeedSortEntry &other) const {
      // put special orbits first.
      if (other.pOrbitType->iType < this->pOrbitType->iType) return true;
      if (this->pOrbitType->iType < other.pOrbitType->iType) return false;
      if (this->p[0] < other.p[0]) return true;
      if (other.p[0] < this->p[0]) return false;
      if (this->p[1] < other.p[1]) return true;
      if (other.p[1] < this->p[1]) return false;
      return this->p[2] < other.p[2];
   }
};


static void PrepareInitialSeeds(FPointArray &SeedPoints, FPointGroup const &PointGroup)
{
   // make sure input points are normalized
   for (size_t iSeed = 0; iSeed != (size_t)SeedPoints.cols(); ++ iSeed) {
      SeedPoints.col(iSeed).normalize();
   }

   // sort them to simplify comparing different integration rules with the same number of points
   // TODO: should probably do this *after* projecting them to their target
   // orbit types...
   std::vector<FSeedSortEntry>
      SeedSortEntries;
   SeedSortEntries.reserve(size_t(SeedPoints.cols()));
   for (ptrdiff_t i = 0; i < SeedPoints.cols(); ++ i) {
      FPoint const &p = SeedPoints.col(i);
      SeedSortEntries.push_back(FSeedSortEntry(p, &PointGroup.GetOrbitType(p)));
   }
   std::stable_sort(SeedSortEntries.begin(), SeedSortEntries.end());
   for (ptrdiff_t i = 0; i < SeedPoints.cols(); ++ i) {
      SeedPoints.col(i) = SeedSortEntries[size_t(i)].p;
   }
}


void FSphereGridOptContext::InitSeedsAndOrbits()
{
   m_nSeeds = m_SeedPoints.cols();

   io.Write("\n** Processing {} grid ({} symmetry-ops) and {} monomials (lmin={}, lmax={}) for {} seed points", m_PointGroup.Name(NAMETYPE_Symmetry), m_SymmetryOps.size(), MakeMonomialList(m_lmin,m_lmax).size(), m_lmin, m_lmax, m_nSeeds);

   io.Write("\n Normalizing and sorting initial seeds...");
   PrepareInitialSeeds(m_SeedPoints, m_PointGroup);


   io.Write(" Classifying orbits and free parameters...");
   m_OrbitTypes = m_PointGroup.GetOrbitTypes(m_SeedPoints);
   // total number of non-redundant free parameters to optimize
   // (these describe only the positions of seed points; the weights are implicit)
   m_nFreeParameters = 0,
   m_nPointsTotal = 0;
   for (size_t iSeed = 0; iSeed != m_nSeeds; ++iSeed) {
      m_nFreeParameters += m_OrbitTypes[iSeed]->nFreeParameters;
      m_nPointsTotal += m_OrbitTypes[iSeed]->Length;
      // project the points to lie exactly onto their target orbit types
      // (in case they come from some outside source, or reading them from file
      // data with finite precision or similar)
      FPoint sq = m_SeedPoints.col(iSeed);
      m_OrbitTypes[iSeed]->ProjectToOrbit(sq, true);
      m_SeedPoints.col(iSeed) = sq;
   }
   if (1) {
      // hack: re-sort (in case the orbit canonicalization moved one far) & re-calc orbit types.
      PrepareInitialSeeds(m_SeedPoints, m_PointGroup);
      m_OrbitTypes = m_PointGroup.GetOrbitTypes(m_SeedPoints);
   }
   io.WriteLine();
   io.WriteCount("Number of grid points", m_nPointsTotal, fmt::format("in {} orbits", m_OrbitTypes.size()));
   io.WriteCount("Number of free parameters", m_nFreeParameters);
}

#ifdef INCLUDE_ABANDONED
// void FSphereGridOptContext::InitTargetFns()
// {
//    // Compute set of target functions which our integration rule is to evaluate
//    // exactly, and compute their overlap matrix and a suitable decomposition of
//    // it to solve conditioning equations with.
//    m_TargetFns = MakeSymmetryUniqueMonomialList(m_lmin, m_lmax, m_PointGroup, m_Options.iPrintLevel, &m_Options);
//    m_nTargetFns = m_TargetFns.size();
//    if (g_UseGroupSymmetrizedFn)
//       io.WriteCount("Number of symmetrized monomials", m_nTargetFns, "symmetry unique");
//    else
//       io.WriteCount("Number of monomials", m_nTargetFns, "symmetry unique");
// 
//    {
//       FDenseMatrix
//          S = MakeOverlapMatrix(m_TargetFns);
//       if (g_UseGroupSymmetrizedFn) {
//          // Group-averaging causes hard singularities in the function space,
//          // because the parts of the functions which lie outside the space
//          // projected to by
//          //
//          //    P_G = \sum_g \sum_j |Rg m_j><Rg m_j|
//          //
//          // are gone. Now, *mathematically* they are *not needed* for the
//          // integration rule optimization (see paper).
//          //
//          // However, their absence also causes some problems. In particular,
//          // the Choleksy solver can't deal with this---this needs a rank-
//          // revealing solver which can explicitly determine the non-singular
//          // space.
//          // So use eigh-solver by default in this case.
//          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Eigh);
// //          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_LDLT);
// //          io.WriteCount("Number of target fn", m_pScd->nEffectiveRank(), "linearly independent");
//       } else {
//          // if not symmetrizing, we will get a wide spectrum of S eigenvalues;
//          // the overlap matrix will have a large condition number for larger L.
//          // However, none of the functions will be strictly linearly dependent,
//          // and so the entire space must be handled accurately. Choleksy
//          // decompositions with triangular Gaussian elimination tend to do this
//          // rather well---better than the Smh version certainly.
//          // So use Scd-solver by default in this case.
//          //
//          // (TODO: hack the LDLT solver to support our use with symmetric
//          // HalfSolve1/HalfSolve2; could totally be done, and it has pivoting,
//          // which could be rather useful.)
// //          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Cholesky);
//          m_pScd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_LDLT);
//       }
//       io.WriteCount("Number of target fn", m_pScd->nEffectiveRank(), "linearly independent");
//    }
//    
//    m_rhs.resize(m_nTargetFns);
// #ifdef INCLUDE_ABANDONED
// //    std::cout << fmt::format("\n** Processing {} symmetry-ops and {} symmetry-unique monomials (lmin={}, lmax={}) for {} seed points", SymmetryOps.size(), nMonomials, lmin, lmax, nSeeds) << std::endl;
// #endif // INCLUDE_ABANDONED
//    for (size_t i = 0; i != m_nTargetFns; ++ i)
//       m_rhs[i] = m_TargetFns[i].CalcUnitSphereIntegral();
//    FScalar
//       fMinRhs = m_rhs.maxCoeff();
//    for (size_t i = 0; i != m_nTargetFns; ++ i)
//       if (m_rhs[i] != 0)
//          fMinRhs = std::min(fMinRhs, m_rhs[i]);
//    std::cout << fmt::format("\n Built rhs vector:  MIN: {:8.2e}  MAX: {:8.2e}", double(fMinRhs), double(m_rhs.maxCoeff())) << std::endl;
// }
#endif // INCLUDE_ABANDONED

void FSphereGridOptContext::InitTargetFns()
{
   io.WriteLine();
   // Compute set of target functions which our integration rule is to evaluate
   // exactly, and compute their overlap matrix and a suitable decomposition of
   // it to solve conditioning equations with.
   m_pTargetFns = 0;
   FTargetFnListOptions
      FnListOptions(PrintQ(PRINT_TargetFn), m_Options.AccuracyPreference, m_Options.NeglectResidualDerivs);
   if (m_Options.TargetFnType == FGridSearchOptions::TARGETFN_TabulatedPolynomials) {
      try {
         m_pTargetFns = new FTargetFnList_PreAveragedPolys(m_lmax, &m_PointGroup, FnListOptions, &m_Options);
      } catch (ETabulatedFnListError const &e) {
         io.Write("!Failed to load tabulated subspace basis:\n {}\nWill resort to monomial span of target space", e.what());
         m_pTargetFns = 0;
      }
   }
   if (!m_pTargetFns) {
      // ^-- if m_pTargetFns == 0 here we either failed to load tabulated
      // polynomials or we were asked to make monomials.
      m_pTargetFns = new FTargetFnList_RawMonomials(m_lmin, m_lmax, &m_PointGroup, FnListOptions, &m_Options);
   }
   m_nTargetFns = m_pTargetFns->size();
//    io.WriteCount(fmt::format("Number of {}", m_pTargetFns->FnTypeName()), m_nTargetFns, m_pTargetFns->FnTypeAnnotation());
   io.WriteCount("Target space basis functions", m_nTargetFns, m_pTargetFns->FnTypeName());
   if (PrintQ(PRINT_TargetFn).MoreQ())
      m_pTargetFns->Print(io);
   if (PrintQ(PRINT_TargetFn).BasicQ())
      m_pTargetFns->WriteInfo(io);

   m_rhs = m_pTargetFns->EvalTargetIntegrals();
   FScalar
      fMinRhs = m_rhs.maxCoeff();
   for (size_t i = 0; i != m_nTargetFns; ++ i)
      if (m_rhs[i] != 0)
         fMinRhs = std::min(fMinRhs, m_rhs[i]);
//    io.Write("\n Built rhs vector:  MIN: {:8.2e}  MAX: {:8.2e}", double(fMinRhs), double(m_rhs.maxCoeff()));
   {
      fmt::MemoryWriter w;
      using std::sqrt;
//       w.write(" Target analytical integrals:");
      w.write("\n Target rhs integrals:       ");
//       char const *pCompFmt =  "{1:8.2e} ({0})";
      char const *pCompFmt =  "{1:8.2e} {0}";
      w.write(pCompFmt, "RMS", double(sqrt(m_rhs.array().abs2().mean())));
      w << " // ";
      w.write(pCompFmt, "MIN", double(fMinRhs));
      w << " // ";
      w.write(pCompFmt, "MAX", double(m_rhs.maxCoeff()));
      io.Write(w.c_str());
//          m_pLog->Write(" {:<30}{:5} (TOT) //{:5} (ALPHA) //{:5} (BETA)", "Number of electrons:", m_WfDecl.nElec(), m_WfDecl.nElecA(), m_WfDecl.nElecB());

   }
}

static void CheckSeedsForCollapse(FPointArray const &SeedPoints, FPointGroup const &PointGroup)
{
   if (1) {
      // check for seeds which ended up at the same place. That would be bad.
      FMappedPointCloud
         PointCloud(PointGroup.GetOps());
      for (size_t iSeed = 0; iSeed != (size_t)SeedPoints.cols(); ++iSeed) {
         FMappedPointCloud::FPointEntry const
            &e = PointCloud.Insert(SeedPoints.col(iSeed), iSeed);
         if (e.iPoint != ptrdiff_t(iSeed))
            throw std::runtime_error(fmt::format("optimization failed -- seeds collapsed to equivalent points :(  ({} maps to {})", iSeed+1, e.iPoint+1));
      }
   }
}


void FSphereGridOptContext::CheckOrbitTypes()
{
   if (1) {
      // recompute orbit types and number of points, in case a point has moved to a higher symmetry orbit.
      // is that supposed to cause any problems? not sure.
      m_OrbitTypes = m_PointGroup.GetOrbitTypes(m_SeedPoints);
      m_nPointsTotal = 0;
      for (size_t iSeed = 0; iSeed != m_nSeeds; ++iSeed) {
         FPoint sq = m_SeedPoints.col(iSeed);
         m_OrbitTypes[iSeed]->ProjectToOrbit(sq, true);
         m_SeedPoints.col(iSeed) = sq;
         m_nPointsTotal += m_OrbitTypes[iSeed]->Length;
      }
   }
   
   if (1) {
      // check if each of the points has stayed on the type of orbit it is
      // supposed to be on.
      for (size_t iSeed = 0; iSeed != m_nSeeds; ++iSeed) {
         size_t nOrbitLength = m_PointGroup.CountOrbitLength(m_SeedPoints.col(iSeed));
         if (nOrbitLength != m_OrbitTypes[iSeed]->Length)
            throw std::runtime_error(fmt::format("optimization failed -- seed #{} moved out of its orbit type (expected orbit length: {}, actual orbit length: {})", iSeed+1, m_OrbitTypes[iSeed]->Length, nOrbitLength));
         if (!m_OrbitTypes[iSeed]->IsOnOrbit(m_SeedPoints.col(iSeed), g_ThrAlmostZero))
            throw std::runtime_error(fmt::format("optimization failed -- seed #{} no longer on orbit classified as '{}', despite having correct length (expected orbit length: {}, actual orbit length: {})", iSeed+1, m_OrbitTypes[iSeed]->pDesc, m_OrbitTypes[iSeed]->Length, nOrbitLength));
      }
//       std::cout << fmt::format("-- check of {} seed points passed.\n", m_nSeeds);
   }
}


void FSphereGridOptContext::RunIterations()
{
   m_FinalResidual = 1e6;
   m_Converged = false;
   m_nExtraDofForWeights = 0;

   // initialize generators of rotation around x,y, and z axis
   m_A[0] << 0,0,0, 0,0,1, 0,-1,0;
   m_A[1] << 0,0,-1, 0,0,0, 1,0,0;
   m_A[2] << 0,1,0, -1,0,0, 0,0,0;

   m_KtkdSolverType = LINSOLVE_PositiveSymmetric_LDLT;
   // ^- may be changed during iterations if this one ends up exploding

   std::vector<FScalar>
      // log of target residuals we achieved in each iteration
      // (for abort check if there is no progress)
      m_fResidualAtItList;

   // next ones are for dynamic step length control, if enabled (Options.MaxStep_DynamicControl)
   FScalar
      // factor adjusted to increase/decrease the current maximum step length,
      // in response to good/bad steps being taken. A value of 1.0 means that
      // the step length generated by the default mechanism is not influenced.
      fMaxStep_DynamicFactor = FScalar(1.0);
   FPointArray
      // last set of seed points which resulted in a decrease of the residual.
      // This is where we will revert to after taking a bad step.
      DynControl_LastGoodSeedPoints = m_SeedPoints;
   FScalar
      // residual at DynControl_LastGoodSeedPoints. New steps will only be
      // accepted if they beat this.
      DynControl_LastGoodResidual = FScalar(1e99);
   bool
      // set (once) after a revert flag is set: it means that the current
      // step being considered is one which has been done with the same
      // trial solution vector before, but is now being re-done with a smaller
      // maximum step width;
      DynControl_StepIsRecalcOfRevertedStep = false;
   bool
      DynControl_Enabled = m_Options.pStepOptions->MaxStep_DynamicControl;
   unsigned
      DynControl_TimesSwitched = 0;
   
   ptrdiff_t
      iIt;
   ct::FTimer
      tMainIterations;
//    io.Write("\n   {:^6} {:^8}   {:8^}     {:8^}  {}  |  WEIGHTS", "ITER.", "RESIDUAL", "UPDATE", "MAX.STEP", "NOTE");
   io.Write("\n  {:>6}  {:^10}  {:^10}  {:^10} {:^3}{:>8} |     WEIGHTS", "ITER.", "RESIDUAL", "UPDATE", "MAX.STEP", "NOTE", "E.TIME");
   
   // this is the main loop of the iterative optimization of a grid.
   for (iIt = 0; iIt < ptrdiff_t(m_Options.pStepOptions->MaxIt); ++ iIt) {
      if (PrintQ(PRINT_Seeds).MostQ()) {
         io.Write("      Checking for seed collapse. Current points:");
         for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
            FPoint sq = m_SeedPoints.col(iSeed);
            size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
            io.Write("        {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
               double(sq[0]), double(sq[1]), double(sq[2]), m_OrbitTypes[iSeed]->pDesc, nOrbitLength);
         }
      }
      CheckSeedsForCollapse(m_SeedPoints, m_PointGroup);
      CheckOrbitTypes();

      FDenseVector
         weights, // dimension: nSeeds (q)
         residual; // dimension: nTargetFn (i)
      FDenseMatrix
         // jacobian: format: i (residual), p_alpha (nVars rotation parameters)
         dK_dp_alpha;
      ComputeWeightsResidualAndJacobian(weights, residual, dK_dp_alpha);
//       ProjectRotationsToOrbits(residual); // not allowed to move pts out of their orbit types. FIXME: should this be included self-consistently in Jacobian determination?


      FDenseVector
         // step vector of rotation angles p_alpha for the seed points.
         // dimension: p_alpha (nVars rotation parameters)
         dx;
      FScalar
         // actual value of trust region size (as result of parameters
         // like residual size, iteration number, etc), which was used
         // for determining 'dx'.
         ActualMaxStepUsed;
      ComputeStep(ActualMaxStepUsed, dx, residual, dK_dp_alpha, iIt, fMaxStep_DynamicFactor);

      FScalar
//          fResidual = Rmsd(m_pScd->HalfSolve1(residual), m_nFreeParameters),
         // ^- compute residual wrt. orthogonalized function set
         // UPDATE: function set now assumed semi-orthogonal as it comes out of FTargetFnList.
         // (so FTargetFnList should do the orthogonalization stuff if required)
         fResidual = Rmsd(residual, m_nFreeParameters),
         fStep = Rmsd(dx, m_nFreeParameters);
      m_FinalResidual = fResidual;
      m_fResidualAtItList.push_back(fResidual);
      
      bool
         AcceptThisStep = true;
      std::string
         StepAnnotation = "";
#ifdef INCLUDE_ABANDONED
//       if (DynControl_Enabled && fStep < FScalar(1e-4) * ActualMaxStepUsed) {
//          // turn off dynamic control once we reach the quadratic region.
//          DynControl_Enabled = false;
//          StepAnnotation = "--";
//       }
#endif // INCLUDE_ABANDONED
      if (DynControl_Enabled) {
//          if (DynControl_StepIsRecalcOfRevertedStep || fResidual <= DynControl_LastGoodResidual || fStep < FScalar(1e-2) * ActualMaxStepUsed) {
         if (DynControl_StepIsRecalcOfRevertedStep || fResidual <= DynControl_LastGoodResidual * 1.0001) {
            AcceptThisStep = true;
         } else {
            AcceptThisStep = false;
            StepAnnotation += "!"; // bad step---this one will be reverted.
         }
         if (DynControl_StepIsRecalcOfRevertedStep)
            StepAnnotation += "D"; // restricted retry on a previously computed seed point set, which was declined (implies R).
//          if (fMaxStep_DynamicFactor < FScalar(1))
//             StepAnnotation += "R"; // step was dynamically restricted
      }
      if (fStep >= 0.99 * ActualMaxStepUsed)
         StepAnnotation += "R"; // step was dynamically restricted
      
//       std::cout << fmt::format("\n NEXT: Accept? {} fResidual = {:8.2e}  LastGood = {:8.2e}\n\n", AcceptThisStep, double(fResidual), double(DynControl_LastGoodResidual));
      
      {
         fmt::MemoryWriter w;
//          w.write(" {:6d}   {:8.2e}   {:8.2e}   {:8.2e}   {:3}  | ", iIt+1, double(fResidual), double(fStep), double(ActualMaxStepUsed), StepAnnotation);
         w.write(" {:6d}  {:10.2e}  {:10.2e}  {:10.2e} {:^3}{:10.2f} | ", iIt+1, double(fResidual), double(fStep), double(ActualMaxStepUsed), StepAnnotation, double(m_tThisGridTotal));
         size_t nSeedsToPrint = std::min(size_t(m_nSeeds), size_t(4));
         for (size_t iSeed = 0; iSeed != nSeedsToPrint; ++ iSeed)
            w.write(" {:21.16f}", double(weights[iSeed]));
         if (nSeedsToPrint < m_nSeeds)
            w.write("   ...");
         io.Write(w.c_str());
      }

      if (DynControl_Enabled) {
         if (!AcceptThisStep) {
            // This step did not reduce the target residual. It should have.
            // Reduce dynamic step length control factor, revert the step, and try again.
            fMaxStep_DynamicFactor *= m_Options.pStepOptions->fMaxStep_DynamicControl_FactorAfterBadStep;
            m_SeedPoints = DynControl_LastGoodSeedPoints;
            DynControl_StepIsRecalcOfRevertedStep = true;
            continue;
         } else {
            if (!DynControl_StepIsRecalcOfRevertedStep) {
               // this step did reduce the residual. Remember the current seed
               // points and new residual, in case we need it as a reference point
               // to revert to later on.
               DynControl_LastGoodResidual = fResidual;
               DynControl_LastGoodSeedPoints = m_SeedPoints;
               // since we got a good step, also gradually increase the control
               // factor for the maximum step length (in case the current step was
               // dynamically restricted)
               fMaxStep_DynamicFactor *= pow(m_Options.pStepOptions->fMaxStep_DynamicControl_FactorAfterBadStep, FScalar(-1)/FScalar(m_Options.pStepOptions->fMaxStep_DynamicControl_NumIterToReset));
               if (fMaxStep_DynamicFactor > FScalar(1.0))
                  fMaxStep_DynamicFactor = FScalar(1.0);
            }
         }
         DynControl_StepIsRecalcOfRevertedStep = false;
      }
      

      if ((fStep < m_Options.fThreshResidual() && iIt >= g_nDampedSigSteps) || fResidual < m_Options.fThreshResidual())
         m_Converged = true;

      if (m_Converged)
         break;
      
      if (m_Options.pStepOptions->nIterProgressCheck > 0) {
         ptrdiff_t
//             nItProgressCheck = 100;
            nItProgressCheck = m_Options.pStepOptions->nIterProgressCheck;
         unsigned
            DynControl_nMaxSwitch = m_Options.pStepOptions->ProgressCheck_DynControl_nMaxSwitch;
         if (iIt > 2*nItProgressCheck && (ptrdiff_t(m_fResidualAtItList.size()) > 2 + nItProgressCheck)) {
            size_t
               iRefIt = m_fResidualAtItList.size() - 1;
            
            FScalar
               fResNow = std::min(m_fResidualAtItList[iRefIt], m_fResidualAtItList[iRefIt-1]),
               fResEarier = std::min(m_fResidualAtItList[iRefIt-nItProgressCheck], m_fResidualAtItList[iRefIt-nItProgressCheck+1]);
//             if (!(fResNow < 0.99 * fResEarier))
//             if (!(abs(fResNow-fResEarier) > 0.01 * fResNow)) {
//             if (!(abs(fResNow-fResEarier) > 1e-6 * fResNow)) {
            if (!(fResNow < FScalar(1 - FScalar(1e-5)) * fResEarier)) {
//                throw std::runtime_error(fmt::format("no progress during last {} iterations: fResNow = {:8.2e}, fResEarier = {:8.2e}. Aborting program.", nItProgressCheck, double(fResNow), double(fResEarier)));
//                m_FinalResidual = 1e20;
               std::string
                  sAction = "Aborting optimization";
               bool
                  AbortOpt = true;
               if (m_Options.pStepOptions->MaxStep_DynamicControl &&
                   DynControl_TimesSwitched < DynControl_nMaxSwitch) {
                  if (DynControl_Enabled) {
                     sAction = "Disabling dynamic step length control";
                  } else {
                     sAction = "Enabling dynamic step length control";
                  }
                  sAction += fmt::format(" (switch {} of {})", DynControl_TimesSwitched+1, DynControl_nMaxSwitch);
                  AbortOpt = false;
               }
               DynControl_TimesSwitched += 1;
               io.Write("\n WARNING: no progress during last {} iterations: fResNow = {:8.2e}, fResEarier = {:8.2e}. {}.\n", nItProgressCheck, double(fResNow), double(fResEarier), sAction);
               if (AbortOpt) {
//                m_FinalResidual = 1e20;
                  m_Converged = false;
                  break;
               } else {
                  DynControl_Enabled = !DynControl_Enabled;
                  fMaxStep_DynamicFactor = 1.0;
                  m_fResidualAtItList.clear();
               }
            }
         }
      }

      UpdateSeedPositions(dx, fResidual, iIt);
   }
   io.WriteLine();
   io.WriteTiming("iteration", double(tMainIterations), iIt);
   io.WriteLine();

   if (m_Converged) {
      io.Write(" Converged. Final seed points & weights (per unique point):");
   } else {
      io.Write(" WARNING: FAILED to converge in {} iterations. Aborted. Final seed points & weights (per unique point):", 1+iIt);
   }
}


#ifdef INCLUDE_ABANDONED
// #if 0
// void FSphereGridOptContext::ComputeWeightsResidualAndJacobian(FDenseVector &weights, FDenseVector &residual, FDenseMatrix &dK_dp_alpha)
// {
//    size_t
//       nVars = 3*m_nSeeds;
//    FDenseMatrix
//       // K[i,q] = \sum_g f_i(Rg r_q): (i: target fn index, q: seed index)
//       K(m_nTargetFns, m_nSeeds),
//       // dKw[i,p] = \sum_{q} dK[i,q]/dp w[q] -- kernel derivative contracted to fixed weights.
//       dKw(m_nTargetFns, nVars),
//       // dKr[q,p] = \sum_{i} dK[i,q]/dp r[i] -- kernel derivative contracted to fixed m_rhs r[i]
//       dKr(m_nSeeds, nVars);
//    bool
//       // In an exact d e/d p_alpha derivative, it is supposed to be there,
//       // according to the formal derivation. However, sometimes these kinds of
//       // methods work better if some derivative contributions are omitted... and
//       // this appears to be one of those cases.
//       Omit_dKr_Term;
//    Omit_dKr_Term = false;
// 
//    FTransformedPointInfo
//       // contains information about seed points transformed by current symmetry operation.
//       // In particular, powers of its x,y,z coordinates.
//       TrafodSeeds(m_SeedPoints, m_lmax);
//    K.setZero();
//    for (size_t iSymOp = 0; iSymOp != m_SymmetryOps.size(); ++ iSymOp) {
//       TrafodSeeds.TransformPoints(m_SymmetryOps[iSymOp]);
//       for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//          for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//             K(iTargetFn, iSeed) += m_TargetFns[iTargetFn].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
//          }
//       }
//    }
// 
// #ifdef INCLUDE_ABANDONED
// //       FDenseMatrix
// //          SmK = Scd.solve(K);
// //       
// //       FDenseMatrix
// //          // [K^T S^{-1} K]_{q_1, q_2}
// //          KtK(m_nSeeds, m_nSeeds);
// //       KtK = SmK.transpose() * K;
// //       Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
// //          KtKd(KtK, Eigen::ComputeThinU | Eigen::ComputeThinV);
// //       // ^- FIXME: try doing this with CD instead of SVD, too?
// //       KtKd.setThreshold(FScalar(1e-16Q));
// //       // ^- note: for REALLY small target ls we might have a KtK which is not
// //       //    positive definite... so we use SVD instead of LLT decomposition
// //       //    (eigenvalues might also do the trick; maybe consider later.)
// //       
// //       // count number of degrees of freedom remaining in space of weigthts (for current seedpoints)
// //       // which are *not* fixed by the demand that the residual is minimal for the current points for the current lmax.
// //       m_nExtraDofForWeights = ptrdiff_t(m_nSeeds);
// //       FDenseVector const
// //          &KtK_sigma = KtKd.singularValues();
// //       for (ptrdiff_t i = 0; i != ptrdiff_t(KtK_sigma.size()); ++ i)
// //          if (KtK_sigma[i] > KtKd.threshold())
// //             m_nExtraDofForWeights -= 1;
// // 
// //       FDenseVector
// //          // dimension: q (seeds)
// //          weights = KtKd.solve(SmK.transpose() * m_rhs),
// //          // dimension: i (monomial)
// //          dm_rhs = Scd.solve(m_rhs - K * weights);
// //       m_LastWeights = weights;
// //       m_nExtraDofForWeights = 
// #endif // INCLUDE_ABANDONED
//       FDenseMatrix
//          SmhK = m_pScd->HalfSolve1(K),
//          // [K^T S^{-1} K]_{q_1, q_2}
//          KtK(m_nSeeds, m_nSeeds);
//       KtK = SmhK.transpose() * SmhK; // <- that's a Syrk
//       FLinearSolverPtr
//          // compute KtK matrix decomposition for linear system solving.
//          pKtKd = MakeLinearSolver(KtK, m_KtkdSolverType);
//       // Notes:
//       // - LDLT is a pivoted Cholesky decomposition.
//       //   using a Cholesky decomp for this is very optimistic in some
//       //   situations... especailly for small L, large rules, and exact
//       //   function symmetrization, the system KtK is often *meant* to be
//       //   singular!
//       //   It is actually quite surprising that Eigen can deal with
//       //   this to this degree on smaller grids. However, for larger grid
//       //   scans, such as ico-13-1, eventually LLT explodes. Interestingly,
//       //   just replacing this by LDLT (pivoted Cholesky) appears to work fine
//       //   as a workaround, although, realistically speaking, the system is
//       //   still supposed to be singular.
//       // - So we next try to see if it exploded (like it should), and if yes,
//       //   replace it by the spectral decomp solver which can isolate the singular
//       //   space if needed.
//       // - TODO: check if this is *really* a good idea. I thought there
//       //   was a reason I put in this LDLT thing... maybe this has something
//       //   to do with some of the large grids looking rather shady?
//       if (1) {
//          bool
//             IsCdSolver = (m_KtkdSolverType == LINSOLVE_PositiveSymmetric_LDLT) ||
//                         (m_KtkdSolverType == LINSOLVE_PositiveSymmetric_Cholesky);
//          if (pKtKd->nEffectiveRank() < size_t(KtK.rows()) && IsCdSolver) {
//             if (m_Options.iPrintLevel >= 1) {
//                std::cout << " NOTE: switching KtK solver from LDLT to Eigh due to singularities." << std::endl;
//             }
//             m_KtkdSolverType = LINSOLVE_PositiveSymmetric_Eigh;
//             pKtKd = MakeLinearSolver(KtK, m_KtkdSolverType);
//          }
//       }
//       
//       if (0) {
//          // count number of degrees of freedom remaining in space of weigthts (for
//          // current seedpoints) which are *not* fixed by the demand that the
//          // residual is minimal for the current points for the current lmax.
//          m_nExtraDofForWeights = m_nSeeds - pKtKd->nEffectiveRank();
//       }
// 
// 
//       // weights: dimension: q (seeds)
// #ifdef INCLUDE_ABANDONED
// //       weights = pKtKd->solve(SmK.transpose() * m_rhs),
// #endif // INCLUDE_ABANDONED
//       weights = pKtKd->Solve(SmhK.transpose() * m_pScd->HalfSolve1(m_rhs)),
//       m_LastWeights = weights;
//       FDenseVector
//          // dm_rhs: dimension: i (target fn)
//          dm_rhs = m_pScd->Solve(m_rhs - K * weights);
// 
//       dKw.setZero();
//       dKr.setZero();
// 
//       for (size_t iSymOp = 0; iSymOp != m_SymmetryOps.size(); ++ iSymOp) {
//          FRotationMatrix const
//             &Rg = m_SymmetryOps[iSymOp],
//             Rg_Ax = Rg * m_A[0],
//             Rg_Ay = Rg * m_A[1],
//             Rg_Az = Rg * m_A[2];
//          FPointArray
//             Rg_Ax_pts = Rg_Ax * m_SeedPoints,
//             Rg_Ay_pts = Rg_Ay * m_SeedPoints,
//             Rg_Az_pts = Rg_Az * m_SeedPoints;
//          TrafodSeeds.TransformPoints(m_SymmetryOps[iSymOp]);
// 
//          // project out rotations which would lead out of the given orbit type.
// #ifdef INCLUDE_ABANDONED
// //          for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed)
// //             m_PointGroup.ProjectRotationToOrbitType(Rg_Ax_pts.col(iSeed), Rg_Ay_pts.col(iSeed), Rg_Az_pts.col(iSeed), TrafodSeeds.TransformedPoints.col(iSeed), *m_OrbitTypes[iSeed]);
// #endif // INCLUDE_ABANDONED
//          m_PointGroup.ProjectRotationToOrbitType(Rg_Ax_pts, Rg_Ay_pts, Rg_Az_pts, TrafodSeeds.TransformedPoints, &m_OrbitTypes[0]);
// 
//          // evaluate contracted derivatives of K for current group operation.
//          for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//             for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//                FPoint
//                   dMdXyz,
//                   Rg_Ax_sq = Rg_Ax_pts.col(iSeed),
//                   Rg_Ay_sq = Rg_Ay_pts.col(iSeed),
//                   Rg_Az_sq = Rg_Az_pts.col(iSeed);
//                m_TargetFns[iTargetFn].EvalPoint(dMdXyz, TrafodSeeds.TransformedPointInfo(iSeed));
// //                m_TargetFns[iTargetFn].EvalPoint(dMdXyz, Rg*m_SeedPoints.col(iSeed).matrix());
//                FPoint
//                   dKiq_dp_alpha(0,0,0);
//                for (size_t beta = 0; beta != 3; ++ beta) {
//                   dKiq_dp_alpha[0] += dMdXyz[beta] * Rg_Ax_sq[beta];
//                   dKiq_dp_alpha[1] += dMdXyz[beta] * Rg_Ay_sq[beta];
//                   dKiq_dp_alpha[2] += dMdXyz[beta] * Rg_Az_sq[beta];
//                }
//                size_t iVar = 3*iSeed;
//                for (size_t beta = 0; beta != 3; ++ beta) {
//                   dKw(iTargetFn, iVar + beta) += dKiq_dp_alpha[beta] * weights[iSeed];
//                   if (!Omit_dKr_Term)
//                      dKr(iSeed, iVar + beta) += dKiq_dp_alpha[beta] * dm_rhs[iTargetFn];
//                }
//             }
//          }
//       }
// 
//       FDenseMatrix
//          // compute S^{-1} K
//          SmK = m_pScd->HalfSolve2(SmhK); // <- TODO: remove this and replace by 2nd HalfSolve2 on the product SmK.transpose() * dKw. SmK not used elsewhere (anymore?)
//          // UPDATE: Note @ fixme: I can do the replacement, but is that really better? dKw is a matrix, too (nTargetFn,nVars),
//          // and SmhK is only (nTargetFn,nSeeds). Seems like doing the transform on SmhK would be the better option
//          // as far as flops are concerned. Stability is another thing, of course...
//       FDenseMatrix
//          dweight_dp_alpha = pKtKd->Solve(dKr - SmK.transpose() * dKw);
// #ifdef INCLUDE_ABANDONED
// //       FDenseMatrix
// //          // compute K^T * S^{-1} dK/dw
// // //          KtSm_dKw = K.transpose() * m_pScd->Solve(dKw);
// // //          KtSm_dKw = m_pScd->HalfSolve2(SmhK).transpose() * dKw;
// //          KtSm_dKw = SmhK.transpose() * m_pScd->HalfSolve1(dKw);
// //       FDenseMatrix
// //          dweight_dp_alpha = pKtKd->Solve(dKr - KtSm_dKw);
// // //          (K^T S^{-1} K)^{-1} K S^{-1} dKw
// #endif // INCLUDE_ABANDONED
//       // this is the final jacobian!
//       // format: i (residual), p_alpha (nVars rotation parameters)
//       dK_dp_alpha = dKw + K * dweight_dp_alpha;
// 
//       residual = K * weights - m_rhs;
//       if (m_Options.iPrintLevel >= 4) {
//          std::cout << " Intermediate residual contributions for current points:\n";
//          for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//             FPoint sq = m_SeedPoints.col(iSeed);
//             size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
//             std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
//                double(sq[0]), double(sq[1]), double(sq[2]), m_OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//             for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//                FMonomial m = m_TargetFns[iTargetFn];
//                std::cout << fmt::format("     |  [{:>6}]  m = {:<20}  K[m,r] = {:24.16e}   D[K(m),p_alpha] = {:24.16e} {:24.16e} {:24.16e}   residual[m] = {:16.8e}\n",
//                   iTargetFn, fmt::format("[{},{},{}]",m.ix, m.iy, m.iz), double(K(iTargetFn,iSeed)), double(dK_dp_alpha(iTargetFn,3*iSeed+0)), double(dK_dp_alpha(iTargetFn,3*iSeed+1)), double(dK_dp_alpha(iTargetFn,3*iSeed+2)), double(residual[iTargetFn]));
//             }
//          }
//       }
// }
// #else
// 
// 
// // this version works, but I couldn't get the Scd thing running on first try.
// // Also, it is not faster, and it does appear as if it is not at all, or not much,
// // more stable than the other version. Hard to say. Only tried 64bit arithmetic,
// // and there both work very similarly. I leave it on for the moment because
// // there are some formal reasons to believe that this *should* be more stable.
// void FSphereGridOptContext::ComputeWeightsResidualAndJacobian(FDenseVector &weights, FDenseVector &residual, FDenseMatrix &dK_dp_alpha)
// {
//    size_t
//       nVars = 3*m_nSeeds;
//    FDenseMatrix
//       // K[i,q] = \sum_g f_i(Rg r_q): (i: target fn index, q: seed index)
//       K(m_nTargetFns, m_nSeeds),
//       // dKw[i,p] = \sum_{q} dK[i,q]/dp w[q] -- kernel derivative contracted to fixed weights.
//       dKw(m_nTargetFns, nVars),
//       // dKe[q,p] = \sum_{i} dK[i,q]/dp e[i] -- kernel derivative contracted to fixed m_rhs e[i]
//       dKe(m_nSeeds, nVars);
// 
//    FTransformedPointInfo
//       // contains information about seed points transformed by current symmetry operation.
//       // In particular, powers of its x,y,z coordinates.
//       TrafodSeeds(m_SeedPoints, m_lmax);
//    K.setZero();
//    for (size_t iSymOp = 0; iSymOp != m_SymmetryOps.size(); ++ iSymOp) {
//       TrafodSeeds.TransformPoints(m_SymmetryOps[iSymOp]);
//       for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//          for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//             K(iTargetFn, iSeed) += m_TargetFns[iTargetFn].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
//          }
//       }
//    }
// 
// //    FLinearSolver
// //       *pScd = m_pScd.get();
//    FDenseMatrix
// //       SmhK = pScd->HalfSolve1(K);
//       SmhK = K;
//    FSvd
//       SvdK(SmhK);
//    SvdK.Compress(g_ThrAlmostZero);
// 
//    // weights: dimension: q (seeds)
// //    weights = SvdK.Solve(pScd->HalfSolve1(m_rhs)),
//    weights = SvdK.Solve(m_rhs),
//    m_LastWeights = weights;
// 
// 
//    FDenseVector
//       // residual: dimension: i (target fn)
// //       err = pScd->Solve(K * weights - m_rhs);
//       err = K * weights - m_rhs;
// 
//    dKw.setZero();
//    dKe.setZero();
// 
//    for (size_t iSymOp = 0; iSymOp != m_SymmetryOps.size(); ++ iSymOp) {
//       FRotationMatrix const
//          &Rg = m_SymmetryOps[iSymOp],
//          Rg_Ax = Rg * m_A[0],
//          Rg_Ay = Rg * m_A[1],
//          Rg_Az = Rg * m_A[2];
//       FPointArray
//          Rg_Ax_pts = Rg_Ax * m_SeedPoints,
//          Rg_Ay_pts = Rg_Ay * m_SeedPoints,
//          Rg_Az_pts = Rg_Az * m_SeedPoints;
//       TrafodSeeds.TransformPoints(m_SymmetryOps[iSymOp]);
// 
//       // project out rotations which would lead out of the given orbit type.
//       m_PointGroup.ProjectRotationToOrbitType(Rg_Ax_pts, Rg_Ay_pts, Rg_Az_pts, TrafodSeeds.TransformedPoints, &m_OrbitTypes[0]);
// 
//       // evaluate contracted derivatives of K for current group operation.
//       for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//          for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//             FPoint
//                dMdXyz,
//                Rg_Ax_sq = Rg_Ax_pts.col(iSeed),
//                Rg_Ay_sq = Rg_Ay_pts.col(iSeed),
//                Rg_Az_sq = Rg_Az_pts.col(iSeed);
//             m_TargetFns[iTargetFn].EvalPoint(dMdXyz, TrafodSeeds.TransformedPointInfo(iSeed));
//             FPoint
//                dKiq_dp_alpha(0,0,0);
//             for (size_t beta = 0; beta != 3; ++ beta) {
//                dKiq_dp_alpha[0] += dMdXyz[beta] * Rg_Ax_sq[beta];
//                dKiq_dp_alpha[1] += dMdXyz[beta] * Rg_Ay_sq[beta];
//                dKiq_dp_alpha[2] += dMdXyz[beta] * Rg_Az_sq[beta];
//             }
//             size_t iVar = 3*iSeed;
//             for (size_t beta = 0; beta != 3; ++ beta) {
//                dKw(iTargetFn, iVar + beta) += dKiq_dp_alpha[beta] * weights[iSeed];
//                dKe(iSeed, iVar + beta) += dKiq_dp_alpha[beta] * err[iTargetFn];
//             }
//          }
//       }
//    }
// 
//    // this is the final jacobian!
//    // format: i (residual), p_alpha (nVars rotation parameters)
// //    FDenseMatrix
// //       Smh_dKw = m_pScd->HalfSolve1(dKw);
// //    dK_dp_alpha = Smh_dKw - SvdK.U * (SvdK.U.transpose() * Smh_dKw).eval() - m_pScd->Solve(SvdK.SolveT(dKe));
// //    dK_dp_alpha = dKw - SvdK.U * (SvdK.U.transpose() * dKw).eval() - m_pScd->HalfSolve1(SvdK.SolveT(dKe));
//    dK_dp_alpha = dKw - SvdK.U * (SvdK.U.transpose() * dKw).eval() - SvdK.SolveT(dKe);
//    // ^-- (1 - U U^T) (dK/dx w) - (K^+)^T (dK/dx res)
// 
//    residual = K * weights - m_rhs;
//    if (m_Options.iPrintLevel >= 4) {
//       std::cout << " Intermediate residual contributions for current points:\n";
//       for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//          FPoint sq = m_SeedPoints.col(iSeed);
//          size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
//          std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
//             double(sq[0]), double(sq[1]), double(sq[2]), m_OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//          for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//             FMonomial m = m_TargetFns[iTargetFn];
//             std::cout << fmt::format("     |  [{:>6}]  m = {:<20}  K[m,r] = {:24.16e}   D[K(m),p_alpha] = {:24.16e} {:24.16e} {:24.16e}   residual[m] = {:16.8e}\n",
//                iTargetFn, fmt::format("[{},{},{}]",m.ix, m.iy, m.iz), double(K(iTargetFn,iSeed)), double(dK_dp_alpha(iTargetFn,3*iSeed+0)), double(dK_dp_alpha(iTargetFn,3*iSeed+1)), double(dK_dp_alpha(iTargetFn,3*iSeed+2)), double(residual[iTargetFn]));
//          }
//       }
//    }
// }
// #endif // 0
#endif // INCLUDE_ABANDONED


void FSphereGridOptContext::ComputeWeightsResidualAndJacobian(FDenseVector &weights, FDenseVector &residual, FDenseMatrix &dK_dp_alpha)
{
   size_t
      nVars = 3*m_nSeeds;
   FDenseMatrix
      // K[i,q] = \sum_g f_i(Rg r_q): (i: target fn index, q: seed index)
      K(m_nTargetFns, m_nSeeds),
      // dKw[i,p] = \sum_{q} dK[i,q]/dp w[q] -- kernel derivative contracted to fixed weights.
      dKw(m_nTargetFns, nVars),
      // dKe[q,p] = \sum_{i} dK[i,q]/dp e[i] -- kernel derivative contracted to fixed m_rhs e[i]
      dKe(m_nSeeds, nVars);

   m_pTargetFns->EvalK(K, m_SeedPoints);
   if (m_Options.AvoidSvds) {
      FDenseMatrix
         // [K^T S^{-1} K]_{q_1, q_2}
         KtK(m_nSeeds, m_nSeeds);
      KtK = K.transpose() * K; // <- that's a Syrk
      FLinearSolverPtr
         // compute KtK matrix decomposition for linear system solving.
         pKtKd = MakeLinearSolver(KtK, m_KtkdSolverType);
      // Notes:
      // - LDLT is a pivoted Cholesky decomposition.
      //   using a Cholesky decomp for this is very optimistic in some
      //   situations... especailly for small L, large rules, and exact
      //   function symmetrization, the system KtK is often *meant* to be
      //   singular!
      //   It is actually quite surprising that Eigen can deal with
      //   this to this degree on smaller grids. However, for larger grid
      //   scans, such as ico-13-1, eventually LLT explodes. Interestingly,
      //   just replacing this by LDLT (pivoted Cholesky) appears to work fine
      //   as a workaround, although, realistically speaking, the system is
      //   still supposed to be singular.
      // - So we next try to see if it exploded (like it should), and if yes,
      //   replace it by the spectral decomp solver which can isolate the singular
      //   space if needed.
      // - TODO: check if this is *really* a good idea. I thought there
      //   was a reason I put in this LDLT thing... maybe this has something
      //   to do with some of the large grids looking rather shady?
      if (1) {
         bool
            IsCdSolver = (m_KtkdSolverType == LINSOLVE_PositiveSymmetric_LDLT) ||
                        (m_KtkdSolverType == LINSOLVE_PositiveSymmetric_Cholesky);
         if (pKtKd->nEffectiveRank() < size_t(KtK.rows()) && IsCdSolver) {
            if (m_Options.iPrintLevel >= 1) {
               std::cout << " NOTE: switching KtK solver from LDLT to Eigh due to singularities." << std::endl;
            }
            m_KtkdSolverType = LINSOLVE_PositiveSymmetric_Eigh;
            pKtKd = MakeLinearSolver(KtK, m_KtkdSolverType);
         }
      }

      // weights: dimension: q (seeds)
      weights = pKtKd->Solve(K.transpose() * m_rhs),
      m_LastWeights = weights;
      // dm_rhs: dimension: i (target fn) ... hm, wait this had the opposite sign?!
      residual = K * weights - m_rhs;
      m_pTargetFns->EvalDerivK(dKw, weights, dKe, residual, m_SeedPoints, m_OrbitTypes, &m_A[0]);

      FDenseMatrix
         dweight_dp_alpha = pKtKd->Solve(dKe + K.transpose() * dKw);
      // this is the final jacobian!
      // format: i (residual), p_alpha (nVars rotation parameters)
      dK_dp_alpha = dKw - K * dweight_dp_alpha;
   } else {
      FSvd
         SvdK(K);
      SvdK.Compress(g_ThrAlmostZero);

      // weights: dimension: q (seeds)
      weights = SvdK.Solve(m_rhs);
      m_LastWeights = weights;

      FDenseVector
         // residual: dimension: i (target fn)
         err = K * weights - m_rhs;

      m_pTargetFns->EvalDerivK(dKw, weights, dKe, err, m_SeedPoints, m_OrbitTypes, &m_A[0]);

      // this is the final jacobian!
      // format: i (residual), p_alpha (nVars rotation parameters)
      dK_dp_alpha = dKw - SvdK.U * (SvdK.U.transpose() * dKw).eval() - SvdK.SolveT(dKe);
      // ^-- (1 - U U^T) (dK/dx w) - (K^+)^T (dK/dx res)

   //    residual = K * weights - m_rhs;
      residual = err;
   }
   
   if (PrintQ(PRINT_Step).AllQ()) {
      std::cout << " Intermediate residual contributions for current points:\n";
      for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
         FPoint sq = m_SeedPoints.col(iSeed);
         size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
         std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
            double(sq[0]), double(sq[1]), double(sq[2]), m_OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
         for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
            std::string sTargetFn = m_pTargetFns->MakeFnDesc(iTargetFn);
            std::cout << fmt::format("     |  [{:>6}]  f = {:<20}  K[m,r] = {:24.16e}   D[K(m),p_alpha] = {:24.16e} {:24.16e} {:24.16e}   residual[m] = {:16.8e}\n",
               iTargetFn, sTargetFn, double(K(iTargetFn,iSeed)), double(dK_dp_alpha(iTargetFn,3*iSeed+0)), double(dK_dp_alpha(iTargetFn,3*iSeed+1)), double(dK_dp_alpha(iTargetFn,3*iSeed+2)), double(residual[iTargetFn]));
         }
      }
   }
#ifdef INCLUDE_ABANDONED
//    if (m_Options.iPrintLevel >= 4) {
//       std::cout << " Intermediate residual contributions for current points:\n";
//       for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
//          FPoint sq = m_SeedPoints.col(iSeed);
//          size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
//          std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
//             double(sq[0]), double(sq[1]), double(sq[2]), m_OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//          for (size_t iTargetFn = 0; iTargetFn != m_nTargetFns; ++ iTargetFn) {
//             FMonomial m = m_TargetFns[iTargetFn];
//             std::cout << fmt::format("     |  [{:>6}]  m = {:<20}  K[m,r] = {:24.16e}   D[K(m),p_alpha] = {:24.16e} {:24.16e} {:24.16e}   residual[m] = {:16.8e}\n",
//                iTargetFn, fmt::format("[{},{},{}]",m.ix, m.iy, m.iz), double(K(iTargetFn,iSeed)), double(dK_dp_alpha(iTargetFn,3*iSeed+0)), double(dK_dp_alpha(iTargetFn,3*iSeed+1)), double(dK_dp_alpha(iTargetFn,3*iSeed+2)), double(residual[iTargetFn]));
//          }
//       }
//    }
#endif // INCLUDE_ABANDONED
}
// #endif // 0



void FSphereGridOptContext::ProjectRotationsToOrbits(FDenseVector &dx)
{
#ifdef INCLUDE_ABANDONED
//    // FIXME: nope... these are not the tangent vectors, are they?
//    // don't I need the anti-symmetric thingy?
//    cx_assert_rt(size_t(dx.size()) == size_t(m_nSeeds*3));
//    for (size_t iSeed = 0; iSeed != m_nSeeds; ++iSeed) {
//       FVector3 v(dx[3*iSeed+0], dx[3*iSeed+1], dx[3*iSeed+2]);
//       FVector3 *pv = &v;
//       m_OrbitTypes[iSeed]->ProjectTangentVectors(&pv, 1, m_SeedPoints.col(iSeed), 0);
//       for (size_t ixyz = 0; ixyz < 3; ++ixyz)
//          dx[3*iSeed + ixyz] = v[ixyz];
//    }
#endif // INCLUDE_ABANDONED
   cx_assert_rt(size_t(dx.size()) == size_t(m_nSeeds*3));
   for (size_t iSeed = 0; iSeed != m_nSeeds; ++iSeed) {
      FVector3 v(dx[3*iSeed+0], dx[3*iSeed+1], dx[3*iSeed+2]);
      m_OrbitTypes[iSeed]->ProjectRotation(v, m_SeedPoints.col(iSeed), 0);
      for (size_t ixyz = 0; ixyz < 3; ++ixyz)
         dx[3*iSeed + ixyz] = v[ixyz];
   }
}


#ifdef INCLUDE_ABANDONED
// void FSphereGridOptContext::ComputeStep(FScalar &CurMaxStep, FDenseVector &dx, FDenseVector const &residual, FDenseMatrix const &Jacobian, size_t iIt, FScalar fMaxStep_DynamicFactor)
// {
//    // note: 'Jacobian' is the (nTargetFn, nVars)-shaped matrix called
//    // dK_dp_alpha at the place it is made and in the paper.
//    
//    // FIXME: this needs lots of code deleting.
//    // We might want to retain both the symmetric and the non-symmetric
//    // step variant (although we implicitly anyway always do the symmetric one
//    // when we compute the implicit weight derivatives...), but certainly
//    // not in 10000123213 different variants, and with lots of broken specials
//    //
//    // Also note the now-common FEighOrSvd class in MathOps.h, which should
//    // be used to merge MakeStep and MakeStepSymmetric.
//    CurMaxStep = 0.0002;
//    if (0) {
//       // MakeStep(dx, Jacobian, residual, CurMaxStep, iIt, fMinRes, m_nFreeParameters, m_Options);
//       throw std::runtime_error("this residual update path is not (anymore) implemented.");
//    } else if (0) {
//       FDenseMatrix
//          JtJ = Jacobian.transpose() * Jacobian;
//       FDenseVector
//          Jtg = Jacobian.transpose() * residual;
//       // FIXME: should I put in something like diag(1/norm(f_i)) before multiplying with Jt?
//       // Or is this already happening implicitly? (program says "Built m_rhs vector:  MIN: 1.00e+00  MAX: 3.00e+00")
//       // Or maybe even explicitly orthogonalize the functions, by multiplying with S^{-1} before J^T?
//       MakeStep(dx, JtJ, Jtg, FEighOrSvd::DECOMP_Eigh, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
//    } else if (0) {
//       std::cout << "   ....ran through J^T S^{-1} J path." << std::endl;
//       FDenseMatrix
//          SmJ = m_pScd->Solve(Jacobian);
//       FDenseMatrix
//          JtSmJ = SmJ.transpose() * Jacobian;
//       FDenseVector
//          JtSmg = SmJ.transpose() * residual;
//       MakeStep(dx, JtSmJ, JtSmg, FEighOrSvd::DECOMP_Eigh, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
//    } else if (1) {
//       // S = L L^T
//       // -> S^{-1} = (L L^T)^{-1} = (L^T)^{-1} L^{-1}
//       FDenseMatrix
// //             SmhJ = Scd.matrixL().solve(Jacobian);
//          SmhJ = m_pScd->HalfSolve1(Jacobian);
// #ifdef INCLUDE_ABANDONED
// //          for (size_t iFn = 0; iFn != m_nTargetFns; ++ iFn)
// //             SmhJ.row(iFn) /= sqrt(FScalar(m_TargetFns[iFn].degeneracy));
// #endif // INCLUDE_ABANDONED
//       FDenseMatrix
//          JtSmJ = SmhJ.transpose() * SmhJ;
//       FDenseVector
// #ifdef INCLUDE_ABANDONED
// //             Smhg = Scd.matrixL().solve(residual),
// #endif // INCLUDE_ABANDONED
//          Smhg = m_pScd->HalfSolve1(residual),
//          JtSmg = SmhJ.transpose() * Smhg;
// #ifdef INCLUDE_ABANDONED
// //          // remove degeneracy factors from column projector
// //          for (size_t iFn = 0; iFn != m_nTargetFns; ++ iFn) {
// // //             JtSmJ.row(iFn) /= FScalar(m_TargetFns[iFn].degeneracy);
// //             JtSmg[iFn] /= sqrt(FScalar(m_TargetFns[iFn].degeneracy));
// //          }
// //          MakeStep(dx, SmhJ, Smhg, CurMaxStep, iIt, fMinRes, m_nFreeParameters, m_Options);
// #endif // INCLUDE_ABANDONED
//       MakeStep(dx, JtSmJ, JtSmg, FEighOrSvd::DECOMP_Eigh, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
// 
//       // ^- this should give the same result as the Scd.solve variant above
//       //
//       // Note: what this does is bypassing the SVD of J by computing the spectral decomposition of J^T J
//       // instead. That should give the squared singular values as eigenvalues.
//       // but isn't the SVD normally numerically better? to get identical results one probably has
//       // to update the sigma+shift vs sigma^2+shift values, though.
//       //
//       // Also, we have many less free parameters than target functions. So the
//       // eigh actually should be quite a bit more efficient
//    } else {
//       // Here: do the equivalent of:
//       //
//       //   dx = diag(h(w)) solve(S^{-1/2} J diag(h(w)), S^{-1/2} e),
//       //
//       // -> so only should use one side of L, not a full Cholesky solve.
//       FDenseMatrix
//          SmhJ = m_pScd->HalfSolve1(Jacobian);
//       FDenseVector
//          Smhg = m_pScd->HalfSolve1(residual);
// #ifdef INCLUDE_ABANDONED
// //          FDenseVector
// //             hw(SmhJ.cols());
// //          assert(size_t(weights.size()) == size_t(m_nSeeds) && hw.size() == weights.size() * 3);
// //          for (int iSeed = 0; iSeed != weights.size(); ++ iSeed) {
// //             FScalar
// //                wi = weights[iSeed],
// //                hwi = FScalar(0);
// // //             if (abs(wi) > 1e-10 / (3*m_nSeeds))
// // // //                hwi = 1/sqrt(abs(wi));
// // //                hwi = 1/abs(wi);
// // //             else
// // //                hwi = 1.;
// //             hwi = sqrt(wi);
// // //             hwi = wi;
// //             for (int alpha = 0; alpha != 3; ++ alpha) {
// //                hw[3*iSeed+alpha] = hwi;
// //             }
// //          }
// // //          SmhJ.array().colwise() *= hw.array();
// //          for (int i = 0; i != SmhJ.cols(); ++ i)
// //             SmhJ.col(i) *= hw[i];
//       // ^- hm... none of those dynamic weight scaling thingies work. I think
//       // scaling with 1/sqrt(wi) *should* be the right variant, but it does not work either.
//       // Guess I'll just leave it alone.
//       //          std::cout << fmt::format("   hw scale: dx.size() = {}  hw.size() = {}  weights.size() = {}  hw[0] = {}\n", dx.size(), hw.size(), weights.size(), hw[0]);
// #endif // INCLUDE_ABANDONED
//       MakeStep(dx, SmhJ, Smhg, FEighOrSvd::DECOMP_Svd, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
// //          dx.array() *= hw.array();
//    }
//    
//    
//    ProjectRotationsToOrbits(dx);
// }
#endif // INCLUDE_ABANDONED

void FSphereGridOptContext::ComputeStep(FScalar &CurMaxStep, FDenseVector &dx, FDenseVector const &residual, FDenseMatrix const &Jacobian, size_t iIt, FScalar fMaxStep_DynamicFactor)
{
   // note: 'Jacobian' is the (nTargetFn, nVars)-shaped matrix called
   // dK_dp_alpha at the place it is made and in the paper.
   
   // note: code now assumes that FTargetFnList does approximate-orthogonality
   // transform by itself (in terms of computing K and/or t). For this reason we
   // do not have any Scd-like stuff in here anymore. Basically, the target
   // functions are assumed to already be more-or-less orthogonal as they are.
   CurMaxStep = 0.0002;
   if (m_Options.AvoidSvds) {
      FDenseMatrix
         JtJ = Jacobian.transpose() * Jacobian;
      FDenseVector
         Jtg = Jacobian.transpose() * residual;
      MakeStep(dx, JtJ, Jtg, FEighOrSvd::DECOMP_Eigh, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
   } else {
      MakeStep(dx, Jacobian, residual, FEighOrSvd::DECOMP_Svd, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
   }
//    if (1) {
//       // S = L L^T
//       // -> S^{-1} = (L L^T)^{-1} = (L^T)^{-1} L^{-1}
//       FDenseMatrix
//          SmhJ = m_pScd->HalfSolve1(Jacobian);
//       FDenseMatrix
//          JtSmJ = SmhJ.transpose() * SmhJ;
//       FDenseVector
//          Smhg = m_pScd->HalfSolve1(residual),
//          JtSmg = SmhJ.transpose() * Smhg;
//       MakeStep(dx, JtSmJ, JtSmg, FEighOrSvd::DECOMP_Eigh, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
//       // ^- this should give the same result as the Scd.solve variant above
//       //
//       // Note: what this does is bypassing the SVD of J by computing the spectral decomposition of J^T J
//       // instead. That should give the squared singular values as eigenvalues.
//       // but isn't the SVD normally numerically better? to get identical results one probably has
//       // to update the sigma+shift vs sigma^2+shift values, though.
//       //
//       // Also, we have many less free parameters than target functions. So the
//       // eigh actually should be quite a bit more efficient
//    } else {
//       // Here: do the equivalent of:
//       //
//       //   dx = diag(h(w)) solve(S^{-1/2} J diag(h(w)), S^{-1/2} e),
//       //
//       // -> so only should use one side of L, not a full Cholesky solve.
//       FDenseMatrix
//          SmhJ = m_pScd->HalfSolve1(Jacobian);
//       FDenseVector
//          Smhg = m_pScd->HalfSolve1(residual);
//       MakeStep(dx, SmhJ, Smhg, FEighOrSvd::DECOMP_Svd, m_nFreeParameters, CurMaxStep, iIt, m_nPointsTotal, fMaxStep_DynamicFactor, m_Options);
//    }
   
   
   ProjectRotationsToOrbits(dx);
}


void FSphereGridOptContext::UpdateSeedPositions(FDenseVector const &dx, FScalar const &fResidual, size_t iIt)
{
   // assemble the actual seed point parameter updates from the step vector
   // 'dx' containing the 3*nVar rotation parameters. I.e., make and apply
   // the corresponding rotation matrices.
   FScalar
      fUpdateFactor = FScalar(1.) - cos(0.2*FScalar(iIt))*m_Options.pStepOptions->Damping;
   if (fResidual < 1e-6)
      fUpdateFactor = 1.;
   FScalar
      // the computed steps should preserve the orbit type. Check & accumulate
      // if this is indeed the case
      // (if not, in the past it was often a result of either some sorts of
      // numerics exploding or of mistakes in the treatment of the point group
      // symmetry; e.g., not having a complete symmetry equivalent set of target
      // functions)
      fOrbitViolation = FScalar(0);
   for (size_t iSeed = 0; iSeed < m_nSeeds; ++ iSeed) {
      // TODO: maybe hack a simple DIIS as option?
      FVector3
         // collect rotation axis & angles p_alpha for current seed point
         // (alpha\in\{x,y,z\})
         dx_sq = FVector3(dx[3*iSeed], dx[3*iSeed+1], dx[3*iSeed+2]);
//       // trace if this rotation would conflict with the seed point's current
//       // position: there should not be an overlap of the rotation axis with the
//       // point's position, as this rotation would not do anything for a point.
//       // So if we got a corresponding component during our optimization,
//       // probably something fishy happened.
//       fOrbitViolation += m_PointGroup.ProjectPointOrStepToOrbitType(
//          0, &dx_sq, m_SeedPoints.col(iSeed), *m_OrbitTypes[iSeed], false);
//       // ^- 'false': do not allow discrete point mappings into fundamental domain.

      if (dx_sq.squaredNorm() == 0)
         continue;

      FRotationGenerator
         // combine seed point #i's rotation angles with corresponding axis-
         // aligned anti-symmetric rotation generators, to make our aligned
         // anti-symmetric rotation generator matrix.
         Ai = fUpdateFactor * (m_A[0] * dx_sq[0] + m_A[1] * dx_sq[1] + m_A[2] * dx_sq[2]);
      FRotationGenerator
         // matrix-exponentiate or CayleyU-it to obtain the corresponding
         // rotation matrix.
         Ui = BuildRotation(Ai, 1.0);
      // ...and apply the resulting rotation to move seed point #i to its
      // new location on the sphere
      m_SeedPoints.col(iSeed) = Ui * m_SeedPoints.col(iSeed);
   }
   if (fOrbitViolation > sqr(g_ThrAlmostZero)) {
      std::cout << fmt::format("      ^- violation of <p,R> = 0 condition: {:8.2e}\n", double(sqrt(fOrbitViolation)));
   }
}


FScalar FSphereGridOptContext::FinishAndExportData()
{
   ct::FTimer
      tFinishAndExport;
   double
      tVerifyRule(0);
   if (1) {
      // map seed points onto positions exactly compatible with orbit type. This
      // time we allow discrete mappings of the points into the fundamental
      // domain (in case they have moved out) in order to obtain a canonical
      // representative of the seed point.
      // FIXME: canonical representative selection not fully implemented yet.
      // It is actually not complicated... generators for limited-dof groups now
      // all work with seed points which are either [x,y,0] or [x,y,y]. We could
      // just iterate over the (full) group orbit to fix the seed into one of those.
      for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
         FPoint sq = m_SeedPoints.col(iSeed);
         m_OrbitTypes[iSeed]->ProjectToOrbit(sq, true);
         m_SeedPoints.col(iSeed) = sq;
      }
   }
   if (1) {
      size_t
         nTotalPoints = 0;
      FScalar
         fMinWt = 1e99,
         fMaxWt = -1e99;
      for (size_t iSeed = 0; iSeed != m_nSeeds; ++ iSeed) {
         FPoint sq = m_SeedPoints.col(iSeed);
         size_t nOrbitLength = m_PointGroup.CountOrbitLength(sq);
         FScalar EffectiveWt = (FScalar(60.)/nOrbitLength) * m_LastWeights[iSeed];
         fMinWt = std::min(fMinWt, EffectiveWt);
         fMaxWt = std::max(fMaxWt, EffectiveWt);
         std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} {:24.16e}   // {}, orbit length: {}",
            double(sq[0]), double(sq[1]), double(sq[2]), double(EffectiveWt), m_OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
         nTotalPoints += nOrbitLength;
      }
      FScalar
         fMaxErrAbs, fMaxErrRel;
      FScalar
         fThrExport =  1e4*m_Options.fThreshResidual();
      fThrExport = std::max(FScalar(1e-16), fThrExport);
      // ^- export rule if either the target accuracy has been reached, or the
      //    rule is numerically exact for double precisoin regardless of whether
      //    or not full convergence has been reached (could happen, e.g., if we
      //    arrived at place very close to full convergence, e.g., 1e-25 res,
      //    which is enough for the target rule, but at which convergence got
      //    stuck.
      tVerifyRule -= double(tFinishAndExport);
      bool brrked = !TestSphereGrid(&fMaxErrAbs, &fMaxErrRel, m_lmax, fThrExport, m_SeedPoints, m_LastWeights, m_PointGroup, m_Options.VerifyRuleOptions);
      tVerifyRule += double(tFinishAndExport);
#ifdef INCLUDE_ABANDONED
//       brrked = false; // FIXME: remove this. large grid test (pre-normalization)
//       std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, edof = {})\n", nTotalPoints, m_lmax, double(fMaxWt/fMinWt), double(m_FinalResidual), m_nExtraDofForWeights);
//       std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, mxe = {:8.2e}, ree = {:8.2e}, edof = {})\n", nTotalPoints, m_lmax, double(fMaxWt/fMinWt), double(m_FinalResidual), double(fMaxErrAbs), double(fMaxErrRel), m_nExtraDofForWeights);
#endif // INCLUDE_ABANDONED
      FScalar
         wtsp = fMaxWt/fMinWt;
      std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, mxe = {:8.2e}, edof = {})\n", nTotalPoints, m_lmax, double(wtsp), double(m_FinalResidual), double(fMaxErrAbs), m_nExtraDofForWeights);
      // replace formal residual by computed maximum absolute error over all monomials from l=0 to lmax.
      m_FinalResidual = fMaxErrAbs;
      if (brrked) {
         std::cout << fmt::format(" WARNING: this rule does NOT integrate l={} exactly (err = {:8.2e}, rel = {:8.2e})\n", m_lmax, double(fMaxErrAbs), double(fMaxErrRel));
      } else {
//          if (m_Converged && m_Options.pExportOptions.get() != 0)
         // ^- converged or not: we got here, and the rule *does* integrate
         // the target l exactly. we should export it.
         if (m_Options.pExportOptions.get() != 0) {
            FExportedGridStats
               ThisGrid,
               &LastGrid = m_LastExportedGridStats;
            ThisGrid.wtsp = wtsp;
            ThisGrid.res = m_FinalResidual;
            ThisGrid.lmax = m_lmax;
            
            if (!ThisGrid.IsWorseThan(LastGrid)) {
               if (ThisGrid.lmax == LastGrid.lmax && !ThisGrid.IsEquivalent(LastGrid)) {
                  std::cout << fmt::format(" NOTE: New grid (lmax={}, wtsp={:.3f}, res = {:8.2e}) vs last exported grid (lmax={}, wtsp={:.3f}, res = {:8.2e}): Accepting new grid, exporting to '{}'\n", ThisGrid.lmax, double(ThisGrid.wtsp), double(ThisGrid.res), LastGrid.lmax, double(LastGrid.wtsp), double(LastGrid.res), m_Options.pExportOptions->FileName);
               }
               ExportGrid(m_SeedPoints, m_LastWeights, m_FinalResidual, m_lmax, *m_Options.pExportOptions, m_Options);
               m_LastExportedGridStats = ThisGrid;
               m_DidExportSomething = true;
            } else {
               if (ThisGrid.lmax == LastGrid.lmax) {
                  std::cout << fmt::format(" NOTE: New grid (lmax={}, wtsp={:.3f}, res = {:8.2e}) vs last exported grid (lmax={}, wtsp={:.3f}, res = {:8.2e}): New grid deemed a regression. Rejecting new grid. Will not overwrite '{}'\n", ThisGrid.lmax, double(ThisGrid.wtsp), double(ThisGrid.res), LastGrid.lmax, double(LastGrid.wtsp), double(LastGrid.res), m_Options.pExportOptions->FileName);
               }
            }
         }
      }
   }
   if (1) {
      io.WriteLine();
      char const *pVerifyType = "mode: unrecognized";
      if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_AllMonomials)
         pVerifyType = "M/all monomials";
      else if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_SymmetryMonomialsOnly)
         pVerifyType = "M/symmetry-unique";
      else if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_MinimalMonomialSetOnly)
         pVerifyType = "M/minimal set";
      io.WriteTiming("rule cross-verification", double(tVerifyRule), pVerifyType);
      io.WriteTiming("optimization step", double(m_tThisGridTotal));
#ifdef INCLUDE_ABANDONED
//       io.WriteTiming("optimization step", double(m_tThisGridTotal));
//       io.WriteTiming("optimization total", double(m_tThisGridTotal), "total");
//       char const *pVerifyType = "mode: unrecognized";
//       if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_AllMonomials)
//          pVerifyType = "mode: test all monomials";
//       else if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_SymmetryMonomialsOnly)
//          pVerifyType = "mode: test symmetry-adapted";
//       else if (m_Options.VerifyRuleOptions == FGridSearchOptions::VERIFY_MinimalMonomialSetOnly)
//          pVerifyType = "mode: test minimal-sufficient-set";
#endif // INCLUDE_ABANDONED
   }
   return m_FinalResidual;
   // if we cannot exactly integrate all monomials for lmax and lmax - 1 (and thereby also all lower ones),
   // we'd probably want to exactly integrate all l with l <= lmax -2 and lstsq the remaining higher order ones...
}


#ifdef INCLUDE_ABANDONED
// FScalar OptimizeSphereGrid(FPointArray &SeedPoints, size_t lmin, size_t lmax, FGridSearchOptions const &Options)
// {
//    FPointGroup const
//       &PointGroup = *Options.pPointGroup;
//    FRotationMatrixList const
//       &SymmetryOps = PointGroup.GetOps();
//    size_t
//       nSeeds = SeedPoints.cols();
//    int
//       iPrintLevel = Options.iPrintLevel;
// 
//    std::cout << fmt::format("\n** Processing {} grid ({} symmetry-ops) and {} monomials (lmin={}, lmax={}) for {} seed points", PointGroup.Name(NAMETYPE_Symmetry), SymmetryOps.size(), MakeMonomialList(lmin,lmax).size(), lmin, lmax, nSeeds) << std::endl;
// 
//    std::cout << "\n Normalizing and sorting initial seeds..." << std::endl;
//    PrepareInitialSeeds(SeedPoints, PointGroup);
// 
// 
//    std::cout << " Classifying orbits and free parameters..." << std::endl;
//    FOrbitTypeRawPtrList
//       OrbitTypes = PointGroup.GetOrbitTypes(SeedPoints);
//    size_t
//       // total number of non-redundant free parameters to optimize
//       // (these describe only the positions of seed points; the weights are implicit)
//       nFreeParameters = 0,
//       nPointsTotal = 0;
//    for (size_t iSeed = 0; iSeed != nSeeds; ++iSeed) {
//       nFreeParameters += OrbitTypes[iSeed]->nFreeParameters;
//       nPointsTotal += OrbitTypes[iSeed]->Length;
//    }
//    io.WriteLine();
//    io.WriteCount("Number of points", nPointsTotal, fmt::format("in {} orbits", OrbitTypes.size()));
//    io.WriteCount("Number of free parameters", nFreeParameters);
// 
//    // Compute set of target functions which our integration rule is to evaluate
//    // exactly, and compute their overlap matrix and a suitable decomposition of
//    // it to solve conditioning equations with.
//    FMonomialList
//       Monomials = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, Options.iPrintLevel, &Options);
//    size_t
//       nMonomials = Monomials.size();
//    io.WriteCount("Number of monomials", nMonomials, "symmetry unique");
// 
//    FLinearSolverPtr
//       pScd;
//    {
//       FDenseMatrix
//          S = MakeOverlapMatrix(Monomials);
//       if (g_UseGroupSymmetrizedFn) {
//          // Group-averaging causes hard singularities in the function space,
//          // because the parts of the functions which lie outside the space
//          // projected to by
//          //
//          //    P_G = \sum_g \sum_j |Rg m_j><Rg m_j|
//          //
//          // are gone. Now, *mathematically* they are *not needed* for the
//          // integration rule optimization (see paper).
//          //
//          // However, their absence also causes some problems. In particular,
//          // the Choleksy solver can't deal with this---this needs a rank-
//          // revealing solver which can explicitly determine the non-singular
//          // space.
//          // So use eigh-solver by default in this case.
//          pScd = FLinearSolverPtr(new FSymmetricPositiveSolverEigh(S));
//          io.WriteCount("Number of target fn", pScd->nEffectiveRank(), "linearly independent");
//       } else {
//          // if not symmetrizing, we will get a wide spectrum of S eigenvalues;
//          // the overlap matrix will have a large condition number for larger L.
//          // However, none of the functions will be strictly linearly dependent,
//          // and so the entire space must be handled accurately. Choleksy
//          // decompositions with triangular Gaussian elimination tend to do this
//          // rather well---better than the Smh version certainly.
//          // So use Scd-solver by default in this case.
//          pScd = FLinearSolverPtr(new FSymmetricPositiveSolverCd(S));
//       }
//    }
//    
// //    std::cout << " Building rhs vector..." << std::endl;
//    FDenseVector
//       // vector of exact sphere integrals over the monomials
//       rhs(nMonomials),
//       last_weights;
// //    std::cout << fmt::format("\n** Processing {} symmetry-ops and {} symmetry-unique monomials (lmin={}, lmax={}) for {} seed points", SymmetryOps.size(), nMonomials, lmin, lmax, nSeeds) << std::endl;
//    for (size_t i = 0; i != nMonomials; ++ i)
//       rhs[i] = Monomials[i].CalcUnitSphereIntegral();
//    FScalar
//       fMinRhs = rhs.maxCoeff();
//    for (size_t i = 0; i != nMonomials; ++ i)
//       if (rhs[i] != 0)
//          fMinRhs = std::min(fMinRhs, rhs[i]);
//    std::cout << fmt::format("\n Built rhs vector:  MIN: {:8.2e}  MAX: {:8.2e}", double(fMinRhs), double(rhs.maxCoeff())) << std::endl;
//    
//    FScalar FinalResidual = 1e6;
// 
//    typedef Eigen::Matrix<FScalar,3,3>
//       FRotationGenerator; // rotation generators.
//    FRotationGenerator
//       // generators of rotation around x,y,z axis
//       A[3];
//    A[0] << 0,0,0, 0,0,1, 0,-1,0;
//    A[1] << 0,0,-1, 0,0,0, 1,0,0;
//    A[2] << 0,1,0, -1,0,0, 0,0,0;
// 
//    bool Converged = false;
//    FScalar fMinRes = 1e99; // minimum residual encountered during the iterations.
//    ptrdiff_t
//       nExtraDofForWeights = 0;
//    ptrdiff_t
//       iIt;
//    std::cout << fmt::format("\n   {:^6} {:^8}   {:8^}     {:8^} |  WEIGHTS", "ITER.", "RESIDUAL", "UPDATE", "MAX.STEP") << std::endl;
//    for (iIt = 0; iIt < ptrdiff_t(Options.MaxIt); ++ iIt) {
//       if (iPrintLevel >= 2) {
//          std::cout << "      Checking for seed collapse. Current points:\n";
//          for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//             FPoint sq = SeedPoints.col(iSeed);
//             size_t nOrbitLength = PointGroup.CountOrbitLength(sq);
//             std::cout << fmt::format("        {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
//                double(sq[0]), double(sq[1]), double(sq[2]), OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//          }
//       }
//       CheckSeedsForCollapse(SeedPoints, PointGroup);
//       size_t
//          nVars = 3*nSeeds;
//       FDenseMatrix
//          // K[i,q]: (i: monomial index, q: seed index)
//          K(nMonomials, nSeeds),
//          // dKw[i,p] = \sum_{q} dK[i,q]/dp w[q] -- kernel derivative contracted to fixed weights.
//          dKw(nMonomials, nVars),
//          // dKr[q,p] = \sum_{i} dK[i,q]/dp r[i] -- kernel derivative contracted to fixed rhs r[i]
//          dKr(nSeeds, nVars);
// 
//       FTransformedPointInfo
//          // contains information about seed points transformed by current symmetry operation.
//          // In particular, powers of its x,y,z coordinates.
//          TrafodSeeds(SeedPoints, lmax);
//       K.setZero();
//       for (size_t iSymOp = 0; iSymOp != SymmetryOps.size(); ++ iSymOp) {
//          TrafodSeeds.TransformPoints(SymmetryOps[iSymOp]);
//          for (size_t iMonomial = 0; iMonomial != nMonomials; ++ iMonomial) {
//             for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//                K(iMonomial, iSeed) += Monomials[iMonomial].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
//             }
//          }
//       }
// 
// //       FDenseMatrix
// //          SmK = Scd.solve(K);
// //       
// //       FDenseMatrix
// //          // [K^T S^{-1} K]_{q_1, q_2}
// //          KtK(nSeeds, nSeeds);
// //       KtK = SmK.transpose() * K;
// //       Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
// //          KtKd(KtK, Eigen::ComputeThinU | Eigen::ComputeThinV);
// //       // ^- FIXME: try doing this with CD instead of SVD, too?
// //       KtKd.setThreshold(FScalar(1e-16Q));
// //       // ^- note: for REALLY small target ls we might have a KtK which is not
// //       //    positive definite... so we use SVD instead of LLT decomposition
// //       //    (eigenvalues might also do the trick; maybe consider later.)
// //       
// //       // count number of degrees of freedom remaining in space of weigthts (for current seedpoints)
// //       // which are *not* fixed by the demand that the residual is minimal for the current points for the current lmax.
// //       nExtraDofForWeights = ptrdiff_t(nSeeds);
// //       FDenseVector const
// //          &KtK_sigma = KtKd.singularValues();
// //       for (ptrdiff_t i = 0; i != ptrdiff_t(KtK_sigma.size()); ++ i)
// //          if (KtK_sigma[i] > KtKd.threshold())
// //             nExtraDofForWeights -= 1;
// // 
// //       FDenseVector
// //          // dimension: q (seeds)
// //          weights = KtKd.solve(SmK.transpose() * rhs),
// //          // dimension: i (monomial)
// //          drhs = Scd.solve(rhs - K * weights);
// //       last_weights = weights;
// 
//       
//       FDenseMatrix
// //          SmhK = Scd.matrixL().solve(K);
//          SmhK = pScd->HalfSolve1(K);
//       FDenseMatrix
// //          SmK = Scd.solve(K);
// //          SmK = Scd.matrixL().transpose().solve(SmhK);
//          SmK = pScd->HalfSolve2(SmhK);
//       
//       FDenseMatrix
//          // [K^T S^{-1} K]_{q_1, q_2}
//          KtK(nSeeds, nSeeds);
//       KtK = SmhK.transpose() * SmhK;
// //       Eigen::LLT<FDenseMatrix>
// //          KtKd(KtK); // KtK decomposition (Cholesky or SVD)
// //       Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
// //          KtKd(KtK, Eigen::ComputeThinU | Eigen::ComputeThinV);
// //       KtKd.setThreshold(FScalar(1e-16Q)); // <- hm.. turning this on made things worse (but changing it to 1e-25Q didn't make things better---gives results identical to LDLT)
//       Eigen::LDLT<FDenseMatrix>
//          KtKd(KtK); // KtK decomposition (Cholesky or SVD; LDLT is a pivoted Cholesky decomposition)
//       // ^- using a Cholesky decomp for this is very optimistic in some
//       //    situations... especailly for small L, large rules, and exact
//       //    function symmetrization, the system KtK is often *meant* to be
//       //    singular. It is actually quite surprising that Eigen can deal with
//       //    this to smaller grids. However, for larger grid scans, such as
//       //    ico-13-1, eventually LLT explodes. Interestingly, just replacing
//       //    this by LDLT (pivoted Cholesky) appears to work fine as a
//       //    workaround, although, realistically speaking, the system is still
//       //    supposed to be singular.
//       //  TODO:
//       //  - If this blows up again, like it should, replace the fixed decomp scheme
//       //    by something which first attempts to pre-screen the Cholesky decomp to
//       //    see if it is positive definite, and otherwise uses a spectral decomp for
//       //    this (should be better than SVD for this kind of thing, KtK is symmetric).
//       //  - Make a separate SymmetricMatrixDecomp or whatever to hold it.
//       nExtraDofForWeights = 0;
// 
// 
//       FDenseVector
//          // dimension: q (seeds)
// //          weights = KtKd.solve(SmhK.transpose() * Scd.matrixL().solve(rhs)),
//          weights = KtKd.solve(SmhK.transpose() * pScd->HalfSolve1(rhs)),
// //          weights = KtKd.solve(SmK.transpose() * rhs),
//          // dimension: i (monomial)
// //          drhs = Scd.solve(rhs - K * weights);
//          drhs = pScd->Solve(rhs - K * weights);
//       last_weights = weights;
// 
//       dKw.setZero();
//       dKr.setZero();
// 
// 
//       for (size_t iSymOp = 0; iSymOp != SymmetryOps.size(); ++ iSymOp) {
//          FRotationMatrix const
//             &Rg = SymmetryOps[iSymOp],
//             Rg_Ax = Rg * A[0],
//             Rg_Ay = Rg * A[1],
//             Rg_Az = Rg * A[2];
//          FPointArray
//             Rg_Ax_pts = Rg_Ax * SeedPoints,
//             Rg_Ay_pts = Rg_Ay * SeedPoints,
//             Rg_Az_pts = Rg_Az * SeedPoints;
//          TrafodSeeds.TransformPoints(SymmetryOps[iSymOp]);
// 
//          // project out rotations which would lead out of the given orbit type.
// //          for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed)
// //             PointGroup.ProjectRotationToOrbitType(Rg_Ax_pts.col(iSeed), Rg_Ay_pts.col(iSeed), Rg_Az_pts.col(iSeed), TrafodSeeds.TransformedPoints.col(iSeed), *OrbitTypes[iSeed]);
//          PointGroup.ProjectRotationToOrbitType(Rg_Ax_pts, Rg_Ay_pts, Rg_Az_pts, TrafodSeeds.TransformedPoints, &OrbitTypes[0]);
// 
//          // evaluate contracted derivatives of K for current group operation.
//          for (size_t iMonomial = 0; iMonomial != nMonomials; ++ iMonomial) {
//             for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//                FPoint
//                   dMdXyz,
//                   Rg_Ax_sq = Rg_Ax_pts.col(iSeed),
//                   Rg_Ay_sq = Rg_Ay_pts.col(iSeed),
//                   Rg_Az_sq = Rg_Az_pts.col(iSeed);
//                Monomials[iMonomial].EvalPoint(dMdXyz, TrafodSeeds.TransformedPointInfo(iSeed));
// //                Monomials[iMonomial].EvalPoint(dMdXyz, Rg*SeedPoints.col(iSeed).matrix());
//                FPoint
//                   dKiq_dp_alpha(0,0,0);
//                for (size_t beta = 0; beta != 3; ++ beta) {
//                   dKiq_dp_alpha[0] += dMdXyz[beta] * Rg_Ax_sq[beta];
//                   dKiq_dp_alpha[1] += dMdXyz[beta] * Rg_Ay_sq[beta];
//                   dKiq_dp_alpha[2] += dMdXyz[beta] * Rg_Az_sq[beta];
//                }
//                size_t iVar = 3*iSeed;
//                for (size_t beta = 0; beta != 3; ++ beta) {
//                   dKw(iMonomial, iVar + beta) += dKiq_dp_alpha[beta] * weights[iSeed];
//                   dKr(iSeed, iVar + beta) += dKiq_dp_alpha[beta] * drhs[iMonomial];
//                }
//             }
//          }
//       }
// 
//       FDenseMatrix
//          dweight_dp_alpha = KtKd.solve(dKr - SmK.transpose() * dKw),
//          // this is the final jacobian!
//          // format: i (residual), p_alpha (nVars rotation parameters)
//          dK_dp_alpha = dKw + K * dweight_dp_alpha;
// 
//       FDenseVector
//          residual = K * weights - rhs;
//       if (iPrintLevel >= 4) {
//          std::cout << " Intermediate residual contributions for current points:\n";
//          for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//             FPoint sq = SeedPoints.col(iSeed);
//             size_t nOrbitLength = PointGroup.CountOrbitLength(sq);
//             std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} // {}, orbit length: {}",
//                double(sq[0]), double(sq[1]), double(sq[2]), OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//             for (size_t iMonomial = 0; iMonomial != nMonomials; ++ iMonomial) {
//                FMonomial m = Monomials[iMonomial];
//                std::cout << fmt::format("     |  [{:>6}]  m = {:<20}  K[m,r] = {:24.16e}   D[K(m),p_alpha] = {:24.16e} {:24.16e} {:24.16e}   residual[m] = {:16.8e}\n",
//                   iMonomial, fmt::format("[{},{},{}]",m.ix, m.iy, m.iz), double(K(iMonomial,iSeed)), double(dK_dp_alpha(iMonomial,3*iSeed+0)), double(dK_dp_alpha(iMonomial,3*iSeed+1)), double(dK_dp_alpha(iMonomial,3*iSeed+2)), double(residual[iMonomial]));
//             }
//          }
//       }
// 
//       FDenseVector
//          dx;
//       FScalar
//          CurMaxStep = 0.0002;
//       if (0) {
// //          MakeStep(dx, dK_dp_alpha, residual, CurMaxStep, iIt, fMinRes, nFreeParameters, Options);
//          throw std::runtime_error("this residual update path is not (anymore) implemented.");
//       } else if (0) {
//          FDenseMatrix
//             JtJ = dK_dp_alpha.transpose() * dK_dp_alpha;
//          FDenseVector
//             Jtg = dK_dp_alpha.transpose() * residual;
//          // FIXME: should I put in something like diag(1/norm(f_i)) before multiplying with Jt?
//          // Or is this already happening implicitly? (program says "Built rhs vector:  MIN: 1.00e+00  MAX: 3.00e+00")
//          // Or maybe even explicitly orthogonalize the functions, by multiplying with S^{-1} before J^T?
//          MakeStepSymmetric(dx, JtJ, Jtg, CurMaxStep, iIt, fMinRes, nFreeParameters, nPointsTotal, Options);
//       } else if (0) {
//          std::cout << "   ....ran through J^T S^{-1} J path." << std::endl;
//          FDenseMatrix
// //             SmJ = Scd.solve(dK_dp_alpha);
//             SmJ = pScd->Solve(dK_dp_alpha);
// //          // remove degeneracy factors from column projector -- not what it is for.
// //          for (size_t iFn = 0; iFn != nMonomials; ++ iFn)
// //             SmJ.row(iFn) /= FScalar(Monomials[iFn].degeneracy);
//          // ^- hm... nope, that is now how it works (and neither is the symmetric variant below).
//          FDenseMatrix
//             JtSmJ = SmJ.transpose() * dK_dp_alpha;
//          FDenseVector
//             JtSmg = SmJ.transpose() * residual;
//          MakeStepSymmetric(dx, JtSmJ, JtSmg, CurMaxStep, iIt, fMinRes, nFreeParameters, nPointsTotal, Options);
//       } else if (1) {
//          // S = L L^T
//          // -> S^{-1} = (L L^T)^{-1} = (L^T)^{-1} L^{-1}
//          FDenseMatrix
// //             SmhJ = Scd.matrixL().solve(dK_dp_alpha);
//             SmhJ = pScd->HalfSolve1(dK_dp_alpha);
// //          for (size_t iFn = 0; iFn != nMonomials; ++ iFn)
// //             SmhJ.row(iFn) /= sqrt(FScalar(Monomials[iFn].degeneracy));
//          FDenseMatrix
//             JtSmJ = SmhJ.transpose() * SmhJ;
//          FDenseVector
// //             Smhg = Scd.matrixL().solve(residual),
//             Smhg = pScd->HalfSolve1(residual),
//             JtSmg = SmhJ.transpose() * Smhg;
// //          // remove degeneracy factors from column projector
// //          for (size_t iFn = 0; iFn != nMonomials; ++ iFn) {
// // //             JtSmJ.row(iFn) /= FScalar(Monomials[iFn].degeneracy);
// //             JtSmg[iFn] /= sqrt(FScalar(Monomials[iFn].degeneracy));
// //          }
// //          MakeStep(dx, SmhJ, Smhg, CurMaxStep, iIt, fMinRes, nFreeParameters, Options);
//          MakeStepSymmetric(dx, JtSmJ, JtSmg, CurMaxStep, iIt, fMinRes, nFreeParameters, nPointsTotal, Options);
//          // ^- this should give the same result as the Scd.solve variant above
//          //
//          // Note: what this does is bypassing the SVD of J by computing the spectral decomposition of J^T J
//          // instead. That should give the squared singular values as eigenvalues.
//          // but isn't the SVD normally numerically better? to get identical results one probably has
//          // to update the sigma+shift vs sigma^2+shift values, though.
//          //
//          // Also, we have many less free parameters than target functions. So the
//          // eigh actually should be quite a bit more efficient
//       } else {
//          // Here: do the equivalent of:
//          //
//          //   dx = diag(h(w)) solve(S^{-1/2} J diag(h(w)), S^{-1/2} e),
//          //
//          // -> so only should use one side of L, not a full Cholesky solve.
//          FDenseMatrix
// //             SmhJ = Scd.matrixL().solve(dK_dp_alpha);
//             SmhJ = pScd->HalfSolve1(dK_dp_alpha);
//          FDenseVector
// //             Smhg = Scd.matrixL().solve(residual);
//             Smhg = pScd->HalfSolve1(residual);
// //          FDenseVector
// //             hw(SmhJ.cols());
// //          assert(size_t(weights.size()) == size_t(nSeeds) && hw.size() == weights.size() * 3);
// //          for (int iSeed = 0; iSeed != weights.size(); ++ iSeed) {
// //             FScalar
// //                wi = weights[iSeed],
// //                hwi = FScalar(0);
// // //             if (abs(wi) > 1e-10 / (3*nSeeds))
// // // //                hwi = 1/sqrt(abs(wi));
// // //                hwi = 1/abs(wi);
// // //             else
// // //                hwi = 1.;
// //             hwi = sqrt(wi);
// // //             hwi = wi;
// //             for (int alpha = 0; alpha != 3; ++ alpha) {
// //                hw[3*iSeed+alpha] = hwi;
// //             }
// //          }
// // //          SmhJ.array().colwise() *= hw.array();
// //          for (int i = 0; i != SmhJ.cols(); ++ i)
// //             SmhJ.col(i) *= hw[i];
//          // ^- hm... none of those dynamic weight scaling thingies work. I think
//          // scaling with 1/sqrt(wi) *should* be the right variant, but it does not work either.
//          // Guess I'll just leave it alone.
//          MakeStep(dx, SmhJ, Smhg, CurMaxStep, iIt, fMinRes, nFreeParameters, nPointsTotal, Options);
// //          std::cout << fmt::format("   hw scale: dx.size() = {}  hw.size() = {}  weights.size() = {}  hw[0] = {}\n", dx.size(), hw.size(), weights.size(), hw[0]);
// //          dx.array() *= hw.array();
//       }
// 
//       FScalar
// //          fResidual = Rmsd(residual, nFreeParameters),
// //          fResidual = Rmsd(Scd.matrixL().solve(residual), nFreeParameters),
//          fResidual = Rmsd(pScd->HalfSolve1(residual), nFreeParameters),
//          // ^- compute residual wrt. orthogonalized function set
//          fStep = Rmsd(dx, nFreeParameters);
//       FinalResidual = fResidual;
//       std::cout << fmt::format(" {:6d}   {:8.2e}   {:8.2e}   {:8.2e} | ", iIt+1, double(fResidual), double(fStep), double(CurMaxStep));
//       for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed)
//          std::cout << fmt::format(" {:22.16f}", double(weights[iSeed]));
//       std::cout << std::endl;
// 
// 
//       if ((fStep < Options.fThreshResidual() && iIt >= g_nDampedSigSteps) || fResidual < Options.fThreshResidual()) {
//          // stop after NEXT iteration.
//          Converged = true;
//       }
// 
//       if (Converged)
//          break;
// 
//       FScalar
//          fUpdateFactor = FScalar(1.) - Options.Damping;
//       if (fResidual < 1e-6)
//          fUpdateFactor = 1.;
//       FScalar
//          fViolation = FScalar(0);
//       for (size_t iSeed = 0; iSeed < nSeeds; ++ iSeed) {
//          // TODO: maybe hack a simple DIIS as option?
//          FVector3
//             dx_sq = FVector3(dx[3*iSeed], dx[3*iSeed+1], dx[3*iSeed+2]);
//          fViolation += PointGroup.ProjectPointOrStepToOrbitType(0, &dx_sq, SeedPoints.col(iSeed), *OrbitTypes[iSeed]);
// 
//          if (dx_sq.squaredNorm() == 0)
//             continue;
// 
//          FRotationGenerator
//             Ai = fUpdateFactor * (A[0] * dx_sq[0] + A[1] * dx_sq[1] + A[2] * dx_sq[2]);
//          FRotationGenerator
//             Ui = BuildRotation(Ai, 1.0);
//          SeedPoints.col(iSeed) = Ui * SeedPoints.col(iSeed);
//       }
//       if (fViolation > 1e-40Q) {
//          std::cout << fmt::format("      ^- violation of <p,R> = 0 condition: {:8.2e}\n", double(sqrt(fViolation)));
//       }
//       
//    }
// 
//    if (Converged) {
//       std::cout << "\n Converged. Final seed points & weights (per unique point):\n";
//    } else {
//       std::cout << fmt::format("\n WARNING: FAILED to converge in {} iterations. Aborted. Final seed points & weights (per unique point):\n", Options.MaxIt);
//    }
//    if (1) {
//       size_t
//          nTotalPoints = 0;
//       FScalar
//          fMinWt = 1e99,
//          fMaxWt = -1e99;
//       for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
//          FPoint sq = SeedPoints.col(iSeed);
//          size_t nOrbitLength = PointGroup.CountOrbitLength(sq);
//          FScalar EffectiveWt = (FScalar(60.)/nOrbitLength) * last_weights[iSeed];
//          fMinWt = std::min(fMinWt, EffectiveWt);
//          fMaxWt = std::max(fMaxWt, EffectiveWt);
//          std::cout << fmt::format("   {:24.16e} {:24.16e} {:24.16e} {:24.16e}   // {}, orbit length: {}",
//             double(sq[0]), double(sq[1]), double(sq[2]), double(EffectiveWt), OrbitTypes[iSeed]->pDesc, nOrbitLength) << std::endl;
//          nTotalPoints += nOrbitLength;
//       }
//       FScalar
//          fMaxErrAbs, fMaxErrRel;
//       bool brrked = !TestSphereGrid(&fMaxErrAbs, &fMaxErrRel, lmax, 1e4*Options.fThreshResidual(), SeedPoints, last_weights, PointGroup);
// //       brrked = false; // FIXME: remove this. large grid test (pre-normalization)
// //       std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, edof = {})\n", nTotalPoints, lmax, double(fMaxWt/fMinWt), double(FinalResidual), nExtraDofForWeights);
// //       std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, mxe = {:8.2e}, ree = {:8.2e}, edof = {})\n", nTotalPoints, lmax, double(fMaxWt/fMinWt), double(FinalResidual), double(fMaxErrAbs), double(fMaxErrRel), nExtraDofForWeights);
//       std::cout << fmt::format(" Total number of points in integration rule: {}  (lmax = {}, wtsp = {:.3f}, res = {:8.2e}, mxe = {:8.2e}, edof = {})\n", nTotalPoints, lmax, double(fMaxWt/fMinWt), double(FinalResidual), double(fMaxErrAbs), nExtraDofForWeights);
//       // replace formal residual by computed maximum absolute error over all monomials from l=0 to lmax.
//       FinalResidual = fMaxErrAbs;
//       if (brrked) {
//          std::cout << fmt::format(" WARNING: this rule does NOT integrate l={} exactly (err = {:8.2e}, rel = {:8.2e})\n", lmax, double(fMaxErrAbs), double(fMaxErrRel));
//       } else {
//          if (Converged && Options.pExportOptions.get() != 0)
//             ExportGrid(SeedPoints, last_weights, FinalResidual, lmax, *Options.pExportOptions, Options);
//       }
//    }
//    return FinalResidual;
// 
//    // if we cannot exactly integrate all monomials for lmax and lmax - 1 (and thereby also all lower ones),
//    // we'd probably want to exactly integrate all l with l <= lmax -2 and lstsq the remaining higher order ones...
// }
#endif


void ExportGrid(FPointArray const &SeedPoints_, FDenseVector const &LastWeights_, FScalar FinalResidual,
   int lmax, FExportOptions const &ExportOptions, FGridSearchOptions const &Options)
{
   FPointGroup const
      &PointGroup = *Options.pPointGroup;
   FPointArray
      SeedPoints = SeedPoints_;
   size_t
      nSeeds = SeedPoints.cols();
   PrepareInitialSeeds(SeedPoints, PointGroup);
   FOrbitTypeRawPtrList
      // make sure this comes AFTER sorting the points!!
      OrbitTypes = PointGroup.GetOrbitTypes(SeedPoints);
   FDenseVector
      Weights = LastWeights_;

   size_t
      nTotalPoints = 0;
   FScalar
      fMinWt = 1e99,
      fMaxWt = -1e99;
   for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
      FPoint sq = SeedPoints.col(iSeed);
      size_t nOrbitLength = PointGroup.CountOrbitLength(sq);
      // input weights are made with respect to a full orbit of the seed point. For
      // special points, this may mean that they occur multiple times in identical positions
      // (e.g., there are only 12 vertices in an icosahedron, in the full size 60 orbit, a seed
      // on a vertex occurs 5 times -> make only one point, and multiply its weight by 5)
      FScalar nRedundancyOfFullOrbit = FScalar(PointGroup.GetOps().size()) / FScalar(nOrbitLength);
      Weights[iSeed] = nRedundancyOfFullOrbit * LastWeights_[iSeed];
      fMinWt = std::min(fMinWt, Weights[iSeed]);
      fMaxWt = std::max(fMaxWt, Weights[iSeed]);
      nTotalPoints += nOrbitLength;
   }

   std::ofstream
      out(ExportOptions.FileName.c_str(), (ExportOptions.AppendToFile? (std::ofstream::app | std::ofstream::ate) : std::ofstream::trunc) | std::ofstream::out);

   out << fmt::format("# {} integration rule: npts = {}, lmax = {}, wtsp = {:.3f}, res = {:8.2e}\n", PointGroup.Name(NAMETYPE_Symmetry), nTotalPoints, lmax, double(fMaxWt/fMinWt), double(FinalResidual));
   if (ExportOptions.StoreVersionInfo || ExportOptions.StoreIdentityInfo) {
      out << "# generated ";
      if (ExportOptions.StoreIdentityInfo) {
         out << MakeCurrentTimeString();
      }
      if (ExportOptions.StoreVersionInfo) {
         if (ExportOptions.StoreIdentityInfo)
            out << " ";
         out << "by " << MakeAiggVersionString();
      }
      out << "\n";
   }
   if (ExportOptions.StoreOptions)
      out << fmt::format("# options: '{}'\n", Options.MakeOptionString());
   
   // emit any extra comment lines we were explicitly asked to include in
   // the exported output file
   for (std::string const &it : ExportOptions.ExtraComments) {
      out << "# " << it << "\n";
   }

   for (size_t iSeed = 0; iSeed != nSeeds; ++ iSeed) {
      if (!out.good())
         throw std::runtime_error(fmt::format("I/O error during export of grid to file '{}'", ExportOptions.FileName.c_str()));
      FPoint
         sq = SeedPoints.col(iSeed);
      FScalar
         wt = Weights[iSeed];
      FPointList
         // total list of points we intend to export.
         pl;
      if (ExportOptions.PointsToExport == FExportOptions::EXPORT_SeedsAndWeightsOnly) {
         pl.push_back(sq);
      } else {
         assert(ExportOptions.PointsToExport == FExportOptions::EXPORT_FullOrbits);
         // make list of unique points.
         FPointCloud
            PointCloud;
         for (size_t iOp = 0; iOp != PointGroup.GetOps().size(); ++ iOp) {
            FPoint
               Rg_sq = PointGroup.GetOps()[iOp] * sq;
            FPointCloud::FPointEntry const
               &e = PointCloud.Insert(Rg_sq, iOp);
            if (e.iPoint != (ptrdiff_t)iOp) {
               // transformed point already exported. Don't do it again.
               // (note: for special orbits)
            } else {
               pl.push_back(Rg_sq);
            }
         }
         if (pl.size() != size_t(OrbitTypes[iSeed]->Length))
            throw std::runtime_error(fmt::format("programming error: orbit lengths inconsistent in export of grid (iseed={}, orbit-length: {}, expected: {}).", iSeed, pl.size(), OrbitTypes[iSeed]->Length));
      }
      for (size_t iImagePt = 0; iImagePt != pl.size(); ++ iImagePt) {
         FPoint const
            &Rg_sq = pl[iImagePt];
         std::string
            AdditionalText;
         if (ExportOptions.PointsToExport == FExportOptions::EXPORT_SeedsAndWeightsOnly)
            AdditionalText = fmt::format(" {:4d}   ({}-{})", OrbitTypes[iSeed]->Length, PointGroup.Name(NAMETYPE_Symmetry), OrbitTypes[iSeed]->pDesc);
//          out << fmt::format("{:24.16e} {:24.16e} {:24.16e} {:24.16e}{}\n",
//                double(Rg_sq[0]), double(Rg_sq[1]), double(Rg_sq[2]), double(wt), AdditionalText);
         int prec = ExportOptions.Precision;
         out << fmt::format("{} {} {} {}{}\n",
               FmtExportFloat(Rg_sq[0], prec), FmtExportFloat(Rg_sq[1], prec), FmtExportFloat(Rg_sq[2], prec), FmtExportFloat(wt, prec), AdditionalText);
      }
   }
   out.flush();
}



bool TestSphereGrid(FScalar *pMaxErrAbs, FScalar *pMaxErrRel, size_t lmax, FScalar ThreshAbs, FPointArray &SeedPoints, FDenseVector const &weights, FPointGroup const &PointGroup, FGridSearchOptions::FVerifyRuleOptions VerifyRuleOptions)
{
   cx_assert_rt(SeedPoints.cols() == weights.size());
   FRotationMatrixList const
      &SymmetryOps = PointGroup.GetOps();
   // collect list of monomials to test over
   FMonomialList
      Monomials;
   {
      size_t lmin = 0;
      if (VerifyRuleOptions == FGridSearchOptions::VERIFY_MinimalMonomialSetOnly) {
         if (lmax > 0)
            lmin = lmax - 1;
      }
      if (VerifyRuleOptions == FGridSearchOptions::VERIFY_AllMonomials) {
         Monomials = MakeMonomialList(lmin, lmax);
         // ^- note: that are *literally all* monomials, from 0 to lmax. Not
         // only symmetry- unique ones, and not only lmax and lmax-1 ones (which
         // should be sufficient). just to be on the absolutely safe side...
      } else if (VerifyRuleOptions == FGridSearchOptions::VERIFY_SymmetryMonomialsOnly || VerifyRuleOptions == FGridSearchOptions::VERIFY_MinimalMonomialSetOnly) {
         // testing for the full set above seemed to work fine in all cases.
         // But for very large rules it becomes rather slow. So I allow testing
         // only smaller sets (which *should* give the same result, unless
         // PointGroup's FPointGroup::CouldMonomialBeATargetFunction() is broken)
         Monomials = MakeSymmetryUniqueMonomialList(lmin, lmax, PointGroup, -1, 0);
      }
   }
   size_t
      nMonomials = Monomials.size();

   FDenseVector
      // results of numerical integration of the monomials using the specified integration rule
      RhsNumerical(nMonomials);
   RhsNumerical.setZero();
      
   FTransformedPointInfo
      // contains information about seed points transformed by current symmetry operation.
      // In particular, powers of its x,y,z coordinates.
      TrafodSeeds(SeedPoints, lmax);
   for (size_t iSymOp = 0; iSymOp != SymmetryOps.size(); ++ iSymOp) {
      TrafodSeeds.TransformPoints(SymmetryOps[iSymOp]);
//       #pragma omp parallel for schedule(dynamic)
      // ^- hmpf... needs a reduction, and I do not have the OmpAccBlock
      //    stuff around...
      //    FIXME: Note: we are doing too many monomials here (to be safe),
      //    but there is also the point that this is actually quite slow.
      //    This might not *only* be the linear algebra...
      //    We should probably OMP this and the corresponding parts of the
      //    main optimization routine properly.
      for (int iMonomial_ = 0; iMonomial_ != int(nMonomials); ++ iMonomial_) {
         size_t iMonomial = size_t(iMonomial_);
//       for (size_t iMonomial = 0; iMonomial != nMonomials; ++ iMonomial) {
         for (size_t iSeed = 0; iSeed != size_t(SeedPoints.cols()); ++ iSeed) {
            RhsNumerical[iMonomial] += weights[iSeed] * Monomials[iMonomial].EvalPoint(TrafodSeeds.TransformedPointInfo(iSeed));
         }
      }
   }

   FScalar
      fMaxErrAbs = 0,
      fMaxErrRel = 0;
//    for (size_t i = 0; i != nMonomials; ++ i) {
//    #pragma omp parallel for reduction(max:fMaxErrAbs,fMaxErrRel)
   // ^- 'error: user defined reduction not found for 'fMaxErrAbs'?!
   #pragma omp parallel for schedule(dynamic)
   for (int i_ = 0; i_ != int(nMonomials); ++ i_) {
      size_t i = size_t(i_);
      FMonomial const
         &mi = Monomials[i];
      FScalar
         // compute (numerically) exact integrals of the monomials over the unit sphere
//          I_Analytical = mi.CalcUnitSphereIntegralRaw(),
         I_Analytical = mi.CalcUnitSphereIntegral(), // <- this one should include the same norm and degeneracy factors as EvalPoint does.
         I_Numerical = RhsNumerical[i];
      FScalar
         fErrAbs = abs(I_Analytical - I_Numerical),
         fErrRel = 0;
      if (mi.nUnevenPow() == 0) {
         // monomials with at least one odd exponent have an analytical
         // sphere integral of zero (testing if the numerical rule reproduces
         // this is still necessary!)
         // So we cannot determine a relative error for those.
         fErrRel = fErrAbs / abs(I_Analytical);
      }
      if (0) {
         std::cout << fmt::format("   verify: m = [{:2},{:2},{:2}]  //  I_anl = {:20.12f}  //  I_num = {:20.12f}  //  erel = {:8.2e}{}", mi.ix, mi.iy, mi.iz, I_Analytical, I_Numerical, fErrRel, (abs(I_Analytical)<1e-10? std::string() : fmt::format("   In/Ia = {:.12f}",I_Numerical/I_Analytical))) << std::endl;
      }

      #pragma omp critical
      {
         if (fErrAbs > fMaxErrAbs)
            fMaxErrAbs = fErrAbs;
         if (fErrRel > fMaxErrRel)
            fMaxErrRel = fErrRel;
      }
   }
   
   if (pMaxErrAbs) *pMaxErrAbs = fMaxErrAbs;
   if (pMaxErrRel) *pMaxErrRel = fMaxErrRel;

   return fMaxErrAbs <= ThreshAbs;
}


void FSphereGridOptContext::SetGridExportReference(FExportedGridStats const &RefStats)
{
   m_LastExportedGridStats = RefStats;
}


FExportedGridStats const *FSphereGridOptContext::pGetLastExportedGridStats() const
{
   if (m_DidExportSomething)
      return &m_LastExportedGridStats;
   else
      return 0;
}


FScalar OptimizeSphereGrid(FPointArray &SeedPoints, size_t lmin, size_t lmax, FExportedGridStats *pGridStats, FGridSearchOptions const &Options)
{
   FSphereGridOptContext
      Context(SeedPoints, lmin, lmax, Options);
   if (pGridStats)
      Context.SetGridExportReference(*pGridStats);
   FScalar
      res = Context.Run();
   if (Context.pGetLastExportedGridStats() && pGridStats != 0) {
      *pGridStats = *Context.pGetLastExportedGridStats();
   }
   return res;
}



void RunGridSearchSchedule(FGridSearchOptions const &Options)
{
   FPointList
      InitialPoints;

   io.WriteLine();
   if (1) {
      std::cout << fmt::format(" Initial seed points -- cartesian coordinates") << std::endl;
      std::cout << "\n   SEED         X-POS            Y-POS          Z-POS" << std::endl;
      FMappedPointCloud
         PointCloud(Options.pPointGroup->GetOps());

      for (size_t iPt = 0; iPt != Options.InitialPoints.size(); ++ iPt) {
         FPoint const &p = Options.InitialPoints[iPt];
         FMappedPointCloud::FPointEntry const
            &e = PointCloud.Insert(p, iPt);
         std::string Comment;
         if (e.iPoint != (ptrdiff_t)iPt) {
//             throw std::runtime_error("Initial point collides with symmetry equivalent other initial point.");
            Comment = fmt::format("    IGNORED -- collides with #{}", e.iPoint);
         } else {
            InitialPoints.push_back(p);
         }
         std::cout << fmt::format("  {:>4}  {:15.8f}  {:15.8f}  {:15.8f}{}", iPt+1, double(p[0]), double(p[1]), double(p[2]), Comment) << std::endl;
      }
   }

   // copy initial points into array format.
   FPointArray
      pts(3, InitialPoints.size());
   for (size_t i = 0; i != InitialPoints.size(); ++ i)
      pts.col(i) = InitialPoints[i];

   FExportedGridStats
      LastExportedGridStats;
   // start with the actual schedule. First, check what we are supposed to be doing:
   // - a single l, or an explicitly given point l?
   // - or a incremental set of ls (upwards or downwards), starting at some given starting point?
   if (Options.lStep == 0) {
      // single point or explicit point list.
      FDegreeSchedule
         lList = Options.lList;
      if (lList.empty()) {
         // no list given. just lStart then. Make a single-entry list to simplify processing.
         if (Options.lStart < 1)
            throw std::runtime_error(fmt::format("RunGridSearchSchedule: no lStep given. Then either lMin or an explicit list of ls must be given."));
         lList.push_back(Options.lStart);
      }
      for (size_t il = 0; il != lList.size(); ++ il) {
         ptrdiff_t lx = lList[il];
         if (lx < 1)
            throw std::runtime_error(fmt::format("RunGridSearchSchedule: explicitly given 'l's of absolute schedule must be non-negative and >= 1. Cannot optimize grid of negative order '{}'", lx));
         
         OptimizeSphereGrid(pts, lx-1, lx, &LastExportedGridStats, Options);
      }
   } else {
      // incremental schedule.
      ptrdiff_t
         lStart = Options.lStart;
      if (lStart < 0) { // not given.
         if (Options.lStep > 0)
            lStart = 1;
         else
            throw std::runtime_error(fmt::format("RunGridSearchSchedule: processing downward optimization schedule (step={}). In this case, a positive starting l must be explicitly given.", Options.lStep));
      }

      ptrdiff_t
         lCur = lStart;
      for ( ; lCur >= 1; lCur += Options.lStep) {
         FScalar
            residual = OptimizeSphereGrid(pts, lCur-1, lCur, &LastExportedGridStats, Options);
         bool
            Exact = residual < Options.fThreshConsiderGridAsExact;
         if (Options.lStep < 0 && Exact) {
            std::cout << "\n-- Found first grid which integrates exactly -- stopping incremental downward search." << std::endl;
            break;
         }
         if (Options.lStep > 0 && !Exact) {
            std::cout << fmt::format("\n-- Found first grid which does not integrate exactly (mxe > {:.2e}) -- stopping incremental upward search.", double(Options.fThreshConsiderGridAsExact)) << std::endl;
            break;
         }
      }

      if (!Options.lList.empty()) {
         std::stringstream str;
         for (size_t il = 0; il != Options.lList.size(); ++ il) {
            if (il != 0) str << ", ";
            str << Options.lList[il];
         }
         std::cout << fmt::format("\n-- Resuming with schedule [{}] relative to current l = {}", str.str(), lCur) << std::endl;
         for (size_t il = 0; il != Options.lList.size(); ++ il) {
            ptrdiff_t lx = lCur + Options.lList[il];
            FScalar residual = OptimizeSphereGrid(pts, lx-1, lx, &LastExportedGridStats, Options);
            bool
               Exact = residual < Options.fThreshConsiderGridAsExact;
            if (Exact) {
               std::cout << "\n-- Found first grid which integrates exactly -- stopping process of explicit point list after incremental search." << std::endl;
               break;
            }
         }

      }
   }
   std::cout << "\n-- Grid optimization completed." << std::endl;
}

