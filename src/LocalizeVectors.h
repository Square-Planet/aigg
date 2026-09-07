#ifndef CX_LOCALIZE_VECTORS_H
#define CX_LOCALIZE_VECTORS_H

#include <stdint.h> // for size_t
#include <ostream>
#include <cmath>
// ^-- note: if using high-precision numbers, make sure they are included before
// this and that they bring their math functions (sqrt, sin, atan2, pow, etc)
// via argument-dependent lookup!


namespace ct_loc_t {

struct FLocalizeOptions
{
   unsigned
      nLocExp; //  2 or 4, exponent p in localization functional: L = sum[iA] n(i,A)^p
   unsigned
      nMaxIt;
   double
      ThrLoc;
   std::ostream
      *pxout; // if non-0, write output here.
//    FLog
//       *pLog; // if non-0, write output here.
   unsigned
      Verbosity;
   FLocalizeOptions();
};
FLocalizeOptions::FLocalizeOptions()
   : nLocExp(2), nMaxIt(2048), ThrLoc(1e-8), pxout(0), Verbosity(1)
{}



/// execute a 2x2 rotation of the vectors pA and pB (both with length nSize)
/// of angle phi. In-place.
template<class FScalar>
void Rot2x2T(FScalar *IR_RP pA, FScalar *IR_RP pB, size_t nSize, FScalar phi)
{
   using std::cos;
   using std::sin;
   FScalar
      cs = cos(phi),
      ss = sin(phi);
   for (size_t i = 0; i < nSize; ++ i) {
      FScalar
         tempA =  cs * pA[i] + ss * pB[i],
         tempB = -ss * pA[i] + cs * pB[i];
      pA[i] = tempA;
      pB[i] = tempB;
   }
}


template<class FScalar> FScalar pow2(FScalar x) { return x*x; }
template<class FScalar> FScalar pow3(FScalar x) { return x*x*x; }
template<class FScalar> FScalar pow4(FScalar x) { FScalar x2 = x*x; return x2*x2; }

// NumCols(CVec)

template<class FScalar, class FDenseMatrix>
void LocalizeVectorsT(FDenseMatrix &CVec, size_t nRows, size_t nCols, size_t const *pGroupOffsets, size_t nGroups, FLocalizeOptions const &Opt)
{
   using std::pow;
   using std::sqrt;
   using std::atan2;
   size_t
      iIt;
   FScalar const
      f025 = (1/FScalar(4));
   FScalar
      L = -1, fVar2 = -1;
   if (Opt.pxout && Opt.Verbosity >= 2) {
      *Opt.pxout << "   ITER.       LOC(Orb)      GRADIENT" << std::endl;
   }
   for (iIt = 0; iIt < Opt.nMaxIt; ++ iIt) {
      // calculate the value of the functional.
      L = 0;
      // loop over groups (normally a ``group'' is all functions on an atom)
      for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp)
         for (size_t iVec = 0; iVec < nCols; ++ iVec) {
            FScalar nA = 0; // population on the atom.
            for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn)
               nA += pow2(CVec(iFn,iVec));
            L += pow(nA, (int)Opt.nLocExp);
         }
      L = pow(L, 1/FScalar(Opt.nLocExp));  // <- easier to compare different powers that way.

      fVar2 = 0;
      // loop over vector pairs.
      for (size_t iVec = 0; iVec < nCols; ++ iVec)
         for (size_t jVec = 0; jVec < iVec; ++ jVec) {
            FScalar
               Aij = 0, // hessian for 2x2 rotation (with variable atan(phi/4) substituted)
               Bij = 0; // gradient for 2x2 rotation (with variable atan(phi/4) substituted)
            // loop over groups (normally atoms)
            for (size_t iGrp = 0; iGrp < nGroups; ++ iGrp) {
               // calculate the charge matrix elements
               //     Cii = <i|A|i>,
               //     Cjj = <j|A|j>,
               // and Cij = <i|A|j>.
               FScalar
                  Cii = 0, Cij = 0, Cjj = 0;
               for (size_t iFn = pGroupOffsets[iGrp]; iFn != pGroupOffsets[iGrp+1]; ++ iFn) {
                  Cii += CVec(iFn,iVec) * CVec(iFn,iVec);
                  Cjj += CVec(iFn,jVec) * CVec(iFn,jVec);
                  Cij += CVec(iFn,iVec) * CVec(iFn,jVec);
               }

               // see Supporting Information for `` Intrinsic atomic orbitals:
               // An unbiased bridge between quantum theory and chemical concepts''
               // for a description of what this does.
               if (Opt.nLocExp == 2 || Opt.nLocExp == 3) {
                  Aij += sqr(Cij) - f025*sqr(Cii - Cjj);
                  Bij += Cij*(Cii - Cjj);
               } else if (Opt.nLocExp == 4) {
                  Aij += -pow4(Cii) - pow4(Cjj) + 6*(pow2(Cii) + pow2(Cjj))*pow2(Cij) + pow3(Cii)*Cjj + Cii*pow3(Cjj);
                  Bij += 4*Cij*(pow3(Cii) - pow3(Cjj));
               }
            }

            // non-degenerate rotation? This condition is not quite right,
            // and in Python it was not required. However, the Fortran atan2
            // sometimes did mysterious things without it.
            if (pow2(Aij) + pow2(Bij) > FScalar(Opt.ThrLoc*1e-2)) {
               FScalar
                  phi = f025 * atan2(Bij,-Aij);
               // 2x2 rotate vectors iVec and jVec.
               Rot2x2T(&CVec(0,iVec), &CVec(0,jVec), nRows, phi);
               fVar2 += pow2(phi);
            }
         }
      fVar2 = sqrt(fVar2 / FScalar(pow2(nCols)));

      if (Opt.pxout && Opt.Verbosity >= 2) {
         *Opt.pxout << fmt::format("{:6}     {:12.6f}      {:8.2e}\n", (1+iIt), double(L), double(fVar2));
      }

      if (fVar2 < Opt.ThrLoc)
         break;
   }

   if (Opt.pxout && Opt.Verbosity >= 1) {
      if (Opt.Verbosity >= 2)
         *Opt.pxout << "\n";
      *Opt.pxout << fmt::format(" Iterative localization: IB/PM, {} iter; Final gradient {:8.2e}", iIt, fVar2);
      Opt.pxout->flush();
   }
}

} // namespace ct_loc

#endif // CX_LOCALIZE_VECTORS_H
