#include "Aigg.h"
#include "MathOps.h"

#include <sstream>
#include <fstream>
#include <typeinfo>
#include <ctime>
#include <stdlib.h> // for atoi
#include <Eigen/Cholesky>
#include <Eigen/SVD>

// #include <unsupported/Eigen/MatrixFunctions> // for matrix exponential

#define EIGEN_SVD_ALGO JacobiSVD // straight up jacobi SVD (more accurate)
// #define EIGEN_SVD_ALGO BDCSVD // bidirectional divide&conquer SVD (faster)


FSvd::FSvd(FDenseMatrix const &A)
{
   Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
      svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
   sigma = svd.singularValues();
   U = svd.matrixU();
   V = svd.matrixV();
}


void FSvd::Compress(FScalar thr) {
   size_t nSig = 0;
   for (size_t i = 0; i < size_t(sigma.size()); ++i) {
      if (sigma[i] >= thr)
         nSig += 1;
   }
   U.conservativeResize(U.rows(), nSig);
   V.conservativeResize(V.rows(), nSig);
   sigma.conservativeResize(nSig);
}


FDenseMatrix FSvd::Solve(FDenseMatrix const &Y) {
   FDenseMatrix
      UtY = U.transpose() * Y;
   for (size_t i = 0; i < size_t(sigma.size()); ++i) {
      UtY.row(i) /= sigma[i];
   }
   return V * UtY;
}


FDenseMatrix FSvd::SolveT(FDenseMatrix const &Y) {
   FDenseMatrix
      VtY = V.transpose() * Y;
   for (size_t i = 0; i < size_t(sigma.size()); ++i) {
      VtY.row(i) /= sigma[i];
   }
   return U * VtY;
}


// - TODO: maybe make a common class for both, for MakeStep/MakeStepSymmetric?
//   and which one would be used would depend on whether the input matrix is
//   symmetric.
// - Note that one could even just compute an actual SVD from an Eigh
//   if the input matrix A is symmetric (even with non-negative sigma if
//   needed: by first setting U=V=Ev, and then flipping the signs of both
//   ew[i] columns U[:,i] of U if ew[i] is negative).
FEigh::FEigh(FDenseMatrix const &H, unsigned Flags)
{
   Eigen::SelfAdjointEigenSolver<FDenseMatrix>
      sd(H, Eigen::ComputeEigenvectors); // compute spectral decomposition
   if (Flags & LargeEwFirst) {
      lambda = sd.eigenvalues().reverse();
      V = sd.eigenvectors().rowwise().reverse();
      // ^- note: A.rowwise().reverse() reverses the order of *COLUMNS* of A,
      //    because that is what you do when you apply the reverse operation
      //    along the data stored in fixed rows (``rowwise'').
   } else {
      lambda = sd.eigenvalues();
      V = sd.eigenvectors(); 
   }
}



struct FEighOrSvdImpl
{
   Eigen::SelfAdjointEigenSolver<FDenseMatrix>
      m_eigh;
   Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
      m_svd;
   FEighOrSvd::FDecompType
      m_DecompType;
};

//    enum FDecompType {
//       DECOMP_Auto, // make spectral decomposition if 'A' is square and symmetric/hermitian, SVD otherwise.
//       DECOMP_Eigh, // make spectral decomposition (requires 'A' to be square and symmetric/hermitian)
//       DECOMP_Svd // make singular value decomposition (required if 'A' is not square or not symmetric)
//    };
   
FEighOrSvd::FEighOrSvd(FDenseMatrix const &A, FDecompType DecompType)
   : p(new FEighOrSvdImpl)
{
   if (DecompType == DECOMP_Auto) {
      // by default use SVD---always works.
      DecompType = DECOMP_Svd;
      // check for eigh special case conditions...
      // ...is this a square matrix?
      if (A.rows() == A.cols()) {
         // is it additionally symmetric?
         FScalar
            SymVio = (A - A.adjoint()).array().maxCoeff();
         if (IsAlmostZero(SymVio))
            // got a symmetric one. Use eigh instead.
            DecompType = DECOMP_Eigh;
      }
   }
   assert(DecompType == DECOMP_Eigh || DecompType == DECOMP_Svd);
   p->m_DecompType = DecompType;
   
   // now that we have decided on what to do, ask eigen
   // to make the actual matrix decomposition.
   if (p->m_DecompType == DECOMP_Eigh) {
      p->m_eigh.compute(A, Eigen::ComputeEigenvectors);
   } else {
      assert(p->m_DecompType == DECOMP_Svd);
      p->m_svd.compute(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
   }
}


FEighOrSvd::~FEighOrSvd()
{
   delete p;
}
   
Eigen::Ref<FDenseMatrix const> FEighOrSvd::U() const
{
   if (p->m_DecompType == DECOMP_Eigh) {
      // note: the 'reverse' is for getting the eigenvalues into large->small order.
      return p->m_eigh.eigenvectors().rowwise().reverse();
   } else {
      assert(p->m_DecompType == DECOMP_Svd);
      return p->m_svd.matrixU();
   }
}

Eigen::Ref<FDenseMatrix const> FEighOrSvd::V() const
{
   if (p->m_DecompType == DECOMP_Eigh) {
      return U();
   } else {
      assert(p->m_DecompType == DECOMP_Svd);
      return p->m_svd.matrixV();
   }
}

Eigen::Ref<FDenseArray const> FEighOrSvd::vals() const
{
   if (p->m_DecompType == DECOMP_Eigh) {
      // note: the 'reverse' is for getting the eigenvalues into large->small order.
      return p->m_eigh.eigenvalues().reverse();
   } else {
      assert(p->m_DecompType == DECOMP_Svd);
      return p->m_svd.singularValues();
   }
}

#ifdef INCLUDE_ABANDONED
// protected:
//    FDecompType
//       m_DecompType;
//    FDenseMatrix
//       m_U, m_V;
//    FDenseVector
//       m_sigma;
#endif // INCLUDE_ABANDONED

FDenseMatrix FLinearSolver::HalfMxm1(Eigen::Ref<FDenseMatrix const> const &A) {
   throw std::runtime_error(fmt::format("linear solver {}::HalfMxm1 not implemented.", typeid(*this).name()));
}

FDenseMatrix FLinearSolver::HalfMxm2(Eigen::Ref<FDenseMatrix const> const &A) {
   throw std::runtime_error(fmt::format("linear solver {}::HalfMxm2 not implemented.", typeid(*this).name()));
}

FDenseMatrix FLinearSolver::Mxm(Eigen::Ref<FDenseMatrix const> const &A) {
   throw std::runtime_error(fmt::format("linear solver {}::Mxm check operation order.", typeid(*this).name()));
   return this->HalfMxm2(this->HalfMxm1(A)); // <-- FIXME: right order?
}


// Implements the symmetric positive solver using a spectral decomposition of S:
//
//      S = Ev diag(ew) Ev.T
//
// where ew is a diagonal matrix of eigenvalues and Ev are the column
// vectors of eigenvectors (Ev is an unitary matrix).
//
// From this we construct an object
//
//     SmhEvT = diag(1/sqrt(ew)) * Ev.T
//
// where we only take the column of Ev and elements of ew which belong
// to ew > Thr (so the object is *not* necessarily square, as it
// only represents the non-singular subspace of S's spectrum).
//
// Inserting this SmhEvT two times then affords:
//
//    A*S^{-1}*B = A * Ev * diag(1/ew) * Ev.T * B
//               = (A * Ev * diag(1/sqrt(ew))) * (diag(1/sqrt(ew)) * Ev.T * B)
//               = (A * SmhEvT.T) * (SmhEvT * B)
//               = (SmhEvT * A.T)^T * (SmhEvT * B)
//
// ...which are operations fulfilling the symmetric half solver requirements.
//
// Unlike the Scd solver below, this version should be capable of gracefully
// dealing with numerical zero eigenvalues (as long as they are reasonably
// separated from the small-but-significant parts of the spectrum), because it
// explicitly constructs a representation of the non-singular subspace.
struct FSymmetricPositiveSolverEigh : public FLinearSolver
{
#ifdef INCLUDE_ABANDONED
//    explicit FSymmetricPositiveSolverEigh(FDenseMatrix const &S)
//    {
//        // compute spectral decomposition
//       Eigen::SelfAdjointEigenSolver<FDenseMatrix>
//          sd(S);
//       FDenseVector
//          ew = sd.eigenvalues();
//       FDenseMatrix const
//          &ev = sd.eigenvectors();
//       PrintVector(ew, fmt::format("TARGET FUNCTION OVERLAP // EIGENVALUES"));
//       PrintVector(S.diagonal(), fmt::format("TARGET FUNCTION OVERLAP // DIAGONAL"));
//       // locate subset of eigenvalues and eigenvectors for non-singular space
//       // (all negative ones are considered singular, too).
//       // Note that eigenvalues come out ordered small to large.
//       size_t
//          iFirstCol = 0;
//       while (iFirstCol < size_t(ew.size()) && ew[iFirstCol] < g_ThrVerySmall)
// //       while (iFirstCol < size_t(ew.size()) && ew[iFirstCol] < 1e-33Q) // FIXME: this needs normalization of the monomials!
//          iFirstCol += 1;
//       size_t
//          nTotal = S.rows(),
//          nRank = nTotal - iFirstCol;
//       // make the diag(ew^{-1/2}) * Ev.T object (we store the non-transposed
//       // object because that allows more more efficient multiplications. In
//       // principle, that is, not that it actually matters with emulated quad
//       // precision arithmetic...)
//       m_SmhEvT.resize(nRank, nTotal);
//       for (size_t i = 0; i < nRank; ++ i) {
//          size_t iew = i + iFirstCol;
//          m_SmhEvT.row(i).array() = (1./sqrt(ew[iew])) * ev.col(iew).array();
//       }
//       m_nEffectiveRank = nRank;
//    }
#endif // INCLUDE_ABANDONED
   explicit FSymmetricPositiveSolverEigh(FDenseMatrix const &S, bool AllowMxm)
      : m_AllowMxm(AllowMxm)
   {
       // compute spectral decomposition
      Eigen::SelfAdjointEigenSolver<FDenseMatrix>
         sd(S);
      FDenseVector
         ew = sd.eigenvalues().reverse();
      Eigen::Ref<FDenseMatrix const>
         ev = sd.eigenvectors().rowwise().reverse();
         // ^- note: this inverts the order of COLUMNS! (because each individual
         //    row contains data for all columns, and that is what is individually
         //    reversed for each row ("rowwise")
      // locate subset of eigenvalues and eigenvectors for non-singular space
      // (all negative ones are considered singular, too).
      // Note that we reversed the data to get eigenvalues in order large to small.
      size_t
         nRank = 0;
      while (nRank < size_t(ew.size()) && ew[nRank] > g_ThrVerySmall)
         nRank += 1;
      size_t
         nTotal = S.rows();
      if (0) {
         PrintVector(ew, fmt::format("TARGET FUNCTION OVERLAP // EIGENVALUES (nRank = {}, ThrEw = {:8.2e}, nTot = {})", nRank, double(g_ThrVerySmall), nTotal));
         PrintVector(S.diagonal(), fmt::format("TARGET FUNCTION OVERLAP // DIAGONAL"));
      }
      // make the diag(ew^{-1/2}) * Ev.T object (we store the non-transposed
      // object because that allows more more efficient multiplications. In
      // principle, that is, not that it actually matters with emulated quad
      // precision arithmetic...)
      m_SmhEvT.resize(nRank, nTotal);
      for (size_t iew = 0; iew < nRank; ++ iew)
         m_SmhEvT.row(iew).array() = (1./sqrt(ew[iew])) * ev.col(iew).array();
      if (m_AllowMxm) {
         m_ShEv.resize(nTotal, nRank);
         for (size_t iew = 0; iew < nRank; ++ iew)
            m_ShEv.col(iew).array() = (sqrt(ew[iew])) * ev.col(iew).array();
      }
      m_nEffectiveRank = nRank;
   }
   
   FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_SmhEvT * A;
   }
   FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_SmhEvT.transpose() * A;
   }
   FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return HalfSolve2(HalfSolve1(A));
   }
   
   FDenseMatrix HalfMxm1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      cx_assert_rt(m_AllowMxm);
      return m_ShEv * A;
   }
   FDenseMatrix HalfMxm2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      cx_assert_rt(m_AllowMxm);
      return m_ShEv.transpose() * A;
   }
protected:
   bool
      m_AllowMxm;
   FDenseMatrix
      // (nRank, nTotal)-shape matrix of ew^{-1/2}-scaled and transposed eigenvectors of S
      m_SmhEvT,
      // (nTotal, nRank)-shape matrix of ew^{+1/2}-scaled and transposed eigenvectors of S,
      // but only if m_AllowMxm is true. (should probably just have kept the eigenvalues.
      // too lazy to fix)
      m_ShEv;
};


// This one implements the symmetric solver using a Choleksy decomposition of S:
//
//    S = L * L.T
//
// where L is a lower triangular matrix (which allows for efficient thread-
// parallel Gaussian elimination and triangular matrix multiplication... though
// not in eigen and not with emulated high-precision arithmetic).
//
// This is faster and more stable than the eigensolver above for cases in which
// one knows that *real* singularities never occur in the matrix S (i.e., S is
// positive definite, not only positive semi-definite), and that even small
// eigenvalues of S, so they occur, actually represent non-negligible
// interactions which *must* be treated.
//
// In practice this occurs for us when using non-symmetrized target functions,
// in which case we do not technically get the degeneracies.
struct FSymmetricPositiveSolverCd : public FLinearSolver
{
   explicit FSymmetricPositiveSolverCd(FDenseMatrix const &S)
   {
      m_Scd.compute(S); // compute Cholesky decomposition (unpivoted! results assummed symmetric in some parts of code)
//       PrintVector(m_Scd.matrixL().diagonal(), fmt::format("TARGET FUNCTION OVERLAP // Scd/L matrix diagonal"));
      // A Choleksy decomposition cannot actually compute the rank. Using this
      // means one already assumes this has full rank.
      m_nEffectiveRank = S.rows();
   }
   
   FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_Scd.matrixL().solve(A);
   }
   FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_Scd.matrixL().transpose().solve(A);
   }
   FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_Scd.solve(A);
   }
   
   FDenseMatrix HalfMxm1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_Scd.matrixL() * A;
   }
   FDenseMatrix HalfMxm2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return m_Scd.matrixL().transpose() * A;
   }
protected:
   Eigen::LLT<FDenseMatrix>
      m_Scd;
   size_t
      m_nEffectiveRank;
};


// This one implements a partially pivoted Cholesky-like decomposition
// with extracted diagonal elements:
//
//    S = P.T L diag(d) L.T P
//
// where P is a permutation matrix.
//
// It is technically not supposed to be rank revealing, but in practice
// appears to work rather stably in many situations in which it really
// is not supposed to.
//
// WARNING:
// - unlike for Eigh and Cd, for this one HalfSolve1
//   and HalfSolve2 are not just transposes of each other; at least
//   not necessarily in a standard implementation.
//
// - we *could* technically implement this with symmetric half-solves,
//   too: we'd just need to split the D matrix into two times
//   sqrt(D), and we could even just eliminate latter the columns of L
//   (and entries of D) which belong to the singular entries.
//
//   I am currently of the impression that for *actual* hard singularities,
//   (i.e., exact zeros), that would actually be equivalent to solving
//   the system on the non-singular subspace. But I admittedly did not
//   think about this very hard. Maybe just try?
//
// - UPDATE: yes... implemented it. sqrt split and HalfSolve / HalfMxm
//   work just fine
#ifdef INCLUDE_ABANDONED
// struct FSymmetricPositiveSolverLDLT : public FLinearSolver
// {
// //    explicit FSymmetricPositiveSolverLDLT(FDenseMatrix const &S)
// //    {
// //       m_Scd.compute(S);
// //       m_DiagVals = m_Scd.vectorD();
// //       m_nEffectiveRank = 0;
// //       for (size_t i = 0; i != size_t(m_DiagVals.size()); ++ i)
// //          if (!IsAlmostZero(m_DiagVals[i]))
// //             // yes... I am aware we are not supposed to use this
// //             // for computing an effective rank.
// //             m_nEffectiveRank += 1;
// //    }
//    explicit FSymmetricPositiveSolverLDLT(FDenseMatrix const &S)
//    {
//       m_Scd.compute(S);
//       FDenseVector const
//          &D = m_Scd.vectorD();
//       m_nEffectiveRank = 0;
//       for (size_t i = 0; i != size_t(D.size()); ++ i)
//          if (!IsAlmostZero(D[i]))
//             // yes... I am aware we are not supposed to use this
//             // for computing an effective rank.
//             m_nEffectiveRank += 1;
//    }
//    
//    FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       // The decomposition of S itself is S = P.T L diag(D) L.T P.
//       // So the decomposition of the inverse translates to
//       //
//       //    S^{-1} = inv(P.T L diag(D) L.T P)
//       //           = inv(P) inv(L.T) inv(diag(D)) inv(L) inv(P.T)
//       //
//       // We will split diag(D) into diag(sqrt(D)) * diag(sqrt(D)) and afterwards
//       // solve for the right half of this:
//       //
//       //    out = inv(diag(sqrt(D))) inv(L) inv(P.T) A
//       //        = diag(1/sqrt(D)) inv(L) P
//       //
//       // (note: for permutation matrices inv(P) = transpose(P), so inv(P.T) = P itself)
//       //
//       // let's give it a try. 
//       FScdType::TranspositionType const
//          &P = m_Scd.transpositionsP();
//       FDenseMatrix
//          out = P * A;
//       m_Scd.matrixL().solveInPlace(out);
//       FScalar
//          Thr = MakeDiagThrForSolve();
//       // now out = inv(L) inv(P.T) A; only diag(1/sqrt(D)) left.
//       for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
//          if (m_DiagVals[i] > Thr)
//             out.row(i) /= sqrt(m_DiagVals[i]);
//          else
//             out.row(i).setZero();
//       }
//       return out;
// //       throw std::runtime_error("HalfSolve1 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       // see comment on HalfSolve1. We're now processing the
//       // remaining terms, proceeding towards the left.
//       FDenseMatrix
//          out = A;
//       FScalar
//          Thr = MakeDiagThrForSolve();
//       for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
//          if (m_DiagVals[i] > Thr)
//             out.row(i) /= sqrt(m_DiagVals[i]);
//          else
//             out.row(i).setZero();
//       }
//       m_Scd.matrixL().transpose().solveInPlace(out);
//       FScdType::TranspositionType const
//          &P = m_Scd.transpositionsP();
//       return P.transpose() * A;
// //       throw std::runtime_error("HalfSolve2 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       return m_Scd.solve(A);
// //       return HalfSolve2(HalfSolve1(A));
//    }
//    
// protected:
//    typedef Eigen::LDLT<FDenseMatrix>
//       FScdType;
//    FScdType
//       m_Scd;
//    FDenseVector
//       m_DiagVals; // entries on diagonal of diag(D) matrix.
//    FScalar MakeDiagThrForSolve() const {
//       return g_ThrVerySmall;
//       FScalar
//          Thr = (std::numeric_limits<FScalar>::min)();
//       // ^- in this one we can deal with *very* small thresholds...
//       //    if it can be represented in a float, we probably can
//       //    deal with it!
//       return Thr;
//    }
// };
#endif // INCLUDE_ABANDONED
struct FSymmetricPositiveSolverLDLT : public FLinearSolver
{
   explicit FSymmetricPositiveSolverLDLT(FDenseMatrix const &S)
   {
      m_Scd.compute(S);
      m_DiagVals = m_Scd.vectorD();
      m_nEffectiveRank = 0;
      FScalar
         Thr = MakeDiagThrForSolve();
      for (size_t i = 0; i != size_t(m_DiagVals.size()); ++ i)
//          if (!IsAlmostZero(m_DiagVals[i]))
         if ((m_DiagVals[i] > Thr) && !IsAlmostZero(sqrt(m_DiagVals[i])))
            // yes... I am aware we are not supposed to use this
            // for computing an effective rank.
            m_nEffectiveRank += 1;
   }
   
   FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      // The decomposition of S itself is S = P.T L diag(D) L.T P.
      // So the decomposition of the inverse translates to
      //
      //    S^{-1} = inv(P.T L diag(D) L.T P)
      //           = inv(P) inv(L.T) inv(diag(D)) inv(L) inv(P.T)
      //
      // We will split diag(D) into diag(sqrt(D)) * diag(sqrt(D)) and afterwards
      // solve for the right half of this:
      //
      //    out = inv(diag(sqrt(D))) inv(L) inv(P.T) A
      //        = diag(1/sqrt(D)) inv(L) P
      //
      // (note: for permutation matrices inv(P) = transpose(P), so inv(P.T) = P itself)
      //
      // let's give it a try. 
      FScdType::TranspositionType const
         &P = m_Scd.transpositionsP();
      FDenseMatrix
         out = P * A;
      m_Scd.matrixL().solveInPlace(out);
      FScalar
         Thr = MakeDiagThrForSolve();
      // now out = inv(L) inv(P.T) A; only diag(1/sqrt(D)) left.
      for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
         if (m_DiagVals[i] > Thr)
            out.row(i) /= sqrt(m_DiagVals[i]);
         else
            out.row(i).setZero();
      }
      return out;
   }
//    FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       FScdType::TranspositionType const
//          &P = m_Scd.transpositionsP();
//       FDenseMatrix
//          out = P * A;
//       m_Scd.matrixL().solveInPlace(out);
//       FScalar
//          Thr = MakeDiagThrForSolve();
//       // now out = inv(L) inv(P.T) A; only diag(1/sqrt(D)) left.
//       FDenseVector
//          InvSqrtD = m_DiagVals;
//       for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
//          if (m_DiagVals[i] > Thr)
//             InvSqrtD[i] /= sqrt(m_DiagVals[i]);
//          else
//             InvSqrtD[i] = 0;
//       }
// //       out.array().rowwise() *= InvSqrtD.array().transpose();
//       return out;
// //       return (out.array().rowwise() * InvSqrtD.array().transpose()).matrix();
// //       throw std::runtime_error("HalfSolve1 not implemented for LDLT solver! (it could be...)");
//    }
   FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      // see comment on HalfSolve1. We're now processing the
      // remaining terms, proceeding towards the left.
      FDenseMatrix
         out = A;
      FScalar
         Thr = MakeDiagThrForSolve();
      for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
         if (m_DiagVals[i] > Thr)
            out.row(i) /= sqrt(m_DiagVals[i]);
         else
            out.row(i).setZero();
      }
      m_Scd.matrixL().transpose().solveInPlace(out);
      FScdType::TranspositionType const
         &P = m_Scd.transpositionsP();
      return P.transpose() * out;
   }
   
   FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
      return HalfSolve2(HalfSolve1(A));
   }
   
   FDenseMatrix HalfMxm1(Eigen::Ref<FDenseMatrix const> const &A) { // override
      // The decomposition of S itself is S = P.T L diag(D) L.T P.
      // For HalfSolve1 we implemented:
      //
      //    out = inv(diag(sqrt(D))) inv(L) inv(P.T) A
      //        = diag(1/sqrt(D)) inv(L) P A
      //
      // so here we want the inverse of that. Which is 
      //
      //    out = P.T L diag(sqrt(D)) A
      //
      // if I'm not very mistaken.
//       FDenseVector
//          SqrtD = sqrt(m_DiagVals.array());
//       FDenseMatrix
//          Adl = m_Scd.matrixL() * ((A.array().rowwise() * SqrtD.array().transpose())).matrix();
//       FScdType::TranspositionType const
//          &P = m_Scd.transpositionsP();
//       return P.transpose() * Adl;
      FDenseMatrix
         Ad = A;
      for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
         Ad.row(i) *= sqrt(m_DiagVals[i]);
      }
      FDenseMatrix
         Adl = m_Scd.matrixL() * Ad;
      FScdType::TranspositionType const
         &P = m_Scd.transpositionsP();
      return P.transpose() * Adl;
   }
   FDenseMatrix HalfMxm2(Eigen::Ref<FDenseMatrix const> const &A) { // override
      FScdType::TranspositionType const
         &P = m_Scd.transpositionsP();
      FDenseMatrix
         PA = P * A,
         out = m_Scd.matrixL().transpose() * PA;
      for (size_t i = 0; i < size_t(m_DiagVals.size()); ++ i) {
         out.row(i) *= sqrt(m_DiagVals[i]);
      }
      return out;
   }
protected:
   typedef Eigen::LDLT<FDenseMatrix>
      FScdType;
   FScdType
      m_Scd;
   FDenseVector
      m_DiagVals; // entries on diagonal of diag(D) matrix.
   FScalar MakeDiagThrForSolve() const {
//       return g_ThrVerySmall;
      FScalar
         Thr = (std::numeric_limits<FScalar>::min)();
      // ^- in this one we can deal with *very* small thresholds...
      //    if it can be represented in a float, we probably can
      //    deal with it!
      return Thr;
   }
};

#ifdef INCLUDE_ABANDONED
// struct FSymmetricPositiveSolverLDLT : public FLinearSolver
// {
//    explicit FSymmetricPositiveSolverLDLT(FDenseMatrix const &S)
//    {
//       m_Scd.compute(S);
//       FDenseVector const
//          &D = m_Scd.vectorD();
//       m_nEffectiveRank = 0;
//       for (size_t i = 0; i != size_t(D.size()); ++ i)
//          if (!IsAlmostZero(D[i]))
//             // yes... I am aware we are not supposed to use this
//             // for computing an effective rank.
//             m_nEffectiveRank += 1;
//    }
//    
//    FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       throw std::runtime_error("HalfSolve1 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       throw std::runtime_error("HalfSolve2 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       return m_Scd.solve(A);
//    }
//    
// protected:
//    Eigen::LDLT<FDenseMatrix>
//       m_Scd;
// };


// struct FSymmetricPositiveSolverLDLT : public FLinearSolver
// {
//    explicit FSymmetricPositiveSolverLDLT(FDenseMatrix const &S)
//    {
//       m_Scd.compute(S);
//       FDenseVector const
//          &D = m_Scd.vectorD();
//       m_nEffectiveRank = 0;
//       for (size_t i = 0; i != size_t(D.size()); ++ i)
//          if (!IsAlmostZero(D[i]))
//             // yes... I am aware we are not supposed to use this
//             // for computing an effective rank.
//             m_nEffectiveRank += 1;
//    }
//    
//    FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       throw std::runtime_error("HalfSolve1 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       throw std::runtime_error("HalfSolve2 not implemented for LDLT solver! (it could be...)");
//    }
//    FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) { // override
//       return m_Scd.solve(A);
//    }
//    
// protected:
//    Eigen::LDLT<FDenseMatrix>
//       m_Scd;
//    size_t
//       m_nEffectiveRank;
// };
#endif // INCLUDE_ABANDONED

FLinearSolverPtr MakeLinearSolver(FDenseMatrix const &S, FLinearSolverType SolverType, bool AllowMxm)
{
   switch (SolverType) {
      case LINSOLVE_PositiveSymmetric_Eigh:
         return FLinearSolverPtr(new FSymmetricPositiveSolverEigh(S, AllowMxm));
      case LINSOLVE_PositiveSymmetric_Cholesky:
         return FLinearSolverPtr(new FSymmetricPositiveSolverCd(S));
      case LINSOLVE_PositiveSymmetric_LDLT:
         return FLinearSolverPtr(new FSymmetricPositiveSolverLDLT(S));
      default: throw std::runtime_error(fmt::format("MakeLinearSolver: solver type index '{}' not recognized", int(SolverType)));
   }
}


// Makes an orthogonal matrix from the anti-symmetric matrix (tau*A)
// via a Cayley transformation.
//
// Notes:
// - This is *NOT* the self-inverse Cayley transform!
// - With tau = 1.0, this function yields the positive Cayley transform:
//
//   CayleyU(A) = la.solve(I - 0.5*A, I + 0.5*A)
//              = 1 + A + (1/2) A^2 + O(A^3)
//
//   This agrees up to second order terms with exp(A).
// - The inverse of this is CayleyA(U) = 2 * la.solve(I + U, U - I) = -2 * la.solve(I + U, I - U)
FMatrixRR CayleyU(FMatrixRR A, FScalar tau)
{
   FMatrixRR
      I;
   I.setIdentity();

//    Eigen::EIGEN_SVD_ALGO<FMatrixRR>
   // "SVDBase: thin U and V are only available when your matrix has a dynamic number of columns."' failed.
   // (does that matter?)
   Eigen::EIGEN_SVD_ALGO<FDenseMatrix>
      svd(I - (tau/FScalar(2))*A, Eigen::ComputeThinU | Eigen::ComputeThinV);
   svd.setThreshold(1e-12);
   return svd.solve(I + (tau/FScalar(2))*A);
}


// form a (orthogonal) rotation matrix R from an input antisymmetric matrix (tau*A),
// similar to R = exp(tau*A)
FMatrixRR BuildRotation(FMatrixRR A, FScalar tau)
{
//    return (tau*A).exp();
   return CayleyU(A, tau);
   // - ...why not use exp(A) directly? Apparently the matrix exponential of eigen
   //   does not work with 128bit floats (at the time of this writing, anyway).
   //   Well, nothing to complain, it IS marked as "unsupported", and probably
   //   hard to get right in general
   //   (and I do not want put either their non-hermitian or hermitian complex
   //   float128 spectral decompositions to the test to build it myself).
}


FVectorR Solve(FMatrixRR const &M, FVectorR const &rhs)
{
   // ...I guess partial pivot LU would also do the trick. But since it's
   // for 3x3 stuff only...
   Eigen::FullPivLU<FMatrixRR>
      lu(M);
   return lu.solve(rhs);
}




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

FScalar Rmsd(FDenseMatrix const &m) {
   return Rmsd(m, size_t(m.rows() * m.cols()));
}

FScalar Rmsd(FDenseMatrix const &m, size_t nDegreesOfFreedom) {
   if (m.size() == 0)
      return 0.;
   if (nDegreesOfFreedom == 0)
      nDegreesOfFreedom = 1;
   return sqrt(FScalar(m.array().abs2().sum()) / FScalar(nDegreesOfFreedom));
}


FScalar RmsdFromIdentity(FDenseMatrix const &m) {
   using std::sqrt;
   cx_assert_rt(m.rows() == m.cols()); // <-- "rmsd from identity" only makes sense for square matrices.
   FScalar fAcc(0);
   size_t N = m.rows();
   for (size_t iCol = 0; iCol < N; ++ iCol) {
      FScalar
         a = m.col(iCol).topRows(iCol).array().abs2().sum(),
         b = m.col(iCol).bottomRows(N-iCol-1).array().abs2().sum();
      fAcc += (a + b) + sqr(m(iCol,iCol) - 1);
   }
   return sqrt(fAcc)/N;
}


#ifdef INCLUDE_ABANDONED
// FScalar RmsdFromIdentity(FDenseMatrix const &m) {
//    using std::sqrt;
//    cx_assert_rt(m.rows() == m.cols()); // <-- "rmsd from identity" only makes sense for square matrices.
//    FScalar fAcc(0);
//    size_t N = m.rows();
//    for (size_t iCol = 0; iCol < N; ++ iCol) {
//       for (size_t iRow = 0; iRow < N; ++ iRow) {
//          fAcc += sqr(m(iRow,iCol) - ((iRow == iCol)? 1 : 0));
//       }
//    }
//    return sqrt(fAcc)/N;
// }



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
#endif // INCLUDE_ABANDONED

void PrintVector(FDenseVector const &ew, std::string Title)
{
   std::cout << fmt::format("      {}", Title) << std::endl;
   for (size_t i = 0; i != size_t(ew.size()); ++ i) {
      std::cout << fmt::format("    {:>8d}  {:24.12f}  (={:8.2e})\n", i+1, double(ew[i]), double(ew[i]));
   }
}

void PrintMatrix(FDenseMatrix const &M, std::string Title)
{
   std::cout << fmt::format("      {}", Title) << std::endl;
   for (size_t iRow = 0; iRow != size_t(M.rows()); ++ iRow) {
      std::cout << fmt::format("    {:>8d}", iRow+1);
      for (size_t iCol = 0; iCol != size_t(M.cols()); ++ iCol)
         std::cout << fmt::format("  {:24.12f}", double(M(iRow,iCol)));
      std::cout << "\n";
   }
}



// Makes and returns (3,3)-shape matrix initialized from the explicitly given
// list of data elements. Data is provided in col-major order (i.e., first three
// arguments to the function describe the components of the column vector in
// m.col(0))
//
// NOTE:
// - This data element order is the *transpose* of what one would get when using
//   Eigen's '<<' operator to initialize the matrix elements. It is also the
//   transpose of the order in which data elements come in the list-of-list
//   constructors of numpy.array and sympy.Matrix.
// - Eigen::Matrix technically does have a braced-init-list constructor
//   ...however, it appears to only be allowed for compile-time matrix types
//   which represent either pure column vector vectors (i.e., have one column
//   only) or pur row vectors (i.e., have one row only), if I understand
//   correctly. When trying this one:
//
//        FRotationMatrix{0,1,2, 3,4,5, 6,7,8};
//
//   I definitely ran into a STATIC_ASSERT which's message I would interpret
//   like that.
FRotationMatrix Matrix3x3(FScalar const &a00, FScalar const &a10, FScalar const &a20,  FScalar const &a01, FScalar const &a11, FScalar const &a21,  FScalar const &a02, FScalar const &a12, FScalar const &a22)
{
   FRotationMatrix a;
   // first column vector
   a(0,0) = a00;
   a(1,0) = a10;
   a(2,0) = a20;

   // second column vector
   a(0,1) = a01;
   a(1,1) = a11;
   a(2,1) = a21;

   // third column vector
   a(0,2) = a02;
   a(1,2) = a12;
   a(2,2) = a22;

   return a;
}
