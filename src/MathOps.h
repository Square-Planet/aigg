#ifndef AIGG_MATH_OPS_H
#define AIGG_MATH_OPS_H

#include "Aigg.h"
#include <string>
#include "CxIntrusivePtr.h"
#include "CxNewtonAux.h" // for trust region solver

using ct_opt::FEwType;
using ct_opt::ComputeTrustRegionStep;
using ct_opt::EWTYPE_Eigh; // trust-region update for symmetric problem with pEw being eigenvalues (U diag(1/(ew_i + lambda)) U^T)
using ct_opt::EWTYPE_Svd; // trust-region update for non-symmetryic problem with pEw being singular values (V diag(sig_i/(sig_i^2 + lambda^2)) U^T)

struct FEighOrSvdImpl;

// Computes and holds the results of EITHER a singular value decomposition (SVD):
//
// A = U * diag(vals) * V.T
//
// OR of a symmetric/hermitian spectral decomposition (EIGH):
//
// A = V * diag(vals) * V.T
//
// In the second case, A must be square and symmetric/hermitian.
//
// Notes:
// - In either case, the diagonal values 'vals' (eigenvalues for EIGH,
//   singular values for SVD) will be sorted from LARGE to SMALL.
// - Note that in the case of a spectral decomposition, these may be NEGATIVE.
// - If a spectral decomposition is used, V will be square and U will be V.T.
//   For a SVD, neither U nor V will generally be square.
struct FEighOrSvd : public ct::FIntrusivePtrDest1
{
   enum FDecompType {
      DECOMP_Auto, // make spectral decomposition if 'A' is square and symmetric/hermitian, SVD otherwise.
      DECOMP_Eigh, // make spectral decomposition (requires 'A' to be square and symmetric/hermitian)
      DECOMP_Svd // make singular value decomposition (required if 'A' is not square or not symmetric)
   };
   
   explicit FEighOrSvd(FDenseMatrix const &A, FDecompType DecompType);
   ~FEighOrSvd(); // override
   
   // returns either left singular vectors U (svd) or eigenvector matrix V (eigh)
   Eigen::Ref<FDenseMatrix const> U() const;
   // returns either right singular vectors V (svd) or (same) eigenvector matrix V (eigh)
   Eigen::Ref<FDenseMatrix const> V() const;
   // returns either singular values (svd) or eigenvalues (eigh). In either
   // case, they will be ordered, and large values will come first. Note that in
   // eigh, the values may be negative (in svd this is absorbed in U column signs).
   Eigen::Ref<FDenseArray const> vals() const;

   // returns transpose of whatever V() returns.
   inline Eigen::Ref<FDenseMatrix const> Vt() const { return V().transpose(); }
   // returns transpose of whatever U() returns.
   inline Eigen::Ref<FDenseMatrix const> Ut() const { return U().transpose(); }
protected:
   FEighOrSvdImpl *p;
private:
   void operator = (FEighOrSvd const &); // not implemented
   FEighOrSvd(FEighOrSvd const &); // not implemented
};
typedef ct::TIntrusivePtr<FEighOrSvd>
   FEighOrSvdPtr;


// Computes and holds the result of a singular value decomposition:
//
// A = U * diag(sigma) * V.T
//
// where U and V's columns are composed orthogonal vectors, and
// the singular values sigma are >= 0 and sorted (reversely) by
// magnitude (largest first).
//
// If A is A (N,M)-matrix, U is a (N,K)-matrix and V a (M,K)-matrix,
// where K <= min(N,M).
struct FSvd : public ct::FIntrusivePtrDest1
{
   FDenseMatrix
      U, V;
   FDenseVector
      sigma;
   explicit FSvd(FDenseMatrix const &A);
   // deletes singular values below thr, and the corresponding columns in U and V.
   void Compress(FScalar thr);
   // solve lstsq(A, X) ... A X = Y
   FDenseMatrix Solve(FDenseMatrix const &Y);
   // solve lstsq(A.T, X)
   FDenseMatrix SolveT(FDenseMatrix const &Y);
};
typedef ct::TIntrusivePtr<FSvd>
   FSvdPtr;


// Computes and holds the result of a spectral decomposition of
// a symmetric/hermitian matrix H:
//
// H = V * diag(lambda) * V.T
//
// if H is a (N,N)-shape symmetric/hermitian matrix, then V is a (N,N)-shape
// orthgonal/unitary matrix (with eigenvectors of H in the columns).
//
// Eigenvalues are ordered; by default from smallest to largest (opposite to
// SVD). And unlike the singular values in SVD, eigenvalues can be negative.
struct FEigh : public ct::FIntrusivePtrDest1
{
   FDenseMatrix
      V;
   FDenseVector
      lambda;
   enum FEighOptions {
      LargeEwFirst = 0x01
   };
   explicit FEigh(FDenseMatrix const &H, unsigned Flags = 0);
};
typedef ct::TIntrusivePtr<FEigh>
   FEighPtr;



// A class for computing and holding the results of matrix decompositions
// which may be used for linear system solutions.
//
// Comments on how it is used concretely in the code:
// - It encapsulates both Cholesky decompositions and spectral decompositions in
//   a common interface. We need that because in various places either
//   decomposition may be called for, depending on circumstances which become
//   clear only at runtime.
//
// - Objects of the class are instanciated via call like
//
//   FLinearSolverPtr pSd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Cholesky)
//   FLinearSolverPtr pSd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_Eigh)
//   FLinearSolverPtr pSd = MakeLinearSolver(S, LINSOLVE_PositiveSymmetric_LDLT)
//
//   etc. These actually create subclasses implementing the virtual interface,
//   which are currently not explicitly accessible outside MathOps.cpp
//
// - In the code, both EIGH and CD are employed with full/half solver routines:
//
//   The defining feature of both, which is used in applications of this, is
//   that one can write a transform
//
//     A S^{-1} B
//
//   as dot(HalfSolve1(A.T).T, HalfSolve1(B)), where both HalfSolve1 calls
//   represent the same physical operation.
//   (- for spectral decomp: it's multiplication with a rank-reduced (U S^{-1/2})
//      with a (N,M)-shape matrix U of orthogonal columns,
//    - for the Cholesky decomp: it's a triangular solve with the L matrix,
//      because with S = L L^T, we have S^{-1} = (L L^T)^{-1} = (L^T)^{-1} L^{-1}).
//   That is used in some places in which, particularly, A = B, so
//   one makes a symmetric result with only one instead of two trafos,
//   and in such a way that the result is guranteed stable symmetric.
//
//   For this reason HalfSolve1 and HalfSolve2 are explicitly exposed in the
//   interface (in other decompositions, e.g., LU, SVD, or QR, these may be
//   doing conceptually different operations, and then this access to their
//   separate parts is likely less useful).
//
// - The class can represent other kinds of linear solvers, but atm
//   we only have uses for symmetric positive definite/semi-definite
//   solvers, and that are all which we have.
struct FLinearSolver : public ct::FIntrusivePtrDest1
{
   virtual FDenseMatrix HalfSolve1(Eigen::Ref<FDenseMatrix const> const &A) = 0;
   virtual FDenseMatrix HalfSolve2(Eigen::Ref<FDenseMatrix const> const &A) = 0;
   virtual FDenseMatrix Solve(Eigen::Ref<FDenseMatrix const> const &A) = 0;

   // should do the inverse of HalfSolve1, such that HalfMxm1(HalfSolve1(x)) == x. Default implementation just crashed with "not implemented"-error
   virtual FDenseMatrix HalfMxm1(Eigen::Ref<FDenseMatrix const> const &A);
   virtual FDenseMatrix HalfMxm2(Eigen::Ref<FDenseMatrix const> const &A);
   virtual FDenseMatrix Mxm(Eigen::Ref<FDenseMatrix const> const &A);
   
   inline size_t nEffectiveRank() const { return m_nEffectiveRank; }
protected:
   size_t
      m_nEffectiveRank;
};
typedef ct::TIntrusivePtr<FLinearSolver>
   FLinearSolverPtr;

// see comments in MathOps.cpp regarding details of what the solver classes do
// and what they are good and/or bad for.
enum FLinearSolverType {
   // make FSymmetricPositiveSolverEigh: based on spectral decomposition
   // S = Ev*diag(ew)*Ev.T (with unitary E)
   LINSOLVE_PositiveSymmetric_Eigh,
   // make FSymmetricPositiveSolverCd: based on regular Cholesky
   // decomposition S = L*L.T (i.e., not explicitly pivoted, and neither
   // with diagonals extracted. *Will* blow up if faced with hard singularities.
   LINSOLVE_PositiveSymmetric_Cholesky,
   // make FSymmetricPositiveSolverLDLT: based on pivoted semi-definite
   // Cholesky-like decomposition S = P.T L diag(d) L.T P (with permutation matrix P)
   LINSOLVE_PositiveSymmetric_LDLT
   
   // Note: in many cases one probably *could* add a rank-revealing QR
   // decomposition instead of the spectral decomposition; that should be
   // enough. But for the moment I prefer the spectral decomposition, because it
   // makes clearer what is going on.
};

// constructs and returns a matrix decomposition of S of the specified type, for
// linear system solutions (technically could also be used for some other stuff,
// but we don't need that here).
FLinearSolverPtr MakeLinearSolver(FDenseMatrix const &S, FLinearSolverType SolverType, bool AllowMxm = false);


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
// - The inverse of this is CayleyA(U) = 2 * la.solve(I + U, U - I)
FMatrixRR CayleyU(FMatrixRR A, FScalar tau);

// form a (orthogonal) rotation matrix R from an input antisymmetric matrix (tau*A),
// similar to R = exp(tau*A)
FMatrixRR BuildRotation(FMatrixRR A, FScalar tau);


// convenience function to solve a guaranteed-nonsingular generic (3,3)-shape linear equation system.
FVectorR Solve(FMatrixRR const &M, FVectorR const &rhs);


void PrintMatrix(FDenseMatrix const &M, std::string Title);
void PrintVector(FDenseVector const &ew, std::string Title);
FScalar Rmsd(FDenseMatrix const &m);
FScalar Rmsd(FDenseMatrix const &m, size_t nDegreesOfFreedom);
FScalar RmsdFromIdentity(FDenseMatrix const &m);
// FScalar Rmsd(FDenseVector const &v);
// FScalar Rmsd(FDenseVector const &v, size_t nDegreesOfFreedom);

#ifdef INCLUDE_ABANDONED
// // this is an auxiliary structure for computing trust-region Newton/Quasi-Newton steps
// // It's c/p stuff from MicroScf's geometry optimizer...
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
#endif // INCLUDE_ABANDONED

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
FRotationMatrix Matrix3x3(FScalar const &a00, FScalar const &a10, FScalar const &a20,  FScalar const &a01, FScalar const &a11, FScalar const &a21,  FScalar const &a02, FScalar const &a12, FScalar const &a22);

#endif // AIGG_MATH_OPS_H
