// This is the main header for the sphere grid optimization program (aigg). It
// contains mostly type definitions and binds for headers & libraries used
// throughout the program.
#ifndef ANGULAR_INTEGRATION_GRID_GENERATOR_H
#define ANGULAR_INTEGRATION_GRID_GENERATOR_H

#include "CxDefs.h" // for assert mainly
#include <stddef.h> // for size_t/ptrdiff_t
#include <iostream>
#include <string>
#include <vector>

// Next one sets up scalar types and default epsilons/numerical
// thresholds.
// - NOTE: it must go before any eigen stuff, because it may declare overloads
//   and/or templates involving the scalar type (FScalar) we will be doing our
//   calculations with.
// - ...and these no longer can be legally introduced at the time of binding the
//   templates using them.
#include "ScalarTypes.h"

#include <math.h>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include "format.h"
#include "CxIo.h"
#include "CxIterTools.h"
#include "CxPrintLevel.h"
using cx::iter::enumerate;
using ct::FPrintLevel;


// Set up global integral constants defining the space dimension and the
// type of phase factors we support.
//
// Comments:
// - Currently this program can only compute grids for the unit sphere of R^3
//   (which is itself a two-dimensional object, and therefore commonly called S^2).
//
// - However, the core algorithm and many components are rather general, and
//   should allow optimizing more general objects, including such of other
//   geometry (e.g., domes instead of full spheres), other manifold spaces
//   (e.g., the unit ball instead of the unit sphere), host space dimension
//   (e.g., R^4 hyperspheres), and even other underlying scalar fields (e.g.,
//   C^2 unit spheres subject to SU(2) symmetry groups).
//
// - For this reasons many program parts are written for general space
//   dimensions, and some contain other generalizations, using these constants
//   to control the actually used specialization. For *actual* unit sphere grid
//   constructions these are, strictly speaking, not needed, but should simplify
//   the adaption to other problems in the future---at the very least they
//   should simplify locating the parts of the code which require adjustments.
enum {
   // number of spatial (Cartesian) axes in our host space R^n.
   AIG_SPACE_AXES = 3,
   // maximum number of different discretized phase angles a scalar matrix
   // element m_{ij} of a symmetry operator can take to create phase factors
   // e^{(2πı/N)*k}. For real numbers, the maximum is 2 (the two possible
   // phase factors being e^0 = +1 and e^{ıπ} = -1)
   AIG_SYMOP_MAX_PHASE_ANGLES = 2
};


// general size matrices/arrays/vectors
typedef Eigen::Matrix<FScalar, Eigen::Dynamic, Eigen::Dynamic>
   FDenseMatrix;
typedef Eigen::Matrix<FScalar, Eigen::Dynamic, 1>
   FDenseVector;
typedef Eigen::Array<FScalar, Eigen::Dynamic, Eigen::Dynamic>
   // almost the same as a matrix, but with somewhat different default operations.
   FDenseArray;
typedef Eigen::Array<FScalar, Eigen::Dynamic, 1>
   FDenseArray1;


// Matrices/arrays/vectors/points for Real-space (3D) coordinates
typedef Eigen::Matrix<FScalar, AIG_SPACE_AXES, 1>
   FPoint;
typedef Eigen::Matrix<double, AIG_SPACE_AXES, 1>
   FPointd;
typedef Eigen::Matrix<FScalar, AIG_SPACE_AXES, 1>
   FVectorR; // R -> "real space" = host space. Will change if AIG_SPACE_AXES (see above) changes.
typedef Eigen::Matrix<FScalar, AIG_SPACE_AXES, AIG_SPACE_AXES>
   FMatrixRR;
typedef Eigen::Matrix<FScalar, 3, 1>
   FVector3;
typedef Eigen::Matrix<FScalar, 3, 3>
   FMatrix33;
// typedef std::vector<FPoint>
//    FPointList;
typedef Eigen::Matrix<FScalar, AIG_SPACE_AXES, Eigen::Dynamic>
   FPointArray;


// 3D-space rotations
typedef Eigen::Quaternion<FScalar>
   FQuaternion;
typedef FQuaternion::RotationMatrixType
   FRotationMatrix;
// ^- Notes:
//    - It is just a typedef for Matrix<Scalar,Dim,Dim> (in RotationBase.h);
//      does not enforce any special properties on the matrix.
//    - In the meantime, we may also use this to represent some symmetry
//      operations with a determinant of -1 instead of +1
//      (the matrices are still members of O(3), but no longer necessarily SO(3)).
//      These are not *technically* rotations, but we keep the term anyway.


typedef std::vector<FQuaternion>
   FQuaternionList;
typedef std::vector<FRotationMatrix>
   FRotationMatrixList;

typedef Eigen::Matrix<FScalar,AIG_SPACE_AXES,AIG_SPACE_AXES>
   // type of SO(3) rotation generators: anti-symmetric (3,3)-shape matrices.
   FRotationGenerator;


template<class T> T sqr(T const &x) { return x*x; }
template<class T> T pow2(T const &x) { return x*x; }
template<class T> T pow3(T const &x) { return x*x*x; }
template<class T> T pow4(T const &x) { T xx = x*x; return xx*xx; }



// makes a string describing this version of aigg---for data export purposes.
std::string MakeAiggVersionString();
// makes a string describing the current time---for data export purposes.
std::string MakeCurrentTimeString();
// concatenate a string `pattern` together for `count` times (like Python `count * pattern` would do);
// if given, append `left` to the left and `right` to the right.
std::string Repeat(size_t count, std::string const &pattern, std::string const &left = std::string(), std::string const right = std::string());


extern ct::FLogStdStream io;

#ifdef INCLUDE_ABANDONED
// namespace io_detail { extern fmt::MemoryWriter g_IoWriteProxy; }
// void IoFlush();
// 
// void Write(fmt::BasicStringRef<Char> format) {
//    io_detail::g_IoWriteProxy << format << '\n'; IoFlush();
// }
// void Write(fmt::BasicStringRef<Char> format, fmt::ArgList args) {
//    fmt::BasicFormatter<Char>(io_detail::g_IoWriteProxy).format(format, args);
//    io_detail::g_IoWriteProxy << '\n';
//    IoFlush();
// }
// 
// void WriteCount(std::string const &Name, ptrdiff_t iValue, fmt::BasicStringRef<char> Annotation="");
// void WriteInfo(std::string const &Name, fmt::BasicStringRef<char> Info, fmt::BasicStringRef<char> Annotation="");
// void WriteInfoExpf(std::string const &Name, FScalar const &f, fmt::BasicStringRef<char> Annotation="");
// void WriteResult(std::string const &Name, FScalar const &f, fmt::BasicStringRef<char> Annotation="");
// void WriteTiming(std::string const &Name, FScalar const &TimeInSec, size_t nCount = 0);
// void WriteLine();
// 
// 
// extern char const
//    *p1ResultFmt,
//    *p1ResultFmtAnnotated,
//    *p1CountFmt,
//    *p1CountFmtAnnotated,
//    *p1InfoFmt,
//    *p1InfoFmtAnnotated,
//    *p1InfoExpfFmt,
//    *p1InfoExpfFmtAnnotated,
//    *p1TimingFmt,
//    *p1TimingFmtAnnotated;
#endif // INCLUDE_ABANDONED




// Upcoming: A simple trick to determine what kind of type the compiler deduces
// for a type/expression. I stole it from here:
// https://diego.assencio.com/?index=05e466a9a1e6dfc958769f80bf986f65.
// 
// To this end, all we do is declare a template class:
//
//     template<class T> class MagicTypeRevealer;
//
// ...but then just leave it at that. That is, we grandly proclaim its
// existence, but never actually define it. So if we do something like this:
//
//     template<class T>
//     struct TMyCoolTemplateType {
//         MagicTypeRevealer<T> _wheeee;
//     }
//
// or this:
//
//     auto x = ...;
//     MagicTypeRevealer<decltype(x)> _wheeee;
//
// we will have a ill-formed program, because for these the compiler would need
// to instantiate the template class MagicTypeRevealer<T> −− which it cannot
// do, because it is undefined. And this is good! Because the compiler will
// verbosely complain about this abhorrently unacceptable state of affairs, and
// in the process produces a nifty error message describing (among other things)
// exactly what type T is.
template<class T> class MagicTypeRevealer;



#endif // ANGULAR_INTEGRATION_GRID_GENERATOR_H
