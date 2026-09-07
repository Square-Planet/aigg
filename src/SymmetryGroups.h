#ifndef SYMMETRY_GROUPS_H
#define SYMMETRY_GROUPS_H

#include "Aigg.h"

#include <memory>
#include <string>
// #include <tuple>

#include "CxIntrusivePtr.h"
#include "PolyXyz.h" // really only for FMonomialN objects for packing/passing parameters on target space queries
#include "AxisPermuflections.h" // for representing symmetry operators which can be efficiently applied to monomials

// typedef std::vector<bool>
//    FBoolArray;

FRotationMatrixList ConvertToMatrixList(FQuaternionList const &ql);

struct FOrbitType : public ct::FIntrusivePtrDest1
{
   size_t
      // internal orbit type id (0: first type, 1: second type, ...).
      //
      // NOTE:
      // - The orbit order is not arbitrary. The classifier will assume
      //   that orbits are ordered from most-general first to most special
      //   last.
      //
      // - This is important for resolving points which lie on multiple
      //   different types of orbits simultaneously.
      //
      //   For example, *every* point r⃗ lies on the ‘general’ orbit, but
      //   the point may also additionally lie on special orbits with one or
      //   zero geometric degrees of freedom. In this case, it will be assigned
      //   the orbit type of the highest `iType` (i.e., most special one) of
      //   which it fulfills the orbit condition.
      iType;
   char const
      // description of the orbit type; e.g., "general", "vertex",
      // "center-of-face", "center-of-edge"
      *pDesc;
   size_t
      // number of elements in a group orbit of this type
      // (e.g., for icosahedral group (I): 60 for general, 30 for center-of-
      // edge, 20 for center-of-face, 12 for vertex)
      Length,
      // number of continuous geometric degrees of freedom which a point on the
      // sphere still has while constrained to lie on an orbit of this type.
      // For three spatial dimensions, this can be 0, 1, or 2 (see below)
      nFreeParameters;

   FOrbitType(size_t iType_, char const *pDesc_, size_t Length_, size_t nFreeParameters_);
   virtual ~FOrbitType();

   // returns whether the point 'p⃗' on the unit sphere lies on this
   // particular type of orbit, up to a given tolerance `fEpsilon`
   virtual bool IsOnOrbit(FPoint const &vPt, FScalar const &fEpsilon) const = 0;
   // project point vPt to lie *exactly* on the type of orbit represented by *this
   //
   // If `AllowDiscreteChange == true`, this may replace a point on the orbit by
   // the canonical representative of this orbit (i.e., perform a non-continuous
   // change to `vPt`). otherwise, it will attempt to move the point to the
   // closest orbit representative of this orbit (may be relevant if using
   // Krylov subspace convergence accelerators/stabilizers like DIIS or similar)
   //
   // Default implementation just normalizes `vPt`, such that it definitely lies
   // on the unit sphere afterwards.
   virtual void ProjectToOrbit(FPoint &vPt, bool AllowDiscreteChange) const;

   // project out subspaces of tangent vectors of vPt which would lead out of
   // the given orbit type.
   // If a `pErrSq` ≠ 0 is given, some representative of the amount of fixes
   // performed will be accumulated to `*pErrSq`.
   //
   // Default implementation projects out components of the tangent vectors
   // which are colinear with `vPt`, and accumulates the corresponding errors.
   // This projection is valid for all points `vPt` constrained to a sphere
   // (since `vPt` lies on the sphere, tangent vectors keeping it on the sphere
   // cannot possibly have a component colinear with it).
   virtual void ProjectTangentVectors(FVector3 **ppvTangents, size_t nTangents, FPoint const &vPt, FScalar *pErrSq) const;
   virtual void ProjectRotation(FVector3 &vAxisAndAngle, FPoint const &vPt, FScalar *pErrSq) const;
private:
   void operator = (FOrbitType const &) const; // not implemented
   FOrbitType(FOrbitType const &); // not implemented
};
typedef ct::TIntrusivePtr<FOrbitType>
   FOrbitTypePtr;
typedef ct::TIntrusivePtr<FOrbitType const>
   FOrbitTypeCptr;

typedef std::vector<FOrbitTypePtr>
   FOrbitTypeList;
typedef std::vector<FOrbitType const*>
   FOrbitTypeRawPtrList;


// general orbit, with two geometric degrees of freedom (see below)
struct FOrbitType_General : public FOrbitType
{
   // Construct a general orbit of the group; points 'r' on this orbit are only
   // constrained by the condition that \norm{\vec r} = 1, and this orbit has
   // the same length as the order of the group.
   FOrbitType_General(size_t iType_, char const *pDesc_, FRotationMatrixList const &SymmetryOps);
   ~FOrbitType_General();

   bool IsOnOrbit(FPoint const &vPt, FScalar const &fEpsilon) const; // override
   void ProjectToOrbit(FPoint &vPt, bool AllowDiscreteChange) const; // override
   void ProjectTangentVectors(FVector3 **ppvTangents, size_t nTangents, FPoint const &vPt, FScalar *pErrSq) const; // override
};


// orbit of points restricted to group images of a given plane, with one geometric degree of freedom (see below)
struct FOrbitType_SpecialPlanes : public FOrbitType
{
   // Construct a special orbit, with points on the orbit still having one
   // continuous geometric degree of freedom. For example, hypothetically, if
   // points were restricted to a special line in the fundamental domain, this
   // would be the case.
   //
   // One can show that all orbits of this type originate from reflection
   // operations (concretely, reflection at some plane through the
   // coordinate origin), and are given by a great circle around the unit
   // sphere, which is in turn defined by the two conditions that:
   // (a) any point r⃗ of the orbit lies on the unit sphere, and
   // (b) any point r⃗ of the orbit fulfills a plane equation ⟨n⃗,r⃗⟩=0,
   // where n⃗ denotes the normal of the reflection plane (through the
   // coordinate origin), or any of its group images.
   FOrbitType_SpecialPlanes(size_t iType_, char const *pDesc_, FVector3 const &vPlaneNormal, FRotationMatrixList const &SymmetryOps);
   ~FOrbitType_SpecialPlanes();

   bool IsOnOrbit(FPoint const &vPt, FScalar const &fEpsilon) const; // override
   void ProjectToOrbit(FPoint &vPt, bool AllowDiscreteChange) const; // override
   void ProjectTangentVectors(FVector3 **ppvTangents, size_t nTangents, FPoint const &vPt, FScalar *pErrSq) const; // override
   void ProjectRotation(FVector3 &vAxisAndAngle, FPoint const &vPt, FScalar *pErrSq) const; // override
protected:
   typedef std::vector<FVector3>
      FNormalList;
   FNormalList
      // images of plane normals under the operations of the group;
      // each point r⃗ on the orbit is constrained to fulfill
      //
      //    ⟨n⃗,r⃗⟩ = 0
      //
      // where n⃗ denotes one of the plane normals in this list.
      m_PlaneNormals;

   // find and return the plane normal n⃗ in `m_PlaneNormals` for which `vPt`
   // has the smallest deviation of the plane equation ⟨n⃗,r⃗⟩ = 0.
   FVector3 const &FindClosestPlane(FPoint const &vPt) const;
};


// orbit of isolated points, with no geometric degrees of freedom (see below)
struct FOrbitType_SpecialPoints : public FOrbitType
{
   // construct a special orbit consisting of `Length` isolated points;
   // a point belongs to this orbit if and only if it is identical to one
   // of the group images of a specified orbit representative. So there
   // are no continuous geometric degrees of freedom.
   FOrbitType_SpecialPoints(size_t iType_, char const *pDesc_, FPoint const &vExamplePoint, FRotationMatrixList const &SymmetryOps);
   ~FOrbitType_SpecialPoints();

   bool IsOnOrbit(FPoint const &vPt, FScalar const &fEpsilon) const; // override
   void ProjectToOrbit(FPoint &vPt, bool AllowDiscreteChange) const; // override
   void ProjectTangentVectors(FVector3 **ppvTangents, size_t nTangents, FPoint const &vPt, FScalar *pErrSq) const; // override
protected:
   typedef std::vector<FPoint>
      FPointList;
   FPointList
      // images of example points under the operations of the group; each point
      // r⃗ on the orbit is constrained to fulfill
      //
      //    r⃗ = p⃗
      //
      // for one of the points p⃗ in this list.
      m_OrbitPoints;

   // find and return the point p⃗ in `m_OrbitPoints` for which `vPt`
   // has the smallest deviation from ‖p⃗ − r⃗‖ = 0.
   FPoint const &FindClosestPoint(FPoint const &vPt) const;
};





struct FGeneratorInfo {
   unsigned
      nOrder; // smallest non-negative integer 'n' such that pow(op,n) == id
   char const
      *pDesc;
   FRotationMatrix
      Op;

   FGeneratorInfo(unsigned nOrder_, char const *pDesc_, FRotationMatrix const &Op_);
};
typedef std::vector<FGeneratorInfo>
   FGeneratorInfoList;




enum FPointGroupNameType {
   NAMETYPE_Symmetry,
   NAMETYPE_Grid,
   NAMETYPE_Group,
   NAMETYPE_Symbol
};


struct FTriangle
{
   enum FBarycentricCoordType {
      // express triangular barycentric coordinates in terms of standard planar geometry
      GEOMTYPE_Planar,
      // express barycentric coordinates in terms of spherical geometry
      // (with all points on the sphere (will be internally normalized)
      // and the areas representing spherical triangle areas)
      GEOMTYPE_Spherical
   };

   FVector3
      v0, v1, v2;

   // return the barycentric representation of the 3D point vCartesianCoords;
   //
   // If `GeomType` is `GEOMTYPE_Planar`, this return [τ₀,τ₁,τ₂] such that
   //
   //     vCartesianCoords = τ₀⋅v⃗₀ + τ₁⋅v⃗₁ + τ₂⋅ v⃗₂.
   //
   // For `GEOMTYPE_Spherical`, it will solve the corresponding spherical geometry
   // equations and return the area fractions [τ₀,τ₁,τ₂] which
   // denote the corresponding spherical triangle areas of the split
   // (in this case `vCartesianCoords`,`v0`,`v1`,`v2` should lie on the sphere; they
   // will be internally normalized before the transformation is applied)
   //
   // Will raise exception if the equation system cannot be solved (e.g., because
   // the input point does not actually lie on the triangle)
   FPoint MakeBarycentric(FPoint const &vCartesianCoords, FBarycentricCoordType GeomType) const;
   FPoint MakeCartesian(FPoint const &vBarycentricCoords, FBarycentricCoordType GeomType) const;

   // returns whether the point `vPt`, given in terms of Cartesian coordinates,
   // lies inside the triangle (up to a numerical accuracy of `Epsilon`).
   //
   // It is assumed that the point given by `vCartesianCoords` lies on the
   // triangle plane.
   bool IsInside_Cartesian(FPoint const &vCartesianCoords, FScalar const &Epsilon, FBarycentricCoordType GeomType) const;
   // returns whether the point `vPt`, given in terms of Barycentric coordinates,
   // lies inside the triangle (up to a numerical accuracy of `Epsilon`).
   bool IsInside_Barycentric(FPoint const &vBarycentricCoords, FScalar const &Epsilon) const;

   bool operator == (FTriangle const &other) const { return v0 == other.v0 && v1 == other.v1 && v2 == other.v2; }
   bool operator != (FTriangle const &other) const { return !(*this == other); }
};


namespace details { struct FGensAndOps; }
// ⬑ `detail` conflicts with `boost:multi_precision::detail`. apparently I have
// a using namspace…. somewhere instead of relying on ADL. Do not want to track
// it atm. It is dangerous if I get it wrong and end up using `double` math instead
// of high-precision floats in some instances.


struct FPointGroup : public ct::FIntrusivePtrDest1
{
//    FQuaternionList const &GetOpsQuatForm() const { return m_SymmetryOpsQuaternionForm; };
   FRotationMatrixList const &GetOps() const { return m_SymmetryOps; };

   virtual size_t CountOrbitLength(FPoint const &p) const;
   // determine the type of orbit the representative p belongs to (this is not a free operation).
   virtual FOrbitType const &GetOrbitType(FPoint const &p) const;

   // array version of GetOrbitType(). Determine orbit type for all orbits.
   FOrbitTypeRawPtrList GetOrbitTypes(FPointArray const &SeedPoints) const;

   // used during optimization progress to remove rotation parameters
   // which would lead out of the prescribed orbit type.
   // Rg_Ax_sq is Rg * Ax * s, where 's' is the seed point, Rg the rotation matrix for the current
   // symmetry operation, and Ax/y/z are anti-symmetric generator matrices for rotations around the x/y/z axes.
   // Rg_sq is the transformed seed point itself.
   virtual FScalar ProjectRotationToOrbitType(FPointArray &Rg_Ax_sq, FPointArray &Rg_Ay_sq, FPointArray &Rg_Az_sq, FPointArray const &Rg_sq, FOrbitType const *const *pOrbitTypes) const;
   // project either point itself (ps != 0) to the given type of orbit, or the update in the point (pdx != 0).
   // sq is the current value of the seed point which is supposed to move.
   virtual FScalar ProjectPointOrStepToOrbitType(FPoint *ps, FVector3 *pdx, FVector3 const &sq, FOrbitType const &OrbitType, bool AllowDiscreteChange) const;

   virtual std::string Name(FPointGroupNameType NameType = NAMETYPE_Group) const = 0;
#ifdef INCLUDE_ABANDONED
//    enum FBarycentricCoordType {
//       // attempt to convert coordinates in terms of 01c coordinates, which
//       // represent a
//       BARYC_01c
//       // attempt to convert coordinates in terms of 012 representing the closest
//       // thing to an equilateral triangle we have on the solid, which is
//       // compatible with the fundamental region. E.g., in an icosahedron, the
//       // three vertices of one of its faces. in a cube, three adjacent face
//       // centers (that would be the coordinates of faces of the dual octahedron
//       // of the cube)
//       BARYC_012,
//    };

//    // input: [0]: tau0, [1]: tau1, [2]: tauc.
//    // For point groups O, T, I, this will return 'true' provided all three coordinates are >= -epsilon.
//    // for points groups with additional symmetry, Td, Th, Oh, Ih, this will take these additional
//    // symmetries into account.
//    virtual bool IsInsideFundamentalRegion(FPoint Barycentric);
//
//
//    // input: [0]: tau0, [1]: tau1, [2]: tauc.
//    virtual FPoint ConvertBarycentricToCartesian(FPoint Barycentric);
#endif // INCLUDE_ABANDONED

   // identifies various triangles used in conjunction with the fundamental region
   // identification and the creation of initial guesses. Some of those
   // may be chosen identical with each other.
   //
   // NOTE: Even if differing from each other, all triangles MUST lie in the
   // same plane! This is required for converting points point in space
   // given in terms of barycentric coordinates of each other.
   enum FTriangleRegionType {
      // the 012 region should represent the closest thing to an equilateral
      // triangle we have on the solid, which is compatible with the fundamental
      // region. E.g., in an icosahedron, octahedron, or tetrahedron, the three
      // vertices of one of its faces. In an octahedral group set up as in a
      // cube, three adjacent face centers (that would be the coordinates of
      // faces of the dual octahedron of the cube)
      //
      // NOTE: region 012 *MUST* include the fundamental region!
      REGION_012,
      // the 01c region should represent vertices of a multiple of the
      // fundamental domain (triangle), which are used for ‘historic’ reasons
      // in conjunction with the Ahrend & Beylkin initial grid construction. In
      // most cases, v⃗₀ v⃗₁ of this would agree with m_Region012, and v⃗_c would be
      // the center of the face m_Region012 ((v⃗₀+v⃗₁+v⃗₂)/3).
      //
      // region 01c may or may not include the fundamental region (normally it does)
      REGION_01c,
      // vertices of fundamental domain (triangle) which exactly tile across
      // the surface of the platonic solid, without repetition.
      // For groups without reflections (T, O, I), this should be set
      // identical to m_Region01c.
      REGION_Fundamental,
      REGION_Count
   };

   FTriangle const &GetTriangle(FTriangleRegionType RegionType_) { assert(RegionType_ < REGION_Count); return m_Regions[RegionType_]; };
   // returns whether the triangle identified by REGION_Fundamental
   // lies completely inside region REGION_01c (or is identical to it).
   // This is the for T, Td, O, Oh, I, Ih, but not Th.
   //
   // Note: region `012` *MUST* include the fundamental region!
   bool IsFundamentialRegionIncludedIn01c() const { return m_Region01cIncludesFundamental; }


   // TODO: delete this. Replace by non-virtual `pMonomialSymmetry->SelectAsCanonialReprQ`
   // used in pre-screening Cartesian monomials `x^ix y^iy z^iz` in target function selection.
   // Possible criteria would be:
   // - return `false` if we know the monomial to have a non-zero group average
   // - return `false` if this is not the monomial we'd like to use as canonical
   //   representative of its group-averaged equivalent set.
   // If in doubt, return `true`. That's what the default implementation does (i.e., no prescreening).
   virtual bool CouldMonomialBeATargetFunction(unsigned ix, unsigned iy, unsigned iz) const;

   FAxisPermuflectionSubgroup const *pAxprs() const { return &*m_pAxprSubgroup; }
   FMonomialSymmetry const *pMonomialSymmetry() const { return &*m_pMonomialSymmetry; };
   enum {
      MaxAxes = AIG_SPACE_AXES, // normally 3 (=3d real space)
      MaxPhases = AIG_SYMOP_MAX_PHASE_ANGLES // normally 2 (=real numbers)
   };
   // number of Cartesian axes to span host space (=cartesian space dimension)
   unsigned nAxes() const { return MaxAxes; }
   // number of different phase factors in scalar field of symmetry ops (2? ⟶ phase ∈ {+1, −1})
   unsigned nPhases() const { return MaxPhases; }


   FPointGroup();
   virtual ~FPointGroup();
protected:
//    FQuaternionList
//       m_SymmetryOpsQuaternionForm;
   FGeneratorInfoList
      m_Generators;
   FRotationMatrixList
      m_SymmetryOps;
   FOrbitTypeList
      m_OrbitTypes;

   FTriangle
      m_Regions[REGION_Count];
   bool
      m_Region01cIncludesFundamental;

   FAxisPermuflectionSubgroupPtr
      m_pAxprSubgroup;
   FMonomialSymmetryPtr
      m_pMonomialSymmetry;

   void _AssignGeneratorsAndOps(FGeneratorInfoList &&Generators_, FRotationMatrixList &&GroupElements_);
   void _AssignGeneratorsAndOps(details::FGensAndOps &&Args_);
//    void _AssignGeneratorsAndOps(std::tuple<FGeneratorInfoList, FRotationMatrixList> &&Args_);
public:
   // add a container-like interface for accessing the controlled symmetry
   // operations as (3,3)-shape rotation matrices.
   typedef FRotationMatrixList::const_iterator const_iterator;
   typedef FRotationMatrixList::value_type value_type;
   size_t size() const { return m_SymmetryOps.size(); }
   bool empty() const { return m_SymmetryOps.empty(); }
   const_iterator begin() const { return m_SymmetryOps.begin(); }
   const_iterator end() const { return m_SymmetryOps.end(); }
   value_type const &operator[] (size_t iSymOp) const { return m_SymmetryOps[iSymOp]; }
};
typedef ct::TIntrusivePtr<FPointGroup>
   FPointGroupPtr;

FPointGroupPtr MakePointGroup(std::string const &Name);


struct FIcosahedralGroup : public FPointGroup
{
   enum FIcosahedralGroupType {
      // 60 elements, proper rotations only
      GROUP_RotationsOnly_I,
      // 120 elements, rotations + inversion at origin
      GROUP_Full_Ih
   };
   typedef FIcosahedralGroupType
      FSubType;

   explicit FIcosahedralGroup(FSubType SubType = GROUP_RotationsOnly_I);
   ~FIcosahedralGroup();
//    virtual FOrbitType const &GetOrbitType(FPoint const &p) const; // override

//    FScalar ProjectRotationToOrbitType(FPointArray &Rg_Ax_sq, FPointArray &Rg_Ay_sq, FPointArray &Rg_Az_sq, FPointArray const &Rg_sq, FOrbitType const **pOrbitTypes) const;
//    FScalar ProjectPointOrStepToOrbitType(FPoint *ps, FVector3 *pdx, FVector3 const &sq, FOrbitType const &OrbitType) const;
   std::string Name(FPointGroupNameType NameType) const;
   bool CouldMonomialBeATargetFunction(unsigned ix, unsigned iy, unsigned iz) const; // override

protected:
   FIcosahedralGroupType
      m_SubType;
};


struct FTetrahedralGroup : public FPointGroup
{
   enum FTetrahedralGroupType {
      // 12 elements, proper rotations only
      GROUP_RotationsOnly_T,
      // 24 elements, rotations + roto-reflections; no inversion center.
      // This is the “normal” tetrahedral-full group
      GROUP_FullAchiral_Td,
      // also 24 elements, rotations + inversion at origin; It's called Th, but it is not
      // really tetrahedral, but semi-dodecahedral. We call it ‘pyritohedral group’.
      GROUP_Chiral_Th
   };
   typedef FTetrahedralGroupType
      FSubType;

   explicit FTetrahedralGroup(FSubType SubType = GROUP_RotationsOnly_T);
   ~FTetrahedralGroup();
//    virtual FOrbitType const &GetOrbitType(FPoint const &p) const; // override

//    FScalar ProjectRotationToOrbitType(FPointArray &Rg_Ax_sq, FPointArray &Rg_Ay_sq, FPointArray &Rg_Az_sq, FPointArray const &Rg_sq, FOrbitType const **pOrbitTypes) const;
//    FScalar ProjectPointOrStepToOrbitType(FPoint *ps, FVector3 *pdx, FVector3 const &sq, FOrbitType const &OrbitType) const;
   std::string Name(FPointGroupNameType NameType) const;
   bool CouldMonomialBeATargetFunction(unsigned ix, unsigned iy, unsigned iz) const; // override
protected:
   FPointArray
      // stored for distinguishing orbit types
      m_TetrahedronVertices;
   FSubType
      m_SubType;
};


struct FOctahedralGroup : public FPointGroup
{
   enum FOctahedralGroupType {
      // 24 elements, proper rotations only
      GROUP_RotationsOnly_O,
      // 48 elements, rotations + inversion at origin
      GROUP_Full_Oh
   };
   typedef FOctahedralGroupType
      FSubType;
   explicit FOctahedralGroup(FSubType SubType = GROUP_RotationsOnly_O, bool SetupAsCube_=false);
   ~FOctahedralGroup();
//    virtual FOrbitType const &GetOrbitType(FPoint const &p) const; // override

//    FScalar ProjectRotationToOrbitType(FPointArray &Rg_Ax_sq, FPointArray &Rg_Ay_sq, FPointArray &Rg_Az_sq, FPointArray const &Rg_sq, FOrbitType const **pOrbitTypes) const;
//    FScalar ProjectPointOrStepToOrbitType(FPoint *ps, FVector3 *pdx, FVector3 const &sq, FOrbitType const &OrbitType) const;
   std::string Name(FPointGroupNameType NameType) const;
   bool CouldMonomialBeATargetFunction(unsigned ix, unsigned iy, unsigned iz) const; // override
protected:
   FSubType
      m_SubType;
   bool
      // if set, set up the fundamental region 01c in terms of a cube.
      // If not, set it up in terms of an octahedron.
      m_FundamentalRegionIsOfCube;
};


#endif // SYMMETRY_GROUPS_H
