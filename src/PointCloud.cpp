#include "Aigg.h"
#include <iostream>
#include <stdexcept>
#include "PointCloud.h"

// static double const s_RandomRotation[9] = {-0.0093011468749436572, -0.431409304997183, 0.9021083639057218, 0.99982454646126095, -0.018681167373731322, 0.0013748742082411725, 0.016259303808730044, 0.90196287370789408, 0.43150736898951497};

typedef ct::TIntrusivePtr<FPointCloudNode>
   FPointCloudNodePtr;


FPointCloud::FPointEntry::FPointEntry()
   : iPoint(0xbadc0de), nRefs(0)
{}


FPointCloud::FPointEntry::FPointEntry(FPoint const &v_, ptrdiff_t i_, ptrdiff_t nRefs_)
   : v(v_), iPoint(i_), nRefs(nRefs_)
{}


struct FNodeStats {
   size_t Depth;
   size_t nMinPoints, nMaxPoints;
   size_t nPointsTotal;
};

struct FPointCloudNode : public ct::FIntrusivePtrDest1 {
   explicit FPointCloudNode(FPointCloud *parent_);
   typedef FPointCloud::FPointEntry
      FPointEntry;
   FPointEntry const &Insert(FPoint const &p, ptrdiff_t iPoint);

   void GetStats(FNodeStats &n) const;
protected:
   typedef std::vector<FPointEntry>
      FPointEntryList;
   FPointCloud
      *m_pParentObj;
   FVector3
      m_vPlane;
   FScalar
      m_fDivider;
   FPointCloudNodePtr
      m_pFrontNode,
      m_pBackNode;
   FPointEntryList
      m_Points;
   void SplitNode();
   FPointEntry const &InsertOnThisNode(FPoint const &p, ptrdiff_t iPoint, bool PointIsOnPlane);
   inline bool HaveLeafs() const;
   inline FScalar CalcPlaneDist(FPoint const &p) const;
};

void FPointCloudNode::GetStats(FNodeStats &ns) const
{
   if (HaveLeafs()) {
      FNodeStats
         f, b;
      m_pFrontNode->GetStats(f);
      m_pBackNode->GetStats(b);
      ns.Depth = 1 + std::max(f.Depth, b.Depth);
      ns.nMaxPoints = std::max(std::max(f.nMaxPoints, b.nMaxPoints), m_Points.size());
      ns.nMinPoints = m_Points.size();
      if (f.nMinPoints != 0) ns.nMinPoints = std::min(std::max(ns.nMinPoints,f.nMinPoints), f.nMinPoints);
      if (b.nMinPoints != 0) ns.nMinPoints = std::min(std::max(ns.nMinPoints,b.nMinPoints), b.nMinPoints);
      ns.nPointsTotal = m_Points.size() + f.nPointsTotal + b.nPointsTotal;
   } else {
      ns.Depth = 1;
      ns.nMinPoints = m_Points.size();
      ns.nMaxPoints = m_Points.size();
      ns.nPointsTotal = m_Points.size();
   }
}



FPointCloud::FPointCloud()
{
   SetEpsilon(1e-10);
   m_nMaxPointsPerNode = 12;
   m_pRootNode = 0;
}

FPointCloud::~FPointCloud()
{
//    if (m_pRootNode) {
//       FNodeStats ns;
//       m_pRootNode->GetStats(ns);
//       std::cout << fmt::format("     ~FPointCloud:  nPts = {}  Depth = {}  MinPt = {}  MaxPt = {}", ns.nPointsTotal, ns.Depth, ns.nMinPoints, ns.nMaxPoints) << std::endl;
//    }
   delete m_pRootNode;
   m_pRootNode = 0;
}


void FPointCloud::SetEpsilon(FScalar fEpsilon_)
{
   m_fEpsilon = fEpsilon_;
   m_fEpsilonSq = m_fEpsilon * m_fEpsilon;
}


bool FPointCloud::IsApproxEqual(FPoint const &a, FPoint const &b) const
{
   return (a - b).squaredNorm() <= m_fEpsilonSq;
}



FPointCloudNode::FPointCloudNode(FPointCloud *parent_)
   : m_pParentObj(parent_)
{
   m_Points.reserve(m_pParentObj->m_nMaxPointsPerNode);
   m_vPlane.setZero();
}


bool FPointCloudNode::HaveLeafs() const
{
   assert((m_pFrontNode.get() != 0) == (m_pBackNode.get() != 0));
   return m_pFrontNode.get() != 0;
}


FScalar FPointCloudNode::CalcPlaneDist(FPoint const &p) const
{
   return this->m_vPlane.dot(p) - m_fDivider;
}


FPointCloud::FPointEntry const &FPointCloudNode::Insert(FPoint const &p, ptrdiff_t iPoint)
{
   // check if we have an entry for this point here already
   for (size_t i = 0; i < m_Points.size(); ++ i) {
      if (m_pParentObj->IsApproxEqual(m_Points[i].v, p)) {
         m_Points[i].nRefs += 1;
         return m_Points[i];
      }
   }


   // do we have child nodes? (if there are no child nodes, there is also no plane)
   if (HaveLeafs()) {
      // classify point according to dividing plane
      FScalar
         f = CalcPlaneDist(p);
      bool
         IsOnPlane = abs(f) < m_pParentObj->m_fEpsilon;
      // almost exactly on the plane? leave it here (not likely to happen).
      if (IsOnPlane) {
         // (note: this may exceed the intended maximum number of points
         // per node -- but very unlikely and not actually that bad)
         return InsertOnThisNode(p, iPoint, true);
      } else if (f > 0)
         return m_pFrontNode->Insert(p, iPoint);
      else {
         assert(f < 0);
         return m_pBackNode->Insert(p, iPoint);
      }
   } else {
      return InsertOnThisNode(p, iPoint, false);
   }
}


FPointCloud::FPointEntry const &FPointCloudNode::InsertOnThisNode(FPoint const &p, ptrdiff_t iPoint, bool PointIsOnPlane)
{
   assert(!HaveLeafs() || PointIsOnPlane);
   // (note: when getting here we already checked if the point lies on the current node.
   // it does not. Add new point entry has to be added, and it has to be added inside the
   // *current* node, as there are no children)

   if (HaveLeafs() || m_Points.size() < m_pParentObj->m_nMaxPointsPerNode) {
      // add point to this node's list.
      m_Points.push_back(FPointEntry(p, iPoint, 1));
      return m_Points.back();
   } else {
      // make leaf nodes
      SplitNode();
      // try to insert recursively again.
      return Insert(p, iPoint);
   }
}


void FPointCloudNode::SplitNode()
{
   // compute center of mass and inertial tensor orientation for
   // points on this plane.
   FPoint
      vCenterOfMass = FPoint::Zero(AIG_SPACE_AXES);
   for (size_t iPoint = 0; iPoint != m_Points.size(); ++ iPoint) {
      vCenterOfMass += m_Points[iPoint].v;
   }
   vCenterOfMass *= FScalar(1)/m_Points.size();
   FMatrixRR
      mInertial;
   mInertial.setZero();
   for (size_t iPoint = 0; iPoint != m_Points.size(); ++ iPoint) {
      FPoint
         vRelPos = m_Points[iPoint].v - vCenterOfMass;
//       FScalar
//          fRsq = vRelPos.squaredNorm();
      for (size_t i = 0; i < AIG_SPACE_AXES; ++ i)
         for (size_t j = 0; j < AIG_SPACE_AXES; ++ j)
//             mInertial(i,j) += FScalar((i==j)? 1 : 0)*Rsq - vRelPos[i]*vRelPos[j];
            // ^- don't really need this here. For just getting the axes, the term ~identity
            //    is not helpful. We'll also invert the other one -- with this orientation
            //    of the outer product term, we want to get the axis which produces the *maximum*
            //    eigenvalue (rather than the minimum eigenvalue, as we'd need in the inertial case)
            mInertial(i,j) += vRelPos[i]*vRelPos[j];
   }
   // find principal axes
   Eigen::SelfAdjointEigenSolver<FMatrixRR>
      es(mInertial, Eigen::ComputeEigenvectors);
   FMatrixRR const
      &mPrincipalAxes = es.eigenvectors();

   // for dividing plane normal, use direction of highest "inertia" (i.e.,
   // largest eigenvalue; it's the last because they come out sorted)
   m_vPlane = mPrincipalAxes.col(mPrincipalAxes.cols() - 1);
   m_fDivider = m_vPlane.dot(vCenterOfMass);

   // sort current set of points on the node into three classes:
   //  - points exactly on the dividing plane -> stay here
   //  - point on the front of the dividing plane -> go to front plane
   //  - points on the back of the dividing plane -> go to back plane
   FPointEntryList
      // points exactly on the dividing plane will stay here.
      // (not likely to happen)
      PointsOnDividingPlane;
   m_pFrontNode = new FPointCloudNode(m_pParentObj);
   m_pBackNode = new FPointCloudNode(m_pParentObj);
   for (size_t iPoint = 0; iPoint != m_Points.size(); ++ iPoint) {
      FPointEntry const
         &p = m_Points[iPoint];
      FScalar
         f = CalcPlaneDist(p.v);
      if (abs(f) < m_pParentObj->m_fEpsilon) {
         PointsOnDividingPlane.push_back(p);
      } else if (f > 0) {
         m_pFrontNode->m_Points.push_back(p);
      } else {
         assert(f < 0);
         m_pBackNode->m_Points.push_back(p);
      }
   }
   m_Points.swap(PointsOnDividingPlane);
}


//FPointArray FPointCloud::GetUniquePoints() const
//{
//   
//}





FPointCloud::FPointEntry const &FPointCloud::Insert(FPoint const &p, ptrdiff_t iPoint)
{
   if (m_pRootNode == 0)
      m_pRootNode = new FPointCloudNode(this);
   return m_pRootNode->Insert(p, iPoint);
}



FMappedPointCloud::FMappedPointCloud(FRotationMatrixList const &TrafoList_)
   : m_TrafoList(TrafoList_)
{}


FMappedPointCloud::FPointEntry const &FMappedPointCloud::Insert(FPoint const &p, ptrdiff_t iPoint)
{
   FRotationMatrixList::const_iterator
      itTrafo;
   FPointEntry const
      *pe = 0;
   ptrdiff_t const
      NOT_SET = -23429037472;
   ptrdiff_t
      iLastPoint = NOT_SET;
   for (itTrafo = m_TrafoList.begin(); itTrafo != m_TrafoList.end(); ++ itTrafo) {
      FPoint
         Rp = (*itTrafo) * p;
      pe = &m_PointCloud.Insert(Rp, iPoint);
      // ^- note: need to return LAST one. point entry objects can be
      // reallocated in subsequent Inert()!
#ifdef _DEBUG
      if (iLastPoint == NOT_SET)
         iLastPoint = pe->iPoint;
      else {
         if (pe->iPoint != iLastPoint)
            throw std::runtime_error("encountered symmetry equivalent points with different point IDs. I assumed this should not happen.");
      }
#else
      if (pe->iPoint != iPoint) {
         assert(iLastPoint == NOT_SET);
         if (iLastPoint != NOT_SET)
            throw std::runtime_error("encountered symmetry point collision AFTER first point...?");
         return *pe;
      }
#endif
   }
   assert(pe != 0);
   return *pe;
}


