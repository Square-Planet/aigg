#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include "Aigg.h" // for scalars and points
#include <stddef.h> // for size_t/ptrdiff_t
#include "CxIntrusivePtr.h"

struct FPointCloudNode;

// stores a, possibly large, set of points. Intended to look up the existence of
// previous points. Data is stored in the form of a semi-balanced BSP tree at
// this moment.
struct FPointCloud
{
   struct FPointEntry {
      FPoint
         v; // position
      ptrdiff_t
         iPoint, // some externally supplied index.
         nRefs; // number of times a corresponding point was inserted.
      FPointEntry();
      FPointEntry(FPoint const &v_, ptrdiff_t i_, ptrdiff_t nRefs_);
   };

   // returns ref to point entry -- either the original, or the one which was there before.
   FPointEntry const &Insert(FPoint const &p, ptrdiff_t iPoint);

   FPointCloud();
   ~FPointCloud();
protected:
   FPointCloudNode
      *m_pRootNode;
   FScalar
      m_fEpsilonSq,
      m_fEpsilon;
   size_t
      m_nMaxPointsPerNode;
   void SetEpsilon(FScalar fEpsilon_);
   bool IsApproxEqual(FPoint const &a, FPoint const &b) const;
//   FPointArray GetUniquePoints() const;
   friend class FPointCloudNode;
private:
   FPointCloud(FPointCloud const&); // not implemented
   void operator = (FPointCloud const&); // not implemented
};


struct FMappedPointCloud
{
   typedef FPointCloud::FPointEntry
      FPointEntry;

   explicit FMappedPointCloud(FRotationMatrixList const &TrafoList);


   FPointEntry const &Insert(FPoint const &p, ptrdiff_t iPoint);
protected:
   FRotationMatrixList const
      &m_TrafoList;
   FPointCloud
      m_PointCloud;
};


#endif // POINT_CLOUD_H
