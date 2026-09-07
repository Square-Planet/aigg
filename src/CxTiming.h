// This file is released under the GNU General Public License ("GPL", version 2)
// as part of the CT8K program. Copying, modification, creating derivative works
// and redistribution of the file is allowed, but _only_ subject to the terms
// of that GPL. You should have received a version of this license along
// with this source code. The program comes "as is", without any kind of
// warranty.
//
// Authors/Copyright holders:
//  - Gerald Knizia, 2006 (tag: cgk, contact: cgk.d@gmx.net)


#ifndef _TIMING_H
#define _TIMING_H

#include <string>
#include <ostream>
#include <stddef.h> // for size_t/ptrdiff_t

namespace ct {

double Second();

struct FTimer
{
   FTimer() : Start(Second()) {};
   operator double () { return Second() - Start; };
   void Reset() { Start = Second(); }
private:
   double Start;
};


enum {
   TIMER_MaxEntries = 0x1000
};

struct FTimerSet
{
   struct FEntry {
      double fTime;
      size_t iInvk;
      char const *pDesc;
   };

   // sets all timers to 0.
   FTimerSet();

   inline void Resume(size_t iClockId, char const *pDesc);
   inline void Pause(size_t iClockId);
   inline void Pause(size_t iClockId, char const *pDesc) { Pause(iClockId); (void)pDesc; } ;

   // increase/decrease the current level of detail timing.
   inline void Enter() { m_CurrentDetailLevel += 1; };
   inline void Leave() { m_CurrentDetailLevel -= 1; };

   inline void Enter(size_t iClockId, char const *pDesc) { Resume(iClockId, pDesc); Enter(); }
   inline void Leave(size_t iClockId, char const *pDesc) { Leave(); Pause(iClockId, pDesc); };
   inline void Leave(size_t iClockId) { Leave(); Pause(iClockId); };

   void Reset();
   void SetLevel(int iNewLevel) { m_TargetDetailLevel = iNewLevel; };
   void PrintReport(std::ostream &out, double ThrPrint = 0.01);
protected:
   FEntry
      Entries[TIMER_MaxEntries];
   ptrdiff_t
      m_CurrentDetailLevel,
      m_TargetDetailLevel; // currently enabled level of detail timing. May be negative to disable.
   size_t nMaxEntries() const { return TIMER_MaxEntries; }
};

void FTimerSet::Resume(size_t iClockId, char const *pDesc)
{
//    _CxAssert(iClockId < TIMER_MaxEntries);
   if (m_CurrentDetailLevel <= m_TargetDetailLevel) {
      FEntry &e = Entries[iClockId];
      e.iInvk += 1;
      e.fTime -= Second();
      e.pDesc = pDesc;
   }
}

void FTimerSet::Pause(size_t iClockId)
{
//    _CxAssert(iClockId < TIMER_MaxEntries);
   if (m_CurrentDetailLevel <= m_TargetDetailLevel) {
      FEntry &e = Entries[iClockId];
      e.fTime += Second();
   }
}

struct FTimeSection
{
   inline FTimeSection(FTimerSet *pTimers_, size_t iClockId_, char const *pDesc_)
      : pTimers(pTimers_), iClockId(iClockId_), pDesc(pDesc_)
//    { if(pTimers) pTimers->Enter(iClockId, pDesc_); }
   { pTimers->Enter(iClockId, pDesc_); }

//    ~FTimeSection() { if (pTimers) pTimers->Leave(iClockId); }
   ~FTimeSection() { pTimers->Leave(iClockId); }
public:
   FTimerSet *pTimers;
   size_t iClockId;
   char const *pDesc;
};

#define TIME_SECTION(a,b,c) FTimeSection IR_UNIQUE_NAME(a,b,c)

// Get time stamp counter (TSC) from current CPU core.
// WARNING: This is fragile and unportable! Use only for developer-side
// performance analysis and debugging!
// Notes:
//   - This counts CPU cycles on the current core
//   - Function should be fast (my guess: ~10-30 cycles), so it can be
//     used for timing small parts of code.
//   - TSCs on different CPU cores may differ. If the calling thread gets
//     re-distributed across cores between two successive calls, the TSC
//     may jump.
//   - In the presence of power management, espect varying cpu rates and
//     large jumps in the counter (suspended mode etc)
unsigned long long dbgGetTimeStampCounter();

} // namespace ct

#endif // _TIMING_H
