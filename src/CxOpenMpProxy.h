// include OpenMP header if available or define inline dummy functions
// for OMP primitives if not.
#ifndef OPENMP_PROXY_H
#define OPENMP_PROXY_H

#ifdef _OPENMP
   #include <omp.h>
#else
   // Current thread id, within the current workteam (that is, within the innermost
   // currently running parallel region). Supposed to return a value ∈ {0, 1, …,
   // (omp_get_num_threads() − 1)}
   inline int omp_get_thread_num() { return 0; }

   // Set number of threads to use in workteams *BELOW THE CURRENT CONTEXT ONLY*.
   // Can be used to control nested parallelism, by calling it from inside running parallel sections.
   // May need to be used together with calls like
   // `omp_set_nested(1);` `omp_set_max_active_levels(2);` `omp_set_dynamic(0);`
   // to be effective.
   //
   // Warning:
   // - not a good idea if you want the program to run on VC...
   // - If I am not very mistaken, this is actually *NOT* the counterpart of
   //   omp_get_num_threads():
   //   + omp_set_num_threads() sets how many threads to open for workteams in
   //     the NEXT parallel region below the current one, while
   //   + omp_get_num_threads() returns how many threads are active in the
   //     current parallel region.
   //
   // see also: https://stackoverflow.com/questions/59434959/openmp-omp-get-num-threads-v-s-omp-get-max-threads
   inline void omp_set_num_threads(int) {}

   // Number of threads which are active in the "current work team" of the current
   // parallel region (will return 1 if in sequential part of the program).
   // If I understand this correctly, this is how many threads the current
   // parallel region is split into, NOT taking into account any outer parallel
   // regions in case of nesting. E.g., if you have something like this:
   // ```
   // omp_set_nested(1);
   // omp_set_max_active_levels(2);
   // omp_set_num_threads(4);
   // #pragma omp parallel for
   // for ( ... ) {
   //    omp_set_num_threads(3);
   //    #pragma omp parallel for
   //    for ( ...) {
   //       // <-- here
   //    }
   // }
   // ```
   // My understanding is that you'd only get `3` as answer for omp_get_num_threads() if called
   // in the inner region.
   inline int omp_get_num_threads() { return 1; }

   // Max number of threads supposed to be running simultaneously (unless
   // entering a parallel section specifying overrides via num_threads() pragma)
   inline int omp_get_max_threads() { return 1; }

   // total number of "virtual" CPU cores (includes hyperthreading cores)
   inline int omp_get_num_procs() { return 1; }

   struct omp_lock_t {};
   inline void omp_destroy_lock(omp_lock_t *){}
   inline void omp_init_lock(omp_lock_t *){}
   inline void omp_set_lock(omp_lock_t *){}
   inline void omp_set_nested(int){}
#endif

#endif // OPENMP_PROXY_H
