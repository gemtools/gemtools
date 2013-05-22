/*
 * PROJECT: GEM-Tools library
 * FILE: gt_profiler.c
 * DATE: 16/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_PROFILER_H_
#define GT_PROFILER_H_

#include "gt_commons.h"
#include "gt_mm.h"

#define GT_PROFILE_COUNTERS_ALLOCATED 1000

// Operative macros
#define GT_TIME_DIFF(start,end) ((end.tv_sec + end.tv_usec/1E6) - (start.tv_sec + start.tv_usec/1E6))

/*
 * Handy macros to display statistics
 */
#define GT_GET_TIME_PERCENTAGE(timer,tota_timer) GT_GET_PERCENTAGE(GT_GET_TIMER(timer),GT_GET_TIMER(tota_timer))
#define GT_GET_RANK_PERCENTAGE(rank_counter,total_rank_counter) GT_GET_PERCENTAGE(GT_GET_RANK_COUNTER(rank_counter),GT_GET_RANK_COUNTER(total_rank_counter))
#define GT_GET_COUNT_PERCENTAGE(counter,total_counter) GT_GET_PERCENTAGE(GT_GET_COUNTER(counter),GT_GET_COUNTER(total_counter))
#define GT_GET_COUNT_DIV(counter1,counter2) ((float)GT_GET_COUNTER(counter1)/(float)GT_GET_COUNTER(counter2))
/* Time per call in milliseconds*/
#define GT_GET_TIME_PER_CALL(timer,counter) \
  ((GT_GET_COUNTER(counter)!=0)?1000.0*(double)GT_GET_TIMER(timer))/((double)GT_GET_COUNTER(counter):0.0)

#ifndef GT_NOPROFILE

  /*
   * Profiling & Statistics counters
   */

  // AC:: Aggregated Counters (TC,RC,EC) [0,400)
  // #define GT_AC_

  // TC:: Time Counters [400,...]
  // #define GT_TC_

  // RC:: Ranks Counters [400,...]
  // #define GT_RC_

  // EC:: Event Counters [400,...]
  // #define GT_EC_

  // TIME counters
  extern struct timeval *gt_prof_begin_timer;
  extern struct timeval *gt_prof_end_timer;
  extern double *gt_prof_timer;
  // TIME functions
  #define GT_START_TIMER(num) gettimeofday((gt_prof_begin_timer+num),NULL)
  #define GT_STOP_TIMER(num) \
    gettimeofday((gt_prof_end_timer+num),NULL); \
    gt_prof_timer[num]+=GT_TIME_DIFF(gt_prof_begin_timer[num],gt_prof_end_timer[num])
  #define GT_GET_TIMER(num) gt_prof_timer[num]
  #define GT_RESET_TIME(num) gt_prof_timer[num]=0.0

  // Rank/LF counters
  extern uint64_t gt_prof_rank_count;
  extern uint64_t *gt_prof_begin_rank_counter;
  extern uint64_t *gt_prof_rank_counter;
  // Rank/LF functions
  #define GT_START_RANK_COUNTER(num) gt_prof_begin_rank_counter[num]=gt_prof_rank_count
  #define GT_STOP_RANK_COUNTER(num) gt_prof_rank_counter[num]+=gt_prof_rank_count-gt_prof_begin_rank_counter[num]
  #define GT_GET_RANK_COUNTER(num) gt_prof_rank_counter[num]

  // General counters & functions
  extern uint64_t *gt_prof_counter;
  #define GT_GET_COUNTER(num) gt_prof_counter[num]
  #define GT_SET_COUNTER(num,val) gt_prof_counter[num]=val
  #define GT_INC_COUNTER(num) gt_prof_counter[num]++
  #define GT_DEC_COUNTER(num) gt_prof_counter[num]--
  #define GT_ADD_COUNTER(num,sum) gt_prof_counter[num]+=sum
  #define GT_RESET_COUNTER(num) gt_prof_counter[num]=0

  // AC operators
  #define GT_START_AC_PROF(AC_COUNTER) \
    GT_INC_COUNTER(AC_COUNTER); \
    GT_START_RANK_COUNTER(AC_COUNTER); \
    GT_START_TIMER(AC_COUNTER)
  #define GT_STOP_AC_PROF(AC_COUNTER) \
    GT_STOP_TIMER(AC_COUNTER); \
    GT_STOP_RANK_COUNTER(AC_COUNTER)

  /*
   * Setup
   */
  void gt_profiler_init();
  void gt_profiler_exit();
#else
  #define GT_GET_TIMER(num) 0.0
  #define GT_RESET_TIME(num)
  // Rank/LF functions
  #define GT_START_RANK_COUNTER(num)
  #define GT_STOP_RANK_COUNTER(num)
  #define GT_GET_RANK_COUNTER(num) 0
  // General counters & functions
  #define GT_GET_COUNTER(num) 0.0
  #define GT_SET_COUNTER(num,val)
  #define GT_INC_COUNTER(num)
  #define GT_DEC_COUNTER(num)
  #define GT_ADD_COUNTER(num,sum)
  #define GT_RESET_COUNTER(num)
  // AC operators
  #define GT_START_AC_PROF(AC_COUNTER)
  #define GT_STOP_AC_PROF(AC_COUNTER)
  // Setup
  #define gt_profiler_init()
  #define gt_profiler_exit()
#endif



#endif /* GT_PROFILER_H_ */
