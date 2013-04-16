/*
 * PROJECT: GEM-Tools library
 * FILE: gt_profiler.c
 * DATE: 16/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_profiler.h"

// TIME counters
struct timeval* gt_prof_begin_timer = NULL;
struct timeval* gt_prof_end_timer = NULL;
double* gt_prof_timer = NULL;
// Rank/LF counters
uint64_t gt_prof_rank_count = 0;
uint64_t *gt_prof_begin_rank_counter = NULL;
uint64_t *gt_prof_rank_counter = NULL;
// General counters & functions
uint64_t *gt_prof_counter = NULL;

#ifndef GT_NOPROFILE

/*
 * Setup
 */
void gt_profiler_init() {
  // Allocate all profiler counters
  gt_prof_begin_timer = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,struct timeval,true);
  gt_prof_end_timer = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,struct timeval,true);
  gt_prof_timer = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,double,true);
  gt_prof_rank_count = 0;
  gt_prof_begin_rank_counter = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,uint64_t,true);
  gt_prof_rank_counter = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,uint64_t,true);
  gt_prof_counter = gt_calloc(GT_PROFILE_COUNTERS_ALLOCATED,uint64_t,true);
}
void gt_profiler_exit() {
  // Release all profiler counters
  gt_free(gt_prof_begin_timer);
  gt_free(gt_prof_end_timer);
  gt_free(gt_prof_timer);
  gt_free(gt_prof_begin_rank_counter);
  gt_free(gt_prof_rank_counter);
  gt_free(gt_prof_counter);
}

#endif
