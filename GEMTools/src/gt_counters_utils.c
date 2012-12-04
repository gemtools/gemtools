/*
 * PROJECT: GEM-Tools library
 * FILE: gt_counters_utils.c
 * DATE: 20/08/2012
 * DESCRIPTION: // TODO
 */

#include "gt_counters_utils.h"

GT_INLINE int64_t gt_counters_get_uniq_degree(gt_vector* const counters) {
  GT_VECTOR_CHECK(counters);
  register bool found_uniq_strata = false;
  register int64_t uniq_degree = 0;
  GT_VECTOR_ITERATE(counters,counter,counter_pos,uint64_t) {
    if (*counter==0) {
      if (found_uniq_strata) ++uniq_degree;
    } else if (*counter==1) {
      if (found_uniq_strata) return uniq_degree;
      found_uniq_strata=true;
    } else if (*counter>1) {
      if (found_uniq_strata) return uniq_degree;
      return GT_NO_STRATA;
    }
  }
  if (found_uniq_strata) return uniq_degree;
  return GT_NO_STRATA;
}
GT_INLINE bool gt_counters_get_next_matching_strata(
    gt_vector* const counters,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps) {
  register const uint64_t num_counters = gt_vector_get_used(counters);
  register uint64_t i;
  for (i=begin_strata;i<num_counters;++i) {
    register const uint64_t counter_val = *gt_vector_get_elm(counters,i,uint64_t);
    if (counter_val!=0) {
      *next_matching_strata = i;
      *num_maps = counter_val;
      return true;
    }
  }
  return false;
}
GT_INLINE int64_t gt_counters_get_min_matching_strata(gt_vector* const counters) {
  GT_VECTOR_ITERATE(counters,counter,counter_pos,uint64_t) {
    if (*counter!=0) return counter_pos+1;
  }
  return GT_NO_STRATA;
}
