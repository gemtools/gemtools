/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_utils.c
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_alignment_utils.h"

/*
 * Alignment high-level insert (Update global state: counters, ...)
 */
GT_INLINE gt_map* gt_alignment_put_map(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map,const bool replace_dup) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  // Handle duplicates
  gt_map* found_map;
  uint64_t found_map_pos;
  if (gt_expect_false(gt_alignment_find_map_fx(gt_map_cmp_fx,alignment,map,&found_map_pos,&found_map))) {
    if (gt_expect_true(replace_dup)) {
      // Remove old map
      gt_alignment_dec_counter(alignment,gt_map_get_global_distance(found_map)+1);
      gt_map_delete(found_map); // TODO: v2. Removal of a mmap member in duplicates resolution of conflict
      // Replace old map
      gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map)+1);
      gt_alignment_set_map(alignment,map,found_map_pos);
      return map;
    } else {
      return found_map;
    }
  } else {
    // Add new map
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map)+1);
    gt_alignment_add_map(alignment,map);
    return map;
  }
}

GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_alignment_insert_map_fx(gt_map_cmp,alignment,map);
}
GT_INLINE void gt_alignment_insert_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  // Prevent duplicates
  gt_alignment_put_map(gt_map_cmp_fx,alignment,map,true);
}
GT_INLINE void gt_alignment_insert_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map_it,map_pos,gt_map*) {
    gt_alignment_insert_map(alignment,*map_it);
  }
}
GT_INLINE void gt_alignment_insert_map_fx_gt_vector(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map_it,map_pos,gt_map*) {
    gt_alignment_insert_map_fx(gt_map_cmp_fx,alignment,*map_it);
  }
}

GT_INLINE void gt_alignment_remove_map(gt_alignment* const alignment,gt_map* const map) {
  // TODO: Scheduled for v2.0
}
GT_INLINE void gt_alignment_remove_map_fx(
      int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  // TODO: Scheduled for v2.0
}

GT_INLINE bool gt_alignment_find_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map,
    uint64_t* const found_map_pos,gt_map** const found_map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  GT_NULL_CHECK(found_map_pos); GT_NULL_CHECK(found_map);
  // Search for the map
  register uint64_t pos = 0;
  GT_ALIGNMENT_ITERATE(alignment,map_it) {
   if (gt_map_cmp_fx(map_it,map)==0) {
     *found_map_pos = pos;
     *found_map = map_it;
     return true;
   }
   ++pos;
  }
  return false;
}
GT_INLINE bool gt_alignment_is_map_contained(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_map* found_map;
  uint64_t found_map_pos;
  return gt_alignment_find_map_fx(gt_map_cmp,alignment,map,&found_map_pos,&found_map);
}
GT_INLINE bool gt_alignment_is_map_contained_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_map* found_map;
  uint64_t found_map_pos;
  return gt_alignment_find_map_fx(gt_map_cmp_fx,alignment,map,&found_map_pos,&found_map);
}

/*
 * Alignment' Counters operators
 */
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_clean(gt_alignment_get_counters_vector(alignment));
  // Recalculate counters
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  register gt_map* map;
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map)+1);
  }
}
GT_INLINE int64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (*counter!=0) return counter_pos+1;
  }
  return GT_NO_STRATA;
}
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const int64_t max_allowed_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (counter_pos>=max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register const bool unique_flag = gt_alignment_get_not_unique_flag(alignment);
  return unique_flag || gt_alignment_is_thresholded_mapped(alignment,UINT64_MAX);
}

/*
 * Alignment' Maps set-operators
 */
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    register gt_map* const map_src_cp = gt_map_copy(map_src);
    gt_alignment_insert_map(alignment_dst,map_src_cp);
  }
}
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    register gt_map* const map_src_cp = gt_map_copy(map_src);
    gt_alignment_insert_map_fx(gt_map_cmp_fx,alignment_dst,map_src_cp);
  }
}

GT_INLINE void gt_alignment_remove_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    gt_alignment_remove_map(alignment_dst,map_src); // TODO: Scheduled for v2.0
  }
}
GT_INLINE void gt_alignment_remove_alignment_maps_fx(
    int64_t (*gt_map_cmp)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_NULL_CHECK(gt_map_cmp);
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    gt_alignment_remove_map_fx(gt_map_cmp,alignment_dst,map_src); // TODO: Scheduled for v2.0
  }
}

GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ZERO_CHECK(num_src_alignments);
  // Create new alignment
  register gt_alignment* const alignment_union = gt_alignment_new();
  gt_alignment_handler_copy(alignment_union,alignment_src);
  // Merge alignment sources into alignment_union
  register uint64_t num_alg_merged = 0;
  register gt_alignment* alignment_target = alignment_src;
  do {
    GT_ALIGNMENT_CHECK(alignment_target);
    gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_union,alignment_target);
    alignment_target = (++num_alg_merged < num_src_alignments) ? va_arg(v_args,gt_alignment*) : NULL;
  } while (gt_expect_true(alignment_target!=NULL));
  return alignment_union;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ZERO_CHECK(num_src_alignments);
  GT_ALIGNMENT_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  register gt_alignment* const alignment_dst =
      gt_alignment_union_alignment_maps_fx_v(gt_map_cmp_fx,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
  return alignment_dst;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_va(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_ZERO_CHECK(num_src_alignments);
  GT_ALIGNMENT_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  register gt_alignment* const alignment_dst =
      gt_alignment_union_alignment_maps_fx_v(gt_map_cmp,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
  return alignment_dst;
}

GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend) {
  GT_ALIGNMENT_CHECK(alignment_minuend);
  GT_ALIGNMENT_CHECK(alignment_subtrahend);
  // Create new alignment
  register gt_alignment* const alignment_difference = gt_alignment_new();
  gt_alignment_handler_copy(alignment_difference,alignment_minuend);
  // Copy not common maps
  GT_ALIGNMENT_ITERATE(alignment_minuend,map_minuend) {
    if (!gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_subtrahend,map_minuend)) { // TODO Improvement: Scheduled for v2.0
      register gt_map* const map_cp = gt_map_copy(map_minuend);
      gt_alignment_insert_map_fx(gt_map_cmp_fx,alignment_difference,map_cp);
    }
  }
  return alignment_difference;
}
GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps(
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend) {
  GT_ALIGNMENT_CHECK(alignment_minuend);
  GT_ALIGNMENT_CHECK(alignment_subtrahend);
  return gt_alignment_subtract_alignment_maps_fx(gt_map_cmp,alignment_minuend,alignment_subtrahend);
}

GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B) {
  GT_ALIGNMENT_CHECK(alignment_src_A);
  GT_ALIGNMENT_CHECK(alignment_src_B);
  // Create new alignment
  register gt_alignment* const alignment_intersection = gt_alignment_new();
  gt_alignment_handler_copy(alignment_intersection,alignment_src_A);
  // Copy common maps
  GT_ALIGNMENT_ITERATE(alignment_src_A,map_A) {
    if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_src_B,map_A)) { // TODO Improvement: Scheduled for v2.0
      register gt_map* const map_cp = gt_map_copy(map_A);
      gt_alignment_insert_map_fx(gt_map_cmp_fx,alignment_intersection,map_cp);
    }
  }
  return alignment_intersection;
}
GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps(
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B) {
  GT_ALIGNMENT_CHECK(alignment_src_A);
  GT_ALIGNMENT_CHECK(alignment_src_B);
  return gt_alignment_intersect_alignment_maps_fx(gt_map_cmp,alignment_src_A,alignment_src_B);
}
