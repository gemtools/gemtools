/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_utils.c
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_alignment_utils.h"

/*
 * Alignment high-level insert (Update global state: counters, ...)
 */
GT_INLINE gt_map* gt_alignment_put_map(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,
    gt_map* const map,const bool replace_duplicated) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  /*
   * Once the @map is given to @gt_alignment_put_map, it belong to the given @alignment.
   * So, in case of duplication, one of the maps will be deleted wrt @replace_duplicated
   */
  // Handle duplicates
  gt_map* found_map;
  uint64_t found_map_pos;
  if (gt_expect_false(gt_alignment_find_map_fx(gt_map_cmp_fx,alignment,map,&found_map_pos,&found_map))) {
    if (gt_expect_true(replace_duplicated && gt_map_get_global_distance(map) < gt_map_get_global_distance(found_map))) {
      /* (gt_map_get_global_bases_aligned(map) <= gt_map_get_global_bases_aligned(found_map) ||
       *  gt_map_get_global_distance(map) < gt_map_get_global_distance(found_map)) */
      // Remove old map
      gt_alignment_dec_counter(alignment,gt_map_get_global_distance(found_map));
      gt_map_delete(found_map);
      // Replace old map
      gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
      gt_alignment_set_map(alignment,map,found_map_pos);
      return map;
    } else {
      gt_map_delete(map);
      return found_map;
    }
  } else {
    // Add new map
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
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
  gt_alignment_put_map(gt_map_cmp_fx,alignment,map,true); /* TODO: Why replace_duplicated?? */
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
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register const bool unique_flag = gt_alignment_get_not_unique_flag(alignment);
  return unique_flag || gt_alignment_is_thresholded_mapped(alignment,UINT64_MAX);
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
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_clear(gt_alignment_get_counters_vector(alignment));
  // Recalculate counters
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  register gt_map* map;
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
  }
}

GT_INLINE int64_t gt_alignment_get_uniq_degree(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_uniq_degree(gt_alignment_get_counters_vector(alignment));
}
GT_INLINE int64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_min_matching_strata(gt_alignment_get_counters_vector(alignment));
}
GT_INLINE bool gt_alignment_get_next_matching_strata(
    gt_alignment* const alignment,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_next_matching_strata(
      gt_alignment_get_counters_vector(alignment),begin_strata,next_matching_strata,num_maps);
}

/*
 * Alignment' Maps set-operators
 */
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  // Initialice new alignment dictionary
  if (alignment_dst->alg_dictionary == NULL) {
    alignment_dst->alg_dictionary = gt_alignment_dictionary_new(alignment_dst);
  }
  gt_ihash_element* ihash_element_b = NULL;
  gt_ihash_element* ihash_element_e = NULL;
  uint64_t found_vector_position;
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    register const uint64_t vector_position = gt_vector_get_used(alignment_dst->maps);
    register const uint64_t begin_position = gt_map_get_position(map_src)-gt_map_get_left_trim_length(map_src);
    register const uint64_t end_position = gt_map_get_position(map_src)+gt_map_get_length(map_src);
    // Try add
    if (gt_alignment_dictionary_try_add(alignment_dst->alg_dictionary,map_src,
        begin_position,end_position,vector_position,&found_vector_position,ihash_element_b,ihash_element_e)) {
      register gt_map* const map_src_cp = gt_map_copy(map_src);
      gt_alignment_inc_counter(alignment_dst,gt_map_get_global_distance(map_src_cp));
      gt_alignment_add_map(alignment_dst,map_src_cp);
    } else {
      // Solve conflict
      register gt_map* map_found = gt_alignment_get_map(alignment_dst,found_vector_position);
      if (gt_map_get_global_distance(map_src) < gt_map_get_global_distance(map_found)) {
        register gt_map* const map_src_cp = gt_map_copy(map_src);
        gt_alignment_dictionary_replace(alignment_dst->alg_dictionary,
            map_src_cp,map_found,begin_position,end_position,found_vector_position,ihash_element_b,ihash_element_e);
      }
    }
  }
  gt_alignment_set_mcs(alignment_dst,GT_MIN(gt_alignment_get_mcs(alignment_dst),gt_alignment_get_mcs(alignment_src)));
}
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    register gt_map* const map_src_cp = gt_map_copy(map_src);
    gt_alignment_put_map(gt_map_cmp_fx,alignment_dst,map_src_cp,true);
  }
  gt_alignment_set_mcs(alignment_dst,GT_MIN(gt_alignment_get_mcs(alignment_dst),gt_alignment_get_mcs(alignment_src)));
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

GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_v(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_ZERO_CHECK(num_src_alignments);
  // Create new alignment
  register gt_alignment* const alignment_union = gt_alignment_copy(alignment_src,false);
  gt_alignment_merge_alignment_maps(alignment_union,alignment_src);
  // Merge alignment sources into alignment_union
  register uint64_t num_alg_merged = 1;
  while (num_alg_merged < num_src_alignments) {
    register gt_alignment* alignment_target = va_arg(v_args,gt_alignment*);
    GT_ALIGNMENT_CHECK(alignment_target);
    gt_alignment_merge_alignment_maps(alignment_union,alignment_target);
    ++num_alg_merged;
  }
  return alignment_union;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_va(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_ZERO_CHECK(num_src_alignments);
  GT_ALIGNMENT_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  register gt_alignment* const alignment_dst =
      gt_alignment_union_alignment_maps_v(num_src_alignments,alignment_src,v_args);
  va_end(v_args);
  return alignment_dst;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ZERO_CHECK(num_src_alignments);
  // Create new alignment
  register gt_alignment* const alignment_union = gt_alignment_copy(alignment_src,false);
  gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_union,alignment_src);
  // Merge alignment sources into alignment_union
  register uint64_t num_alg_merged = 1;
  while (num_alg_merged < num_src_alignments) {
    register gt_alignment* alignment_target = va_arg(v_args,gt_alignment*);
    GT_ALIGNMENT_CHECK(alignment_target);
    gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_union,alignment_target);
    ++num_alg_merged;
  }
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

GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend) {
  GT_ALIGNMENT_CHECK(alignment_minuend);
  GT_ALIGNMENT_CHECK(alignment_subtrahend);
  // Create new alignment
  register gt_alignment* const alignment_difference = gt_alignment_copy(alignment_minuend,false);
  // Copy not common maps
  GT_ALIGNMENT_ITERATE(alignment_minuend,map_minuend) {
    if (!gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_subtrahend,map_minuend)) { // TODO Improvement: Scheduled for v2.0
      register gt_map* const map_cp = gt_map_copy(map_minuend);
      gt_alignment_put_map(gt_map_cmp_fx,alignment_difference,map_cp,false);
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
  register gt_alignment* const alignment_intersection = gt_alignment_copy(alignment_src_A,false);
  // Copy common maps
  GT_ALIGNMENT_ITERATE(alignment_src_A,map_A) {
    if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_src_B,map_A)) { // TODO Improvement: Scheduled for v2.0
      register gt_map* const map_cp = gt_map_copy(map_A);
      gt_alignment_put_map(gt_map_cmp_fx,alignment_intersection,map_cp,false);
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

/*
 * Alignment realignment
 */
GT_INLINE void gt_alignment_realign_hamming(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_hamming_sa(map,alignment->read,sequence_archive);
  }
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_realign_levenshtein(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_levenshtein_sa(map,alignment->read,sequence_archive);
  }
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_realign_weighted(
    gt_alignment* const alignment,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_NULL_CHECK(gt_weigh_fx);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_weighted_sa(map,alignment->read,sequence_archive,gt_weigh_fx);
  }
  gt_alignment_recalculate_counters(alignment);
}
