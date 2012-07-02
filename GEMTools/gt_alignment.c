/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_alignment.h"

#define GT_ALIGNMENT_NUM_INITIAL_MAPS 20
#define GT_ALIGNMENT_INITIAL_LENGTH_COUNTERS 5

/*
 * Setup
 */
GT_INLINE gt_alignment* gt_alignment_new() {
  gt_alignment* alignment = malloc(sizeof(gt_alignment));
  gt_cond_fatal_error(!alignment,MEM_HANDLER);
  alignment->tag = NULL;
  alignment->tag_length = 0;
  alignment->read = NULL;
  alignment->read_length = 0;
  alignment->qualities = NULL;
  alignment->counters = gt_vector_new(GT_ALIGNMENT_INITIAL_LENGTH_COUNTERS,sizeof(uint64_t));
  alignment->max_complete_strata = 0;
  alignment->maps = gt_vector_new(GT_ALIGNMENT_NUM_INITIAL_MAPS,sizeof(gt_map));
  alignment->maps_txt = NULL;
  return alignment;
}
GT_INLINE void gt_alignment_clear(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  alignment->tag = NULL;
  alignment->tag_length = 0;
  alignment->read = NULL;
  alignment->read_length = 0;
  alignment->qualities = NULL;
  gt_vector_clean(alignment->counters);
  alignment->max_complete_strata = 0;
  gt_alignment_clear_maps(alignment);
  alignment->maps_txt = NULL;
}
GT_INLINE void gt_alignment_delete(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_delete(alignment->counters);
  gt_alignment_clear_maps(alignment);
  gt_vector_delete(alignment->maps);
  free(alignment);
}

/*
 * Accessors
 */
GT_INLINE char* gt_alignment_get_tag(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->tag;
}
GT_INLINE void gt_alignment_set_tag(gt_alignment* const alignment,char* const tag) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(tag);
  alignment->tag = tag;
  alignment->tag_length = strlen(tag);
}
GT_INLINE char* gt_alignment_get_read(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->read;
}
GT_INLINE void gt_alignment_set_read(gt_alignment* const alignment,char* const read) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(read);
  alignment->read =  read;
  alignment->read_length = strlen(read);
  gt_fatal_check(alignment->qualities!=NULL &&
      strlen(alignment->qualities)!=alignment->read_length,ALIGN_READ_QUAL_LENGTH);
}
GT_INLINE char* gt_alignment_get_qualities(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->qualities;
}
GT_INLINE void gt_alignment_set_qualities(gt_alignment* alignment,char* qualities) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(qualities);
  gt_fatal_check(alignment->read!=NULL &&
      strlen(alignment->qualities)!=alignment->read_length,ALIGN_READ_QUAL_LENGTH);
  alignment->qualities = qualities;
}

GT_INLINE gt_vector* gt_alignment_get_counters_vector(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->counters;
}
GT_INLINE void gt_alignment_set_counters_vector(gt_alignment* const alignment,gt_vector* const counters) {
  GT_ALIGNMENT_CHECK(alignment);
  alignment->counters = counters;
}
GT_INLINE uint64_t gt_alignment_get_num_counters(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_vector_get_used(alignment->counters);
}
GT_INLINE uint64_t gt_alignment_get_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  return *gt_vector_get_elm(alignment->counters,stratum-1,uint64_t);
}
GT_INLINE void gt_alignment_set_counter(gt_alignment* const alignment,const uint64_t stratum,const uint64_t value) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(alignment->counters,stratum,true);
  if (gt_vector_get_used(alignment->counters)<stratum) {
    gt_vector_set_used(alignment->counters,stratum);
  }
  *gt_vector_get_elm(alignment->counters,stratum-1,uint64_t) = value;
}
GT_INLINE void gt_alignment_dec_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(alignment->counters,stratum,true);
  if (gt_vector_get_used(alignment->counters)<stratum) {
    gt_vector_set_used(alignment->counters,stratum);
  }
  --(*gt_vector_get_elm(alignment->counters,stratum-1,uint64_t));
}
GT_INLINE void gt_alignment_inc_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(alignment->counters,stratum,true);
  if (gt_vector_get_used(alignment->counters)<stratum) {
    gt_vector_set_used(alignment->counters,stratum);
  }
  ++(*gt_vector_get_elm(alignment->counters,stratum-1,uint64_t));
}
GT_INLINE uint64_t gt_alignment_get_mcs(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->max_complete_strata;
}
GT_INLINE void gt_alignment_set_mcs(gt_alignment* const alignment,const uint64_t max_complete_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  alignment->max_complete_strata = max_complete_strata;
}

/*
 * Maps Handlers
 */
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(map);
  gt_vector_insert(alignment->maps,map,gt_map*);
}
GT_INLINE void gt_alignment_clear_maps(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  register gt_map* map;
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    gt_map_delete(map);
  }
  gt_vector_clean(alignment->maps);
}
GT_INLINE gt_map* gt_alignment_get_map(gt_alignment* const alignment,const uint64_t position) {
  GT_ALIGNMENT_EDITABLE_CHECK(alignment);
  return *gt_vector_get_elm(alignment->maps,position,gt_map*);
}
GT_INLINE uint64_t gt_alignment_get_num_maps(gt_alignment* const alignment) {
  GT_ALIGNMENT_EDITABLE_CHECK(alignment);
  return gt_vector_get_used(alignment->maps);
}

/*
 * Higher-level Methods
 */
GT_INLINE void gt_alignment_is_contained(gt_alignment* const alignment,gt_map* const map,gt_map** const dup_map) {
  gt_alignment_map_iterator alignment_map_iterator;
  gt_alignment_new_map_iterator(alignment,&alignment_map_iterator);
  register gt_map* aux_map;
  while ((aux_map=gt_alignment_next_map(&alignment_map_iterator))) {
    // TODO
  }
}
GT_INLINE void gt_alignment_insert_map__check_dup(gt_alignment* const alignment,gt_map** const map) {
//  // TODO
//  GT_ALIGNMENT_CHECK(alignment);
//  GT_NULL_CHECK(*map);
//  gt_map* dup_map; // TODO: Solve Duplication
//  if (gt_alignment_is_contained(alignment,*map,&dup_map)) {
//    *map = dup_map;
//  } else {
//    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
//    gt_alignment_add_map(alignment,map);
//  }
}
GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(map);
  gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map)+1);
  gt_alignment_add_map(alignment,map);
}
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment) {
  // TODO
}
GT_INLINE uint64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (*counter!=0) return counter_pos+1;
  }
  return UINT64_MAX;
}
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const uint64_t max_allowed_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if ((counter_pos+1)>=max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_alignment_is_thresholded_mapped(alignment,UINT64_MAX);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment* alignment_cpy = malloc(sizeof(gt_alignment));
  gt_cond_fatal_error(!alignment_cpy,MEM_HANDLER);
  alignment_cpy->tag = alignment->tag;
  alignment_cpy->tag_length = alignment->tag_length;
  alignment_cpy->read = alignment->read;
  alignment_cpy->read_length = alignment->read_length;
  alignment_cpy->qualities = alignment->qualities;
  alignment_cpy->counters = gt_vector_new(gt_vector_get_used(alignment->counters),sizeof(uint64_t));
  gt_vector_copy(alignment_cpy->counters,alignment->counters);
  alignment_cpy->max_complete_strata = alignment->max_complete_strata;
  alignment_cpy->maps = gt_vector_new(gt_vector_get_used(alignment->maps),sizeof(gt_map*));
  gt_vector_copy(alignment_cpy->maps,alignment->maps);
  alignment_cpy->maps_txt = alignment->maps_txt;
  return alignment_cpy;
}
GT_INLINE gt_alignment* gt_alignment_deep_copy(gt_alignment* const alignment) {
  // TODO
  return NULL;
}

/*
 * Alignment's Maps iterator
 */
GT_INLINE void gt_alignment_new_map_iterator(gt_alignment* const alignment,gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_EDITABLE_CHECK(alignment);
  alignment_map_iterator->alignment = alignment;
  alignment_map_iterator->total_pos = gt_vector_get_used(alignment->maps);
  alignment_map_iterator->next_map = gt_vector_get_mem(alignment->maps,gt_map*);
  alignment_map_iterator->next_pos = 0;
}
GT_INLINE gt_map* gt_alignment_next_map(gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_EDITABLE_CHECK(alignment_map_iterator->alignment);
  if (gt_expect_true(alignment_map_iterator->next_pos<alignment_map_iterator->total_pos)) {
    register gt_map* const map = *(alignment_map_iterator->next_map);
    ++alignment_map_iterator->next_map;
    ++alignment_map_iterator->next_pos;
    return map;
  } else {
    return NULL;
  }
}
GT_INLINE gt_map* gt_alignment_dinamic_next_map(gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_EDITABLE_CHECK(alignment_map_iterator->alignment);
  register const gt_alignment* const alignment = alignment_map_iterator->alignment;
  if (gt_expect_true(alignment_map_iterator->next_pos<gt_vector_get_used(alignment->maps))) {
    register gt_map* const map = *gt_vector_get_elm(alignment->maps,alignment_map_iterator->next_pos,gt_map*);
    ++alignment_map_iterator->next_pos;
    return map;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_alignment_next_map_pos(gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_EDITABLE_CHECK(alignment_map_iterator->alignment);
  return alignment_map_iterator->next_pos;
}
