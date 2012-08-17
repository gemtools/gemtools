/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_alignment.h"

#define GT_ALIGNMENT_NUM_INITIAL_MAPS 20
#define GT_ALIGNMENT_INITIAL_LENGTH_COUNTERS 5

typedef struct {
  uint64_t dummy;
} gt_alignment_dictionary;

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
  alignment->max_complete_strata = UINT64_MAX;
  alignment->maps = gt_vector_new(GT_ALIGNMENT_NUM_INITIAL_MAPS,sizeof(gt_map));
  alignment->maps_txt = NULL;
  alignment->maps_dictionary = gt_shash_new();
  return alignment;
}
GT_INLINE void gt_alignment_clear(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_clear_maps(alignment);
  gt_shash_clean(alignment->maps_dictionary,true,true);
  gt_alignment_clear_handler(alignment);
}
GT_INLINE void gt_alignment_clear_handler(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  alignment->tag = NULL;
  alignment->tag_length = 0;
  alignment->read = NULL;
  alignment->read_length = 0;
  alignment->qualities = NULL;
  gt_vector_clean(alignment->counters);
  alignment->max_complete_strata = UINT64_MAX;
  gt_vector_clean(alignment->maps);
  alignment->maps_txt = NULL;
}
GT_INLINE void gt_alignment_delete(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_clear_maps(alignment);
  gt_shash_delete(alignment->maps_dictionary,true,true);
  gt_alignment_delete_handler(alignment);
}
GT_INLINE void gt_alignment_delete_handler(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_delete(alignment->counters);
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
      strlen(qualities)!=alignment->read_length,ALIGN_READ_QUAL_LENGTH);
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

GT_INLINE bool gt_alignment_has_qualities(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return (gt_alignment_get_qualities(alignment)!=NULL);
}
GT_INLINE bool gt_alignment_get_not_unique_flag(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return alignment->not_unique_flag;
}

/*
 * Maps Handlers
 */
GT_INLINE char* gt_alignment_record_seq_name(gt_alignment* const alignment,char* const seq_name) {
  register char* key;
  if ((key=gt_shash_get_key(alignment->maps_dictionary,seq_name))) {
    return key;
  } else {
    key = gt_string_cpy(seq_name,strlen(seq_name));
    gt_shash_insert(alignment->maps_dictionary,key,malloc(sizeof(gt_alignment_dictionary)));
    return key;
  }
}
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(map);
  // Alias all map keys into a common dictionary
  map->seq_name = gt_alignment_record_seq_name(alignment,map->seq_name);
  // Insert the map
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
 * Miscellaneous
 */
GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment* alignment_cpy = gt_alignment_new();
  gt_cond_fatal_error(!alignment_cpy,MEM_HANDLER);
  gt_alignment_dup(alignment_cpy,alignment);
  return alignment_cpy;
}
GT_INLINE gt_alignment* gt_alignment_deep_copy(gt_alignment* const alignment) {
  // TODO
  return NULL;
}
GT_INLINE void gt_alignment_handler_dup(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_src);
  alignment_dst->tag = alignment_src->tag;
  alignment_dst->tag_length = alignment_src->tag_length;
  alignment_dst->read = alignment_src->read;
  alignment_dst->read_length = alignment_src->read_length;
  alignment_dst->qualities = alignment_src->qualities;
//  alignment_dst->counters = gt_vector_new(gt_vector_get_used(alignment_src->counters),sizeof(uint64_t));
  gt_vector_copy(alignment_dst->counters,alignment_src->counters);
  alignment_dst->not_unique_flag = alignment_src->not_unique_flag; // FIXME
  alignment_dst->max_complete_strata = alignment_src->max_complete_strata;
}
GT_INLINE void gt_alignment_dup(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_src);
  gt_alignment_handler_dup(alignment_dst,alignment_src);
//  alignment_dst->maps = gt_vector_new(gt_vector_get_used(alignment_src->maps),sizeof(gt_map*));
  gt_vector_copy(alignment_dst->maps,alignment_src->maps);
  alignment_dst->maps_dictionary = alignment_src->maps_dictionary;
  alignment_dst->maps_txt = alignment_src->maps_txt;
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
