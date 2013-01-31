/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_alignment.h"

#define GT_ALIGNMENT_TAG_INITIAL_LENGTH 100
#define GT_ALIGNMENT_READ_INITIAL_LENGTH 150
#define GT_ALIGNMENT_NUM_INITIAL_MAPS 20
#define GT_ALIGNMENT_NUM_INITIAL_COUNTERS 5

typedef struct {
  uint64_t dummy; // TODO: Extensible dictionary element
} gt_maps_dictionary_element;

/*
 * Setup
 */
GT_INLINE gt_alignment* gt_alignment_new() {
  gt_alignment* alignment = malloc(sizeof(gt_alignment));
  gt_cond_fatal_error(!alignment,MEM_HANDLER);
  alignment->alignment_id = UINT32_MAX;
  alignment->in_block_id = UINT32_MAX;
  alignment->tag = gt_string_new(GT_ALIGNMENT_TAG_INITIAL_LENGTH);
  alignment->read = gt_string_new(GT_ALIGNMENT_READ_INITIAL_LENGTH);
  alignment->qualities = gt_string_new(GT_ALIGNMENT_READ_INITIAL_LENGTH);
  alignment->counters = gt_vector_new(GT_ALIGNMENT_NUM_INITIAL_COUNTERS,sizeof(uint64_t));
  alignment->maps = gt_vector_new(GT_ALIGNMENT_NUM_INITIAL_MAPS,sizeof(gt_map));
  alignment->maps_txt = NULL;
  alignment->attributes = gt_attribute_new();
  return alignment;
}
GT_INLINE void gt_alignment_clear_handler(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_string_clear(alignment->tag);
  gt_string_clear(alignment->read);
  gt_string_clear(alignment->qualities);
  alignment->maps_txt = NULL;
  gt_attribute_clear(alignment->attributes);
}
GT_INLINE void gt_alignment_clear(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_clear_maps(alignment);
  gt_vector_clear(alignment->counters);
  gt_alignment_clear_handler(alignment);
}
GT_INLINE void gt_alignment_delete(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_clear_maps(alignment);
  gt_string_delete(alignment->tag);
  gt_string_delete(alignment->read);
  gt_string_delete(alignment->qualities);
  gt_vector_delete(alignment->counters);
  gt_vector_delete(alignment->maps);
  gt_attribute_delete(alignment->attributes);
  free(alignment);
}

/*
 * Accessors
 */
GT_INLINE char* gt_alignment_get_tag(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_string_get_string(alignment->tag);
}
GT_INLINE void gt_alignment_set_tag(gt_alignment* const alignment,char* const tag,const uint64_t length) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(tag);
  gt_string_set_nstring(alignment->tag,tag,length);
}
GT_INLINE uint64_t gt_alignment_get_tag_length(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_string_get_length(alignment->tag);
}
GT_INLINE char* gt_alignment_get_read(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_string_get_string(alignment->read);
}
GT_INLINE void gt_alignment_set_read(gt_alignment* const alignment,char* const read,const uint64_t length) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(read);
  gt_string_set_nstring(alignment->read,read,length);
  gt_fatal_check(!gt_string_is_null(alignment->qualities) &&
      gt_string_get_length(alignment->qualities)!=gt_string_get_length(alignment->read),ALIGNMENT_READ_QUAL_LENGTH);
}
GT_INLINE uint64_t gt_alignment_get_read_length(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_string_get_length(alignment->read);
}
GT_INLINE char* gt_alignment_get_qualities(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_string_get_string(alignment->qualities);
}
GT_INLINE void gt_alignment_set_qualities(gt_alignment* const alignment,char* const qualities,const uint64_t length) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(qualities);
  gt_string_set_nstring(alignment->qualities,qualities,length);
  gt_fatal_check(!gt_string_is_null(alignment->read) &&
      gt_string_get_length(alignment->qualities)!=gt_string_get_length(alignment->read),ALIGNMENT_READ_QUAL_LENGTH);
}
GT_INLINE bool gt_alignment_has_qualities(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return !gt_string_is_null(alignment->qualities);
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
  return *gt_vector_get_elm(alignment->counters,stratum,uint64_t);
}
GT_INLINE void gt_alignment_dynamically_allocate_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  register const uint64_t used_strata = stratum+1;
  gt_vector_reserve(alignment->counters,used_strata,true);
  if (gt_vector_get_used(alignment->counters) < used_strata) {
    gt_vector_set_used(alignment->counters,used_strata);
  }
}
GT_INLINE void gt_alignment_set_counter(gt_alignment* const alignment,const uint64_t stratum,const uint64_t value) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_dynamically_allocate_counter(alignment,stratum);
  *gt_vector_get_elm(alignment->counters,stratum,uint64_t) = value;
}
GT_INLINE void gt_alignment_dec_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_dynamically_allocate_counter(alignment,stratum);
  --(*gt_vector_get_elm(alignment->counters,stratum,uint64_t));
}
GT_INLINE void gt_alignment_inc_counter(gt_alignment* const alignment,const uint64_t stratum) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment_dynamically_allocate_counter(alignment,stratum);
  ++(*gt_vector_get_elm(alignment->counters,stratum,uint64_t));
}

/*
 * Predefined attributes
 */
GT_INLINE uint64_t gt_alignment_get_mcs(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register uint64_t* const mcs = (uint64_t*)gt_attribute_get(alignment->attributes,GT_ATTR_MAX_COMPLETE_STRATA);
  if (mcs == NULL) return UINT64_MAX;
  return *mcs;
}
GT_INLINE void gt_alignment_set_mcs(gt_alignment* const alignment,uint64_t max_complete_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_attribute_set(alignment->attributes,GT_ATTR_MAX_COMPLETE_STRATA,&max_complete_strata,sizeof(uint64_t));
}
GT_INLINE bool gt_alignment_get_not_unique_flag(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register bool* const not_unique_flag = (bool*)gt_attribute_get(alignment->attributes,GT_ATTR_NOT_UNIQUE);
  if (not_unique_flag==NULL) return false;
  return *not_unique_flag;
}
GT_INLINE void gt_alignment_set_not_unique_flag(gt_alignment* const alignment,bool is_not_unique) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_attribute_set(alignment->attributes,GT_ATTR_NOT_UNIQUE,&is_not_unique,sizeof(bool));
}

/*
 * Maps Handlers
 */
GT_INLINE uint64_t gt_alignment_get_num_maps(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_vector_get_used(alignment->maps);
}
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(map);
  // Insert the map
  gt_vector_insert(alignment->maps,map,gt_map*);
}
GT_INLINE void gt_alignment_add_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map_it,map_count,gt_map*) {
    GT_MAP_CHECK(*map_it);
    gt_alignment_add_map(alignment,*map_it);
  }
}
GT_INLINE gt_map* gt_alignment_get_map(gt_alignment* const alignment,const uint64_t position) {
  GT_ALIGNMENT_CHECK(alignment);
  return *gt_vector_get_elm(alignment->maps,position,gt_map*);
}
GT_INLINE void gt_alignment_set_map(gt_alignment* const alignment,gt_map* const map,const uint64_t position) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_MAP_CHECK(map);
  // Insert the map
  *gt_vector_get_elm(alignment->maps,position,gt_map*) = map;
}
GT_INLINE void gt_alignment_clear_maps(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_ITERATE(alignment->maps,alg_map,alg_map_pos,gt_map*) {
    gt_map_delete(*alg_map);
  }
  gt_vector_clear(alignment->maps);
}
GT_INLINE bool gt_alignment_locate_map_reference(gt_alignment* const alignment,gt_map* const map,uint64_t* const position) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_MAP_CHECK(map);
  GT_VECTOR_ITERATE(alignment->maps,alg_map,alg_map_pos,gt_map*) {
    if (*alg_map==map) { /* Cmp references */
      *position = alg_map_pos;
      return true;
    }
  }
  return false;
}

/*
 * Miscellaneous
 */
GT_INLINE void gt_alignment_handler_copy(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_src);
  alignment_dst->alignment_id = alignment_src->alignment_id;
  alignment_dst->in_block_id = alignment_src->in_block_id;
  gt_string_copy(alignment_dst->tag,alignment_src->tag);
  gt_string_copy(alignment_dst->read,alignment_src->read);
  gt_string_copy(alignment_dst->qualities,alignment_src->qualities);
  alignment_dst->maps_txt = alignment_src->maps_txt;
  // Copy attributes
  gt_attribute_delete(alignment_dst->attributes); // FIXME
  alignment_dst->attributes = gt_shash_deep_copy(alignment_src->attributes);
}

GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment,const bool copy_maps) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_alignment* alignment_cp = gt_alignment_new();
  gt_cond_fatal_error(!alignment_cp,MEM_HANDLER);
  // Copy handler
  gt_alignment_handler_copy(alignment_cp,alignment);
  // Copy maps
  if (copy_maps) {
    // Copy map related fields (deep copy) {MAPS,MAPS_DICCTIONARY,COUNTERS,ATTRIBUTES}
    gt_vector_copy(alignment_cp->counters,alignment->counters);
    GT_VECTOR_ITERATE(alignment->maps,alg_map,alg_map_pos,gt_map*) {
      gt_alignment_add_map(alignment_cp,gt_map_copy(*alg_map));
    }
  }
  return alignment_cp;
}

/*
 * Alignment's Maps iterator
 */
GT_INLINE void gt_alignment_new_map_iterator(gt_alignment* const alignment,gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_CHECK(alignment);
  alignment_map_iterator->alignment = alignment;
  alignment_map_iterator->next_pos = 0;
}
GT_INLINE gt_map* gt_alignment_next_map(gt_alignment_map_iterator* const alignment_map_iterator) {
  GT_NULL_CHECK(alignment_map_iterator);
  GT_ALIGNMENT_CHECK(alignment_map_iterator->alignment);
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
  GT_ALIGNMENT_CHECK(alignment_map_iterator->alignment);
  return alignment_map_iterator->next_pos;
}
