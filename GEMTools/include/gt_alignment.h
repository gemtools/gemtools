/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_ALIGNMENT_H_
#define GT_ALIGNMENT_H_

#include "gt_commons.h"
#include "gt_map.h"

typedef struct {
  uint32_t alignment_id;
  uint32_t in_block_id;
  gt_string* tag;
  gt_string* read;
  gt_string* qualities;
  gt_vector* counters;
  gt_vector* maps; /* (gt_map*) */
  char* maps_txt;
  gt_shash* maps_dictionary;
  gt_shash* attributes;
} gt_alignment;

typedef struct {
  gt_alignment* alignment;
  uint64_t next_pos;
} gt_alignment_map_iterator;

/*
 * Checkers
 */
#define GT_ALIGNMENT_CHECK(alignment) \
  gt_fatal_check(alignment==NULL,NULL_HANDLER); \
  GT_STRING_CHECK(alignment->tag); \
  GT_STRING_CHECK(alignment->read); \
  GT_STRING_CHECK(alignment->qualities); \
  GT_VECTOR_CHECK(alignment->counters); \
  GT_VECTOR_CHECK(alignment->maps); \
  GT_HASH_CHECK(alignment->maps_dictionary); \
  GT_HASH_CHECK(alignment->attributes)
/*
 *  TODO: Scheduled for v2.0 (all lazy parsing)
 *  TODO: Check
 *  #define GT_ALIGNMENT_EDITABLE_CHECK(alignment) \
 *    GT_ALIGNMENT_CHECK(alignment); \
 *    gt_fatal_check(alignment->maps_txt!=NULL,ALIGN_MAPS_NOT_PARSED)
 */

/*
 * Setup
 */
GT_INLINE gt_alignment* gt_alignment_new(void);
GT_INLINE void gt_alignment_clear(gt_alignment* const alignment);
GT_INLINE void gt_alignment_delete(gt_alignment* const alignment);

/*
 * Accessors
 */
GT_INLINE char* gt_alignment_get_tag(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_tag(gt_alignment* const alignment,char* const tag,const uint64_t length);
GT_INLINE uint64_t gt_alignment_get_tag_length(gt_alignment* const alignment);

GT_INLINE char* gt_alignment_get_read(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_read(gt_alignment* const alignment,char* const read,const uint64_t length);
GT_INLINE uint64_t gt_alignment_get_read_length(gt_alignment* const alignment);

GT_INLINE char* gt_alignment_get_qualities(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_qualities(gt_alignment* const alignment,char* const qualities,const uint64_t length);
GT_INLINE bool gt_alignment_has_qualities(gt_alignment* const alignment);

GT_INLINE gt_vector* gt_alignment_get_counters_vector(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_counters_vector(gt_alignment* const alignment,gt_vector* const counters);
GT_INLINE uint64_t gt_alignment_get_num_counters(gt_alignment* const alignment);
GT_INLINE uint64_t gt_alignment_get_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_set_counter(gt_alignment* const alignment,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_alignment_dec_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_inc_counter(gt_alignment* const alignment,const uint64_t stratum);

/*
 * Attribute accessors
 */
// Predefined attributes
#define GT_ATTR_MAX_COMPLETE_STRATA "MCS"
#define GT_ATTR_NOT_UNIQUE "NOT-UNIQUE"
GT_INLINE void* gt_alignment_get_attr(
    gt_alignment* const alignment,char* const attribute_id);
GT_INLINE void gt_alignment_set_attr(
    gt_alignment* const alignment,char* const attribute_id,
    void* const attribute,const size_t element_size);
#define gt_alignment_get_attribute(alignment,attribute_id,type) ((type*)gt_alignment_get_attr(alignment,attribute_id))
#define gt_alignment_set_attribute(alignment,attribute_id,attribute,type) gt_alignment_set_attr(alignment,attribute_id,attribute,sizeof(type))

GT_INLINE uint64_t gt_alignment_get_mcs(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_mcs(gt_alignment* const alignment,uint64_t max_complete_strata);
GT_INLINE void gt_alignment_set_not_unique_flag(gt_alignment* const alignment,bool is_not_unique);
GT_INLINE bool gt_alignment_get_not_unique_flag(gt_alignment* const alignment);

/*
 * Maps Handlers
 */
GT_INLINE uint64_t gt_alignment_get_num_maps(gt_alignment* const alignment);
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_add_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector);
GT_INLINE gt_map* gt_alignment_get_map(gt_alignment* const alignment,const uint64_t position);
GT_INLINE void gt_alignment_set_map(gt_alignment* const alignment,gt_map* const map,const uint64_t position);
GT_INLINE void gt_alignment_clear_maps(gt_alignment* const alignment);

GT_INLINE bool gt_alignment_locate_map_reference(gt_alignment* const alignment,gt_map* const map,uint64_t* const position);

/*
 * Miscellaneous
 */
GT_INLINE void gt_alignment_handler_copy(gt_alignment* const alignment_dst,gt_alignment* const alignment_src);
GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_alignment_deep_copy(gt_alignment* const alignment);

// Alignment's Maps iterator
GT_INLINE void gt_alignment_new_map_iterator(gt_alignment* const alignment,gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE gt_map* gt_alignment_next_map(gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE uint64_t gt_alignment_next_map_pos(gt_alignment_map_iterator* const alignment_map_iterator);

/*
 * Iterate over the map of an alignment
 *   Alignment = {(map1),(map2.block1,map2.block2),(map3),(map4)}
 *   GT_ALIGNMENT_ITERATE := {(map1),(map2.block1,map2.block2),(map3),(map4)}
 */
#define GT_ALIGNMENT_ITERATE(alignment,map) \
  /* Alignment. Iterate over all maps */ \
  gt_alignment_map_iterator __##map##_iterator; \
  register gt_map* map; \
  gt_alignment_new_map_iterator(alignment,&(__##map##_iterator)); \
  while ((map=gt_alignment_next_map(&(__##map##_iterator))))

#endif /* GT_ALIGNMENT_H_ */
