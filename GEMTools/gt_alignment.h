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
  char* tag;
  uint64_t tag_length;
  char* read;
  uint64_t read_length;
  char* qualities;
  gt_vector* counters;
  uint64_t max_complete_strata;
  gt_vector* maps; /* (gt_map*) */
  char* maps_txt;
} gt_alignment;

typedef struct {
  gt_alignment* alignment;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_map** next_map;
} gt_alignment_map_iterator;

/*
 * Checkers
 */
#define GT_ALIGNMENT_CHECK(alignment) gt_fatal_check( \
    alignment==NULL||alignment->counters==NULL||alignment->maps==NULL,NULL_HANDLER)
#define GT_ALIGNMENT_EDITABLE_CHECK(alignment) \
  GT_ALIGNMENT_CHECK(alignment); \
  gt_fatal_check(alignment->maps_txt!=NULL,ALIGN_MAPS_NOT_PARSED)

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
GT_INLINE void gt_alignment_set_tag(gt_alignment* const alignment,char* const tag);

GT_INLINE char* gt_alignment_get_read(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_read(gt_alignment* const alignment,char* const read);

GT_INLINE char* gt_alignment_get_qualities(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_qualities(gt_alignment* const alignment,char* const qualities);

GT_INLINE gt_vector* gt_alignment_get_counters_vector(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_counters_vector(gt_alignment* const alignment,gt_vector* const counters);

GT_INLINE uint64_t gt_alignment_get_num_counters(gt_alignment* const alignment);
GT_INLINE uint64_t gt_alignment_get_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_set_counter(gt_alignment* const alignment,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_alignment_dec_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_inc_counter(gt_alignment* const alignment,const uint64_t stratum);

GT_INLINE uint64_t gt_alignment_get_mcs(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_mcs(gt_alignment* const alignment,const uint64_t max_complete_strata);

/*
 * Maps Handlers
 */
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_clear_maps(gt_alignment* const alignment);
GT_INLINE gt_map* gt_alignment_get_map(gt_alignment* const alignment,const uint64_t position);
GT_INLINE uint64_t gt_alignment_get_num_maps(gt_alignment* const alignment);

/*
 * Higher-level Methods
 *   (Update global state: counters, ...)
 */
GT_INLINE void gt_alignment_insert_map__check_dup(gt_alignment* const alignment,gt_map** const map);
GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment);
GT_INLINE uint64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment);
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const uint64_t max_allowed_strata);
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment);

/*
 * Miscellaneous
 */
GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_alignment_deep_copy(gt_alignment* const alignment);

// Alignment's Maps iterator
GT_INLINE void gt_alignment_new_map_iterator(gt_alignment* const alignment,gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE gt_map* gt_alignment_next_map(gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE gt_map* gt_alignment_dinamic_next_map(gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE uint64_t gt_alignment_next_map_pos(gt_alignment_map_iterator* const alignment_map_iterator);

#endif /* GT_ALIGNMENT_H_ */
