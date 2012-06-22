/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_TEMPLATE_H_
#define GT_TEMPLATE_H_

#include "gt_commons.h"
#include "gt_alignment.h"

// Codes gt_status
#define GT_TEMPLATE_OK 1
#define GT_TEMPLATE_FAIL 0

typedef struct {
  uint32_t template_id;
  uint32_t in_block_id;
  char* tag;
  uint64_t tag_length;
  uint64_t max_complete_strata;
  gt_vector* blocks; /* (gt_alignment*) */ /* paired::blocks->used=2 */
  gt_vector* counters;
  gt_vector* maps_relation; /* ( (gt_map*) ) */
  char* maps_txt;
} gt_template;

typedef struct {
  gt_template* template;
  uint64_t num_blocks;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_map** next_map_array;
} gt_template_iterator;

/*
 * Checkers
 */
#define GT_TEMPLATE_CHECK(template) gt_fatal_check( \
  template==NULL||template->blocks==NULL||  \
  template->counters==NULL||template->maps_relation==NULL,NULL_HANDLER)
#define GT_TEMPLATE_CONSISTENCY_CHECK(template) \
  GT_TEMPLATE_CHECK(template); \
  gt_check(gt_vector_get_used(template->blocks)==0,TEMPLATE_ZERO_BLOCKS)
#define GT_TEMPLATE_EDITABLE_CHECK(template) \
  GT_TEMPLATE_CONSISTENCY_CHECK(template); \
  gt_fatal_check(template->maps_txt!=NULL,TEMPLATE_MAPS_NOT_PARSED)

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new();
GT_INLINE void gt_template_delete(gt_template* const template);
GT_INLINE void gt_template_clear(gt_template* const template);

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template);
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag);

GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template);
GT_INLINE void gt_template_set_mcs(gt_template* const template,const uint64_t max_complete_strata);

GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,uint64_t const position);
GT_INLINE void gt_template_clear_blocks(gt_template* const template);
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template);

GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum);

/*
 * Template's matches handlers (Map relation)
 */
// TODO varg funcs()
GT_INLINE void gt_template_add_match_gtvector(gt_template* const template,gt_vector* maps/*(gt_map*)*/);
GT_INLINE uint64_t gt_template_get_num_matches(gt_template* const template);
GT_INLINE void gt_template_get_match(gt_template* const template,const uint64_t position,gt_map** map_array);
GT_INLINE void gt_template_clear_matches(gt_template* const template);

/*
 * Higher-level Accessors
 *   (update global state: counters, ...)
 */
#define GT_TEMPLATE_CALCULATE_DISTANCE UINT64_MAX
// TODO varg funcs()
GT_INLINE void gt_template_insert_match_gtvector(gt_template* const template,const uint64_t total_distance,/*(gt_map*)*/gt_vector* maps);
GT_INLINE void gt_template_recalculate_counters(gt_template* const template);

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_copy(gt_template* const template);
GT_INLINE gt_template* gt_template_deep_copy(gt_template* const template);

// Template's Maps iterator
GT_INLINE void gt_template_iterator_new(gt_template* const template,gt_template_iterator* const template_iterator);
GT_INLINE gt_status gt_template_next_map(gt_template_iterator* const template_iterator,gt_map** const map_array);
GT_INLINE gt_status gt_template_dinamic_next_map(gt_template_iterator* const template_iterator,gt_map** const map_array);
GT_INLINE uint64_t gt_template_next_map_pos(gt_template_iterator* const template_iterator);

#endif /* GT_TEMPLATE_H_ */
