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
#include "gt_iterators.h"

// Codes gt_status
#define GT_TEMPLATE_OK 1
#define GT_TEMPLATE_FAIL 0

typedef struct {
  uint32_t template_id;
  uint32_t in_block_id;
  char* tag;
  uint64_t tag_length;
  uint64_t max_complete_strata;
  bool not_unique_flag;
  gt_vector* blocks; /* (gt_alignment*) */ /* paired::blocks->used=2 */
  gt_vector* counters;
  gt_vector* mmaps; /* ( (gt_map*) ) */
  gt_vector* mmaps_attributes; /* ( (gt_mmap_attributes) ) */
  char* maps_txt;
} gt_template;
typedef struct {
  uint64_t distance;
  uint64_t score;
} gt_mmap_attributes;

// Iterators
typedef struct {
  gt_template* template;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_alignment** next_alignment;
} gt_template_alignment_iterator;
typedef struct {
  gt_template* template;
  uint64_t num_blocks;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_map** next_map_array;
} gt_template_maps_iterator;

/*
 * Checkers
 */
#define GT_TEMPLATE_CHECK(template) gt_fatal_check( \
  template==NULL||template->blocks==NULL|| \
  template->counters==NULL||template->mmaps==NULL|| \
  template->mmaps_attributes==NULL,NULL_HANDLER)
#define GT_TEMPLATE_CONSISTENCY_CHECK(template) \
  GT_TEMPLATE_CHECK(template); \
  gt_check(gt_vector_get_used(template->blocks)==0,TEMPLATE_ZERO_BLOCKS)
#define GT_TEMPLATE_EDITABLE_CHECK(template) \
  GT_TEMPLATE_CONSISTENCY_CHECK(template); \
  gt_fatal_check(template->maps_txt!=NULL,TEMPLATE_MAPS_NOT_PARSED)

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new(void);
GT_INLINE void gt_template_delete(gt_template* const template);
GT_INLINE void gt_template_clear(gt_template* const template);
GT_INLINE void gt_template_clear_mmap_attributes(gt_mmap_attributes* const mmap_attr);

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template);
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag);

GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position);
GT_INLINE gt_alignment* gt_template_dyn_get_block(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_clear_blocks(gt_template* const template);
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template);

GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template);
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters);

GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum);

GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template);
GT_INLINE void gt_template_set_mcs(gt_template* const template,const uint64_t max_complete_strata);

GT_INLINE bool gt_template_has_qualities(gt_template* const template);
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template);

/*
 * Template's multimaps handlers (Map relation)
 */
GT_INLINE uint64_t gt_template_get_num_mmap(gt_template* const template);
GT_INLINE gt_mmap_attributes* gt_template_get_mmap_attr(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_clear_mmap(gt_template* const template);
/* */
GT_INLINE void gt_template_add_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE gt_map** gt_template_get_mmap(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_set_mmap(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
/* */
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const maps,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const maps,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_set_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const maps,gt_mmap_attributes* const mmap_attr);
/* */
GT_INLINE void gt_template_add_mmap_v(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args);
GT_INLINE void gt_template_set_mmap_v(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr,va_list v_args);
/* */
GT_INLINE void gt_template_add_mmap_va(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,...);
GT_INLINE void gt_template_set_mmap_va(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr,...);

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_copy(gt_template* const template);
GT_INLINE gt_template* gt_template_deep_copy(gt_template* const template);
// Template's Alignments iterator (end1,end2, ... )
GT_INLINE void gt_template_new_alignment_iterator(gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE gt_alignment* gt_template_dinamic_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE uint64_t gt_template_next_alignment_pos(gt_template_alignment_iterator* const template_alignment_iterator);
// Template's Maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
GT_INLINE void gt_template_new_maps_iterator(gt_template* const template,gt_template_maps_iterator* const template_maps_iterator);
GT_INLINE gt_status gt_template_next_maps(gt_template_maps_iterator* const template_maps_iterator,gt_map*** const map_array);
GT_INLINE gt_status gt_template_dinamic_next_maps(gt_template_maps_iterator* const template_maps_iterator,gt_map*** const map_array);
GT_INLINE uint64_t gt_template_next_maps_pos(gt_template_maps_iterator* const template_maps_iterator);

#endif /* GT_TEMPLATE_H_ */
