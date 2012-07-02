/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_template.h"

#define GT_TEMPLATE_NUM_INITIAL_BLOCKS 2
#define GT_TEMPLATE_INITIAL_LENGTH_COUNTERS 10
#define GT_TEMPLATE_NUM_INITIAL_MMAPS 20

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new() {
  gt_template* template = malloc(sizeof(gt_template));
  gt_cond_fatal_error(!template,MEM_HANDLER);
  template->template_id = UINT32_MAX;
  template->in_block_id = UINT32_MAX;
  template->tag = NULL;
  template->tag_length = 0;
  template->max_complete_strata=0;
  template->blocks = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_BLOCKS,sizeof(gt_alignment*));
  template->counters = gt_vector_new(GT_TEMPLATE_INITIAL_LENGTH_COUNTERS,sizeof(uint64_t));
  template->mmaps = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_map*));
  template->mmaps_attributes = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_mmap_attributes));
  return template;
}
GT_INLINE void gt_template_delete(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_ITERATE(template->blocks,alg_it,alg_pos,gt_alignment*) {
    gt_alignment_delete(*alg_it);
  }
  gt_vector_delete(template->blocks);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->mmaps);
  gt_vector_delete(template->mmaps_attributes);
}
GT_INLINE void gt_template_clear(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  template->tag = NULL;
  template->tag_length = 0;
  template->max_complete_strata=0;
  gt_template_clear_blocks(template);
  gt_vector_clean(template->blocks);
  gt_vector_clean(template->counters);
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
}
GT_INLINE void gt_template_clear_mmap_attributes(gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(mmap_attr);
  mmap_attr->score = 0;
  mmap_attr->distance = 0;
}

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return template->tag;
}
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag) {
  GT_TEMPLATE_CHECK(template);
  template->tag = tag;
  template->tag_length = strlen(tag);
}

GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  gt_vector_insert(template->blocks,alignment,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  gt_cond_fatal_error(position >= gt_template_get_num_blocks(template),
      POSITION_OUT_OF_RANGE_INFO,(uint64_t)position,0ul,gt_template_get_num_blocks(template)-1);
  return *gt_vector_get_elm(template->blocks,position,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_dyn_get_block(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  if (position < num_blocks) {
    return *gt_vector_get_elm(template->blocks,position,gt_alignment*);
  } else if (position == num_blocks) { // Dynamically allocate a new block (the next one)
    register gt_alignment* dynamic_allocated_alignment = gt_alignment_new();
    gt_template_add_block(template,dynamic_allocated_alignment);
    return dynamic_allocated_alignment;
  } else {
    gt_fatal_error(POSITION_OUT_OF_RANGE_INFO,(uint64_t)position,0ul,num_blocks-1);
  }
}
GT_INLINE void gt_template_clear_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_ITERATE(template->blocks,alg_it,alg_pos,gt_alignment*) {
    gt_alignment_delete(*alg_it);
  }
}
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return gt_vector_get_used(template->blocks);
}

GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return template->counters;
}
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters) {
  GT_TEMPLATE_CHECK(template);
  template->counters = counters;
}
GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return gt_vector_get_used(template->counters);
}
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  return *gt_vector_get_elm(template->counters,stratum-1,uint64_t);
}
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(template->counters,stratum,true);
  ++(*gt_vector_get_elm(template->counters,stratum-1,uint64_t));
}
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(template->counters,stratum,true);
  --(*gt_vector_get_elm(template->counters,stratum-1,uint64_t));
}
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(template->counters,stratum,true);
  *gt_vector_get_elm(template->counters,stratum-1,uint64_t) = value;
}

GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return template->max_complete_strata;
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,const uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  template->max_complete_strata = max_complete_strata;
}

/*
 * Template's matches handlers
 */
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const maps,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_check(gt_template_get_num_blocks(template)!=gt_vector_get_used(maps),TEMPLATE_ADD_BAD_NUM_BLOCKS);
  GT_VECTOR_ITERATE(maps,map_ptr,map_pos,gt_map*) {
    register gt_map* map = *map_ptr;
    GT_MAP_CHECK(map);
    gt_vector_insert(template->mmaps,map,gt_map*);
  }
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const maps,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(maps); GT_NULL_CHECK(mmap_attr);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  // Reset output vector
  gt_vector_prepare(maps,gt_map*,num_blocks);
  // Retrieve the maps from the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** mmap = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  register uint64_t i;
  for(i=0;i<num_blocks;i++) {
    gt_vector_insert(maps,mmap[i],gt_map*);
  }
  *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
}
GT_INLINE void gt_template_set_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const maps,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(maps); GT_NULL_CHECK(mmap_attr);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  gt_fatal_check(gt_vector_get_used(maps)!=num_blocks,TEMPLATE_ADD_BAD_NUM_BLOCKS);
  // Store the maps into the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** mmap = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  register uint64_t i = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    mmap[i++] = *map;
  }
  *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes) = *mmap_attr;
}
GT_INLINE uint64_t gt_template_get_num_mmap(gt_template* const template) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t num_mmaps = gt_vector_get_used(template->mmaps);
  gt_fatal_check(num_mmaps % num_blocks != 0,TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION);
  return num_mmaps/num_blocks;
}
GT_INLINE void gt_template_clear_mmap(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
}
//GT_INLINE void gt_template_add_match_v(gt_template* const template,va_list vl_maps/*(gt_map*)*/) {
//GT_INLINE void gt_template_add_match_va(gt_template* const template,.../*(gt_map*)*/) {
//GT_INLINE void gt_template_get_mmap(  // FIXME: Not public yet
//    gt_template* const template,const uint64_t position,gt_map** mmap_array,gt_mmap_attributes* mmap_attr) {
//  GT_TEMPLATE_EDITABLE_CHECK(template);
//  gt_fatal_check(mmap_array==NULL,NULL_HANDLER);
//  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/gt_template_get_num_blocks(template)),POSITION_OUT_OF_RANGE);
//  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
//  register const uint64_t begin_map_pos = num_blocks*position;
//  register gt_map** map_relation_at_pos = gt_vector_get_elm(template->mmaps,begin_map_pos,gt_map*);
//  register uint64_t i;
//  for(i=0;i<num_blocks;i++) {
//    mmap_array[i] = map_relation_at_pos[i];
//  }
//  mmap_attr = gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
//}

/*
 * Higher-level Procedures
 */
GT_INLINE void gt_template_insert_match_gtvector(gt_template* const template,gt_vector* const maps,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(maps); GT_NULL_CHECK(mmap_attr);
  gt_template_inc_counter(template,mmap_attr->distance);
  gt_template_add_mmap_gtvector(template,maps,mmap_attr);
}
GT_INLINE void gt_template_recalculate_counters(gt_template* const template) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  // TODO
}
GT_INLINE uint64_t gt_template_get_min_matching_strata(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  register gt_vector* vector = gt_template_get_counters_vector(template);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (*counter!=0) return counter_pos+1;
  }
  return UINT64_MAX;
}
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  register gt_vector* vector = gt_template_get_counters_vector(template);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if ((counter_pos+1)>=max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE bool gt_template_is_mapped(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  return gt_template_is_thresholded_mapped(template,UINT64_MAX);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_copy(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_template* template_cpy = malloc(sizeof(gt_template));
  gt_cond_fatal_error(!template_cpy,MEM_HANDLER);
  template_cpy->template_id = template->template_id;
  template_cpy->in_block_id = template->in_block_id;
  template_cpy->blocks = gt_vector_new(gt_vector_get_used(template->blocks),sizeof(gt_alignment*));
  gt_vector_copy(template_cpy->blocks,template->blocks);
  template_cpy->counters = gt_vector_new(gt_vector_get_used(template->counters),sizeof(uint64_t));
  gt_vector_copy(template_cpy->counters,template->counters);
  template_cpy->mmaps = gt_vector_new(gt_vector_get_used(template->mmaps),sizeof(gt_map*));
  gt_vector_copy(template_cpy->mmaps,template->mmaps);
  return template_cpy;
}
GT_INLINE gt_template* gt_template_deep_copy(gt_template* const template) {
  return NULL; // TODO
}

/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(template_alignment_iterator);
  template_alignment_iterator->template = template;
  template_alignment_iterator->next_pos = 0;
  template_alignment_iterator->total_pos = gt_vector_get_used(template->blocks);
  template_alignment_iterator->next_alignment = gt_vector_get_mem(template->blocks,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  if (gt_expect_true(template_alignment_iterator->next_pos<template_alignment_iterator->total_pos)) {
    register gt_alignment* const alignment = *template_alignment_iterator->next_alignment;
    ++template_alignment_iterator->next_alignment;
    ++template_alignment_iterator->next_pos;
    return alignment;
  } else {
    return NULL;
  }
}
GT_INLINE gt_alignment* gt_template_dinamic_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  register gt_template* const template = template_alignment_iterator->template;
  if (gt_expect_true(template_alignment_iterator->next_pos<gt_vector_get_used(template->blocks))) {
    register gt_alignment* const alignment =
        *gt_vector_get_elm(template->blocks,template_alignment_iterator->next_pos,gt_alignment*);
    ++template_alignment_iterator->next_pos;
    return alignment;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_template_next_alignment_pos(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  return template_alignment_iterator->next_pos;
}

/*
 * Template's Maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
 */
GT_INLINE void gt_template_new_maps_iterator(
    gt_template* const template,gt_template_maps_iterator* const template_maps_iterator) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(template_maps_iterator);
  template_maps_iterator->template = template;
  template_maps_iterator->num_blocks = gt_vector_get_used(template->blocks);
  template_maps_iterator->next_pos = 0;
  template_maps_iterator->total_pos = gt_vector_get_used(template->mmaps);
  template_maps_iterator->next_map_array = gt_vector_get_mem(template->mmaps,gt_map*);
}
GT_INLINE gt_status gt_template_next_maps(
    gt_template_maps_iterator* const template_maps_iterator,gt_map*** const map_array) {
  GT_NULL_CHECK(template_maps_iterator); GT_NULL_CHECK(map_array);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_maps_iterator->template);
  register const uint64_t num_blocks = template_maps_iterator->num_blocks;
  if (gt_expect_true(template_maps_iterator->next_pos+num_blocks<=template_maps_iterator->total_pos)) {
    *map_array = template_maps_iterator->next_map_array;
    template_maps_iterator->next_map_array+=num_blocks;
    template_maps_iterator->next_pos+=num_blocks;
    return GT_TEMPLATE_OK;
  } else {
    *map_array = NULL;
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE gt_status gt_template_dinamic_next_maps(
    gt_template_maps_iterator* const template_maps_iterator,gt_map*** const map_array) {
  GT_NULL_CHECK(template_maps_iterator); GT_NULL_CHECK(map_array);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_maps_iterator->template);
  register const uint64_t num_blocks = template_maps_iterator->num_blocks;
  register const gt_template* const template = template_maps_iterator->template;
  if (gt_expect_true(template_maps_iterator->next_pos+num_blocks<=gt_vector_get_used(template->blocks))) {
    *map_array = gt_vector_get_elm(template->mmaps,template_maps_iterator->next_pos,gt_map*);
    template_maps_iterator->next_pos+=num_blocks;
    return GT_TEMPLATE_OK;
  } else {
    *map_array = NULL;
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE uint64_t gt_template_next_maps_pos(gt_template_maps_iterator* const template_maps_iterator) {
  GT_NULL_CHECK(template_maps_iterator);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_maps_iterator->template);
  return template_maps_iterator->next_pos/template_maps_iterator->num_blocks;
}
