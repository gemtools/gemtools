/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_template.h"
#include "gt_iterators.h"

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
  gt_template_delete_handler(template);
}
GT_INLINE void gt_template_delete_handler(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_vector_delete(template->blocks);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->mmaps);
  gt_vector_delete(template->mmaps_attributes);
  free(template);
}
GT_INLINE void gt_template_clear(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_template_clear_blocks(template);
  gt_template_clear_handler(template);
}
GT_INLINE void gt_template_clear_handler(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  template->tag = NULL;
  template->tag_length = 0;
  template->max_complete_strata=0;
  gt_vector_clean(template->blocks);
  gt_vector_clean(template->counters);
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
}

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_ALINGMENT(template) {
    gt_alignment_get_tag(GT_TEMPLATE_REDUCED_ALINGMENT(template));
  }
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
  return template->max_complete_strata; // FIXME: Not defined.......
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,const uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  template->max_complete_strata = max_complete_strata;
}

GT_INLINE bool gt_template_has_qualities(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  if (num_blocks==0) {
    return false;
  } else {
    register uint64_t i;
    for (i=0;i<num_blocks;++i) {
      register gt_alignment* alignment = gt_template_get_block(template,i);
      if (gt_alignment_get_qualities(alignment)==NULL) {
        return false;
      }
    }
    return true;
  }
}
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return template->not_unique_flag;
}

/*
 * Template's matches handlers
 */
GT_INLINE uint64_t gt_template_get_num_mmap(gt_template* const template) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_ALINGMENT(template) { // Alignment Reduction
    return gt_alignment_get_num_maps(gt_template_get_block(template,0));
  }
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t num_mmaps = gt_vector_get_used(template->mmaps);
  gt_fatal_check(num_mmaps % num_blocks != 0,TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION);
  return num_mmaps/num_blocks;
}
GT_INLINE gt_mmap_attributes* gt_template_get_mmap_attr(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  return gt_expect_true(gt_vector_get_used(template->mmaps_attributes)>0) ?
      gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes) : NULL;
}
GT_INLINE void gt_template_clear_mmap(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
}
GT_INLINE void gt_template_clear_mmap_attributes(gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(mmap_attr);
  mmap_attr->score = 0;
  mmap_attr->distance = 0;
}
/* */
GT_INLINE void gt_template_add_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    GT_MAP_CHECK(mmap[i]);
    gt_vector_insert(template->mmaps,mmap[i],gt_map*);
  }
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE gt_map** gt_template_get_mmap(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  // Retrieve the maps from the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** mmap = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  if (mmap_attr) *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
  return mmap;
}
GT_INLINE void gt_template_set_mmap(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  // Store the maps into the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** mmap_template = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  register uint64_t i = 0;
  for (i=0;i<num_blocks;++i) {
    mmap_template[i] = mmap[i];
  }
  gt_vector_set_elm(template->mmaps_attributes,position,gt_mmap_attributes,*mmap_attr);
}
/* */
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
  if (mmap_attr) *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
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
  gt_vector_set_elm(template->mmaps_attributes,position,gt_mmap_attributes,*mmap_attr);
}
/* */
GT_INLINE void gt_template_add_mmap_v(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args) {
  GT_TEMPLATE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    register gt_map* const map = va_arg(v_args,gt_map*);
    GT_MAP_CHECK(map);
    gt_vector_insert(template->mmaps,map,gt_map*);
  }
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE void gt_template_set_mmap_v(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr,va_list v_args) {
  GT_TEMPLATE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  // Store the maps into the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** mmap_template = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  register uint64_t i = 0;
  for (i=0;i<num_blocks;++i) {
    register gt_map* const map = va_arg(v_args,gt_map*);
    GT_MAP_CHECK(map);
    mmap_template[i] = map;
  }
  gt_vector_set_elm(template->mmaps_attributes,position,gt_mmap_attributes,*mmap_attr);
}
/* */
GT_INLINE void gt_template_add_mmap_va(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,...) {
  GT_TEMPLATE_CHECK(template);
  va_list v_args;
  va_start(v_args,mmap_attr);
  gt_template_add_mmap_v(template,mmap_attr,v_args);
}
GT_INLINE void gt_template_set_mmap_va(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr,...) {
  GT_TEMPLATE_CHECK(template);
  va_list v_args;
  va_start(v_args,mmap_attr);
  gt_template_set_mmap_v(template,position,mmap_attr,v_args);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_copy(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_template* template_cpy = gt_template_new();
  gt_cond_fatal_error(!template_cpy,MEM_HANDLER);
  gt_template_dup(template_cpy,template);
  return template_cpy;
}
GT_INLINE gt_template* gt_template_deep_copy(gt_template* const template) {
  return NULL; // TODO
}
GT_INLINE void gt_template_handler_dup(gt_template* template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_src);
  template_dst->tag = template_src->tag;
  template_dst->tag_length = template_src->tag_length;
  template_dst->template_id = template_src->template_id;
  template_dst->in_block_id = template_src->in_block_id;
//  template_dst->blocks = gt_vector_new(gt_vector_get_used(template_src->blocks),sizeof(gt_alignment*));
  gt_vector_copy(template_dst->blocks,template_src->blocks);
//  template_dst->counters = gt_vector_new(gt_vector_get_used(template_src->counters),sizeof(uint64_t));
  gt_vector_copy(template_dst->counters,template_src->counters);
  template_dst->not_unique_flag = template_src->not_unique_flag;
  template_dst->max_complete_strata = template_src->max_complete_strata;
}
GT_INLINE void gt_template_dup(gt_template* template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_src);
  gt_template_handler_dup(template_dst,template_src);
//  template_dst->mmaps = gt_vector_new(gt_vector_get_used(template_src->mmaps),sizeof(gt_map*));
  gt_vector_copy(template_dst->mmaps,template_src->mmaps);
//  template_dst->mmaps_attributes = gt_vector_new(gt_vector_get_used(template_src->mmaps_attributes),sizeof(gt_mmap_attributes));
  gt_vector_copy(template_dst->mmaps_attributes,template_src->mmaps_attributes);
  template_dst->maps_txt = template_src->maps_txt;
}

/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_TEMPLATE_CHECK(template);
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
