/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_template.h"

#define GT_TEMPLATE_NUM_INITIAL_BLOCKS 2
#define GT_TEMPLATE_INITIAL_LENGTH_COUNTERS 10
#define GT_TEMPLATE_NUM_INITIAL_MAPS_RELATION 20

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
  template->maps_relation = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MAPS_RELATION,sizeof(gt_map*));
  return template;
}
GT_INLINE void gt_template_delete(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_ITERATE(template->blocks,alg_it,alg_pos,gt_alignment*) {
    gt_alignment_delete(*alg_it);
  }
  gt_vector_delete(template->blocks);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->maps_relation);
}
GT_INLINE void gt_template_clear(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  template->tag = NULL;
  template->tag_length = 0;
  template->max_complete_strata=0;
  gt_template_clear_blocks(template);
  gt_vector_clean(template->blocks);
  gt_vector_clean(template->counters);
  gt_vector_clean(template->maps_relation);
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

GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return template->max_complete_strata;
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,const uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  template->max_complete_strata = max_complete_strata;
}

GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  gt_vector_insert(template->blocks,alignment,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,uint64_t const position) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  return *gt_vector_get_elm(template->blocks,position,gt_alignment*);
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
  (*gt_vector_get_elm(template->counters,stratum-1,uint64_t))++;
}
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  gt_vector_reserve(template->counters,stratum,true);
  (*gt_vector_get_elm(template->counters,stratum-1,uint64_t))--;
}
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  // Check distant counter
  register const uint64_t counts_used = gt_vector_get_used(template->counters);
  if (stratum > counts_used) { // Clear intermediate counters
    gt_vector_reserve(template->counters,stratum,false);
    memset(gt_vector_get_mem(template->counters,uint64_t)+counts_used,0,
        sizeof(uint64_t)*(stratum-counts_used));
  }
  *gt_vector_get_elm(template->counters,stratum-1,uint64_t) = value;
}

/*
 * Template's matches handlers
 */
//GT_INLINE void gt_template_add_match_v(gt_template* const template,va_list vl_maps/*(gt_map*)*/) {
//  GT_TEMPLATE_EDITABLE_CHECK(template);
//  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
//  register uint64_t i;
//  va_start(vl_maps,num_blocks);
//  for(i=0;i<num_blocks;i++){
//    register gt_map* map = va_arg(vl_maps,gt_map*);
//    GT_MAP_CHECK(map);
//    gt_vector_insert(template->maps_relation,map,gt_map*);
//  }
//  va_end(vl_maps);
//}
//GT_INLINE void gt_template_add_match_va(gt_template* const template,.../*(gt_map*)*/) {
//  GT_TEMPLATE_EDITABLE_CHECK(template);
//  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
//  va_list vl_maps;
//  va_start(vl_maps,num_blocks);
//  gt_template_add_match_v(template,vl_maps);
//  va_end(vl_maps);
//}
GT_INLINE void gt_template_add_match_gtvector(gt_template* const template,gt_vector* maps/*(gt_map*)*/) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  GT_VECTOR_ITERATE(maps,map_ptr,map_pos,gt_map*) {
    register gt_map* map = *map_ptr;
    GT_MAP_CHECK(map);
    gt_vector_insert(template->maps_relation,map,gt_map*);
  }
}
GT_INLINE uint64_t gt_template_get_num_matches(gt_template* const template) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t num_maps_relation = gt_vector_get_used(template->maps_relation);
  gt_fatal_check(num_maps_relation % num_blocks != 0,TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION);
  return num_maps_relation/num_blocks;
}
GT_INLINE void gt_template_get_match(gt_template* const template,const uint64_t position,gt_map** map_array) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  gt_fatal_check(map_array==NULL,NULL_HANDLER);
  gt_fatal_check(position>=(gt_vector_get_used(template->maps_relation)/gt_template_get_num_blocks(template)),
      POSITION_OUT_OF_RANGE);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t begin_map_pos = num_blocks*position;
  register gt_map** map_relation_at_pos = gt_vector_get_elm(template->maps_relation,begin_map_pos,gt_map*);
  register uint64_t i;
  for(i=0;i<num_blocks;i++) {
    map_array[i] = map_relation_at_pos[i];
  }
}
GT_INLINE void gt_template_clear_matches(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_vector_clean(template->maps_relation);
}

/*
 * Higher-level Accessors
 *   (update global state: counters, ...)
 */
GT_INLINE void gt_template_insert_match_gtvector(gt_template* const template,uint64_t total_distance,/*(gt_map*)*/gt_vector* maps) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  if (total_distance==GT_TEMPLATE_CALCULATE_DISTANCE) {
    total_distance = 0;
    GT_VECTOR_ITERATE(maps,map_ptr,map_pos,gt_map*) {
      GT_MAP_CHECK((*map_ptr));
      total_distance += (*map_ptr)->distance;
    }
  }
  gt_template_inc_counter(template,total_distance);
  gt_template_add_match_gtvector(template,maps);
}
/*
 * GT_INLINE void gt_template_insert_match_v(gt_template* const template,const uint64_t total_distance,va_list vl_maps);
 * GT_INLINE void gt_template_insert_match_va(gt_template* const template,const uint64_t total_distance,gt_map* map,...);
 */
//GT_INLINE void gt_template_insert_match_v(gt_template* const template,const uint64_t total_distance,va_list vl_maps/*(gt_map*)*/) {
//  GT_TEMPLATE_EDITABLE_CHECK(template);
//  if (total_distance==GT_TEMPLATE_CALCULATE_DISTANCE) {
//    register const uint64_t num_blocks = gt_template_get_num_blocks(template);
//    register uint64_t i;
//    va_list vl_maps;
//    va_start(vl_maps,num_blocks);
//    for(i=0;i<num_blocks;i++){
//      register gt_map* map = va_arg(vl_maps,gt_map*);
//      GT_MAP_CHECK(map);
//      gt_vector_insert(template->maps_relation,map,gt_map*);
//    }
//    va_end(vl_maps);
//  }
//
//}
//GT_INLINE void gt_template_insert_match_va(gt_template* const template,const uint64_t total_distance,gt_map* map,...) {
//  GT_TEMPLATE_EDITABLE_CHECK(template);
//  va_list vl_maps;
//  va_start(vl_maps,map);
//  gt_template_insert_match_v(template,total_distance,vl_maps);
//  va_end(vl_maps);
//}
GT_INLINE void gt_template_recalculate_counters(gt_template* const template) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  // TODO
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
  template_cpy->maps_relation = gt_vector_new(gt_vector_get_used(template->maps_relation),sizeof(gt_map*));
  gt_vector_copy(template_cpy->maps_relation,template->maps_relation);
  return template_cpy;
}
GT_INLINE gt_template* gt_template_deep_copy(gt_template* const template) {
  return NULL; // TODO
}

/*
 * Template's Maps iterator
 */
GT_INLINE void gt_template_iterator_new(
    gt_template* const template,gt_template_iterator* const template_iterator) {
  GT_TEMPLATE_EDITABLE_CHECK(template);
  gt_fatal_check(template_iterator==NULL,NULL_HANDLER);
  template_iterator->template = template;
  template_iterator->num_blocks = gt_vector_get_used(template->blocks);
  template_iterator->next_pos = 0;
  template_iterator->total_pos = gt_vector_get_used(template->maps_relation);
  template_iterator->next_map_array = gt_vector_get_mem(template->maps_relation,gt_map*);
}
GT_INLINE gt_status gt_template_next_map(
    gt_template_iterator* const template_iterator,gt_map** const map_array) {
  gt_fatal_check(template_iterator==NULL||map_array==NULL,NULL_HANDLER);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_iterator->template);
  register const uint64_t num_blocks = template_iterator->num_blocks;
  if (gt_expect_true(template_iterator->next_pos+num_blocks<=template_iterator->total_pos)) {
    register gt_map** const map_mem = template_iterator->next_map_array;
    register uint64_t i;
    for (i=0;i<num_blocks;++i) {
      map_array[i]=map_mem[i];
    }
    template_iterator->next_map_array+=num_blocks;
    template_iterator->next_pos+=num_blocks;
    return GT_TEMPLATE_OK;
  } else {
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE gt_status gt_template_dinamic_next_map(
    gt_template_iterator* const template_iterator,gt_map** const map_array) {
  gt_fatal_check(template_iterator==NULL||map_array==NULL,NULL_HANDLER);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_iterator->template);
  register const uint64_t num_blocks = template_iterator->num_blocks;
  register const gt_template* const template = template_iterator->template;
  if (gt_expect_true(template_iterator->next_pos+num_blocks<=gt_vector_get_used(template->blocks))) {
    register gt_map** const map_mem =
        gt_vector_get_elm(template->maps_relation,template_iterator->next_pos,gt_map*);
    register uint64_t i;
    for (i=0;i<num_blocks;++i) {
      map_array[i]=map_mem[i];
    }
    template_iterator->next_map_array+=num_blocks;
    template_iterator->next_pos+=num_blocks;
    return GT_TEMPLATE_OK;
  } else {
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE uint64_t gt_template_next_map_pos(gt_template_iterator* const template_iterator) {
  gt_fatal_check(template_iterator==NULL,NULL_HANDLER);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_iterator->template);
  return (template_iterator->next_map_array -
      gt_vector_get_mem(template_iterator->template->maps_relation,gt_map*)) / template_iterator->num_blocks;
}
