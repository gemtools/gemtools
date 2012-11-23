/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * DESCRIPTION: Data structure modeling sequences' templates.
 *   That is, set of alignments and relationships between their maps.
 */

#include "gt_template.h"

#define GT_TEMPLATE_TAG_INITIAL_LENGTH 100
#define GT_TEMPLATE_NUM_INITIAL_COUNTERS 10
#define GT_TEMPLATE_NUM_INITIAL_BLOCKS 2
#define GT_TEMPLATE_NUM_INITIAL_MMAPS 20

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new() {
  gt_template* template = malloc(sizeof(gt_template));
  gt_cond_fatal_error(!template,MEM_HANDLER);
  template->template_id = UINT32_MAX;
  template->in_block_id = UINT32_MAX;
  template->tag = gt_string_new(GT_TEMPLATE_TAG_INITIAL_LENGTH);
  template->blocks = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_BLOCKS,sizeof(gt_alignment*));
  template->counters = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_COUNTERS,sizeof(uint64_t));
  template->mmaps = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_map*));
  template->mmaps_attributes = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_mmap_attributes));
  template->maps_txt = NULL;
  template->attributes = gt_shash_new();
  return template;
}
GT_INLINE void gt_template_clear_handler(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_string_clear(template->tag);
  template->maps_txt = NULL;
  gt_shash_clean(template->attributes,true,true);
}
GT_INLINE void gt_template_clear(gt_template* const template,const bool delete_alignments,const bool delete_maps) {
  GT_TEMPLATE_CHECK(template);
  if (delete_alignments) gt_template_delete_blocks(template,delete_maps);
  gt_vector_clean(template->counters);
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
  gt_template_clear_handler(template);
}
GT_INLINE void gt_template_delete(gt_template* const template,const bool delete_alignments,const bool delete_maps) {
  GT_TEMPLATE_CHECK(template);
  gt_string_delete(template->tag);
  if (delete_alignments) gt_template_delete_blocks(template,delete_maps);
  gt_vector_delete(template->blocks);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->mmaps);
  gt_vector_delete(template->mmaps_attributes);
  gt_shash_delete(template->attributes,true,true);
  free(template);
}

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return gt_string_get_string(template->tag);
}
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag,const uint64_t length) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(tag);
  gt_string_set_nstring(template->tag,tag,length);
}
GT_INLINE uint64_t gt_template_get_total_length(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_read_length(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  register uint64_t total_length = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    total_length += gt_alignment_get_read_length(alignment);
  }
  return total_length;
}
/* Blocks (single alignments) */
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return gt_vector_get_used(template->blocks);
}
GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  gt_vector_insert(template->blocks,alignment,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  return *gt_vector_get_elm(template->blocks,position,gt_alignment*);
}
GT_INLINE gt_alignment* gt_template_get_block_dyn(gt_template* const template,const uint64_t position) {
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
  gt_vector_clean(template->blocks);
}
GT_INLINE void gt_template_delete_blocks(gt_template* const template,const bool delete_maps) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_ITERATE(template->blocks,alignment_it,alignment_pos,gt_alignment*) {
    gt_alignment_delete(*alignment_it,delete_maps);
  }
  gt_vector_clean(template->blocks);
}
/* Counters */
GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_counters_vector(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return template->counters;
}
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(counters);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_counters_vector(alignment,counters);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  template->counters = counters;
}
GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_num_counters(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_vector_get_used(template->counters);
}
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION;
  return *gt_vector_get_elm(template->counters,stratum-1,uint64_t);
}
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_inc_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_reserve(template->counters,stratum,true);
  ++(*gt_vector_get_elm(template->counters,stratum-1,uint64_t));
}
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_dec_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_reserve(template->counters,stratum,true);
  --(*gt_vector_get_elm(template->counters,stratum-1,uint64_t));
}
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value) {
  GT_TEMPLATE_CHECK(template);
  gt_fatal_check(stratum==0,COUNTERS_POS_STRATUM);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_counter(alignment,stratum,value);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_reserve(template->counters,stratum,true);
  *gt_vector_get_elm(template->counters,stratum-1,uint64_t) = value;
}

/*
 * Attribute accessors
 */
GT_INLINE void* gt_template_get_attr(
    gt_template* const template,char* const attribute_id) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(attribute_id);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_attr(alignment,attribute_id);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_attribute_get(template->attributes,attribute_id);
}
GT_INLINE void gt_template_set_attr(
    gt_template* const template,char* const attribute_id,
    void* const attribute,const size_t element_size) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_ZERO_CHECK(element_size);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_attr(alignment,attribute_id,attribute,element_size);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Insert attribute
  gt_attribute_set(template->attributes,attribute_id,attribute,element_size);
}
/*
 * Predefined attributes
 */
GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template) { // FIXME: reduce
  GT_TEMPLATE_CHECK(template);
  register uint64_t* mcs = gt_template_get_attribute(template,GT_ATTR_MAX_COMPLETE_STRATA,uint64_t);
  if (mcs == NULL) return UINT64_MAX;
  return *mcs;
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  gt_template_set_attribute(template,GT_ATTR_MAX_COMPLETE_STRATA,&max_complete_strata,uint64_t);
}
GT_INLINE bool gt_template_has_qualities(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    register gt_alignment* alignment = gt_template_get_block(template,i);
    if (!gt_alignment_has_qualities(alignment)) return false;
  }
  return true;
}
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  register bool* const not_unique_flag = gt_template_get_attribute(template,GT_ATTR_NOT_UNIQUE,bool);
  if (not_unique_flag==NULL) return false;
  return *not_unique_flag;
}
GT_INLINE void gt_template_set_not_unique_flag(gt_template* const template,bool is_not_unique) {
  GT_TEMPLATE_CHECK(template);
  gt_template_set_attribute(template,GT_ATTR_NOT_UNIQUE,&is_not_unique,bool);
}

/*
 * Multi-maps handlers (Map's relation)
 */
GT_INLINE void gt_template_mmap_attr_new(gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(mmap_attr);
  mmap_attr->score = 0;
  mmap_attr->distance = 0;
}
GT_INLINE gt_mmap_attributes* gt_template_get_mmap_attr(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  return gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
}
GT_INLINE void gt_template_set_mmap_attr(gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes) = *mmap_attr;
}
/* */
GT_INLINE uint64_t gt_template_get_num_mmaps(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_num_maps(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_vector_get_used(template->mmaps)/gt_template_get_num_blocks(template);
}
GT_INLINE void gt_template_clear_mmaps(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_clear_maps(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_clean(template->mmaps);
  gt_vector_clean(template->mmaps_attributes);
}
/* */
GT_INLINE void gt_template_add_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_add_map(alignment,*mmap);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
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
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_vector_get_elm(alignment->maps,position,gt_map*);
  } GT_TEMPLATE_END_REDUCTION;
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
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_map(alignment,*mmap,position);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register gt_map** const template_mmap = gt_vector_get_elm(template->mmaps,num_blocks*position,gt_map*);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    GT_MAP_CHECK(mmap[i]);
    template_mmap[i] = mmap[i];
  }
  gt_vector_set_elm(template->mmaps_attributes,position,gt_mmap_attributes,*mmap_attr);
}
/* */
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_check(gt_template_get_num_blocks(template)!=gt_vector_get_used(mmap),TEMPLATE_ADD_BAD_NUM_BLOCKS);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_add_map_gt_vector(alignment,mmap);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Add mmaps
  GT_VECTOR_ITERATE(mmap,map_ptr,map_pos,gt_map*) {
    register gt_map* const map = *map_ptr;
    GT_MAP_CHECK(map);
    gt_vector_insert(template->mmaps,map,gt_map*);
  }
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  gt_fatal_check(position>=(gt_vector_get_used(template->mmaps)/num_blocks),POSITION_OUT_OF_RANGE);
  // Reset output vector
  gt_vector_prepare(mmap,gt_map*,num_blocks);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_vector_insert(mmap,gt_alignment_get_map(alignment,position),gt_map*);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Retrieve the maps from the mmap vector
  register const uint64_t init_map_pos = num_blocks*position;
  register gt_map** template_mmap = gt_vector_get_elm(template->mmaps,init_map_pos,gt_map*);
  register uint64_t i;
  for(i=0;i<num_blocks;i++) {
    gt_vector_insert(mmap,template_mmap[i],gt_map*);
  }
  if (mmap_attr) *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
}
/* */
GT_INLINE void gt_template_add_mmap_v(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    register gt_map* const map = va_arg(v_args,gt_map*);
    GT_MAP_CHECK(map);
    gt_alignment_add_map(alignment,map);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    register gt_map* const map = va_arg(v_args,gt_map*);
    GT_MAP_CHECK(map);
    gt_vector_insert(template->mmaps,map,gt_map*);
  }
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
/* */
GT_INLINE void gt_template_add_mmap_va(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,...) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  va_list v_args;
  va_start(v_args,mmap_attr);
  gt_template_add_mmap_v(template,mmap_attr,v_args);
  va_end(v_args);
}

/*
 * Miscellaneous
 */
GT_INLINE void gt_template_copy_handler(gt_template* template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  template_dst->template_id = template_src->template_id;
  template_dst->in_block_id = template_src->in_block_id;
  gt_string_copy(template_dst->tag,template_src->tag);
  template_dst->maps_txt = template_src->maps_txt;
  // Copy templates' attributes
  template_dst->attributes = gt_shash_deep_copy(template_src->attributes);
}
GT_INLINE void gt_template_copy_blocks(gt_template* template_dst,gt_template* const template_src,const bool copy_maps,const bool deep_copy_maps) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  gt_vector_clean(template_dst->blocks);
  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    register gt_alignment* alg_copy =
        gt_alignment_copy(gt_template_get_block(template_src,i),copy_maps,deep_copy_maps);
    gt_template_add_block(template_dst,alg_copy);
  }
}
GT_INLINE gt_template* gt_template_copy(
    gt_template* const template,const bool copy_blocks,const bool deep_copy_blocks,
    const bool deep_copy_maps,const bool copy_mmaps) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_template* template_cpy = gt_template_new();
  gt_cond_fatal_error(!template_cpy,MEM_HANDLER);
  // Copy handler
  gt_template_copy_handler(template_cpy,template);
  if (!copy_blocks) return template_cpy;
  // Copy blocks
  if (deep_copy_blocks) {
    gt_template_copy_blocks(template_cpy,template,true,deep_copy_maps);
  } else {
    gt_vector_copy(template_cpy->blocks,template->blocks);
  }
  // Copy mmaps
  if (copy_mmaps) {
    // Copy counters
    gt_vector_copy(template_cpy->counters,template->counters);
    if (deep_copy_maps) {
      // Copy mmaps & mmaps_attributes
      gt_map** mmap;
      gt_mmap_attributes* mmap_attr = NULL;
      gt_template_maps_iterator template_maps_iterator;
      gt_template_new_mmap_iterator(template,&template_maps_iterator);
      register const uint64_t num_blocks = gt_template_get_num_blocks(template);
      while (gt_expect_true(gt_template_next_mmap(&template_maps_iterator,&mmap,&mmap_attr))) {
        uint64_t i, map_pos;
        for (i=0;i<num_blocks;++i) {
          gt_cond_fatal_error(!gt_alignment_locate_map_reference(
              gt_template_get_block(template,i),mmap[i],&map_pos),TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT);
          register gt_map* homologe_map = gt_alignment_get_map(gt_template_get_block(template_cpy,i),map_pos);
          gt_vector_insert(template->mmaps,homologe_map,gt_map*);
        }
        gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
      }
    } else {
      gt_vector_copy(template_cpy->mmaps,template->mmaps);
      gt_vector_copy(template_cpy->mmaps_attributes,template->mmaps_attributes);
    }
  }
  return template_cpy;
}

/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(
    gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(template_alignment_iterator);
  template_alignment_iterator->template = template;
  template_alignment_iterator->next_pos = 0;
}
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_alignment_iterator->template);
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
  GT_TEMPLATE_CONSISTENCY_CHECK(template_alignment_iterator->template);
  return template_alignment_iterator->next_pos;
}

/*
 * Template's Maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
 */
GT_INLINE void gt_template_new_mmap_iterator(
    gt_template* const template,gt_template_maps_iterator* const template_maps_iterator) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(template_maps_iterator);
  template_maps_iterator->template = template;
  template_maps_iterator->next_pos = 0;
}
GT_INLINE gt_status gt_template_next_mmap(
    gt_template_maps_iterator* const template_maps_iterator,
    gt_map*** const mmap_ref,gt_mmap_attributes** const mmap_attr) {
  GT_NULL_CHECK(template_maps_iterator); GT_NULL_CHECK(mmap_ref);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_maps_iterator->template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_maps_iterator->template,alignment) {
    if (gt_expect_true(template_maps_iterator->next_pos<gt_alignment_get_num_maps(alignment))) {
      *mmap_ref = gt_vector_get_elm(alignment->maps,template_maps_iterator->next_pos,gt_map*);
      template_maps_iterator->next_pos++;
      return GT_TEMPLATE_OK;
    } else {
      *mmap_ref = NULL;
      return GT_TEMPLATE_FAIL;
    }
  } GT_TEMPLATE_END_REDUCTION;
  register gt_template* const template = template_maps_iterator->template;
  register const uint64_t num_blocks = gt_vector_get_used(template->blocks);
  if (gt_expect_true(template_maps_iterator->next_pos+num_blocks<=gt_vector_get_used(template->mmaps))) {
    *mmap_ref = gt_vector_get_elm(template->mmaps,template_maps_iterator->next_pos,gt_map*);
    if (mmap_attr) {
      register const uint64_t attr_position = template_maps_iterator->next_pos/num_blocks;
      *mmap_attr = gt_template_get_mmap_attr(template,attr_position);
    }
    template_maps_iterator->next_pos+=num_blocks;
    return GT_TEMPLATE_OK;
  } else {
    *mmap_ref = NULL;
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE uint64_t gt_template_next_mmap_pos(gt_template_maps_iterator* const template_maps_iterator) {
  GT_NULL_CHECK(template_maps_iterator);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_maps_iterator->template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT_(template_maps_iterator->template) {
    return template_maps_iterator->next_pos;
  } GT_TEMPLATE_END_REDUCTION;
  return template_maps_iterator->next_pos/gt_vector_get_used(template_maps_iterator->template->blocks);
}
