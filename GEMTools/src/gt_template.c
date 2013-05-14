/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure modeling sequences' templates.
 *   That is, set of alignments and relationships between their maps.
 */

#include "gt_template.h"
#include "gt_sam_attributes.h"

#define GT_TEMPLATE_TAG_INITIAL_LENGTH 100
#define GT_TEMPLATE_NUM_INITIAL_COUNTERS 10
#define GT_TEMPLATE_NUM_INITIAL_BLOCKS 2
#define GT_TEMPLATE_NUM_INITIAL_MMAPS 20

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new() {
  gt_template* template = gt_alloc(gt_template);
  template->template_id = UINT32_MAX;
  template->in_block_id = UINT32_MAX;
  template->tag = gt_string_new(GT_TEMPLATE_TAG_INITIAL_LENGTH);
  template->alignment_end1=NULL;
  template->alignment_end2=NULL;
  template->counters = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_COUNTERS,sizeof(uint64_t));
  template->mmaps = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_map*));
  template->mmaps_attributes = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_mmap_attributes));
  template->attributes = gt_attributes_new();
  return template;
}
GT_INLINE void gt_template_clear_handler(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_string_clear(template->tag);
  gt_attributes_clear(template->attributes);
}
GT_INLINE void gt_template_clear(gt_template* const template,const bool delete_alignments) {
  GT_TEMPLATE_CHECK(template);
  if (delete_alignments) gt_template_delete_blocks(template);
  gt_vector_clear(template->counters);
  gt_vector_clear(template->mmaps);
  gt_vector_clear(template->mmaps_attributes);
  gt_template_clear_handler(template);
}
GT_INLINE void gt_template_delete(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_string_delete(template->tag);
  gt_template_delete_blocks(template);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->mmaps);
  gt_vector_delete(template->mmaps_attributes);
  gt_attributes_delete(template->attributes);
  gt_free(template);
}

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_tag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_string_get_string(template->tag);
}
GT_INLINE gt_string* gt_template_get_string_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_string_tag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return template->tag;
}
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag,const uint64_t length) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(tag);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_set_tag(alignment,tag,length);
  } GT_TEMPLATE_END_REDUCTION;
  gt_string_set_nstring(template->tag,tag,length);
}
GT_INLINE uint64_t gt_template_get_total_length(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_read_length(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  uint64_t total_length = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    total_length += gt_alignment_get_read_length(alignment);
  }
  return total_length;
}
GT_INLINE int64_t gt_template_get_pair(gt_template* const template){
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_pair(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return *((int64_t*)gt_attributes_get(template->attributes,GT_ATTR_ID_TAG_PAIR));
}
/* Blocks (single alignments) */
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return (gt_expect_true(template->alignment_end1!=NULL) ? (template->alignment_end2!=NULL ? 2 : 1 ) : 0 );
}
GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  if (template->alignment_end1==NULL) {
    template->alignment_end1 = alignment;
  } else if (template->alignment_end2==NULL) {
    template->alignment_end2 = alignment;
  } else {
    gt_fatal_error(TEMPLATE_TOO_MANY_BLOCKS);
  }
}
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  gt_cond_fatal_error(position>1,TEMPLATE_BLOCKS_EXCESS);
  return (position==0 ? template->alignment_end1 : template->alignment_end2);
}
GT_INLINE gt_alignment* gt_template_get_block_dyn(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  if (position==0) {
    if (template->alignment_end1==NULL) template->alignment_end1 = gt_alignment_new();
    return template->alignment_end1;
  } else if (position==1) {
    if (template->alignment_end1==NULL) template->alignment_end1 = gt_alignment_new();
    if (template->alignment_end2==NULL) template->alignment_end2 = gt_alignment_new();
    return template->alignment_end2;
  } else {
    gt_fatal_error(TEMPLATE_BLOCKS_EXCESS);
  }
}
GT_INLINE void gt_template_delete_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  if (template->alignment_end1!=NULL) {
    gt_alignment_delete(template->alignment_end1);
    template->alignment_end1=NULL;
  }
  if (template->alignment_end2!=NULL) {
    gt_alignment_delete(template->alignment_end2);
    template->alignment_end2=NULL;
  }
}
/* SingleEnd/PairedEnd API (just adaptors) */
GT_INLINE void gt_template_set_end1(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  GT_ALIGNMENT_CHECK(alignment);
  template->alignment_end1=alignment;
}
GT_INLINE void gt_template_set_end2(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  GT_ALIGNMENT_CHECK(alignment);
  template->alignment_end2=alignment;
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
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION;
  return *gt_vector_get_elm(template->counters,stratum,uint64_t);
}
GT_INLINE void gt_template_dynamically_allocate_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  const uint64_t used_strata = stratum+1;
  gt_vector_reserve(template->counters,used_strata,true);
  if (gt_vector_get_used(template->counters)<used_strata) {
    gt_vector_set_used(template->counters,used_strata);
  }
}
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_inc_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  ++(*gt_vector_get_elm(template->counters,stratum,uint64_t));
}
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_dec_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  --(*gt_vector_get_elm(template->counters,stratum,uint64_t));
}
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_counter(alignment,stratum,value);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  *gt_vector_get_elm(template->counters,stratum,uint64_t) = value;
}

/*
 * Predefined attributes
 */
GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_mcs(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  uint64_t* mcs = gt_attributes_get(template->attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA);
  if (mcs == NULL) return UINT64_MAX;
  return *mcs;
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_mcs(alignment,max_complete_strata);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_attributes_add(template->attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA,&max_complete_strata,uint64_t);
}
GT_INLINE bool gt_template_has_qualities(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_has_qualities(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    gt_alignment* alignment = gt_template_get_block(template,i);
    if (!gt_alignment_has_qualities(alignment)) return false;
  }
  return true;
}
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_not_unique_flag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  bool* const not_unique_flag = gt_attributes_get(template->attributes,GT_ATTR_ID_NOT_UNIQUE);
  if (not_unique_flag==NULL) return false;
  return *not_unique_flag;
}
GT_INLINE void gt_template_set_not_unique_flag(gt_template* const template,bool is_not_unique) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_not_unique_flag(alignment,is_not_unique);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_attributes_add(template->attributes,GT_ATTR_ID_NOT_UNIQUE,&is_not_unique,bool);
}
GT_INLINE gt_map** gt_template_get_mmap_primary(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  // NOTE: No reduction performed
  return gt_attributes_get(template->attributes,GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT);
}
GT_INLINE void gt_template_set_mmap_primary(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  // NOTE: No reduction performed
  gt_attributes_add(template->attributes,GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT,mmap,gt_map**);
}

/*
 * Multi-maps handlers (Map's relation)
 */
GT_INLINE void gt_template_mmap_attr_new(gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(mmap_attr);
  mmap_attr->distance = 0;
  mmap_attr->gt_score = GT_MAP_NO_GT_SCORE;
  mmap_attr->phred_score = GT_MAP_NO_PHRED_SCORE;
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
  return gt_vector_get_used(template->mmaps)/2;
}
GT_INLINE void gt_template_clear_mmaps(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_clear_maps(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_clear(template->mmaps);
  gt_vector_clear(template->mmaps_attributes);
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
  GT_MMAP_CHECK(mmap);
  gt_vector_insert(template->mmaps,mmap[0],gt_map*);
  gt_vector_insert(template->mmaps,mmap[1],gt_map*);
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE gt_map** gt_template_get_mmap(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_vector_get_elm(alignment->maps,position,gt_map*);
  } GT_TEMPLATE_END_REDUCTION;
  // Retrieve the maps from the mmap vector
  GT_TEMPLATE_CHECK_MMAP_POSITION(template,position);
  if (mmap_attr) *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
  return gt_vector_get_elm(template->mmaps,2*position,gt_map*);
}
GT_INLINE void gt_template_set_mmap(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_map(alignment,*mmap,position);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  GT_MMAP_CHECK(mmap);
  gt_map** const template_mmap = gt_vector_get_elm(template->mmaps,2*position,gt_map*);
  template_mmap[0] = mmap[0];
  template_mmap[1] = mmap[1];
  gt_vector_set_elm(template->mmaps_attributes,position,gt_mmap_attributes,*mmap_attr);
}
/* */
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_check(gt_vector_get_used(mmap)!=1,TEMPLATE_ADD_BAD_NUM_BLOCKS);
    gt_alignment_add_map_gt_vector(alignment,mmap);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Add mmaps
  gt_check(gt_vector_get_used(mmap)!=2,TEMPLATE_ADD_BAD_NUM_BLOCKS);
  gt_map** const _mmap = (gt_map**)((mmap)->memory);
  GT_MMAP_CHECK(_mmap);
  gt_vector_insert(template->mmaps,_mmap[0],gt_map*);
  gt_vector_insert(template->mmaps,_mmap[1],gt_map*);
  gt_vector_insert(template->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
}
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_vector_prepare(mmap,gt_map*,1); // Reset output vector
    gt_vector_insert(mmap,gt_alignment_get_map(alignment,position),gt_map*);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Retrieve the maps from the mmap vector
  GT_TEMPLATE_CHECK_MMAP_POSITION(template,position);
  gt_vector_prepare(mmap,gt_map*,2); // Reset output vector
  gt_map** const template_mmap = gt_vector_get_elm(template->mmaps,2*position,gt_map*);
  gt_vector_insert(mmap,template_mmap[0],gt_map*);
  gt_vector_insert(mmap,template_mmap[1],gt_map*);
  if (mmap_attr) *mmap_attr = *gt_vector_get_elm(template->mmaps_attributes,position,gt_mmap_attributes);
}
/* */
GT_INLINE void gt_template_add_mmap_v(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_map* const map = va_arg(v_args,gt_map*);
    GT_MAP_CHECK(map);
    gt_alignment_add_map(alignment,map);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Add maps
  bool is_mmap_null = true;
  uint64_t i;
  for (i=0;i<2;++i) {
    gt_map* const map = va_arg(v_args,gt_map*);
    is_mmap_null = is_mmap_null && (map==NULL);
    gt_vector_insert(template->mmaps,map,gt_map*);
  }
  gt_cond_fatal_error(is_mmap_null,TEMPLATE_MMAP_NULL);
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
  gt_attributes_copy(template_dst->attributes,template_src->attributes); // Copy templates' attributes
}
GT_INLINE void gt_template_copy_blocks(gt_template* template_dst,gt_template* const template_src,const bool copy_maps) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  gt_template_delete_blocks(template_dst);
  const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    gt_alignment* const alignment_copy = gt_alignment_copy(gt_template_get_block(template_src,i),copy_maps);
    gt_template_add_block(template_dst,alignment_copy);
  }
}
GT_INLINE gt_template* gt_template_copy(gt_template* const template,const bool copy_maps,const bool copy_mmaps) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  gt_template* template_cpy = gt_template_new();
  gt_cond_fatal_error(!template_cpy,MEM_HANDLER);
  // Copy handler
  gt_template_copy_handler(template_cpy,template);
  // Copy blocks
  gt_template_copy_blocks(template_cpy,template,copy_maps);
  // Copy mmaps
  if (copy_maps && copy_mmaps && gt_template_get_num_blocks(template)>1) {
    // Copy counters
    gt_vector_copy(template_cpy->counters,template->counters);
    // Copy mmaps & mmaps_attributes
    GT_TEMPLATE_ITERATE_MMAP__ATTR(template,mmap,mmap_attr) {
      uint64_t map_pos;
      GT_MMAP_ITERATE(mmap,map,end_pos) {
        gt_cond_fatal_error(!gt_alignment_locate_map_reference(
            gt_template_get_block(template,end_pos),map,&map_pos),TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT);
        gt_map* const homologe_map = gt_alignment_get_map(gt_template_get_block(template_cpy,end_pos),map_pos);
        gt_vector_insert(template_cpy->mmaps,homologe_map,gt_map*);
      }
      gt_vector_insert(template_cpy->mmaps_attributes,*mmap_attr,gt_mmap_attributes);
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
  gt_template* const template = template_alignment_iterator->template;
  if (gt_expect_true(template_alignment_iterator->next_pos < gt_template_get_num_blocks(template))) {
    gt_alignment* const alignment = gt_template_get_block(template,template_alignment_iterator->next_pos);
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
  gt_template* const template = template_maps_iterator->template;
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  if (gt_expect_true(template_maps_iterator->next_pos+num_blocks<=gt_vector_get_used(template->mmaps))) {
    *mmap_ref = gt_vector_get_elm(template->mmaps,template_maps_iterator->next_pos,gt_map*);
    if (mmap_attr) {
      const uint64_t attr_position = template_maps_iterator->next_pos/num_blocks;
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
  return template_maps_iterator->next_pos/gt_template_get_num_blocks(template_maps_iterator->template);
}
