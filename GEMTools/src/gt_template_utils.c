/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_utils.h
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_template_utils.h"

/*
 * Template useful functions
 */
GT_INLINE void gt_template_deduce_alignments_tags(gt_template* const template) {
  /*
   * Tag Block Conventions::
   *  (1) line1: NAME => line1&2: NAME
   *      line2: NAME
   *  (2) line1: NAME<SEP>SOMETHING_ELSE_1 => line1&2: NAME<SEP>SOMETHING_ELSE_1|NAME<SEP>SOMETHING_ELSE_2
   *      line2: NAME<SEP>SOMETHING_ELSE_2
   *  (3) line1: NAME<SEP>1 => line1&2: NAME
   *      line2: NAME<SEP>2
   * Tag Decomposition::
   *  (<TAG_BLOCK><SEP>){num_blocks}
   *    -> tag.num_blocks == read.num_blocks
   *    -> <SEP> := ' ' | '/' | '#' | "|"  // TODO
   */
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    register gt_alignment* const alignment = gt_template_get_block(template,i);
    gt_sprintf(alignment->tag,"%s/%"PRIu64,gt_string_get_string(template->tag),i+1);
  }
}
GT_INLINE void gt_template_deduce_template_tag(gt_template* const template,gt_alignment* const alignment) {
  register uint64_t tag_length = gt_alignment_get_tag_length(alignment);
  register gt_string* const tag = alignment->tag;
  if (tag_length>2 && *gt_string_char_at(tag,tag_length-2)==SLASH) {
    tag_length-=2;
  }
  if (tag_length>0) {
    gt_string_set_nstring(template->tag,gt_string_get_string(tag),tag_length);
  }
}

/*
 * Template's MMaps operators (Update global state: counters, ...)
 */
GT_INLINE void gt_template_alias_dup_mmap_members(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_map** const uniq_mmaps) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(uniq_mmaps);
  // Resolve mmap
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    GT_NULL_CHECK(mmap[i]);
    uniq_mmaps[i] = gt_alignment_put_map(gt_map_cmp_fx,gt_template_get_block(template,i),mmap[i],false);
  }
}
GT_INLINE gt_map** gt_template_raw_put_mmap(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  // Resolve mmap aliasing/insertion
  register gt_map** uniq_mmaps = malloc(gt_template_get_num_blocks(template)*sizeof(gt_map*));
  gt_template_alias_dup_mmap_members(gt_map_cmp_fx,template,mmap,uniq_mmaps);
  // Raw mmap insertion
  gt_template_add_mmap(template,uniq_mmaps,mmap_attr);
  register gt_map** template_mmap=gt_template_get_mmap(template,gt_template_get_num_mmaps(template)-1,NULL);
  free(uniq_mmaps); // Free auxiliary vector
  return template_mmap;
}
GT_INLINE gt_map** gt_template_put_mmap(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,const bool replace_duplicated) {
  GT_NULL_CHECK(gt_mmap_cmp_fx); GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  // Check mmap duplicates
  gt_map** found_mmap;
  gt_mmap_attributes found_mmap_attr={0,0};
  uint64_t found_mmap_pos=0;
  register bool is_duplicated = gt_expect_false(gt_template_find_mmap_fx(gt_mmap_cmp_fx,
      template,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attr));
  register gt_map** template_mmap;
  if (!is_duplicated || replace_duplicated) { // FIXME: Chose which to replace (like alignment)
    // Resolve mmap aliasing/insertion
    register gt_map** uniq_mmaps = malloc(gt_template_get_num_blocks(template)*sizeof(gt_map*));
    gt_template_alias_dup_mmap_members(gt_map_cmp_fx,template,mmap,uniq_mmaps);
    // Insert mmap
    if (!is_duplicated) { // Add new mmap
      gt_template_inc_counter(template,mmap_attr->distance);
      gt_template_add_mmap(template,uniq_mmaps,mmap_attr);
      template_mmap=gt_template_get_mmap(template,gt_template_get_num_mmaps(template)-1,NULL);
    } else { // Replace mmap
      gt_template_dec_counter(template,found_mmap_attr.distance); // Remove old mmap
      gt_template_set_mmap(template,found_mmap_pos,uniq_mmaps,mmap_attr); // Replace old mmap
      template_mmap=gt_template_get_mmap(template,found_mmap_pos,NULL);
    }
    free(uniq_mmaps); // Free auxiliary vector
  } else {
    // Delete mmap
    GT_MULTIMAP_ITERATE_BLOCKS(mmap,gt_template_get_num_blocks(template),map,end_pos) {
      gt_map_delete(map);
    }
    template_mmap=gt_template_get_mmap(template,found_mmap_pos,NULL);
  }
  return template_mmap;
}

GT_INLINE void gt_template_insert_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_template_insert_mmap_fx(gt_mmap_cmp,template,mmap,mmap_attr);
}
GT_INLINE void gt_template_insert_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp,template,mmap,mmap_attr,true); // TODO: Why replace, why?
}
GT_INLINE void gt_template_insert_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_template_insert_mmap_gtvector_fx(gt_mmap_cmp,template,mmap,mmap_attr);
}
GT_INLINE void gt_template_insert_mmap_gtvector_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_VECTOR_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_check(gt_vector_get_used(mmap)!=gt_template_get_num_blocks(template),TEMPLATE_ADD_BAD_NUM_BLOCKS);
  gt_template_insert_mmap_fx(gt_mmap_cmp_fx,template,gt_vector_get_mem(mmap,gt_map*),mmap_attr);
}

//// TODO: Scheduled for v2.0
//GT_INLINE void gt_template_remove_mmap(
//    gt_template* const template,gt_map** const mmap) {
//  // TODO
//}
//GT_INLINE void gt_template_remove_mmap_fx(
//    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
//    gt_template* const template,gt_map** const mmap) {
//  // TODO
//}


GT_INLINE bool gt_template_find_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,
    uint64_t* const found_mmap_pos,gt_map*** const found_mmap,gt_mmap_attributes* const found_mmap_attr) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template); GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(found_mmap_pos); GT_NULL_CHECK(found_mmap);
  // Search for the mmap
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t pos = 0;
  GT_TEMPLATE_ITERATE_(template,template_mmap) {
    if (gt_mmap_cmp_fx(template_mmap,mmap,num_blocks)==0) {
      *found_mmap_pos = pos;
      *found_mmap = template_mmap;
      if (found_mmap_attr) *found_mmap_attr = *gt_template_get_mmap_attr(template,pos);
      return true;
    }
    ++pos;
  }
  return false;
}
GT_INLINE bool gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template); GT_NULL_CHECK(mmap);
  gt_map** found_mmap;
  uint64_t found_mmap_pos;
  return gt_template_find_mmap_fx(gt_mmap_cmp,template,mmap,
      &found_mmap_pos,&found_mmap,NULL);
}
GT_INLINE bool gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CONSISTENCY_CHECK(template); GT_NULL_CHECK(mmap);
  gt_map** found_mmap;
  uint64_t found_mmap_pos;
  return gt_template_find_mmap_fx(gt_mmap_cmp_fx,template,mmap,
      &found_mmap_pos,&found_mmap,NULL);
}

/* Insert size should be and estimate of the original fragment size that was sequenced.  We can get this from the following:
 *
 * position of rightmost read + no. bases in rightmost read - position of leftmost read
 *
 * If the read has been trimmed from the start of the read then we can't get the original size, but this is a relatively
 * rare event.  Trimming from the end of the read does not effect the calculation as the rightmost read will be on the
 * negative strand so trimming x bases from the end of the read will shift the mapped position of the read by x.
 *
 * Split mappings require special handling in that we need to consider the number of bases read + the distance between the
 * last block in each mapping as follows:
 *
 * Position of last block of rightmost read + no. bases in rightmost read - (position of last block of leftmost read + no. bases)
 * in all other blocks of leftmost read
 */
GT_INLINE int64_t gt_template_get_insert_size(gt_map** const mmap,uint64_t *gt_error) {
  // Get last block of each map
  gt_map *block[2]={0,0};
  uint64_t length[2]={0,0};
  int64_t x=0;
  *gt_error=GT_TEMPLATE_INSERT_SIZE_OK;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(mmap[0],map_it) {
    block[0]=map_it;
    length[0]+=gt_map_get_base_length(block[0]);
  } GT_END_MAP_BLOCKS_ITERATOR;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(mmap[1],map_it) {
    block[1]=map_it;
    length[1]+=gt_map_get_base_length(block[1]);
  } GT_END_MAP_BLOCKS_ITERATOR;
  if(gt_string_equals(block[0]->seq_name,block[1]->seq_name)) {
    if(block[0]->strand!=block[1]->strand) {
      if(block[0]->strand==FORWARD) x=1+block[1]->position+length[1]-(block[0]->position+length[0]-gt_map_get_base_length(block[0]));
      else x=1+block[0]->position+length[0]-(block[1]->position+length[1]-gt_map_get_base_length(block[1]));
    } else {
      *gt_error=GT_TEMPLATE_INSERT_SIZE_SAME_STRAND;
    }
  } else {
    *gt_error=GT_TEMPLATE_INSERT_SIZE_DIFFERENT_CONTIGS;
  }
  /*  if (gt_map_get_position(mmap[1])>gt_map_get_position(mmap[0])) {
    return ((int64_t)gt_map_get_position(mmap[1])-(int64_t)(gt_map_get_position(mmap[0])+gt_map_get_length(mmap[0])));
  } else {
    return ((int64_t)gt_map_get_position(mmap[0])-(int64_t)(gt_map_get_position(mmap[1])+gt_map_get_length(mmap[1])));
    } */
  return x;
}

/*
 * Template's Counters operators
 */
GT_INLINE bool gt_template_is_mapped(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  return gt_template_is_thresholded_mapped(template,UINT64_MAX);
}
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_is_thresholded_mapped(alignment,max_allowed_strata);
  } GT_TEMPLATE_END_REDUCTION;
  register gt_vector* vector = gt_template_get_counters_vector(template);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if ((counter_pos+1)>=max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE void gt_template_recalculate_counters(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_recalculate_counters(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Clear previous counters
  gt_vector_clear(gt_template_get_counters_vector(template));
  // Recalculate counters
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  GT_TEMPLATE__ATTR_ITERATE_(template,mmap,mmap_attr) {
    register uint64_t i, total_distance = 0;
    for (i=0;i<num_blocks;++i) {
      total_distance+=gt_map_get_global_distance(mmap[i]);
    }
    mmap_attr->distance=total_distance;
    gt_template_inc_counter(template,total_distance);
  }
}

GT_INLINE int64_t gt_template_get_min_matching_strata(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_min_matching_strata(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_min_matching_strata(gt_template_get_counters_vector(template));
}
GT_INLINE int64_t gt_template_get_uniq_degree(gt_template* const template) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_uniq_degree(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_uniq_degree(gt_template_get_counters_vector(template));
}
GT_INLINE bool gt_template_get_next_matching_strata(
    gt_template* const template,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_next_matching_strata(alignment,begin_strata,next_matching_strata,num_maps);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_next_matching_strata(gt_template_get_counters_vector(template),
      begin_strata,next_matching_strata,num_maps);
}

/*
 * Template Set operators
 */
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_dst);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  gt_template_merge_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_dst,template_src);
}
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_dst,gt_template* const template_src) {
  GT_NULL_CHECK(gt_mmap_cmp_fx); GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_dst,template_src);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_dst,alignment_dst) {
      gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_dst,alignment_src);
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Merge mmaps
  GT_TEMPLATE_CONSISTENCY_CHECK(template_dst);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  GT_TEMPLATE__ATTR_ITERATE(template_src,mmap,mmap_attr) {
    register gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
    gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_dst,mmap_copy,mmap_attr,true);
    free(mmap_copy); // Free array handler
  }
  gt_template_set_mcs(template_dst,GT_MIN(gt_template_get_mcs(template_dst),gt_template_get_mcs(template_src)));
}

//// TODO: Scheduled for v2.0
//GT_INLINE void gt_template_remove_template_mmaps(
//    gt_template* const template_dst,gt_template* const template_src) {
//  GT_TEMPLATE_CONSISTENCY_CHECK(template_dst);
//  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
//  GT_TEMPLATE_ITERATE_(template_src,mmap) {
//    gt_template_remove_mmap(template_dst,mmap);
//  }
//}
//GT_INLINE void gt_template_remove_template_mmaps_fx(
//    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
//    gt_template* const template_dst,gt_template* const template_src) {
//  GT_NULL_CHECK(gt_mmap_cmp_fx);
//  GT_TEMPLATE_CONSISTENCY_CHECK(template_dst);
//  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
//  GT_TEMPLATE_ITERATE_(template_src,mmap) {
//    gt_template_remove_mmap_fx(gt_mmap_cmp_fx,template_dst,mmap);
//  }
//}

GT_INLINE gt_template* gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_ZERO_CHECK(num_src_templates);
  // Create new template
  register gt_template* const template_union = gt_template_copy(template_src,false,false);
  gt_template_merge_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_union,template_src);
  // Merge template sources into template_union
  register uint64_t num_tmp_merged = 1;
  while (num_tmp_merged < num_src_templates) {
    register gt_template* template_target = va_arg(v_args,gt_template*);
    GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_union,template_target);
    GT_TEMPLATE_CONSISTENCY_CHECK(template_target);
    gt_template_merge_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_union,template_target);
    ++num_tmp_merged;
  }
  return template_union;
}
GT_INLINE gt_template* gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_ZERO_CHECK(num_src_templates);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  register gt_template* const template_union =
      gt_template_union_template_mmaps_fx_v(gt_mmap_cmp_fx,gt_map_cmp_fx,num_src_templates,template_src,v_args);
  va_end(v_args);
  return template_union;
}
GT_INLINE gt_template* gt_template_union_template_mmaps_va(
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_ZERO_CHECK(num_src_templates);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  register gt_template* const template_union =
      gt_template_union_template_mmaps_fx_v(gt_mmap_cmp,gt_map_cmp,num_src_templates,template_src,v_args);
  va_end(v_args);
  return template_union;
}

GT_INLINE gt_template* gt_template_subtract_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_minuend,gt_template* const template_subtrahend) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_minuend,template_subtrahend);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_minuend,alignment_minuend) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_subtrahend,alignment_subtrahend) {
      register gt_alignment* const alignment_difference =
          gt_alignment_subtract_alignment_maps_fx(gt_map_cmp_fx,alignment_minuend,alignment_subtrahend);
      register gt_template* const template_difference = gt_template_new();
      gt_template_add_block(template_difference,alignment_difference);
      return template_difference;
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION;
  // Subtract
  GT_TEMPLATE_CONSISTENCY_CHECK(template_minuend);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_subtrahend);
  register gt_template* const template_difference = gt_template_copy(template_minuend,false,false);
  uint64_t found_mmap_pos;
  gt_map** found_mmap;
  gt_mmap_attributes found_mmap_attr;
  GT_TEMPLATE__ATTR_ITERATE(template_minuend,mmap,mmap_attr) {
    if (!gt_template_find_mmap_fx(gt_mmap_cmp_fx,template_subtrahend,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attr)) {
      register gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
      gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_difference,mmap_copy,mmap_attr,false);
      free(mmap_copy);
    }
  }
  return template_difference;
}
GT_INLINE gt_template* gt_template_subtract_template_mmaps(
    gt_template* const template_minuend,gt_template* const template_subtrahend) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_minuend);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_subtrahend);
  return gt_template_subtract_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_minuend,template_subtrahend);
}

GT_INLINE gt_template* gt_template_intersect_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_A,gt_template* const template_B) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_A,template_B);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_A,alignment_A) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_B,alignment_B) {
      register gt_alignment* const alignment_intersection =
          gt_alignment_intersect_alignment_maps_fx(gt_map_cmp_fx,alignment_A,alignment_B);
      register gt_template* const template_intersection = gt_template_new();
      gt_template_add_block(template_intersection,alignment_intersection);
      return template_intersection;
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION;
  // Intersect
  GT_TEMPLATE_CONSISTENCY_CHECK(template_A);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_B);
  register gt_template* const template_intersection = gt_template_copy(template_A,false,false);
  uint64_t found_mmap_pos;
  gt_map** found_mmap;
  gt_mmap_attributes found_mmap_attr;
  GT_TEMPLATE__ATTR_ITERATE(template_A,mmap,mmap_attr) {
    if (gt_template_find_mmap_fx(gt_mmap_cmp_fx,template_B,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attr)) {
      register gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
      gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_intersection,mmap_copy,mmap_attr,false);
      free(mmap_copy);
    }
  }
  return template_intersection;
}
GT_INLINE gt_template* gt_template_intersect_template_mmaps(
    gt_template* const template_A,gt_template* const template_B) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template_A);
  GT_TEMPLATE_CONSISTENCY_CHECK(template_B);
  return gt_template_intersect_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_A,template_B);
}

/*
 * Template realignment
 */
GT_INLINE void gt_template_realign_hamming(gt_template* const template,gt_sequence_archive* const sequence_archive) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    gt_alignment_realign_hamming(alignment,sequence_archive);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_realign_levenshtein(gt_template* const template,gt_sequence_archive* const sequence_archive) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    gt_alignment_realign_levenshtein(alignment,sequence_archive);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_realign_weighted(
    gt_template* const template,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_NULL_CHECK(gt_weigh_fx);
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    gt_alignment_realign_weighted(alignment,sequence_archive,gt_weigh_fx);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}






