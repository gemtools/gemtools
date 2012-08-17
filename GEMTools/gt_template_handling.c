/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_handling.h
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_template_handling.h"

// FIXME: Duplicates function
GT_INLINE void gt_template_insert_mmap(gt_template* const template,gt_map** mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_template_inc_counter(template,mmap_attr->distance);
  gt_template_add_mmap(template,mmap,mmap_attr);
}
GT_INLINE void gt_template_insert_match_gtvector(gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(mmap_attr);
  gt_template_inc_counter(template,mmap_attr->distance);
  gt_template_add_mmap_gtvector(template,mmap,mmap_attr);
}
GT_INLINE void gt_template_insert_mmap_v(gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  gt_template_inc_counter(template,mmap_attr->distance);
  gt_template_add_mmap_v(template,mmap_attr,v_args);
}
GT_INLINE void gt_template_insert_mmap_va(gt_template* const template,gt_mmap_attributes* const mmap_attr,.../* gt_map* maps */) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap_attr);
  va_list v_args;
  va_start(v_args,mmap_attr);
  gt_template_insert_mmap_v(template,mmap_attr,v_args);
}

/*
 * Template's Counters operators
 */
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
  GT_TEMPLATE_IF_REDUCES_ALINGMENT(template) {
    return gt_alignment_is_mapped(GT_TEMPLATE_REDUCED_ALINGMENT(template));
  }
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
 * Template's Maps operators
 */
GT_INLINE gt_map** gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_TEMPLATE_ITERATE(template,map_array) {
    if (gt_mmap_cmp(map_array,mmap,__map_array_num_blocks)==0) return map_array;
  }
  return NULL;
}
GT_INLINE gt_map** gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(gt_mmap_cmp_fx); GT_NULL_CHECK(mmap);
  GT_TEMPLATE_ITERATE(template,map_array) {
    if (gt_mmap_cmp_fx(map_array,mmap,__map_array_num_blocks)==0) return map_array;
  }
  return NULL;
}
GT_INLINE void gt_template_add_template_mmaps(gt_template* const template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  register uint64_t position = 0;
  _GT_TEMPLATE_ITERATE(template_src,map_array) {
    gt_template_add_mmap(template_dst,map_array,
        gt_template_get_mmap_attr(template_src,position));
    ++position;
  }
}
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  register uint64_t position = 0;
  _GT_TEMPLATE_ITERATE(template_src,map_array) {
    if (gt_template_is_mmap_contained(template_dst,map_array)==NULL) {
      gt_template_add_mmap(template_dst,map_array,
          gt_template_get_mmap_attr(template_src,position));
    }
    ++position;
  }
}
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  register uint64_t position = 0;
  _GT_TEMPLATE_ITERATE(template_src,map_array) {
    if (gt_template_is_mmap_contained_fx(gt_mmap_cmp_fx,template_dst,map_array)==NULL) {
      gt_template_add_mmap(template_dst,map_array,  // FIXME Insert, not add
          gt_template_get_mmap_attr(template_src,position));
    }
    ++position;
  }
}


/*
 * Template's Maps set-operators
 */
GT_INLINE void gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  // TODO
}
GT_INLINE void gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_union_template_mmaps_fx_v(gt_mmap_cmp_fx,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_template_union_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_union_template_mmaps_fx_v(gt_mmap_cmp,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}

GT_INLINE void gt_template_intersect_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  register uint64_t i;
  for (i=0;i<num_src_templates;++i) {
    register gt_template* aux_template = (i==0) ? template_src : va_arg(v_args,gt_template*);
    register uint64_t position = 0;
    _GT_TEMPLATE_ITERATE(aux_template,mmap) {
      if (gt_template_is_mmap_contained_fx(gt_mmap_cmp_fx,template_base,mmap)!=NULL &&
          gt_template_is_mmap_contained_fx(gt_mmap_cmp_fx,template_dst,mmap)==NULL) {
        gt_template_add_mmap(template_dst,mmap,
            gt_template_get_mmap_attr(aux_template,position)); // FIXME
      }
      ++position;
    }
  }
}
GT_INLINE void gt_template_intersect_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_intersect_template_mmaps_fx_v(gt_mmap_cmp_fx,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_template_intersect_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_intersect_template_mmaps_fx_v(gt_mmap_cmp,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}

GT_INLINE void gt_template_subtract_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  // Store source templates from the arguments
  gt_template** aux_template = malloc(sizeof(gt_template*)*num_src_templates);
  register uint64_t i;
  for (i=0;i<num_src_templates;++i) {
    aux_template[i] = (i==0) ? template_src : va_arg(v_args,gt_template*);
  }
  register uint64_t position = 0;
  _GT_TEMPLATE_ITERATE(template_base,mmap) {
    register bool subtracted = false;
    for (i=0;i<num_src_templates;++i) {
      if (gt_template_is_mmap_contained_fx(gt_mmap_cmp_fx,aux_template[i],mmap)!=NULL) {
        subtracted = true; break;
      }
    }
    if (!subtracted) {
      gt_template_add_mmap(template_dst,mmap,
          gt_template_get_mmap_attr(template_base,position)); // FIXME
    }
    ++position;
  }
  free(aux_template);
}
GT_INLINE void gt_template_subtract_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_subtract_template_mmaps_fx_v(gt_mmap_cmp_fx,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_template_subtract_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(template_dst); GT_NULL_CHECK(template_base);
  GT_ZERO_CHECK(num_src_templates); GT_NULL_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template_subtract_template_mmaps_fx_v(gt_mmap_cmp,template_dst,
      template_base,num_src_templates,template_src,v_args);
  va_end(v_args);
}
