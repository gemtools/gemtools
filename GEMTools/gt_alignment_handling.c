/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_handling.c
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_alignment_handling.h"
#include "gt_iterators.h"

/*
 * Alignment high-level insert (Update global state: counters, ...)
 */
GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment); // FIXME: Avoid duplications => Solve Duplication
  GT_MAP_CHECK(map);
  gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map)+1);
  gt_alignment_add_map(alignment,map);
}
GT_INLINE void gt_alignment_insert_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map,map_pos,gt_map*) {
    GT_MAP_CHECK(*map);
    gt_alignment_insert_map(alignment,*map);
  }
}

/*
 * Alignment' Counters operators
 */
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment) {
  // TODO
}
GT_INLINE uint64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (*counter!=0) return counter_pos+1;
  }
  return UINT64_MAX;
}
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const uint64_t max_allowed_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if ((counter_pos+1)>=max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment); // FIXME: Check unique flag
  return gt_alignment_is_thresholded_mapped(alignment,UINT64_MAX);
}

/*
 * Alignment' Maps operators
 */
GT_INLINE gt_map* gt_alignment_is_map_contained(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_MAP_CHECK(map);
  GT_MAPS_ITERATE(alignment,map_it) {
   if (gt_map_cmp(map_it,map)==0) return map_it;
  }
  return NULL;
}
GT_INLINE gt_map* gt_alignment_is_map_contained_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_MAP_CHECK(map);
  GT_MAPS_ITERATE(alignment,map_it) {
   if (gt_map_cmp_fx(map_it,map)==0) return map_it;
  }
  return NULL;
}
GT_INLINE void gt_alignment_add_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_MAPS_ITERATE(alignment_src,map_src) {
    gt_alignment_insert_map(alignment_dst,map_src);
  }
}
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_MAPS_ITERATE(alignment_src,map_src) {
    if (gt_alignment_is_map_contained(alignment_dst,map_src)==NULL) {
      gt_alignment_insert_map(alignment_dst,map_src);
    }
  }
}
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_MAPS_ITERATE(alignment_src,map_src) {
    if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_dst,map_src)==NULL) {
      gt_alignment_insert_map(alignment_dst,map_src);
    }
  }
}

/*
 * Alignment' Maps set-operators
 */
GT_INLINE void gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  // TODO
}
GT_INLINE void gt_alignment_union_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_union_alignment_maps_fx_v(gt_map_cmp_fx,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_alignment_union_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_union_alignment_maps_fx_v(gt_map_cmp,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}

GT_INLINE void gt_alignment_intersect_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  register uint64_t i;
  for (i=0;i<num_src_alignments;++i) {
    register gt_alignment* aux_alignment = (i==0) ? alignment_src : va_arg(v_args,gt_alignment*);
    GT_MAPS_ITERATE(aux_alignment,map) {
      if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_base,map)!=NULL &&
          gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_dst,map)==NULL) {
        gt_alignment_insert_map(alignment_dst,map);
      }
    }
  }
}
GT_INLINE void gt_alignment_intersect_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_intersect_alignment_maps_fx_v(gt_map_cmp_fx,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_alignment_intersect_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_intersect_alignment_maps_fx_v(gt_map_cmp,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}

GT_INLINE void gt_alignment_subtract_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  // Store source alignments from the arguments
  gt_alignment** aux_alignment = malloc(sizeof(gt_alignment*)*num_src_alignments);
  register uint64_t i;
  for (i=0;i<num_src_alignments;++i) {
    aux_alignment[i] = (i==0) ? alignment_src : va_arg(v_args,gt_alignment*);
  }
  GT_MAPS_ITERATE(alignment_base,map) {
    register bool subtracted = false;
    for (i=0;i<num_src_alignments;++i) {
      if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,aux_alignment[i],map)!=NULL) {
        subtracted = true; break;
      }
    }
    if (!subtracted) {
      gt_alignment_insert_map(alignment_dst,map);
    }
  }
  free(aux_alignment);
}
GT_INLINE void gt_alignment_subtract_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_subtract_alignment_maps_fx_v(gt_map_cmp_fx,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}
GT_INLINE void gt_alignment_subtract_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(alignment_dst); GT_NULL_CHECK(alignment_base);
  GT_ZERO_CHECK(num_src_alignments); GT_NULL_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment_subtract_alignment_maps_fx_v(gt_map_cmp,alignment_dst,alignment_base,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
}

