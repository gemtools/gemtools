/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_handling.h
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */


#ifndef GT_ALIGNMENT_HANDLING_H_
#define GT_ALIGNMENT_HANDLING_H_

#include "gt_commons.h"
#include "gt_alignment.h"

/*
 * Alignment high-level insert (Update global state: counters, ...)
 */
GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_insert_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector);

/*
 * Alignment's Counters operators
 */
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment);
GT_INLINE uint64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment);
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const uint64_t max_allowed_strata);
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment);

/*
 * Alignment's Maps operators
 */
GT_INLINE gt_map* gt_alignment_is_map_contained(gt_alignment* const alignment,gt_map* const map);
GT_INLINE gt_map* gt_alignment_is_map_contained_fx(
    int64_t (*gt_map_cmp)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_add_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src);
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src);
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_ori);

/*
 * Alignment's Maps set-operators
 */
GT_INLINE void gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args);
GT_INLINE void gt_alignment_union_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
GT_INLINE void gt_alignment_union_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
#define gt_alignment_union_alignment_maps(alignment_dst,alignment_base,alignment_src) \
        gt_alignment_union_alignment_maps_fx_va(gt_map_cmp,alignment_dst,alignment_base,1,alignment_src)
#define gt_alignment_union_alignment_maps_fx(gt_map_cmp_fx,alignment_dst,alignment_base,alignment_src) \
        gt_alignment_union_alignment_maps_fx_va(gt_map_cmp_fx,alignment_dst,alignment_base,1,alignment_src)

GT_INLINE void gt_alignment_intersect_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args);
GT_INLINE void gt_alignment_intersect_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
GT_INLINE void gt_alignment_intersect_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
#define gt_alignment_intersect_alignment_maps(alignment_dst,alignment_base,alignment_src) \
        gt_alignment_intersect_alignment_maps_fx_va(gt_map_cmp,alignment_dst,alignment_base,1,alignment_src)
#define gt_alignment_intersect_alignment_maps_fx(gt_map_cmp_fx,alignment_dst,alignment_base,alignment_src) \
        gt_alignment_intersect_alignment_maps_fx_va(gt_map_cmp_fx,alignment_dst,alignment_base,1,alignment_src)

GT_INLINE void gt_alignment_subtract_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args);
GT_INLINE void gt_alignment_subtract_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
GT_INLINE void gt_alignment_subtract_alignment_maps_va(
    gt_alignment* const alignment_dst,gt_alignment* const alignment_base,
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
#define gt_alignment_subtract_alignment_maps(alignment_dst,alignment_base,alignment_src) \
        gt_alignment_subtract_alignment_maps_fx_va(gt_map_cmp,alignment_dst,alignment_base,1,alignment_src)
#define gt_alignment_subtract_alignment_maps_fx(gt_map_cmp_fx,alignment_dst,alignment_base,alignment_src) \
        gt_alignment_subtract_alignment_maps_fx_va(gt_map_cmp_fx,alignment_dst,alignment_base,1,alignment_src)

#endif /* GT_ALIGNMENT_HANDLING_H_ */
