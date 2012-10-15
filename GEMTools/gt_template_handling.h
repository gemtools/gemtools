/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_handling.h
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_TEMPLATE_HANDLING_H_
#define GT_TEMPLATE_HANDLING_H_

#include "gt_commons.h"
#include "gt_template.h"

/*
 * Template high-level insert (Update global state: counters, ...)
 */
// FIXME: Avoid Duplicates
GT_INLINE void gt_template_insert_mmap(gt_template* const template,gt_map** mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_insert_mmap_gtvector(gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_insert_mmap_v(gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args);
GT_INLINE void gt_template_insert_mmap_va(gt_template* const template,gt_mmap_attributes* const mmap_attr,...); /* gt_map* maps */

/*
 * Template's Counters operators
 */
GT_INLINE void gt_template_recalculate_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_min_matching_strata(gt_template* const template);
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata);
GT_INLINE bool gt_template_is_mapped(gt_template* const template);

/*
 * Template's Maps operators
 */
GT_INLINE gt_map** gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap);
GT_INLINE gt_map** gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap);
GT_INLINE void gt_template_add_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_src);

/*
 * Template's Maps set-operators
 */
GT_INLINE void gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE void gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE void gt_template_union_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_union_template_mmaps(template_dst,template_base,template_src) \
        gt_template_union_template_mmaps_fx_va(gt_map_cmp,template_dst,template_base,1,template_src)
#define gt_template_union_template_mmaps_fx(gt_mmap_cmp_fx,template_dst,template_base,template_src) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp_fx,template_dst,template_base,1,template_src)

GT_INLINE void gt_template_intersect_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE void gt_template_intersect_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE void gt_template_intersect_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_intersect_template_mmaps(template_dst,template_base,template_src) \
        gt_template_intersect_template_mmaps_fx_va(gt_map_cmp,template_dst,template_base,1,template_src)
#define gt_template_intersect_template_mmaps_fx(gt_mmap_cmp_fx,template_dst,template_base,template_src) \
        gt_template_intersect_template_mmaps_fx_va(gt_mmap_cmp_fx,template_dst,template_base,1,template_src)

GT_INLINE void gt_template_subtract_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE void gt_template_subtract_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE void gt_template_subtract_template_mmaps_va(
    gt_template* const template_dst,gt_template* const template_base,
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_subtract_template_mmaps(template_dst,template_base,template_src) \
        gt_template_subtract_template_mmaps_fx_va(gt_map_cmp,template_dst,template_base,1,template_src)
#define gt_template_subtract_template_mmaps_fx(gt_mmap_cmp_fx,template_dst,template_base,template_src) \
        gt_template_subtract_template_mmaps_fx_va(gt_mmap_cmp_fx,template_dst,template_base,1,template_src)

#endif /* GT_TEMPLATE_HANDLING_H_ */
