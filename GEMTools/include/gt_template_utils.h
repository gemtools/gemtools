/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_utils.h
 * DATE: 19/07/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_TEMPLATE_UTILS_H_
#define GT_TEMPLATE_UTILS_H_

#include "gt_commons.h"
#include "gt_alignment.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"

/*
 * Template basic tools
 */
GT_INLINE gt_status gt_template_deduce_alignments_tags(gt_template* const template);

/*
 * Template's MMaps high-level insertion (basic building block)
 */
GT_INLINE gt_map** gt_template_raw_put_mmap(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,
    const bool alignment_insertion,const bool delete_unused_mmap_members);
GT_INLINE gt_map** gt_template_put_mmap(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,
    const bool replace_dup,const bool alignment_insertion,const bool delete_unused_mmap_members);

/*
 * Template's MMaps high-level insertion operators (Update global state: counters, ...)
 */
GT_INLINE void gt_template_insert_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,const bool alignment_insertion);
GT_INLINE void gt_template_insert_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,const bool alignment_insertion);
GT_INLINE void gt_template_insert_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr,const bool alignment_insertion);
GT_INLINE void gt_template_insert_mmap_gtvector_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr,const bool alignment_insertion);

// TODO: Scheduled for v2.0
GT_INLINE void gt_template_remove_mmap(
    gt_template* const template,gt_map** const mmap);
GT_INLINE void gt_template_remove_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap);

GT_INLINE bool gt_template_find_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,
    uint64_t* const found_mmap_pos,gt_map*** const found_mmap,gt_mmap_attributes* const found_mmap_attr);
GT_INLINE bool gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap);
GT_INLINE bool gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap);

/*
 * Template's Counters operators
 */
GT_INLINE void gt_template_recalculate_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_min_matching_strata(gt_template* const template);
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata);
GT_INLINE bool gt_template_is_mapped(gt_template* const template);

/*
 * Template Set operators
 */
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_src);

// TODO: Scheduled for v2.0
GT_INLINE void gt_template_remove_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_remove_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_dst,gt_template* const template_src);

GT_INLINE gt_template* gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE gt_template* gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE gt_template* gt_template_union_template_mmaps_va(
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_union_template_mmaps(template_src_A,template_src_B) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp,2,template_src_A,template_src_B)
#define gt_template_union_template_mmaps_fx(gt_mmap_cmp_fx,template_src_A,template_src_B) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp_fx,2,template_src_A,template_src_B)

GT_INLINE gt_template* gt_template_subtract_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_minuend,gt_template* const template_subtrahend);
GT_INLINE gt_template* gt_template_subtract_template_mmaps(
    gt_template* const template_minuend,gt_template* const template_subtrahend);

GT_INLINE gt_template* gt_template_intersect_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template_src_A,gt_template* const template_src_B);
GT_INLINE gt_template* gt_template_intersect_template_mmaps(
    gt_template* const template_src_A,gt_template* const template_src_B);

#endif /* GT_TEMPLATE_UTILS_H_ */
