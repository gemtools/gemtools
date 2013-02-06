/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_utils.h
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
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
GT_INLINE void gt_template_deduce_alignments_tags(gt_template* const template);
GT_INLINE void gt_template_deduce_template_tag(gt_template* const template,gt_alignment* const alignment);

/*
 * Template's MMaps high-level insertion (basic building block)
 */
GT_INLINE gt_map** gt_template_raw_put_mmap(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE gt_map** gt_template_put_mmap(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,const bool replace_duplicated);

/*
 * Template's MMaps high-level insertion operators (Update global state: counters, ...)
 */
GT_INLINE void gt_template_insert_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_insert_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_insert_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_insert_mmap_gtvector_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr);

GT_INLINE bool gt_template_find_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,
    uint64_t* const found_mmap_pos,gt_map*** const found_mmap,gt_mmap_attributes* const found_mmap_attr);
GT_INLINE bool gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap);
GT_INLINE bool gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap);

#define GT_TEMPLATE_INSERT_SIZE_OK 0
#define GT_TEMPLATE_INSERT_SIZE_DIFFERENT_CONTIGS 1
#define GT_TEMPLATE_INSERT_SIZE_SAME_STRAND 2

GT_INLINE int64_t gt_template_get_insert_size(gt_map** const mmap,uint64_t* gt_error);

/*
 * Template's Counters operators
 */
GT_INLINE bool gt_template_is_mapped(gt_template* const template);
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata);
GT_INLINE void gt_template_recalculate_counters(gt_template* const template);

GT_INLINE int64_t gt_template_get_min_matching_strata(gt_template* const template);
GT_INLINE int64_t gt_template_get_uniq_degree(gt_template* const template);
GT_INLINE bool gt_template_get_next_matching_strata(
    gt_template* const template,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps);

/*
 * Template Set operators
 */
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_dst,gt_template* const template_src);

GT_INLINE gt_template* gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE gt_template* gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE gt_template* gt_template_union_template_mmaps_va(
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_union_template_mmaps(template_src_A,template_src_B) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp,gt_map_cmp,2,template_src_A,template_src_B)
#define gt_template_union_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_src_A,template_src_B) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp_fx,gt_map_cmp_fx,2,template_src_A,template_src_B)

GT_INLINE gt_template* gt_template_subtract_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_minuend,gt_template* const template_subtrahend);
GT_INLINE gt_template* gt_template_subtract_template_mmaps(
    gt_template* const template_minuend,gt_template* const template_subtrahend);

GT_INLINE gt_template* gt_template_intersect_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_A,gt_template* const template_B);
GT_INLINE gt_template* gt_template_intersect_template_mmaps(
    gt_template* const template_A,gt_template* const template_B);

/*
 * Template realignment
 */
GT_INLINE void gt_template_realign_hamming(gt_template* const template,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_template_realign_levenshtein(gt_template* const template,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_template_realign_weighted(
    gt_template* const template,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*));

#endif /* GT_TEMPLATE_UTILS_H_ */
