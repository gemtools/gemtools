/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_ALIGN_H_
#define GT_MAP_ALIGN_H_

#include "gt_commons.h"
#include "gt_map.h"
#include "gt_sequence_archive.h"

// Compact Dynamic Programming Pattern (used in Myers' Fast Bit-Vector algorithm)
typedef struct {
  // TODO
} gt_cdp_pattern;
// Compact Vector Pattern (used in Hamming-ASM Bit-Vector algorithm)
typedef struct {
  // TODO
} gt_cv_pattern;

/*
 * Map (Re)alignment operators
 */
GT_INLINE bool gt_map_check_alignment(gt_map* const map,char* const pattern,char* const sequence);
GT_INLINE bool gt_map_check_alignment_sa(gt_map* const map,char* const pattern,gt_sequence_archive* const sequence_archive);

GT_INLINE void gt_map_realign_hamming(gt_map* const map,char* const pattern,char* const sequence);
GT_INLINE void gt_map_realign_hamming_sa(gt_map* const map,char* const pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_realign_levenshtein(gt_map* const map,char* const pattern,char* const sequence);
GT_INLINE void gt_map_realign_levenshtein_sa(gt_map* const map,char* const pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_realign_weighted(
    gt_map* const map,char* const pattern,char* const sequence,
    int32_t (*gt_weigh_fx)(char*,char*));
GT_INLINE void gt_map_realign_weighted_sa(
    gt_map* const map,char* const pattern,
    gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*));

GT_INLINE void gt_map_search_global_alignment_hamming(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,const uint64_t hamming_distance);
GT_INLINE void gt_map_search_global_alignment_levenshtein(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,const uint64_t levenshtein_distance);
GT_INLINE void gt_map_search_global_alignment_weighted(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,
    int32_t (*gt_weigh_fx)(char*,char*),const uint64_t score_threshold);

/*
 * Bit-compressed (Re)alignment operators (Levenshtein)
 */
GT_INLINE gt_cv_pattern* gt_map_compile_cv_pattern(char* const pattern);
GT_INLINE void gt_map_delete_cv_pattern(gt_cv_pattern* const cv_pattern);

GT_INLINE void gt_map_cv_realign(gt_map* const map,gt_cv_pattern* const cv_pattern,char* const sequence);
GT_INLINE void gt_map_cv_realign_sa(gt_map* const map,gt_cv_pattern* const cv_pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_cv_search_global_alignment(
    gt_map* const map,gt_cv_pattern* const cv_pattern,
    char* const sequence,const uint64_t max_scope,const uint64_t levenshtein_distance);

GT_INLINE gt_cdp_pattern* gt_map_compile_cdp_pattern(char* const pattern);
GT_INLINE void gt_map_delete_cdp_pattern(gt_cdp_pattern* const cdp_pattern);

GT_INLINE void gt_map_cdp_realign(gt_map* const map,gt_cdp_pattern* const cdp_pattern,char* const sequence);
GT_INLINE void gt_map_cdp_realign_sa(gt_map* const map,gt_cdp_pattern* const cdp_pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_cdp_search_global_alignment(
    gt_map* const map,gt_cdp_pattern* const cdp_pattern,
    char* const sequence,const uint64_t max_scope,const uint64_t levenshtein_distance);

#endif /* GT_MAP_ALIGN_H_ */
