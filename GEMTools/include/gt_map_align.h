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

/*
 * Error Codes
 */
#define GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ 5
#define GT_MAP_CHECK_ALG_NO_MISMS 10
#define GT_MAP_CHECK_ALG_BAD_MISMS 11
#define GT_MAP_CHECK_ALG_MISMATCH 12
#define GT_MAP_CHECK_ALG_MISMS_OUT_OF_SEQ 13
#define GT_MAP_CHECK_ALG_INS_OUT_OF_SEQ 20
#define GT_MAP_CHECK_ALG_DEL_OUT_OF_SEQ 30

// Compact Dynamic Programming Pattern (used in Myers' Fast Bit-Vector algorithm)
typedef struct {
  // TODO
} gt_cdp_pattern;
// Compact Vector Pattern (used in Hamming-ASM Bit-Vector algorithm)
typedef struct {
  // TODO
} gt_cv_pattern;

/*
 * Map check/recover operators
 */
GT_INLINE gt_status gt_map_block_check_alignment(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length);
GT_INLINE gt_status gt_map_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive);

GT_INLINE gt_status gt_map_block_recover_mismatches(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length);
GT_INLINE gt_status gt_map_recover_mismatches_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive);

/*
 * Map (Re)alignment operators
 */
GT_INLINE gt_status gt_map_block_realign_hamming(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t length);
GT_INLINE gt_status gt_map_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive);

GT_INLINE gt_status gt_map_block_realign_levenshtein(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,const bool ends_free);
GT_INLINE gt_status gt_map_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive);

GT_INLINE gt_status gt_map_block_realign_weighted(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,int32_t (*gt_weigh_fx)(char*,char*));
GT_INLINE gt_status gt_map_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*));

GT_INLINE void gt_map_search_global_alignment_hamming(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,const uint64_t max_hamming_distance);
GT_INLINE void gt_map_search_global_alignment_levenshtein(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,const uint64_t max_levenshtein_distance);
GT_INLINE void gt_map_search_global_alignment_weighted(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t max_scope,
    int32_t (*gt_weigh_fx)(char*,char*),const uint64_t score_threshold);

/*
 * Bit-compressed (Re)alignment operators (Hamming)
 */
GT_INLINE gt_cv_pattern* gt_map_compile_cv_pattern(char* const pattern);
GT_INLINE void gt_map_delete_cv_pattern(gt_cv_pattern* const cv_pattern);

GT_INLINE void gt_map_cv_realign(gt_map* const map,gt_cv_pattern* const cv_pattern,char* const sequence);
GT_INLINE void gt_map_cv_realign_sa(gt_map* const map,gt_cv_pattern* const cv_pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_cv_search_global_alignment(
    gt_map* const map,gt_cv_pattern* const cv_pattern,
    char* const sequence,const uint64_t max_scope,const uint64_t levenshtein_distance);

/*
 * Bit-compressed (Re)alignment operators (Levenshtein)
 */
GT_INLINE gt_cdp_pattern* gt_map_compile_cdp_pattern(char* const pattern);
GT_INLINE void gt_map_delete_cdp_pattern(gt_cdp_pattern* const cdp_pattern);

GT_INLINE void gt_map_cdp_realign(gt_map* const map,gt_cdp_pattern* const cdp_pattern,char* const sequence);
GT_INLINE void gt_map_cdp_realign_sa(gt_map* const map,gt_cdp_pattern* const cdp_pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_map_cdp_search_global_alignment(
    gt_map* const map,gt_cdp_pattern* const cdp_pattern,
    char* const sequence,const uint64_t max_scope,const uint64_t levenshtein_distance);

/*
 * Filters
 */


#endif /* GT_MAP_ALIGN_H_ */
