/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats.h
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 *
 * SCHEDULED TODO: min,avg,max length of mapped reads
 *                 ACGT counting of reads
 *                 INDEL TRANSITIONS => A --> AA
 *                                      A --> AAA
 *                                      A --> AAAA
 *                                      A --> AAAAA..A
 *                 CHECK Total bases aligned
 *                 SELECTIVE TEST SWITCH
 */

#ifndef GT_STATS_H_
#define GT_STATS_H_

#include "gt_commons.h"
#include "gt_compact_dna_string.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

/*
 * Range Definition
 */
#define GT_STATS_MMAP_RANGE_1 0
#define GT_STATS_MMAP_RANGE_5 1
#define GT_STATS_MMAP_RANGE_10 2
#define GT_STATS_MMAP_RANGE_50 3
#define GT_STATS_MMAP_RANGE_100 4
#define GT_STATS_MMAP_RANGE_500 5
#define GT_STATS_MMAP_RANGE_1000 6
#define GT_STATS_MMAP_RANGE_BEHOND 7
#define GT_STATS_MMAP_RANGE 8

#define GT_STATS_INSS_RANGE_INF_V  (-500)
#define GT_STATS_INSS_RANGE_SUP_V   10000
#define GT_STATS_INSS_RANGE_BUCKETS 10000

#define GT_STATS_MISMS_RANGE_0 0
#define GT_STATS_MISMS_RANGE_1 1
#define GT_STATS_MISMS_RANGE_2 2
#define GT_STATS_MISMS_RANGE_3 3
#define GT_STATS_MISMS_RANGE_4 4
#define GT_STATS_MISMS_RANGE_5 5
#define GT_STATS_MISMS_RANGE_6 6
#define GT_STATS_MISMS_RANGE_7 7
#define GT_STATS_MISMS_RANGE_8 8
#define GT_STATS_MISMS_RANGE_9 9
#define GT_STATS_MISMS_RANGE_10 10
#define GT_STATS_MISMS_RANGE_20 11
#define GT_STATS_MISMS_RANGE_50 12
#define GT_STATS_MISMS_RANGE_BEHOND 13
#define GT_STATS_MISMS_RANGE 14

#define GT_STATS_UNIQ_RANGE_0 0
#define GT_STATS_UNIQ_RANGE_1 1
#define GT_STATS_UNIQ_RANGE_2 2
#define GT_STATS_UNIQ_RANGE_3 3
#define GT_STATS_UNIQ_RANGE_10 4
#define GT_STATS_UNIQ_RANGE_50 5
#define GT_STATS_UNIQ_RANGE_100 6
#define GT_STATS_UNIQ_RANGE_500 7
#define GT_STATS_UNIQ_RANGE_BEHOND 8
#define GT_STATS_UNIQ_RANGE_X 9
#define GT_STATS_UNIQ_RANGE 10

#define GT_STATS_LARGE_READ_POS_RANGE 1000
#define GT_STATS_SHORT_READ_POS_RANGE 100

#define GT_STATS_MISMS_BASE_A 0
#define GT_STATS_MISMS_BASE_C 1
#define GT_STATS_MISMS_BASE_G 2
#define GT_STATS_MISMS_BASE_T 3
#define GT_STATS_MISMS_BASE_N 4
#define GT_STATS_MISMS_BASE_RANGE 5

#define GT_STATS_NT_BASE_A 0
#define GT_STATS_NT_BASE_C 1
#define GT_STATS_NT_BASE_G 2
#define GT_STATS_NT_BASE_T 3
#define GT_STATS_NT_BASE_N 4
#define GT_STATS_NT_BASE_RANGE 5

#define GT_STATS_QUAL_SCORE_RANGE 256

#define GT_STATS_NUM_JUNCTION_1      0
#define GT_STATS_NUM_JUNCTION_2      1
#define GT_STATS_NUM_JUNCTION_3      2
#define GT_STATS_NUM_JUNCTION_BEHOND 3
#define GT_STATS_NUM_JUNCTION_RANGE  4

#define GT_STATS_LEN_JUNCTION_100    0
#define GT_STATS_LEN_JUNCTION_1000   1
#define GT_STATS_LEN_JUNCTION_5000   2
#define GT_STATS_LEN_JUNCTION_10000  3
#define GT_STATS_LEN_JUNCTION_50000  4
#define GT_STATS_LEN_JUNCTION_BEHOND 5
#define GT_STATS_LEN_JUNCTION_RANGE  6

typedef struct {
  // Mismatch/Indel Distribution
  uint64_t *mismatches;    /* MISMS_RANGE */
  uint64_t *indel_length;  /* MISMS_RANGE */
  uint64_t *errors_events; /* MISMS_RANGE */
  uint64_t total_mismatches;
  uint64_t total_levenshtein;
  uint64_t total_errors_events;
  uint64_t total_indel_length;
  uint64_t *error_position; /* READ_POS_RANGE */
  // Trim/Mapping stats
  uint64_t total_bases;
  uint64_t total_bases_matching;
  uint64_t total_bases_trimmed;
  // Nucleotide stats (wrt to the maps=read+errors)
  int64_t nt_stats_error;
  // Insert Size Distribution
  uint64_t *inss; /* INSS_RANGE */
  // Mismatch/Errors bases
  uint64_t *misms_transition;  /* MISMS_BASE_RANGE*MISMS_BASE_RANGE */
  uint64_t *indel_transition;  /*  */
  uint64_t *qual_score_misms;  /* QUAL_SCORE_RANGE */
  uint64_t *qual_score_errors; /* QUAL_SCORE_RANGE */
} gt_maps_error_profile;

typedef struct {
  uint64_t num_mapped_with_splitmaps;
  uint64_t num_mapped_only_splitmaps;
  uint64_t total_splitmaps;
  uint64_t total_junctions;
  uint64_t *num_junctions; /* NUM_JUNCTION_RANGE */
  uint64_t *length_junctions; /* LEN_JUNCTION_RANGE */
  uint64_t *junction_position; /* READ_POS_RANGE */
  /* Paired SM combinations */
  uint64_t pe_sm_sm;
  uint64_t pe_sm_rm;
  uint64_t pe_rm_rm;
} gt_splitmaps_profile;

typedef struct {
  // Length stats
  uint64_t min_length;
  uint64_t max_length;
  uint64_t total_bases;
  uint64_t total_bases_aligned;
  uint64_t mapped_min_length;
  uint64_t mapped_max_length;
  // Mapped/Maps
  uint64_t num_blocks;
  uint64_t num_alignments;
  uint64_t num_maps;
  uint64_t num_mapped;
  // MMap Distribution
  uint64_t *mmap; /* MMAP_RANGE */
  // Uniq Distribution
  uint64_t *uniq; /* UNIQ_RANGE */
  // Error Profile
  gt_maps_error_profile *maps_error_profile;
  // Split maps info
  gt_splitmaps_profile* splitmaps_profile;
} gt_stats;


typedef struct {
  bool best_map;
  bool nucleotide_stats;
  bool error_profile;
  bool split_map_stats;
} gt_stats_analysis;

/*
 * STATS Profile
 */
GT_INLINE gt_stats* gt_stats_new();
GT_INLINE void gt_stats_clear(gt_stats *stats);
GT_INLINE void gt_stats_delete(gt_stats *stats);

/*
 * STATS Merge
 */
void gt_stats_merge(gt_stats** const stats,const uint64_t stats_array_size);

/*
 * Calculate stats
 */
GT_INLINE void gt_stats_calculate_alignment_stats(
    gt_stats* const stats,gt_alignment* const alignment,gt_stats_analysis* const stats_analysis);
GT_INLINE void gt_stats_calculate_template_stats(
    gt_stats* const stats,gt_template* const template,gt_stats_analysis* const stats_analysis);

#endif /* GT_STATS_H_ */
