/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats.h
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
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

#define GT_STATS_INSS_RANGE_NEG  0
#define GT_STATS_INSS_RANGE_OVER 1
#define GT_STATS_INSS_RANGE_100  2
#define GT_STATS_INSS_RANGE_200  3
#define GT_STATS_INSS_RANGE_300  4
#define GT_STATS_INSS_RANGE_400  5
#define GT_STATS_INSS_RANGE_500  6
#define GT_STATS_INSS_RANGE_600  7
#define GT_STATS_INSS_RANGE_700  8
#define GT_STATS_INSS_RANGE_800  9
#define GT_STATS_INSS_RANGE_900  10
#define GT_STATS_INSS_RANGE_1000 11
#define GT_STATS_INSS_RANGE_2000 12
#define GT_STATS_INSS_RANGE_5000 13
#define GT_STATS_INSS_RANGE_10000 14
#define GT_STATS_INSS_RANGE_BEHOND 15
#define GT_STATS_INSS_RANGE 16

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

#define GT_STATS_MISMS_1_CONTEXT_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE)*GT_STATS_MISMS_RANGE)
#define GT_STATS_MISMS_2_CONTEXT_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE)*GT_STATS_MISMS_RANGE)
#define GT_STATS_INDEL_TRANSITION_1_RANGE ((GT_STATS_MISMS_BASE_RANGE))
#define GT_STATS_INDEL_TRANSITION_2_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE))
#define GT_STATS_INDEL_TRANSITION_3_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE))
#define GT_STATS_INDEL_TRANSITION_4_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE))
#define GT_STATS_INDEL_1_CONTEXT (2*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_RANGE)
#define GT_STATS_INDEL_2_CONTEXT (2*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_RANGE)

/*
 * Handy Functions
 */
#define GT_STATS_GET_PERCENTAGE(AMOUNT,TOTAL) ((TOTAL)?100.0*(float)(AMOUNT)/(float)(TOTAL):0.0)
#define GT_STATS_DIV(NUMERATOR,DENOMINATOR) ((DENOMINATOR)?(NUMERATOR)/(DENOMINATOR):(0))

/*
 * Stats Data Structures
 */
typedef struct {
  // Mismatch/Indel Profile
  uint64_t *mismatches;        /* GT_STATS_MISMS_RANGE */
  uint64_t *levenshtein;       /* GT_STATS_MISMS_RANGE */
  uint64_t *insertion_length;  /* GT_STATS_MISMS_RANGE */
  uint64_t *deletion_length;   /* GT_STATS_MISMS_RANGE */
  uint64_t *errors_events;     /* GT_STATS_MISMS_RANGE */
  uint64_t total_mismatches;
  uint64_t total_levenshtein;
  uint64_t total_indel_length;
  uint64_t total_errors_events;
  uint64_t *error_position;    /* GT_STATS_LARGE_READ_POS_RANGE */
  // Trim/Mapping stats
  uint64_t total_bases;
  uint64_t total_bases_matching;
  uint64_t total_bases_trimmed;
  // Insert Size Distribution
  uint64_t *inss;               /* GT_STATS_INSS_RANGE */
  // Mismatch/Errors bases
  uint64_t *misms_transition;   /* GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_RANGE */
  uint64_t *qual_score_misms;   /* GT_STATS_QUAL_SCORE_RANGE */
  uint64_t *misms_1context;     /* GT_STATS_MISMS_1_CONTEXT_RANGE */
  uint64_t *misms_2context;     /* GT_STATS_MISMS_2_CONTEXT_RANGE */
  uint64_t *indel_transition_1; /* GT_STATS_INDEL_TRANSITION_1_RANGE */
  uint64_t *indel_transition_2; /* GT_STATS_INDEL_TRANSITION_2_RANGE */
  uint64_t *indel_transition_3; /* GT_STATS_INDEL_TRANSITION_3_RANGE */
  uint64_t *indel_transition_4; /* GT_STATS_INDEL_TRANSITION_4_RANGE */
  uint64_t *indel_1context;     /* GT_STATS_INDEL_1_CONTEXT */
  uint64_t *indel_2context;     /* GT_STATS_INDEL_2_CONTEXT */
  uint64_t *qual_score_errors;  /* GT_STATS_QUAL_SCORE_RANGE */
} gt_maps_profile; // TODO Called map profile

typedef struct {
  // General SM
  /*
   * @num_mapped_with_splitmaps, @num_mapped_only_splitmaps
   *   At the level of alignments... not matter if SE or PE or how many maps
   */
  uint64_t num_mapped_with_splitmaps;
  uint64_t num_mapped_only_splitmaps;
  /*
   * @total_splitmaps :: Number of blocks with SM.
   *   Eg 'chr1:+:12345:10>20*90::chr1:+:12445:100' => 1
   *   Eg 'chr1:+:12345:10>20*90::chr1:+:12445:10>20*90' => 2
   *   Eg 'chr1:+:12345:10>20*10>20*80::chr1:+:12445:10>20*90' => 2
   */
  uint64_t total_splitmaps;
  uint64_t total_junctions;
  uint64_t *num_junctions; /* GT_STATS_NUM_JUNCTION_RANGE */
  uint64_t *length_junctions; /* GT_STATS_LEN_JUNCTION_RANGE */
  uint64_t *junction_position; /* GT_STATS_SHORT_READ_POS_RANGE */
  // Paired SM combinations
  uint64_t pe_sm_sm;
  uint64_t pe_sm_rm;
  uint64_t pe_rm_rm;
} gt_splitmaps_profile;

typedef struct {
  // Length stats
  uint64_t min_length;
  uint64_t max_length;
  uint64_t total_bases; // WRT to the read
  uint64_t total_bases_aligned; // WRT to reads mapped
  uint64_t mapped_min_length;
  uint64_t mapped_max_length;
  // Nucleotide counting (wrt to the maps=read+errors)
  uint64_t *nt_counting; /* GT_STATS_MISMS_BASE_RANGE */
  // Mapped/Maps
  uint64_t num_blocks; // SE => 1 block. PE => 2 blocks
  uint64_t num_alignments;
  uint64_t num_maps;
  uint64_t num_mapped;
  // MMap Distribution
  uint64_t *mmap; /* GT_STATS_MMAP_RANGE */
  // Uniq Distribution
  uint64_t *uniq; /* GT_STATS_UNIQ_RANGE */
  // Error Profile
  gt_maps_profile *maps_profile;
  // Split maps info
  gt_splitmaps_profile* splitmaps_profile;
} gt_stats;

typedef struct {
  bool best_map;
  bool nucleotide_stats;
  bool maps_profile;
  bool indel_profile; // TODO
  bool split_map_stats;
} gt_stats_analysis;

#define GT_STATS_ANALYSIS_DEFAULT() \
  { \
    .best_map=false, \
    .nucleotide_stats=true, \
    .maps_profile=true, \
    .indel_profile=false, \
    .split_map_stats=true \
  }

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
 *   BOTE: @seq_archive==NULL if no indel_profile is requested (default)
 */
GT_INLINE void gt_stats_calculate_template_stats(
    gt_stats* const stats,gt_template* const template,gt_sequence_archive* seq_archive,gt_stats_analysis* const stats_analysis);

#endif /* GT_STATS_H_ */
