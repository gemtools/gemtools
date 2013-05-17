/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats.c
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple utility to output statistics from {MAP,SAM,FASTQ} files
 */

#include "gt_stats.h"

/*
 * Handy Macros
 */
#define GT_STATS_GET_PERCENTAGE_ERROR(percentage,one_per_cent) (uint64_t)((double)percentage*one_per_cent)
#define GT_STATS_VECTOR_ADD(VECTOR_A,VECTOR_B,RANGE) { \
  uint64_t iii; \
  for (iii=0;iii<RANGE;++iii) {VECTOR_A[iii] += VECTOR_B[iii];} \
}

/*
 * MAPS Error Profile
 */
GT_INLINE gt_maps_profile* gt_maps_profile_new() {
  // Allocate handler
  gt_maps_profile* maps_profile = gt_alloc(gt_maps_profile);
  /*
   * Init
   */
  // Mismatch/Indel Profile
  maps_profile->mismatches = gt_calloc(GT_STATS_MISMS_RANGE,uint64_t,true);
  maps_profile->levenshtein = gt_calloc(GT_STATS_MISMS_RANGE,uint64_t,true);
  maps_profile->insertion_length = gt_calloc(GT_STATS_MISMS_RANGE,uint64_t,true);
  maps_profile->deletion_length = gt_calloc(GT_STATS_MISMS_RANGE,uint64_t,true);
  maps_profile->errors_events = gt_calloc(GT_STATS_MISMS_RANGE,uint64_t,true);
  // Mismatch/Indel Distribution
  maps_profile->total_mismatches=0;
  maps_profile->total_levenshtein=0;
  maps_profile->total_indel_length=0;
  maps_profile->total_errors_events=0;
  maps_profile->error_position = gt_calloc(GT_STATS_LARGE_READ_POS_RANGE,uint64_t,true);
  // Trim/Mapping stats
  maps_profile->total_bases=0;
  maps_profile->total_bases_matching=0;
  maps_profile->total_bases_trimmed=0;
  // Strandness combinations
  maps_profile->pair_strand_rf=0;
  maps_profile->pair_strand_fr=0;
  maps_profile->pair_strand_ff=0;
  maps_profile->pair_strand_rr=0;
  // Insert Size Distribution
  maps_profile->inss = gt_calloc(GT_STATS_INSS_RANGE,uint64_t,true);
  maps_profile->inss_fine_grain = gt_calloc(GT_STATS_INSS_FG_RANGE,uint64_t,true);
  // Mismatch/Errors bases
  maps_profile->misms_transition = gt_calloc(GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE,uint64_t,true);
  maps_profile->qual_score_misms = gt_calloc(GT_STATS_QUAL_SCORE_RANGE,uint64_t,true);
  maps_profile->misms_1context = gt_calloc(GT_STATS_MISMS_1_CONTEXT_RANGE,uint64_t,true);
  maps_profile->misms_2context = gt_calloc(GT_STATS_MISMS_2_CONTEXT_RANGE,uint64_t,true);
  maps_profile->indel_transition_1 = gt_calloc(GT_STATS_INDEL_TRANSITION_1_RANGE,uint64_t,true);
  maps_profile->indel_transition_2 = gt_calloc(GT_STATS_INDEL_TRANSITION_2_RANGE,uint64_t,true);
  maps_profile->indel_transition_3 = gt_calloc(GT_STATS_INDEL_TRANSITION_3_RANGE,uint64_t,true);
  maps_profile->indel_transition_4 = gt_calloc(GT_STATS_INDEL_TRANSITION_4_RANGE,uint64_t,true);
  maps_profile->indel_1context = gt_calloc(GT_STATS_INDEL_1_CONTEXT,uint64_t,true);
  maps_profile->indel_2context = gt_calloc(GT_STATS_INDEL_2_CONTEXT,uint64_t,true);
  maps_profile->qual_score_errors = gt_calloc(GT_STATS_QUAL_SCORE_RANGE,uint64_t,true);
  return maps_profile;
}
GT_INLINE void gt_maps_profile_clear(gt_maps_profile* const maps_profile) {
  /*
   * Init
   */
  // Mismatch/Indel Profile
  memset(maps_profile->mismatches,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_profile->levenshtein,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_profile->insertion_length,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_profile->deletion_length,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_profile->errors_events,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  // Mismatch/Indel Distribution
  maps_profile->total_mismatches=0;
  maps_profile->total_levenshtein=0;
  maps_profile->total_indel_length=0;
  maps_profile->total_errors_events=0;
  memset(maps_profile->error_position,0,GT_STATS_LARGE_READ_POS_RANGE*sizeof(uint64_t));
  // Trim/Mapping stats
  maps_profile->total_bases=0;
  maps_profile->total_bases_matching=0;
  maps_profile->total_bases_trimmed=0;
  // Strandness combinations
  maps_profile->pair_strand_rf=0;
  maps_profile->pair_strand_fr=0;
  maps_profile->pair_strand_ff=0;
  maps_profile->pair_strand_rr=0;
  // Insert Size Distribution
  memset(maps_profile->inss,0,GT_STATS_INSS_RANGE*sizeof(uint64_t));
  memset(maps_profile->inss_fine_grain,0,GT_STATS_INSS_FG_RANGE*sizeof(uint64_t));
  // Mismatch/Errors bases
  memset(maps_profile->misms_transition,0,GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*sizeof(uint64_t));
  memset(maps_profile->qual_score_misms,0,GT_STATS_QUAL_SCORE_RANGE*sizeof(uint64_t));
  memset(maps_profile->misms_1context,0,GT_STATS_MISMS_1_CONTEXT_RANGE*sizeof(uint64_t));
  memset(maps_profile->misms_2context,0,GT_STATS_MISMS_2_CONTEXT_RANGE*sizeof(uint64_t));
  memset(maps_profile->indel_transition_1,0,GT_STATS_INDEL_TRANSITION_1_RANGE*sizeof(uint64_t));
  memset(maps_profile->indel_transition_2,0,GT_STATS_INDEL_TRANSITION_2_RANGE*sizeof(uint64_t));
  memset(maps_profile->indel_transition_3,0,GT_STATS_INDEL_TRANSITION_3_RANGE*sizeof(uint64_t));
  memset(maps_profile->indel_transition_4,0,GT_STATS_INDEL_TRANSITION_4_RANGE*sizeof(uint64_t));
  memset(maps_profile->indel_1context,0,GT_STATS_INDEL_1_CONTEXT*sizeof(uint64_t));
  memset(maps_profile->indel_2context,0,GT_STATS_INDEL_2_CONTEXT*sizeof(uint64_t));
  memset(maps_profile->qual_score_errors,0,GT_STATS_QUAL_SCORE_RANGE*sizeof(uint64_t));

}
GT_INLINE void gt_maps_profile_delete(gt_maps_profile* const maps_profile) {
  // Mismatch/Indel Profile
  gt_free(maps_profile->mismatches);
  gt_free(maps_profile->levenshtein);
  gt_free(maps_profile->insertion_length);
  gt_free(maps_profile->deletion_length);
  gt_free(maps_profile->errors_events);
  // Mismatch/Indel Distribution
  gt_free(maps_profile->error_position);
  // Insert Size Distribution
  gt_free(maps_profile->inss);
  gt_free(maps_profile->inss_fine_grain);
  // Mismatch/Errors bases
  gt_free(maps_profile->misms_transition);
  gt_free(maps_profile->qual_score_misms);
  gt_free(maps_profile->misms_1context);
  gt_free(maps_profile->misms_2context);
  gt_free(maps_profile->indel_transition_1);
  gt_free(maps_profile->indel_transition_2);
  gt_free(maps_profile->indel_transition_3);
  gt_free(maps_profile->indel_transition_4);
  gt_free(maps_profile->indel_1context);
  gt_free(maps_profile->indel_2context);
  gt_free(maps_profile->qual_score_errors);
  gt_free(maps_profile);
}
GT_INLINE void gt_maps_profile_merge(
    gt_maps_profile* const maps_profile_dst,gt_maps_profile* const maps_profile_src) {
  // Mismatch/Indel Profile
  GT_STATS_VECTOR_ADD(maps_profile_dst->mismatches,maps_profile_src->mismatches,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->levenshtein,maps_profile_src->levenshtein,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->insertion_length,maps_profile_src->insertion_length,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->deletion_length,maps_profile_src->deletion_length,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->errors_events,maps_profile_src->errors_events,GT_STATS_MISMS_RANGE);
  // Mismatch/Indel Distribution
  maps_profile_dst->total_mismatches+=maps_profile_src->total_mismatches;
  maps_profile_dst->total_levenshtein+=maps_profile_src->total_levenshtein;
  maps_profile_dst->total_indel_length+=maps_profile_src->total_indel_length;
  maps_profile_dst->total_errors_events+=maps_profile_src->total_errors_events;
  GT_STATS_VECTOR_ADD(maps_profile_dst->error_position,maps_profile_src->error_position,GT_STATS_LARGE_READ_POS_RANGE);
  // Trim/Mapping stats
  maps_profile_dst->total_bases+=maps_profile_src->total_bases;
  maps_profile_dst->total_bases_matching+=maps_profile_src->total_bases_matching;
  maps_profile_dst->total_bases_trimmed+=maps_profile_src->total_bases_trimmed;
  // Strandness combinations
  maps_profile_dst->pair_strand_rf+=maps_profile_src->pair_strand_rf;
  maps_profile_dst->pair_strand_fr+=maps_profile_src->pair_strand_fr;
  maps_profile_dst->pair_strand_ff+=maps_profile_src->pair_strand_ff;
  maps_profile_dst->pair_strand_rr+=maps_profile_src->pair_strand_rr;
  // Insert Size Distribution
  GT_STATS_VECTOR_ADD(maps_profile_dst->inss,maps_profile_src->inss,GT_STATS_INSS_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->inss_fine_grain,maps_profile_src->inss_fine_grain,GT_STATS_INSS_FG_RANGE);
  // Mismatch/Errors bases
  GT_STATS_VECTOR_ADD(maps_profile_dst->misms_transition,maps_profile_src->misms_transition,GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->qual_score_misms,maps_profile_src->qual_score_misms,GT_STATS_QUAL_SCORE_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->misms_1context,maps_profile_src->misms_1context,GT_STATS_MISMS_1_CONTEXT_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->misms_2context,maps_profile_src->misms_2context,GT_STATS_MISMS_2_CONTEXT_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_transition_1,maps_profile_src->indel_transition_1,GT_STATS_INDEL_TRANSITION_1_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_transition_2,maps_profile_src->indel_transition_2,GT_STATS_INDEL_TRANSITION_2_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_transition_3,maps_profile_src->indel_transition_3,GT_STATS_INDEL_TRANSITION_3_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_transition_4,maps_profile_src->indel_transition_4,GT_STATS_INDEL_TRANSITION_4_RANGE);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_1context,maps_profile_src->indel_1context,GT_STATS_INDEL_1_CONTEXT);
  GT_STATS_VECTOR_ADD(maps_profile_dst->indel_2context,maps_profile_src->indel_2context,GT_STATS_INDEL_2_CONTEXT);
  GT_STATS_VECTOR_ADD(maps_profile_dst->qual_score_errors,maps_profile_src->qual_score_errors,GT_STATS_QUAL_SCORE_RANGE);
}

/*
 * SPLITMAPS Profile
 */
GT_INLINE gt_splitmaps_profile* gt_splitmaps_profile_new() {
  // Allocate handler
  gt_splitmaps_profile* splitmaps_profile = gt_alloc(gt_splitmaps_profile);
  /*
   * Init
   */
  // General SM
  splitmaps_profile->num_mapped_with_splitmaps = 0;
  splitmaps_profile->num_mapped_only_splitmaps = 0;
  splitmaps_profile->total_splitmaps = 0;
  splitmaps_profile->total_junctions = 0;
  splitmaps_profile->num_junctions = gt_calloc(GT_STATS_NUM_JUNCTION_RANGE,uint64_t,true);
  splitmaps_profile->length_junctions = gt_calloc(GT_STATS_LEN_JUNCTION_RANGE,uint64_t,true);
  splitmaps_profile->junction_position = gt_calloc(GT_STATS_SHORT_READ_POS_RANGE,uint64_t,true);
  // Paired SM combinations
  splitmaps_profile->pe_sm_sm = 0;
  splitmaps_profile->pe_sm_rm = 0;
  splitmaps_profile->pe_rm_rm = 0;
  return splitmaps_profile;
}
GT_INLINE void gt_splitmaps_profile_clear(gt_splitmaps_profile* const splitmaps_profile) {
  /*
   * Init
   */
  // General SM
  splitmaps_profile->num_mapped_with_splitmaps = 0;
  splitmaps_profile->num_mapped_only_splitmaps = 0;
  splitmaps_profile->total_splitmaps = 0;
  splitmaps_profile->total_junctions = 0;
  memset(splitmaps_profile->num_junctions,0,GT_STATS_NUM_JUNCTION_RANGE*sizeof(uint64_t));
  memset(splitmaps_profile->length_junctions,0,GT_STATS_LEN_JUNCTION_RANGE*sizeof(uint64_t));
  memset(splitmaps_profile->junction_position,0,GT_STATS_SHORT_READ_POS_RANGE*sizeof(uint64_t));
  // Paired SM combinations
  splitmaps_profile->pe_sm_sm = 0;
  splitmaps_profile->pe_sm_rm = 0;
  splitmaps_profile->pe_rm_rm = 0;
}
GT_INLINE void gt_splitmaps_profile_delete(gt_splitmaps_profile* const splitmaps_profile) {
  gt_free(splitmaps_profile->num_junctions);
  gt_free(splitmaps_profile->length_junctions);
  gt_free(splitmaps_profile->junction_position);
  gt_free(splitmaps_profile);
}
GT_INLINE void gt_splitmaps_profile_merge(
    gt_splitmaps_profile* const splitmaps_profile_dst,gt_splitmaps_profile* const splitmaps_profile_src) {
  // General SM
  splitmaps_profile_dst->num_mapped_with_splitmaps += splitmaps_profile_src->num_mapped_with_splitmaps;
  splitmaps_profile_dst->num_mapped_only_splitmaps += splitmaps_profile_src->num_mapped_only_splitmaps;
  splitmaps_profile_dst->total_splitmaps += splitmaps_profile_src->total_splitmaps;
  splitmaps_profile_dst->total_junctions += splitmaps_profile_src->total_junctions;
  GT_STATS_VECTOR_ADD(splitmaps_profile_dst->num_junctions,splitmaps_profile_src->num_junctions,GT_STATS_NUM_JUNCTION_RANGE);
  GT_STATS_VECTOR_ADD(splitmaps_profile_dst->length_junctions,splitmaps_profile_src->length_junctions,GT_STATS_LEN_JUNCTION_RANGE);
  GT_STATS_VECTOR_ADD(splitmaps_profile_dst->junction_position,splitmaps_profile_src->junction_position,GT_STATS_SHORT_READ_POS_RANGE);
  // Paired SM combinations
  splitmaps_profile_dst->pe_sm_sm += splitmaps_profile_src->pe_sm_sm;
  splitmaps_profile_dst->pe_sm_rm += splitmaps_profile_src->pe_sm_rm;
  splitmaps_profile_dst->pe_rm_rm += splitmaps_profile_src->pe_rm_rm;
}

/*
 * Calculate Distributions
 */
GT_INLINE void gt_stats_get_misms(uint64_t* const misms_array,const uint64_t misms_num,const uint64_t read_length) {
  const double one_per_cent = read_length/100.0;
  if (misms_num==0) {
    misms_array[GT_STATS_MISMS_RANGE_0]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(1,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_1]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(2,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_2]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(3,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_3]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(4,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_4]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(5,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_5]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(6,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_6]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(7,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_7]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(8,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_8]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(9,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_9]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(10,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_10]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(20,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_20]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(50,one_per_cent)) {
    misms_array[GT_STATS_MISMS_RANGE_50]++;
  } else {
    misms_array[GT_STATS_MISMS_RANGE_BEHOND]++;
  }
}
GT_INLINE void gt_stats_get_mmap_distribution(uint64_t* const mmap,const int64_t num_maps) {
  if (num_maps>0) {
    if (num_maps<=1) {
      ++mmap[GT_STATS_MMAP_RANGE_1];
    } else if (num_maps<=5) {
      ++mmap[GT_STATS_MMAP_RANGE_5];
    } else if (num_maps<=10) {
      ++mmap[GT_STATS_MMAP_RANGE_10];
    } else if (num_maps<=50) {
      ++mmap[GT_STATS_MMAP_RANGE_50];
    } else if (num_maps<=100) {
      ++mmap[GT_STATS_MMAP_RANGE_100];
    } else if (num_maps<=500) {
      ++mmap[GT_STATS_MMAP_RANGE_500];
    } else if (num_maps<=1000) {
      ++mmap[GT_STATS_MMAP_RANGE_1000];
    } else {
      ++mmap[GT_STATS_MMAP_RANGE_BEHOND];
    }
  }
}
GT_INLINE void gt_stats_get_uniq_distribution(uint64_t* const uniq,const int64_t uniq_degree) {
  if (uniq_degree<=GT_NO_STRATA) {
    ++uniq[GT_STATS_UNIQ_RANGE_X];
  } else if (uniq_degree==0) {
    ++uniq[GT_STATS_UNIQ_RANGE_0];
  } else if (uniq_degree<=1) {
    ++uniq[GT_STATS_UNIQ_RANGE_1];
  } else if (uniq_degree<=2) {
    ++uniq[GT_STATS_UNIQ_RANGE_2];
  } else if (uniq_degree<=3) {
    ++uniq[GT_STATS_UNIQ_RANGE_3];
  } else if (uniq_degree<=10) {
    ++uniq[GT_STATS_UNIQ_RANGE_10];
  } else if (uniq_degree<=50) {
    ++uniq[GT_STATS_UNIQ_RANGE_50];
  } else if (uniq_degree<=100) {
    ++uniq[GT_STATS_UNIQ_RANGE_100];
  } else if (uniq_degree<=500) {
    ++uniq[GT_STATS_UNIQ_RANGE_500];
  } else {
    ++uniq[GT_STATS_UNIQ_RANGE_BEHOND];
  }
}
GT_INLINE void gt_stats_get_inss_distribution(uint64_t* const inss,const int64_t insert_size) {
  if (insert_size<=(-100)) {
    ++inss[GT_STATS_INSS_RANGE_NEG];
  } else if (insert_size<=0) {
    ++inss[GT_STATS_INSS_RANGE_OVER];
  } else if (insert_size<=100) {
    ++inss[GT_STATS_INSS_RANGE_100];
  } else if (insert_size<=200) {
    ++inss[GT_STATS_INSS_RANGE_200];
  } else if (insert_size<=300) {
    ++inss[GT_STATS_INSS_RANGE_300];
  } else if (insert_size<=400) {
    ++inss[GT_STATS_INSS_RANGE_400];
  } else if (insert_size<=500) {
    ++inss[GT_STATS_INSS_RANGE_500];
  } else if (insert_size<=600) {
    ++inss[GT_STATS_INSS_RANGE_600];
  } else if (insert_size<=700) {
    ++inss[GT_STATS_INSS_RANGE_700];
  } else if (insert_size<=800) {
    ++inss[GT_STATS_INSS_RANGE_800];
  } else if (insert_size<=900) {
    ++inss[GT_STATS_INSS_RANGE_900];
  } else if (insert_size<=1000) {
    ++inss[GT_STATS_INSS_RANGE_1000];
  } else if (insert_size<=2000) {
    ++inss[GT_STATS_INSS_RANGE_2000];
  } else if (insert_size<=5000) {
    ++inss[GT_STATS_INSS_RANGE_5000];
  } else if (insert_size<=10000) {
    ++inss[GT_STATS_INSS_RANGE_10000];
  } else {
    ++inss[GT_STATS_INSS_RANGE_BEHOND];
  }
}
// Fine-grain version
GT_INLINE void gt_stats_get_inss_fg_distribution(uint64_t* const inss_fg,const int64_t insert_size) {
  // Check boundaries
  if (insert_size<GT_STATS_INSS_FG_MIN) {
    ++inss_fg[0];
  } else if (insert_size>GT_STATS_INSS_FG_MAX) {
    ++inss_fg[GT_STATS_INSS_FG_RANGE-1];
  } else { // Fit insert_size into the buckets
    ++inss_fg[GT_STATS_INSS_FG_GET_BUCKET(insert_size)];
  }
}

GT_INLINE void gt_stats_get_juntions_distribution(uint64_t* const junctions,const uint64_t num_junctions) {
  //if (num_junctions>0) {
    if (num_junctions<=1) {
      ++junctions[GT_STATS_NUM_JUNCTION_1];
    } else if (num_junctions==2) {
      ++junctions[GT_STATS_NUM_JUNCTION_2];
    } else if (num_junctions==3) {
      ++junctions[GT_STATS_NUM_JUNCTION_3];
    } else {
      ++junctions[GT_STATS_NUM_JUNCTION_BEHOND];
    }
  //}
}
GT_INLINE void gt_stats_get_juntions_length_distribution(uint64_t* const junctions_length,const uint64_t length) {
  //if (length>0) {
    if (length<=100) {
      ++junctions_length[GT_STATS_LEN_JUNCTION_100];
    } else if (length<=1000) {
      ++junctions_length[GT_STATS_LEN_JUNCTION_1000];
    } else if (length<=5000) {
      ++junctions_length[GT_STATS_LEN_JUNCTION_5000];
    } else if (length<=10000) {
      ++junctions_length[GT_STATS_LEN_JUNCTION_10000];
    } else if (length<=50000) {
      ++junctions_length[GT_STATS_LEN_JUNCTION_50000];
    } else {
      ++junctions_length[GT_STATS_LEN_JUNCTION_BEHOND];
    }
  //}
}
GT_INLINE void gt_stats_get_nucleotide_stats(uint64_t* const nt_counting,gt_string* const read) {
  GT_STRING_ITERATE(read,read_string,pos) {
    ++nt_counting[gt_cdna_encode[(uint8_t)read_string[pos]]];
  }
}

/*
 * STATS Profile
 */
GT_INLINE gt_stats* gt_stats_new() {
  // Allocate handler
  gt_stats* stats = gt_alloc(gt_stats);
  // Length
  stats->min_length=UINT64_MAX;
  stats->max_length=0;
  stats->total_bases=0;
  stats->total_bases_aligned=0;
  stats->mapped_min_length=UINT64_MAX;
  stats->mapped_max_length=0;
  // Nucleotide counting (wrt to the maps=read+errors)
  stats->nt_counting = gt_calloc(GT_STATS_MISMS_BASE_RANGE,uint64_t,true);
  // Mapped/Maps/MMaps/Uniq...
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  stats->mmap = gt_calloc(GT_STATS_MMAP_RANGE,uint64_t,true); // MMaps
  stats->uniq = gt_calloc(GT_STATS_UNIQ_RANGE,uint64_t,true); // Uniq
  // Maps Error Profile
  stats->maps_profile = gt_maps_profile_new();
  // Split maps Profile
  stats->splitmaps_profile = gt_splitmaps_profile_new();
  return stats;
}
GT_INLINE void gt_stats_clear(gt_stats* const stats) {
  // Length
  stats->min_length=UINT64_MAX;
  stats->max_length=0;
  stats->total_bases=0;
  stats->total_bases_aligned=0;
  stats->mapped_min_length=UINT64_MAX;
  stats->mapped_max_length=0;
  // Nucleotide counting (wrt to the maps=read+errors)
  memset(stats->nt_counting,0,GT_STATS_MISMS_BASE_RANGE*sizeof(uint64_t)); // MMaps
  // Mapped/Maps/MMaps/Uniq...
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  // MMap Distribution
  memset(stats->mmap,0,GT_STATS_MMAP_RANGE*sizeof(uint64_t)); // MMaps
  // Uniq Distribution
  memset(stats->uniq,0,GT_STATS_UNIQ_RANGE*sizeof(uint64_t)); // Uniq
  // Maps Error Profile
  gt_maps_profile_clear(stats->maps_profile);
  // Split maps Profile
  gt_splitmaps_profile_clear(stats->splitmaps_profile);
}
GT_INLINE void gt_stats_delete(gt_stats* const stats) {
  gt_free(stats->mmap);
  gt_free(stats->uniq);
  gt_maps_profile_delete(stats->maps_profile);
  gt_splitmaps_profile_delete(stats->splitmaps_profile);
  gt_free(stats);
}

/*
 * STATS Merge
 */
GT_INLINE void gt_stats_merge(gt_stats** const stats,const uint64_t stats_array_size) {
  uint64_t i;
  for (i=1;i<stats_array_size;++i) {
    // Length
    stats[0]->min_length = GT_MIN(stats[0]->min_length,stats[i]->min_length);
    stats[0]->max_length = GT_MAX(stats[0]->max_length,stats[i]->max_length);
    stats[0]->total_bases += stats[i]->total_bases;
    stats[0]->total_bases_aligned += stats[i]->total_bases_aligned;
    stats[0]->mapped_min_length = GT_MIN(stats[0]->mapped_min_length,stats[i]->mapped_min_length);
    stats[0]->mapped_max_length = GT_MAX(stats[0]->mapped_max_length,stats[i]->mapped_max_length);
    // Nucleotide counting (wrt to the maps=read+errors)
    GT_STATS_VECTOR_ADD(stats[0]->nt_counting,stats[i]->nt_counting,GT_STATS_MISMS_BASE_RANGE);
    // Mapped/Maps
    stats[0]->num_blocks += stats[i]->num_blocks;
    stats[0]->num_alignments += stats[i]->num_alignments;
    stats[0]->num_maps += stats[i]->num_maps;
    stats[0]->num_mapped += stats[i]->num_mapped;
    // MMap Distribution
    GT_STATS_VECTOR_ADD(stats[0]->mmap,stats[i]->mmap,GT_STATS_MMAP_RANGE);
    // Uniq Distribution
    GT_STATS_VECTOR_ADD(stats[0]->uniq,stats[i]->uniq,GT_STATS_UNIQ_RANGE);
    // Merge Maps Error Profile
    gt_maps_profile_merge(stats[0]->maps_profile,stats[i]->maps_profile);
    // Merge SplitMaps Profile
    gt_splitmaps_profile_merge(stats[0]->splitmaps_profile,stats[i]->splitmaps_profile);
    // Free Handlers
    gt_stats_delete(stats[i]);
  }
}

/*
 * Calculate stats
 */
GT_INLINE void gt_stats_make_indel_profile(
    gt_maps_profile *maps_error_profile,
    gt_template* const template,const uint64_t alignment_total_length) {

}
GT_INLINE void gt_stats_make_maps_error_profile(
    gt_maps_profile *maps_error_profile,gt_template* const template,
    const uint64_t alignment_total_length,gt_map** const mmap) {
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t total_mismatches=0;
  uint64_t total_levenshtein=0;
  uint64_t total_ins_length=0, total_del_length=0;
  uint64_t total_errors_events=0;
  uint64_t total_bases=0, total_bases_trimmed=0, total_bases_not_matching=0;
  // Iterate all MMAPS
  GT_MMAP_ITERATE_ENDS(mmap,num_blocks,map,end_pos) { // Iterate over ends (/1,/2)
    gt_alignment* const alignment = gt_template_get_block(template,end_pos);
    gt_string* const read = alignment->read;
    gt_string* const quals = alignment->qualities;
    const bool has_qualities = gt_alignment_has_qualities(alignment);
    char quality_misms = 0;
    GT_MAP_ITERATE(map,map_block) { // Iterate over splits (map blocks)
      total_bases += gt_map_get_base_length(map_block);
      GT_MISMS_ITERATE(map_block,misms) { // Iterate over mismatches/indels
        ++total_errors_events;
        // Records position of misms/indel
        if (misms->position < GT_STATS_LARGE_READ_POS_RANGE) {
          maps_error_profile->error_position[misms->position]++;
        }
        // Record quality of misms/indel
        if (has_qualities) {
          quality_misms = gt_string_get_string(quals)[misms->position];
          maps_error_profile->qual_score_errors[(uint8_t)quality_misms]++;
        }
        switch (misms->misms_type) {
          case MISMS:
            ++total_mismatches;
            ++total_levenshtein;
            ++total_bases_not_matching;
            // Record transition
            gt_check(misms->base==gt_string_get_string(read)[misms->position],MISMS_TRANSITION);
            uint64_t idx = 0;
            idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
            idx *= GT_STATS_MISMS_BASE_RANGE;
            idx += gt_cdna_encode[(uint8_t)misms->base];
            maps_error_profile->misms_transition[idx]++;
            if (misms->position>0 && misms->position<gt_map_get_base_length(map_block)-1) { // 1-context
              idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-1]];
              idx *= GT_STATS_MISMS_BASE_RANGE;
              idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
              idx *= GT_STATS_MISMS_BASE_RANGE;
              idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+1]];
              idx *= GT_STATS_MISMS_BASE_RANGE;
              idx += gt_cdna_encode[(uint8_t)misms->base];
              maps_error_profile->misms_1context[idx]++;
//              // TODO: Think about it... overkilling
//              if (misms->position>1 && misms->position<map->base_length-2) { // 2-context
//                idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-2]];
//                idx *= GT_STATS_MISMS_BASE_RANGE;
//                idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-1]];
//                idx *= GT_STATS_MISMS_BASE_RANGE;
//                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
//                idx *= GT_STATS_MISMS_BASE_RANGE;
//                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+1]];
//                idx *= GT_STATS_MISMS_BASE_RANGE;
//                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+2]];
//                idx *= GT_STATS_MISMS_BASE_RANGE;
//                idx += gt_cdna_encode[(uint8_t)misms->base];
//                maps_error_profile->misms_2context[idx]++;
//              }
            }
            // Record quality of misms
            if (has_qualities) {
              maps_error_profile->qual_score_misms[(uint8_t)quality_misms]++;
            }
            break;
          case DEL:
            total_levenshtein += misms->size;
            total_del_length += misms->size;
            total_bases_not_matching += misms->size;
            // Record trim
            if (misms->position==0 || misms->position+misms->size==gt_map_get_base_length(map_block)) {
              total_bases_trimmed += misms->size;
            }
            break;
          case INS:
            total_levenshtein += misms->size;
            total_ins_length += misms->size;
            break;
        }
      }
    }
  }
  // Record general error stats
  maps_error_profile->total_bases += total_bases;
  maps_error_profile->total_bases_matching += total_bases-total_bases_not_matching;
  maps_error_profile->total_bases_trimmed += total_bases_trimmed;
  maps_error_profile->total_mismatches += total_mismatches;
  maps_error_profile->total_levenshtein += total_levenshtein;
  maps_error_profile->total_indel_length += total_ins_length+total_del_length;
  maps_error_profile->total_errors_events += total_errors_events;
  // Record the distribution of the error stats
  gt_stats_get_misms(maps_error_profile->mismatches,total_mismatches,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->levenshtein,total_levenshtein,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->insertion_length,total_ins_length,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->deletion_length,total_del_length,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->errors_events,total_errors_events,alignment_total_length);
}

GT_INLINE void gt_stats_make_mmaps_profile(
    gt_stats* const stats,gt_template* const template,
    const uint64_t alignment_total_length,gt_stats_analysis* const stats_analysis) {
  // Check not null maps
  if (gt_template_get_num_mmaps(template)==0) return;
  const uint64_t num_blocks_template = gt_template_get_num_blocks(template);
  const uint64_t paired_map = (num_blocks_template==2);
  // Iterate over all/best maps
  gt_maps_profile* const maps_error_profile = stats->maps_profile;
  gt_splitmaps_profile* const splitmaps_profile = stats->splitmaps_profile;
  bool has_splitsmaps = false;
  bool only_splitsmaps = true;
  GT_TEMPLATE_ITERATE(template,mmap) {
    /*
     * Insert Size Distribution
     */
    if (paired_map) {
      gt_status gt_err;
      int64_t ins_size = gt_template_get_insert_size(mmap,&gt_err);
      if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK) {
        gt_stats_get_inss_distribution(maps_error_profile->inss,ins_size);
        gt_stats_get_inss_fg_distribution(maps_error_profile->inss_fine_grain,ins_size);
      }
      if (mmap[0]->strand==FORWARD && mmap[1]->strand==REVERSE) { // F+R
        ++maps_error_profile->pair_strand_fr;
      } else if (mmap[0]->strand==REVERSE && mmap[1]->strand==FORWARD) { // R+F
        ++maps_error_profile->pair_strand_rf;
      } else if (mmap[0]->strand==FORWARD && mmap[1]->strand==FORWARD) { // F+F
        ++maps_error_profile->pair_strand_ff;
      } else { // R+R
        ++maps_error_profile->pair_strand_rr;
      }
    }
    /*
     * Error Profile
     */
    if (stats_analysis->maps_profile) { // TODO: bypass this to the upper level
      gt_stats_make_maps_error_profile(maps_error_profile,template,alignment_total_length,mmap);
    }
    /*
     * SplitMap Stats
     */
    if (stats_analysis->split_map_stats) {
      // SM block stats
      bool has_sm[2] = {true, true};
      GT_MMAP_ITERATE(mmap,map,end_pos) {
        const uint64_t num_blocks = gt_map_get_num_blocks(map);
        // Calculate total_junctions & total_splitmaps
        if (num_blocks > 1) {
          has_splitsmaps = true;
          splitmaps_profile->total_splitmaps++;
          gt_stats_get_juntions_distribution(splitmaps_profile->num_junctions,num_blocks-1);
        } else {
          has_sm[end_pos] = false;
        }
        // Calculate junction_position & length_junctions
        GT_MAP_ITERATE(map,map_block) {
          if (gt_map_has_next_block(map_block)) {
            splitmaps_profile->total_junctions++;
            gt_stats_get_juntions_length_distribution(splitmaps_profile->length_junctions,gt_map_get_junction_size(map_block));
            const uint64_t juntion_position = gt_map_get_base_length(map_block);
            if (juntion_position < GT_STATS_SHORT_READ_POS_RANGE) splitmaps_profile->junction_position[juntion_position]++;
          }
        }
      }
      // {SM,RM} Blocks combinations
      if (paired_map) {
        if (has_sm[0] && has_sm[1]) {
          splitmaps_profile->pe_sm_sm++;
        } else if (has_sm[0] || has_sm[1]) {
          only_splitsmaps = false;
          splitmaps_profile->pe_sm_rm++;
        } else {
          only_splitsmaps = false;
          splitmaps_profile->pe_rm_rm++;
        }
      } else {
        if (!has_sm[0]) {
          only_splitsmaps = false;
        }
      }
    }
    /*
     * BEST-MAP :: Break if we just one the best map
     */
    if (stats_analysis->best_map) break;
  } // GT_TEMPLATE_ITERATE_END;
  /*
   * Global SM stats about the whole alignment
   */
  if (stats_analysis->split_map_stats && has_splitsmaps) {
    splitmaps_profile->num_mapped_with_splitmaps++;
    if (only_splitsmaps) splitmaps_profile->num_mapped_only_splitmaps++;
  }
}
GT_INLINE void gt_stats_calculate_template_stats(
    gt_stats* const stats,gt_template* const template,gt_sequence_archive* seq_archive,gt_stats_analysis* const stats_analysis) {
  // Basic stats
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t num_maps = gt_template_get_num_mmaps(template);
  if (stats_analysis->best_map && num_maps>0) num_maps = 1; // NumMaps correction
  const bool is_mapped = (num_maps>0 || gt_template_is_mapped(template));
  /*
   * Length STATS
   */
  uint64_t alignment_total_length = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    // Read length stats
    uint64_t read_length =
        (gt_alignment_get_read_length(alignment)==0 && gt_alignment_get_num_maps(alignment)>0) ?
         gt_map_get_base_length(gt_alignment_get_map(alignment,0)) : gt_alignment_get_read_length(alignment);
    alignment_total_length += read_length;
    stats->min_length = GT_MIN(stats->min_length,read_length);
    stats->max_length = GT_MAX(stats->max_length,read_length);
    if (is_mapped) {
      stats->mapped_min_length = GT_MIN(stats->mapped_min_length,read_length);
      stats->mapped_max_length = GT_MAX(stats->mapped_max_length,read_length);
    }
    // Nucleotide stats
    gt_stats_get_nucleotide_stats(stats->nt_counting,alignment->read);
  }
  stats->total_bases += alignment_total_length;
  /*
   * Uniq Distribution
   */
  gt_stats_get_uniq_distribution(stats->uniq,gt_template_get_uniq_degree(template));
  /*
   * Maps/MMaps STATS
   */
  ++stats->num_alignments;
  stats->num_blocks += num_blocks;
  if (is_mapped){ // Is mapped?
    ++stats->num_mapped; // Mapped/Maps
    stats->num_maps += num_maps;
    stats->total_bases_aligned += alignment_total_length;
    gt_stats_get_mmap_distribution(stats->mmap,num_maps); // MMap Distribution
  }
  /*
   * MMaps Profile {Insert Size Distribution, Error Profile, SM Profile, ...}
   */
  if (num_maps > 0) {
    gt_stats_make_mmaps_profile(stats,template,alignment_total_length,stats_analysis);
    if (stats_analysis->indel_profile) {
      gt_stats_make_indel_profile(stats->maps_profile,template,alignment_total_length);
    }
  }
}

/*
 * STATS Report Output Printers
 */
GT_INLINE uint64_t gt_stats_sum_misms_pos(uint64_t* const pos_error,uint64_t const begin,uint64_t const end) {
  uint64_t i, accum = 0;
  for (i=begin;i<end;++i) {
    accum += pos_error[i];
  }
  return accum;
}


GT_INLINE void gt_stats_print_mmap_distribution(FILE* stream,uint64_t* const mmap,const uint64_t num_alignments,const uint64_t num_mapped) {
#define GT_STATS_PRINT_MMAP_FORMAT "%6" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_MMAP(RANGE) mmap[RANGE],100.0*(float)mmap[RANGE]/(float)num_alignments
  if(!num_alignments) return;
  fprintf(stream,"MMap.ranges\n");
  fprintf(stream,"  -->        [0] \t=> "GT_STATS_PRINT_MMAP_FORMAT,(num_alignments-num_mapped),100.0*(float)(num_alignments-num_mapped)/(float)num_alignments);
  fprintf(stream,"  -->        [1] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_1));
  fprintf(stream,"  -->      (1,5] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_5));
  fprintf(stream,"  -->     (5,10] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_10));
  fprintf(stream,"  -->    (10,50] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_50));
  fprintf(stream,"  -->   (50,100] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_100));
  fprintf(stream,"  -->  (100,500] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_500));
  fprintf(stream,"  --> (500,1000] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_1000));
  fprintf(stream,"  --> (1000,inf) \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(GT_STATS_MMAP_RANGE_BEHOND));
}
GT_INLINE void gt_stats_print_uniq_distribution(FILE* stream,uint64_t* const uniq,const uint64_t num_alignments) {
#define GT_STATS_PRINT_UNIQ_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_UNIQ(RANGE) uniq[RANGE],100.0*(float)uniq[RANGE]/(float)num_alignments
  if(!num_alignments) return;
  fprintf(stream,"Uniq.ranges\n");
  fprintf(stream,"  -->        [X] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_X));
  fprintf(stream,"  -->        [0] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_0));
  fprintf(stream,"  -->        [1] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_1));
  fprintf(stream,"  -->        [2] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_2));
  fprintf(stream,"  -->        [3] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_3));
  fprintf(stream,"  -->     (3,10] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_10));
  fprintf(stream,"  -->    (10,50] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(GT_STATS_UNIQ_RANGE_50));
  const uint64_t accum = uniq[GT_STATS_UNIQ_RANGE_100]+uniq[GT_STATS_UNIQ_RANGE_500]+uniq[GT_STATS_UNIQ_RANGE_BEHOND];
  fprintf(stream,"  -->   (50,inf) \t=> "GT_STATS_PRINT_UNIQ_FORMAT,accum,100.0*(float)(accum)/(float)num_alignments);
}
GT_INLINE void gt_stats_print_inss_distribution(FILE* stream,uint64_t* const inss,const uint64_t num_maps) {
#define GT_STATS_PRINT_INSS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_INSS(RANGE) inss[RANGE],100.0*(float)inss[RANGE]/(float)num_maps
  if(!num_maps) return;
  fprintf(stream,"InsS.ranges\n");
  fprintf(stream,"  -->   (-inf,-100] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_NEG));
  fprintf(stream,"  -->      (-100,0] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_OVER));
  fprintf(stream,"  -->       (0,100] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_100));
  fprintf(stream,"  -->     (100,200] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_200));
  fprintf(stream,"  -->     (200,300] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_300));
  fprintf(stream,"  -->     (300,400] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_400));
  fprintf(stream,"  -->     (400,500] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_500));
  fprintf(stream,"  -->     (500,600] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_600));
  fprintf(stream,"  -->     (600,700] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_700));
  fprintf(stream,"  -->     (700,800] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_800));
  fprintf(stream,"  -->     (800,900] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_900));
  fprintf(stream,"  -->    (900,1000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_1000));
  fprintf(stream,"  -->   (1000,2000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_2000));
  fprintf(stream,"  -->   (2000,5000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_5000));
  fprintf(stream,"  -->  (5000,10000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_10000));
  fprintf(stream,"  -->   (10000,inf) \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(GT_STATS_INSS_RANGE_BEHOND));
}
GT_INLINE void gt_stats_print_inss_fg_distribution(FILE* stream,uint64_t* const inss,const uint64_t num_maps) {
#define GT_STATS_PRINT_INSS_FG_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_INSS_FG(RANGE) inss[RANGE],100.0*(float)inss[RANGE]/(float)num_maps
  if(!num_maps) return;
  fprintf(stream,"InsS.ranges\n");
  int64_t current_bucket, current_inf=GT_STATS_INSS_FG_MIN;
  for (current_bucket=0;current_bucket<GT_STATS_INSS_FG_RANGE;++current_bucket) {
    if (current_inf==GT_STATS_INSS_FG_MIN) {
      fprintf(stream,"  -->   (-inf,%5"PRIu64") \t=> "GT_STATS_PRINT_INSS_FG_FORMAT,
          current_inf+GT_STATS_INSS_FG_STEP,GT_STATS_PRINT_INSS_FG(0));
    } else if (current_bucket==GT_STATS_INSS_FG_RANGE-1) {
      fprintf(stream,"  -->   [%5"PRIu64",+inf) \t=> "GT_STATS_PRINT_INSS_FG_FORMAT,
          current_inf,GT_STATS_PRINT_INSS_FG(current_bucket));
    } else {
      fprintf(stream,"  -->   [%5"PRIu64",%5"PRIu64") \t=> "GT_STATS_PRINT_INSS_FG_FORMAT,
          current_inf,current_inf+GT_STATS_INSS_FG_STEP,GT_STATS_PRINT_INSS_FG(current_bucket));
    }
    current_inf+=GT_STATS_INSS_FG_STEP;
  }
}
GT_INLINE void gt_stats_print_error_event_distribution(FILE* stream,uint64_t* const error,const uint64_t num_maps) {
#define GT_STATS_PRINT_MISMS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_MISMS(RANGE) error[RANGE],100.0*(float)error[RANGE]/(float)num_maps
#define GT_STATS_PRINT_MISMS_VAL(misms_val) misms_val,100.0*(float)misms_val/(float)num_maps
  if(!num_maps) return;
  fprintf(stream,"  -->         [0]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_0));
  fprintf(stream,"  -->         [1]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_1));
  fprintf(stream,"  -->         [2]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_2));
  const uint64_t misms_accum = error[GT_STATS_MISMS_RANGE_3]+error[GT_STATS_MISMS_RANGE_4]+
      error[GT_STATS_MISMS_RANGE_5]+error[GT_STATS_MISMS_RANGE_6]+error[GT_STATS_MISMS_RANGE_7]+
      error[GT_STATS_MISMS_RANGE_8]+error[GT_STATS_MISMS_RANGE_9]+error[GT_STATS_MISMS_RANGE_10];
  fprintf(stream,"  -->      (2,10]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS_VAL(misms_accum));
  fprintf(stream,"  -->     (10,20]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_20));
  fprintf(stream,"  -->     (20,50]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_50));
  fprintf(stream,"  -->    (50,100]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(GT_STATS_MISMS_RANGE_BEHOND));
}
GT_INLINE void gt_stats_print_read_event_positions(
    FILE* stream,uint64_t* const pos_error,uint64_t const num_errors,uint64_t const max_length) {
  const uint64_t max_length_stats = GT_MIN(max_length,GT_STATS_LARGE_READ_POS_RANGE);
  const uint64_t ten_per_cent = max_length_stats/10;
  uint64_t i, centinel, error_sum;
  if(!num_errors) return;
  fprintf(stream,"Error.position [0,%" PRIu64 "]nt\n",max_length_stats);
  for (i=0,centinel=0;i<100;i+=10,centinel+=ten_per_cent) {
    const uint64_t top = (i==90)?max_length_stats:centinel+ten_per_cent;
    error_sum = gt_stats_sum_misms_pos(pos_error,centinel,top);
    fprintf(stream,"  -->   [%3" PRIu64 " - %3" PRIu64 ")nt \t=> %1.3f%%\n",
        centinel,top,100.0*(float)error_sum/(float)num_errors);
  }
}
GT_INLINE void gt_stats_print_num_junctions_distribution(FILE* stream,uint64_t* const num_junctions,uint64_t const total) {
#define GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_NUM_JUNCTIONS(RANGE) num_junctions[RANGE],100.0*(float)num_junctions[RANGE]/(float)total
  fprintf(stream,"  -->           [1] \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(GT_STATS_NUM_JUNCTION_1));
  fprintf(stream,"  -->           [2] \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(GT_STATS_NUM_JUNCTION_2));
  fprintf(stream,"  -->           [3] \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(GT_STATS_NUM_JUNCTION_3));
  fprintf(stream,"  -->       (3,inf) \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(GT_STATS_NUM_JUNCTION_BEHOND));
}
GT_INLINE void gt_stats_print_length_junctions_distribution(FILE* stream,uint64_t* const length_junctions,uint64_t const total_junctions) {
#define GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_LENGTH_JUNCTION(RANGE) length_junctions[RANGE],100.0*(float)length_junctions[RANGE]/(float)total_junctions
  fprintf(stream,"SM.Length.Junctions\n");
  fprintf(stream,"  -->       [0,100] \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_100));
  fprintf(stream,"  -->    (100,1000] \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_1000));
  fprintf(stream,"  -->   (1000,5000] \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_5000));
  fprintf(stream,"  -->  (5000,10000] \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_10000));
  fprintf(stream,"  --> (10000,50000] \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_50000));
  fprintf(stream,"  -->   (50000,inf) \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(GT_STATS_LEN_JUNCTION_BEHOND));
}
GT_INLINE void gt_stats_print_junction_position_distribution(FILE* stream,uint64_t* const junction_position,uint64_t const total_junctions,uint64_t const max_length) {
  const uint64_t max_length_stats = GT_MIN(max_length,GT_STATS_SHORT_READ_POS_RANGE);
  const uint64_t ten_per_cent = max_length_stats/10;
  uint64_t i, centinel;
  if(!total_junctions) return;
  fprintf(stream,"Juntion.position [0,%" PRIu64 "]nt\n",max_length_stats);
  for (i=0,centinel=0;i<100;i+=10,centinel+=ten_per_cent) {
    const uint64_t top = (i==90)?max_length_stats:centinel+ten_per_cent;
    const uint64_t num_junctions = gt_stats_sum_misms_pos(junction_position,centinel,top);
    fprintf(stream,"  -->   [%3" PRIu64 " - %3" PRIu64")nt \t=> %1.3f%%\n",centinel,top,100.0*(float)num_junctions/(float)total_junctions);
  }
}
GT_INLINE void gt_stats_print_qualities_error_distribution(FILE* stream,uint64_t* const qualities_error,uint64_t const total_error) {
  uint64_t i;
  uint64_t j, qual;
  if(!total_error) return;
  // Print Header
  fprintf(stream,"    ");
  for (j=0;j<16;++j) fprintf(stream,"[%4"PRIu64"]",j);
  fprintf(stream,"\n");
  // Print Values
  for (i=0,qual=0;i<16;++i) {
    fprintf(stream,"%3"PRIu64" ",qual);
    for (j=0;j<16;++j,++qual) {
      fprintf(stream,"[%4.1f]",100.0*((double)qualities_error[qual]/(double)total_error));
    }
    fprintf(stream,"\n");
  }
}
GT_INLINE void gt_stats_print_misms_transition_table(FILE* stream,uint64_t* const misms_trans,uint64_t const total_misms) {
  uint64_t i, pos=0;
  if(!total_misms) return;
  // Print Header
  fprintf(stream,"    [  A  ][  C  ][  G  ][  T  ][  N  ]\n");
  fprintf(stream,"[A] ");
  for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stream,"\n");
  fprintf(stream,"[C] ");
  for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stream,"\n");
  fprintf(stream,"[G] ");
  for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stream,"\n");
  fprintf(stream,"[T] ");
  for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stream,"\n");
  fprintf(stream,"[N] ");
  for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stream,"\n");
}
#define GT_STATS_GET_IXD_TRANSITION_1_CONTEXT(a,b,c,i) ((((a*GT_STATS_MISMS_BASE_RANGE+b)*GT_STATS_MISMS_BASE_RANGE)+c)*GT_STATS_MISMS_BASE_RANGE+i)
GT_INLINE void gt_stats_print_misms_transition_table_1context(FILE* stream,uint64_t* const misms_trans,uint64_t const total_misms) {
  if(!total_misms) return;
  // Print Header
  fprintf(stream,"      [  A  ][  C  ][  G  ][  T  ][  N  ]\n");
  char bases[] = {'A','C','G','T','N'};
  uint64_t a,b,c;
  for (b=0;b<4;++b) {
    for (a=0;a<4;++a) {
      for (c=0;c<4;++c) {
        fprintf(stream,"[%c%c%c] ",bases[a],bases[b],bases[c]);
        uint64_t i;
        for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i) {
          fprintf(stream,"[%5.2f]",100.0*(double)misms_trans[GT_STATS_GET_IXD_TRANSITION_1_CONTEXT(a,b,c,i)]/(double)total_misms);
        }
        fprintf(stream,"\n");
      }
    }
  }
}

GT_INLINE void gt_stats_print_split_maps_stats(FILE* stream,gt_stats* const stats,const bool paired_end) {
  gt_splitmaps_profile* const splitmap_stats = stats->splitmaps_profile;
  if (splitmap_stats->total_splitmaps==0) {
    fprintf(stream,"SM.Total \t 0\n");
  } else {
    fprintf(stream,"SM.Num.mapped.withSM \t %" PRIu64 " (%2.3f%%)\n",splitmap_stats->num_mapped_with_splitmaps,
        GT_GET_PERCENTAGE(splitmap_stats->num_mapped_with_splitmaps,stats->num_mapped));
    fprintf(stream,"SM.Num.mapped.onlyBySM \t %" PRIu64 " (%2.3f%%)\n",splitmap_stats->num_mapped_only_splitmaps,
        GT_GET_PERCENTAGE(splitmap_stats->num_mapped_only_splitmaps,stats->num_mapped));
    fprintf(stream,"SM.Num.Splitted.Segments \t %" PRIu64 " (%2.3f SplitMaps/Maps)\n",splitmap_stats->total_splitmaps,
        GT_GET_PERCENTAGE(splitmap_stats->total_splitmaps,stats->num_blocks));
    fprintf(stream,"SM.Num.Junctions \t %" PRIu64 " (%2.3f Juntions/SplitMaps)\n",splitmap_stats->total_junctions,
        GT_DIV_F(splitmap_stats->total_junctions,splitmap_stats->total_splitmaps));
    gt_stats_print_num_junctions_distribution(stream,splitmap_stats->num_junctions,splitmap_stats->total_splitmaps);
    gt_stats_print_length_junctions_distribution(stream,splitmap_stats->length_junctions,splitmap_stats->total_junctions);
    gt_stats_print_junction_position_distribution(stream,splitmap_stats->junction_position,splitmap_stats->total_junctions,stats->max_length);
    if (paired_end) {
      fprintf(stream,"SM.PE.Combinations\n");
      fprintf(stream,"  --> SM+SM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_sm_sm,100.0*(float)splitmap_stats->pe_sm_sm/(float)stats->num_maps);
      fprintf(stream,"  --> SM+RM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_sm_rm,100.0*(float)splitmap_stats->pe_sm_rm/(float)stats->num_maps);
      fprintf(stream,"  --> RM+RM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_rm_rm,100.0*(float)splitmap_stats->pe_rm_rm/(float)stats->num_maps);
    }
  }
}
GT_INLINE void gt_stats_print_maps_stats(FILE* stream, gt_stats* const stats,const uint64_t num_reads,const bool paired_end) {
  const uint64_t num_templates = paired_end ? num_reads>>1 : num_reads; // SE => 1 template. PE => 1 template
  // Total bases (aligned/trimmed/unaligned)
  const gt_maps_profile* const maps_profile = stats->maps_profile;
  fprintf(stream,"  --> Total.Bases %" PRIu64 " (%2.3f per map) \n",
      maps_profile->total_bases,GT_DIV_F(maps_profile->total_bases,stats->num_maps));
  fprintf(stream,"    --> Bases.Matching %" PRIu64 " (%2.3f) \n",
      maps_profile->total_bases_matching,GT_GET_PERCENTAGE(maps_profile->total_bases_matching,maps_profile->total_bases));
  fprintf(stream,"    --> Bases.Trimmed %" PRIu64 " (%2.3f) \n",
      maps_profile->total_bases_trimmed,GT_GET_PERCENTAGE(maps_profile->total_bases_trimmed,maps_profile->total_bases));
  gt_stats_print_mmap_distribution(stream,stats->mmap,num_templates,stats->num_mapped);
  gt_stats_print_uniq_distribution(stream,stats->uniq,num_templates);
  if (paired_end) {
    fprintf(stream,"  --> Strand.combinations \n");
    fprintf(stream,"    --> F+R %" PRIu64 " (%2.3f%%) \n",maps_profile->pair_strand_fr,GT_GET_PERCENTAGE(maps_profile->pair_strand_rf,stats->num_maps));
    fprintf(stream,"    --> R+F %" PRIu64 " (%2.3f%%) \n",maps_profile->pair_strand_rf,GT_GET_PERCENTAGE(maps_profile->pair_strand_fr,stats->num_maps));
    fprintf(stream,"    --> F+F %" PRIu64 " (%2.3f%%) \n",maps_profile->pair_strand_ff,GT_GET_PERCENTAGE(maps_profile->pair_strand_ff,stats->num_maps));
    fprintf(stream,"    --> R+R %" PRIu64 " (%2.3f%%) \n",maps_profile->pair_strand_rr,GT_GET_PERCENTAGE(maps_profile->pair_strand_rr,stats->num_maps));
    gt_stats_print_inss_fg_distribution(stream,stats->maps_profile->inss_fine_grain,stats->num_maps);
  }
  fprintf(stream,"[ERROR.PROFILE]\n");
  fprintf(stream,"  --> Total.Mismatches %" PRIu64 " (%2.3f per map) \n",
      maps_profile->total_mismatches,GT_DIV_F(maps_profile->total_mismatches,stats->num_maps));
  fprintf(stream,"  --> Total.Errors %" PRIu64 " (%2.3f per map) \n",
      maps_profile->total_errors_events,GT_DIV_F(maps_profile->total_errors_events,stats->num_maps));
  fprintf(stream,"  --> Total.Indels.Length %" PRIu64 " (%2.3f per map) \n",
      maps_profile->total_indel_length,GT_DIV_F(maps_profile->total_indel_length,stats->num_maps));
  fprintf(stream,"  --> Total.Levenshtein %" PRIu64 " (%2.3f per map) \n",
      maps_profile->total_levenshtein,GT_DIV_F(maps_profile->total_levenshtein,stats->num_maps));
  if (stats->num_maps > 0) {
    fprintf(stream,"Mismatches.Distribution\n");
    gt_stats_print_error_event_distribution(stream,stats->maps_profile->mismatches,stats->num_maps);
    fprintf(stream,"Ins.Length.Distribution\n");
    gt_stats_print_error_event_distribution(stream,stats->maps_profile->insertion_length,stats->num_maps);
    fprintf(stream,"Del.Length.Distribution\n");
    gt_stats_print_error_event_distribution(stream,stats->maps_profile->deletion_length,stats->num_maps);
    fprintf(stream,"Error.Events.Distribution\n");
    gt_stats_print_error_event_distribution(stream,stats->maps_profile->errors_events,stats->num_maps);
  }
  gt_stats_print_read_event_positions(stream,stats->maps_profile->error_position,stats->maps_profile->total_errors_events,stats->max_length);
}
GT_INLINE void gt_stats_print_general_stats(FILE* stream,gt_stats* const stats,const uint64_t num_reads,const bool paired_end) {
  const uint64_t num_templates = paired_end ? num_reads>>1 : num_reads; // SE => 1 template. PE => 1 template
  // For the case of zero input lines
  if(stats->min_length>stats->max_length) stats->min_length=0;
  if(stats->mapped_min_length>stats->mapped_max_length) stats->mapped_min_length=0;
  // Print
  fprintf(stream,"  --> Num.Reads %" PRIu64 "\n",num_reads);
  fprintf(stream,"    --> Reads.Length (min,avg,max) (%" PRIu64 ",%" PRIu64 ",%" PRIu64 ")\n",
      stats->min_length,GT_DIV(stats->total_bases,stats->num_blocks),stats->max_length);
  fprintf(stream,"  --> Reads.Mapped %" PRIu64 " (%2.3f%%)\n",stats->num_mapped,
      num_templates?100.0*(float)stats->num_mapped/(float)num_templates:0.0);
  fprintf(stream,"    --> Reads.Mapped.Length (min,avg,max) (%" PRIu64 ",%" PRIu64 ",%" PRIu64 ")\n",
      stats->mapped_min_length,GT_DIV(stats->total_bases_aligned,(paired_end)?stats->num_mapped*2:stats->num_mapped),stats->mapped_max_length);
  fprintf(stream,"  --> Num.Bases %" PRIu64 "\n",stats->total_bases);
  fprintf(stream,"    --> Bases.Prop [A]=%2.3f%% [C]=%2.3f%% [G]=%2.3f%% [T]=%2.3f%% [N]=%2.3f%%\n",
      100.0*(float)stats->nt_counting[0]/(float)stats->total_bases,100.0*(float)stats->nt_counting[1]/(float)stats->total_bases,
      100.0*(float)stats->nt_counting[2]/(float)stats->total_bases,100.0*(float)stats->nt_counting[3]/(float)stats->total_bases,
      100.0*(float)stats->nt_counting[4]/(float)stats->total_bases);
  fprintf(stream,"    --> Bases.Aligned %" PRIu64 " (%2.3f%%)\n",stats->total_bases_aligned,
      GT_GET_PERCENTAGE(stats->total_bases_aligned,stats->total_bases));
  fprintf(stream,"  --> Num.alignments %" PRIu64 "\n",stats->num_alignments);
  fprintf(stream,"  --> Num.maps %" PRIu64 " (%2.3f map/alg)\n",stats->num_maps,GT_DIV_F(stats->num_maps,stats->num_mapped));
}
