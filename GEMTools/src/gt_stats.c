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
  register uint64_t iii; \
  for (iii=0;iii<RANGE;++iii) {VECTOR_A[iii] += VECTOR_B[iii];} \
}

/*
 * MAPS Error Profile
 */
GT_INLINE gt_maps_error_profile* gt_maps_error_profile_new() {
  // Allocate handler
  gt_maps_error_profile* maps_ep = malloc(sizeof(gt_maps_error_profile));
  gt_cond_fatal_error(!maps_ep,MEM_HANDLER);
  /*
   * Init
   */
  // Mismatch/Indel Profile
  maps_ep->mismatches = calloc(GT_STATS_MISMS_RANGE,sizeof(uint64_t));
  maps_ep->levenshtein = calloc(GT_STATS_MISMS_RANGE,sizeof(uint64_t));
  maps_ep->indel_length = calloc(2*GT_STATS_MISMS_RANGE,sizeof(uint64_t));
  maps_ep->errors_events = calloc(GT_STATS_MISMS_RANGE,sizeof(uint64_t));
  // Mismatch/Indel Distribution
  maps_ep->total_mismatches=0;
  maps_ep->total_levenshtein=0;
  maps_ep->total_indel_length=0;
  maps_ep->total_errors_events=0;
  maps_ep->error_position = calloc(GT_STATS_LARGE_READ_POS_RANGE,sizeof(uint64_t));
  // Trim/Mapping stats
  maps_ep->total_bases=0;
  maps_ep->total_bases_matching=0;
  maps_ep->total_bases_trimmed=0;
  // Insert Size Distribution
  maps_ep->inss = calloc(GT_STATS_INSS_RANGE,sizeof(uint64_t));
  // Mismatch/Errors bases
  maps_ep->misms_transition = calloc(GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE,sizeof(uint64_t));
  maps_ep->qual_score_misms = calloc(GT_STATS_QUAL_SCORE_RANGE,sizeof(uint64_t));
  maps_ep->misms_1context = calloc(GT_STATS_MISMS_1_CONTEXT_RANGE,sizeof(uint64_t));
  maps_ep->misms_2context = calloc(GT_STATS_MISMS_2_CONTEXT_RANGE,sizeof(uint64_t));
  maps_ep->indel_transition_1 = calloc(GT_STATS_INDEL_TRANSITION_1_RANGE,sizeof(uint64_t));
  maps_ep->indel_transition_2 = calloc(GT_STATS_INDEL_TRANSITION_2_RANGE,sizeof(uint64_t));
  maps_ep->indel_transition_3 = calloc(GT_STATS_INDEL_TRANSITION_3_RANGE,sizeof(uint64_t));
  maps_ep->indel_transition_4 = calloc(GT_STATS_INDEL_TRANSITION_4_RANGE,sizeof(uint64_t));
  maps_ep->indel_1context = calloc(GT_STATS_INDEL_1_CONTEXT,sizeof(uint64_t));
  maps_ep->indel_2context = calloc(GT_STATS_INDEL_2_CONTEXT,sizeof(uint64_t));
  maps_ep->qual_score_errors = calloc(GT_STATS_QUAL_SCORE_RANGE,sizeof(uint64_t));
  return maps_ep;
}
GT_INLINE void gt_maps_error_profile_clear(gt_maps_error_profile* const maps_ep) {
  /*
   * Init
   */
  // Mismatch/Indel Profile
  memset(maps_ep->mismatches,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->levenshtein,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_length,0,2*GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->errors_events,0,GT_STATS_MISMS_RANGE*sizeof(uint64_t));
  // Mismatch/Indel Distribution
  maps_ep->total_mismatches=0;
  maps_ep->total_levenshtein=0;
  maps_ep->total_indel_length=0;
  maps_ep->total_errors_events=0;
  memset(maps_ep->error_position,0,GT_STATS_LARGE_READ_POS_RANGE*sizeof(uint64_t));
  // Trim/Mapping stats
  maps_ep->total_bases=0;
  maps_ep->total_bases_matching=0;
  maps_ep->total_bases_trimmed=0;
  // Insert Size Distribution
  memset(maps_ep->inss,0,GT_STATS_INSS_RANGE*sizeof(uint64_t));
  // Mismatch/Errors bases
  memset(maps_ep->misms_transition,0,GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*sizeof(uint64_t));
  memset(maps_ep->qual_score_misms,0,GT_STATS_QUAL_SCORE_RANGE*sizeof(uint64_t));
  memset(maps_ep->misms_1context,0,GT_STATS_MISMS_1_CONTEXT_RANGE*sizeof(uint64_t));
  memset(maps_ep->misms_2context,0,GT_STATS_MISMS_2_CONTEXT_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_transition_1,0,GT_STATS_INDEL_TRANSITION_1_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_transition_2,0,GT_STATS_INDEL_TRANSITION_2_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_transition_3,0,GT_STATS_INDEL_TRANSITION_3_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_transition_4,0,GT_STATS_INDEL_TRANSITION_4_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_1context,0,GT_STATS_INDEL_1_CONTEXT*sizeof(uint64_t));
  memset(maps_ep->indel_2context,0,GT_STATS_INDEL_2_CONTEXT*sizeof(uint64_t));
  memset(maps_ep->qual_score_errors,0,GT_STATS_QUAL_SCORE_RANGE*sizeof(uint64_t));

}
GT_INLINE void gt_maps_error_profile_delete(gt_maps_error_profile* const maps_ep) {
  // Mismatch/Indel Profile
  free(maps_ep->mismatches);
  free(maps_ep->levenshtein);
  free(maps_ep->indel_length);
  free(maps_ep->errors_events);
  // Mismatch/Indel Distribution
  free(maps_ep->error_position);
  // Insert Size Distribution
  free(maps_ep->inss);
  // Mismatch/Errors bases
  free(maps_ep->misms_transition);
  free(maps_ep->qual_score_misms);
  free(maps_ep->misms_1context);
  free(maps_ep->misms_2context);
  free(maps_ep->indel_transition_1);
  free(maps_ep->indel_transition_2);
  free(maps_ep->indel_transition_3);
  free(maps_ep->indel_transition_4);
  free(maps_ep->indel_1context);
  free(maps_ep->indel_2context);
  free(maps_ep->qual_score_errors);
  free(maps_ep);
}
GT_INLINE void gt_maps_error_profile_merge(
    gt_maps_error_profile* const maps_ep_dst,gt_maps_error_profile* const maps_ep_src) {
  // GT_STATS_VECTOR_ADD(maps_ep_dst->,maps_ep_src->,);
  // Mismatch/Indel Profile
  GT_STATS_VECTOR_ADD(maps_ep_dst->mismatches,maps_ep_src->mismatches,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->levenshtein,maps_ep_src->levenshtein,GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_length,maps_ep_src->indel_length,2*GT_STATS_MISMS_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->errors_events,maps_ep_src->errors_events,GT_STATS_MISMS_RANGE);
  // Mismatch/Indel Distribution
  maps_ep_dst->total_mismatches+=maps_ep_src->total_mismatches;
  maps_ep_dst->total_levenshtein+=maps_ep_src->total_levenshtein;
  maps_ep_dst->total_indel_length+=maps_ep_src->total_indel_length;
  maps_ep_dst->total_errors_events+=maps_ep_src->total_errors_events;
  GT_STATS_VECTOR_ADD(maps_ep_dst->error_position,maps_ep_src->error_position,GT_STATS_LARGE_READ_POS_RANGE);
  // Trim/Mapping stats
  maps_ep_dst->total_bases+=maps_ep_src->total_bases;
  maps_ep_dst->total_bases_matching+=maps_ep_src->total_bases_matching;
  maps_ep_dst->total_bases_trimmed+=maps_ep_src->total_bases_trimmed;
  // Insert Size Distribution
  GT_STATS_VECTOR_ADD(maps_ep_dst->inss,maps_ep_src->inss,GT_STATS_INSS_RANGE);
  // Mismatch/Errors bases
  GT_STATS_VECTOR_ADD(maps_ep_dst->misms_transition,maps_ep_src->misms_transition,GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->qual_score_misms,maps_ep_src->qual_score_misms,GT_STATS_QUAL_SCORE_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->misms_1context,maps_ep_src->misms_1context,GT_STATS_MISMS_1_CONTEXT_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->misms_2context,maps_ep_src->misms_2context,GT_STATS_MISMS_2_CONTEXT_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_transition_1,maps_ep_src->indel_transition_1,GT_STATS_INDEL_TRANSITION_1_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_transition_2,maps_ep_src->indel_transition_2,GT_STATS_INDEL_TRANSITION_2_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_transition_3,maps_ep_src->indel_transition_3,GT_STATS_INDEL_TRANSITION_3_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_transition_4,maps_ep_src->indel_transition_4,GT_STATS_INDEL_TRANSITION_4_RANGE);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_1context,maps_ep_src->indel_1context,GT_STATS_INDEL_1_CONTEXT);
  GT_STATS_VECTOR_ADD(maps_ep_dst->indel_2context,maps_ep_src->indel_2context,GT_STATS_INDEL_2_CONTEXT);
  GT_STATS_VECTOR_ADD(maps_ep_dst->qual_score_errors,maps_ep_src->qual_score_errors,GT_STATS_QUAL_SCORE_RANGE);
}

/*
 * SPLITMAPS Profile
 */
GT_INLINE gt_splitmaps_profile* gt_splitmaps_profile_new() {
  // Allocate handler
  gt_splitmaps_profile* splitmaps_profile = malloc(sizeof(gt_splitmaps_profile));
  gt_cond_fatal_error(!splitmaps_profile,MEM_HANDLER);
  /*
   * Init
   */
  // General SM
  splitmaps_profile->num_mapped_with_splitmaps = 0;
  splitmaps_profile->num_mapped_only_splitmaps = 0;
  splitmaps_profile->total_splitmaps = 0;
  splitmaps_profile->total_junctions = 0;
  splitmaps_profile->num_junctions = calloc(GT_STATS_NUM_JUNCTION_RANGE,sizeof(uint64_t));
  splitmaps_profile->length_junctions = calloc(GT_STATS_LEN_JUNCTION_RANGE,sizeof(uint64_t));
  splitmaps_profile->junction_position = calloc(GT_STATS_SHORT_READ_POS_RANGE,sizeof(uint64_t));
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
  free(splitmaps_profile->num_junctions);
  free(splitmaps_profile->length_junctions);
  free(splitmaps_profile->junction_position);
  free(splitmaps_profile);
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
  register const double one_per_cent = read_length/100.0;
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
  } else {
    ++inss[GT_STATS_INSS_RANGE_BEHOND];
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
  gt_stats* stats = malloc(sizeof(gt_stats));
  gt_cond_fatal_error(!stats,MEM_HANDLER);
  // Length
  stats->min_length=UINT64_MAX;
  stats->max_length=0;
  stats->total_bases=0;
  stats->total_bases_aligned=0;
  stats->mapped_min_length=UINT64_MAX;
  stats->mapped_max_length=0;
  // Nucleotide counting (wrt to the maps=read+errors)
  stats->nt_counting = calloc(GT_STATS_MISMS_BASE_RANGE,sizeof(uint64_t));
  // Mapped/Maps/MMaps/Uniq...
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  stats->mmap = calloc(GT_STATS_MMAP_RANGE,sizeof(uint64_t)); // MMaps
  stats->uniq = calloc(GT_STATS_UNIQ_RANGE,sizeof(uint64_t)); // Uniq
  // Maps Error Profile
  stats->maps_error_profile = gt_maps_error_profile_new();
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
  gt_maps_error_profile_clear(stats->maps_error_profile);
  // Split maps Profile
  gt_splitmaps_profile_clear(stats->splitmaps_profile);
}
GT_INLINE void gt_stats_delete(gt_stats* const stats) {
  free(stats->mmap);
  free(stats->uniq);
  gt_maps_error_profile_delete(stats->maps_error_profile);
  gt_splitmaps_profile_delete(stats->splitmaps_profile);
  free(stats);
}

/*
 * STATS Merge
 */
GT_INLINE void gt_stats_merge(gt_stats** const stats,const uint64_t stats_array_size) {
  register uint64_t i;
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
    gt_maps_error_profile_merge(stats[0]->maps_error_profile,stats[i]->maps_error_profile);
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
    gt_maps_error_profile *maps_error_profile,
    gt_template* const template,const uint64_t alignment_total_length) {

}
GT_INLINE void gt_stats_make_maps_error_profile(
    gt_maps_error_profile *maps_error_profile,gt_template* const template,
    const uint64_t alignment_total_length,gt_map** const mmap) {
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t total_mismatches=0;
  uint64_t total_levenshtein=0;
  uint64_t total_ins_length=0, total_del_length=0;
  uint64_t total_errors_events=0;
  uint64_t total_bases=0, total_bases_trimmed=0, total_bases_not_matching=0;
  // Iterate all MMAPS
  GT_MULTIMAP_ITERATE_BLOCKS(mmap,num_blocks,map,end_pos) {
    register gt_alignment* const alignment = gt_template_get_block(template,end_pos);
    register gt_string* const read = alignment->read;
    register gt_string* const quals = alignment->qualities;
    register const bool has_qualities = gt_alignment_has_qualities(alignment);
    register char quality_misms = 0;
    GT_MAP_ITERATE(map,map_block) {
      total_bases += map->base_length;
      GT_MISMS_ITERATE(map_block,misms) {
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
            register uint64_t idx = 0;
            idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
            idx *= GT_STATS_MISMS_RANGE;
            idx += gt_cdna_encode[(uint8_t)misms->base];
            maps_error_profile->misms_transition[idx]++;
            if (misms->position>0 && misms->position<map->base_length-1) { // 1-context
              idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-1]];
              idx *= GT_STATS_MISMS_RANGE;
              idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
              idx *= GT_STATS_MISMS_RANGE;
              idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+1]];
              idx *= GT_STATS_MISMS_RANGE;
              idx += gt_cdna_encode[(uint8_t)misms->base];
              maps_error_profile->misms_1context[idx]++;
              if (misms->position>1 && misms->position<map->base_length-2) { // 2-context
                idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-2]];
                idx *= GT_STATS_MISMS_RANGE;
                idx  = gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position-1]];
                idx *= GT_STATS_MISMS_RANGE;
                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]];
                idx *= GT_STATS_MISMS_RANGE;
                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+1]];
                idx *= GT_STATS_MISMS_RANGE;
                idx += gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position+2]];
                idx *= GT_STATS_MISMS_RANGE;
                idx += gt_cdna_encode[(uint8_t)misms->base];
                maps_error_profile->misms_2context[idx]++;
              }
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
            if (misms->position==0 || misms->position+misms->size==map_block->base_length) {
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
  gt_stats_get_misms(maps_error_profile->indel_length,total_ins_length,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->indel_length+GT_STATS_MISMS_RANGE,total_del_length,alignment_total_length);
  gt_stats_get_misms(maps_error_profile->errors_events,total_errors_events,alignment_total_length);
}

GT_INLINE void gt_stats_make_mmaps_profile(
    gt_stats* const stats,gt_template* const template,
    const uint64_t alignment_total_length,gt_stats_analysis* const stats_analysis) {
  // Check not null maps
  if (gt_template_get_num_mmaps(template)==0) return;
  register const uint64_t num_blocks_template = gt_template_get_num_blocks(template);
  register const uint64_t paired_map = (num_blocks_template==2);
  // Iterate over all/best maps
  register gt_maps_error_profile* const maps_error_profile = stats->maps_error_profile;
  register gt_splitmaps_profile* const splitmaps_profile = stats->splitmaps_profile;
  register bool has_splitsmaps = false;
  register bool only_splitsmaps = true;
  GT_TEMPLATE_ITERATE(template,mmap) {
    /*
     * Insert Size Distribution
     */
    if (paired_map) {
      uint64_t gt_err;
      int64_t ins_size = gt_template_get_insert_size(mmap,&gt_err);
      if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK) gt_stats_get_inss_distribution(maps_error_profile->inss,ins_size);
    }
    /*
     * Error Profile
     */
    if (stats_analysis->error_profile) {
      gt_stats_make_maps_error_profile(maps_error_profile,template,alignment_total_length,mmap);
    }
    /*
     * SplitMap Stats
     */
    if (stats_analysis->split_map_stats) {
      // SM block stats
      bool has_sm[2] = {true, true};
      GT_MULTIMAP_ITERATE(mmap,map,end_pos) {
        register const uint64_t num_blocks = gt_map_get_num_blocks(map);
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
            gt_stats_get_juntions_length_distribution(splitmaps_profile->length_junctions,gt_map_get_junction_distance(map_block));
            register const uint64_t juntion_position = gt_map_get_base_length(map_block);
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
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t num_maps = gt_template_get_num_mmaps(template);
  if (stats_analysis->best_map && num_maps>0) num_maps = 1; // NumMaps correction
  register const bool is_mapped = (num_maps>0 || gt_template_is_mapped(template));
  /*
   * Length STATS
   */
  register uint64_t alignment_total_length = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    // Read length stats
    register uint64_t read_length =
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
      gt_stats_make_indel_profile(stats->maps_error_profile,template,alignment_total_length);
    }
  }
}
