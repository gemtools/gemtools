/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats.c
 * DATE: 10/12/2012
 * DESCRIPTION: // TODO
 */

#include "gt_stats.h"

/*
 * Handy Macros
 */
#define GT_STATS_GET_PERCENTAGE_ERROR(percentage,one_per_cent) (uint64_t)((double)percentage*one_per_cent)

// FIXME: Erase to cdna
const uint8_t gt_cdna_encode[256] = {
    [0 ... 255] = MISMS_BASE_N,
    ['A'] = MISMS_BASE_A,['C'] = MISMS_BASE_C,['G'] = MISMS_BASE_G,['T'] = MISMS_BASE_T,
    ['a'] = MISMS_BASE_A,['c'] = MISMS_BASE_C,['g'] = MISMS_BASE_G,['t'] = MISMS_BASE_T,
};

/*
 * MAPS Error Profile
 */
GT_INLINE gt_maps_error_profile* gt_maps_error_profile_new() {
  // Allocate handler
  gt_maps_error_profile* maps_ep = malloc(sizeof(gt_maps_error_profile));
  gt_cond_fatal_error(!maps_ep,MEM_HANDLER);
  // Init
  maps_ep->error_position = calloc(LARGE_READ_POS_RANGE,sizeof(uint64_t));
  maps_ep->mismatches = calloc(MISMS_RANGE,sizeof(uint64_t));
  maps_ep->indel_length = calloc(MISMS_RANGE,sizeof(uint64_t));
  maps_ep->errors_events = calloc(MISMS_RANGE,sizeof(uint64_t));
  maps_ep->inss = calloc(INSS_RANGE,sizeof(uint64_t));
  // Mismatch/Indel Distribution
  maps_ep->total_mismatches=0;
  maps_ep->total_indel_length=0;
  maps_ep->total_errors_events=0;
  // Trim stats
  maps_ep->total_bases_aligned=0;
  maps_ep->total_bases_trimmed=0;
  // Mismatch/Errors bases
  maps_ep->misms_transition = calloc(MISMS_BASE_RANGE*MISMS_BASE_RANGE,sizeof(uint64_t));
  maps_ep->qual_score_misms = calloc(QUAL_SCORE_RANGE,sizeof(uint64_t));
  maps_ep->qual_score_errors = calloc(QUAL_SCORE_RANGE,sizeof(uint64_t));
  return maps_ep;
}
GT_INLINE void gt_maps_error_profile_clear(gt_maps_error_profile* const maps_ep) {
  // Init
  memset(maps_ep->error_position,0,LARGE_READ_POS_RANGE*sizeof(uint64_t));
  memset(maps_ep->mismatches,0,MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->indel_length,0,MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->errors_events,0,MISMS_RANGE*sizeof(uint64_t));
  memset(maps_ep->inss,0,INSS_RANGE*sizeof(uint64_t));
  // Mismatch/Indel Distribution
  maps_ep->total_mismatches=0;
  maps_ep->total_indel_length=0;
  maps_ep->total_errors_events=0;
  // Trim stats
  maps_ep->total_bases_aligned=0;
  maps_ep->total_bases_trimmed=0;
  // Mismatch/Errors bases
  memset(maps_ep->misms_transition,0,MISMS_BASE_RANGE*MISMS_BASE_RANGE*sizeof(uint64_t));
  memset(maps_ep->qual_score_misms,0,QUAL_SCORE_RANGE*sizeof(uint64_t));
  memset(maps_ep->qual_score_errors,0,QUAL_SCORE_RANGE*sizeof(uint64_t));
}
GT_INLINE void gt_maps_error_profile_delete(gt_maps_error_profile* const maps_ep) {
  free(maps_ep->error_position);
  free(maps_ep->mismatches);
  free(maps_ep->indel_length);
  free(maps_ep->errors_events);
  free(maps_ep->inss);
  free(maps_ep->misms_transition);
  free(maps_ep->qual_score_misms);
  free(maps_ep->qual_score_errors);
}
GT_INLINE void gt_maps_error_profile_merge(
    gt_maps_error_profile* const maps_ep_dst,gt_maps_error_profile* const maps_ep_src) {
  register uint64_t j;
  // Mismatch/Indel Distribution
  for (j=0;j<MISMS_RANGE;++j) {
    maps_ep_dst->mismatches[j] += maps_ep_src->mismatches[j];
    maps_ep_dst->indel_length[j] += maps_ep_src->indel_length[j];
    maps_ep_dst->errors_events[j] += maps_ep_src->errors_events[j];
  }
  maps_ep_dst->total_mismatches += maps_ep_src->total_mismatches;
  maps_ep_dst->total_indel_length += maps_ep_src->total_indel_length;
  maps_ep_dst->total_errors_events += maps_ep_src->total_errors_events;
  // Inss Range
  for (j=0;j<INSS_RANGE;++j) maps_ep_dst->inss[j] += maps_ep_src->inss[j];
  // Error Positions
  for (j=0;j<LARGE_READ_POS_RANGE;++j) maps_ep_dst->error_position[j] += maps_ep_src->error_position[j];
  // Trim stats
  maps_ep_dst->total_bases_aligned += maps_ep_src->total_bases_aligned;
  maps_ep_dst->total_bases_trimmed += maps_ep_src->total_bases_trimmed;
  // Mismatch/Errors bases
  for (j=0;j<MISMS_BASE_RANGE*MISMS_BASE_RANGE;++j) {
    maps_ep_dst->misms_transition[j] += maps_ep_src->misms_transition[j];
  }
  for (j=0;j<QUAL_SCORE_RANGE;++j) {
    maps_ep_dst->qual_score_misms[j] += maps_ep_src->qual_score_misms[j];
    maps_ep_dst->qual_score_errors[j] += maps_ep_src->qual_score_errors[j];
  }
}

/*
 * SPLITMAPS Profile
 */
GT_INLINE gt_splitmaps_profile* gt_splitmaps_profile_new() {
  // Allocate handler
  gt_splitmaps_profile* splitmaps_profile = malloc(sizeof(gt_splitmaps_profile));
  gt_cond_fatal_error(!splitmaps_profile,MEM_HANDLER);
  // Init
  splitmaps_profile->num_mapped_with_splitmaps = 0;
  splitmaps_profile->num_mapped_only_splitmaps = 0;
  splitmaps_profile->total_splitmaps = 0;
  splitmaps_profile->total_junctions = 0;
  splitmaps_profile->num_junctions = calloc(NUM_JUNCTION_RANGE,sizeof(uint64_t));
  splitmaps_profile->length_junctions = calloc(LEN_JUNCTION_RANGE,sizeof(uint64_t));
  splitmaps_profile->junction_position = calloc(SHORT_READ_POS_RANGE,sizeof(uint64_t));
  splitmaps_profile->pe_sm_sm = 0;
  splitmaps_profile->pe_sm_rm = 0;
  splitmaps_profile->pe_rm_rm = 0;
  return splitmaps_profile;
}
GT_INLINE void gt_splitmaps_profile_clear(gt_splitmaps_profile* const splitmaps_profile) {
  // Init
  splitmaps_profile->num_mapped_with_splitmaps = 0;
  splitmaps_profile->num_mapped_only_splitmaps = 0;
  splitmaps_profile->total_splitmaps = 0;
  splitmaps_profile->total_junctions = 0;
  memset(splitmaps_profile->num_junctions,0,NUM_JUNCTION_RANGE*sizeof(uint64_t));
  memset(splitmaps_profile->length_junctions,0,LEN_JUNCTION_RANGE*sizeof(uint64_t));
  memset(splitmaps_profile->junction_position,0,SHORT_READ_POS_RANGE*sizeof(uint64_t));
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
  register uint64_t j;
  splitmaps_profile_dst->num_mapped_with_splitmaps += splitmaps_profile_src->num_mapped_with_splitmaps;
  splitmaps_profile_dst->num_mapped_only_splitmaps += splitmaps_profile_src->num_mapped_only_splitmaps;
  splitmaps_profile_dst->total_splitmaps += splitmaps_profile_src->total_splitmaps;
  splitmaps_profile_dst->total_junctions += splitmaps_profile_src->total_junctions;
  for (j=0;j<NUM_JUNCTION_RANGE;++j) splitmaps_profile_dst->num_junctions[j] += splitmaps_profile_src->num_junctions[j];
  for (j=0;j<LEN_JUNCTION_RANGE;++j) splitmaps_profile_dst->length_junctions[j] += splitmaps_profile_src->length_junctions[j];
  for (j=0;j<SHORT_READ_POS_RANGE;++j) splitmaps_profile_dst->junction_position[j] += splitmaps_profile_src->junction_position[j];
  splitmaps_profile_dst->pe_sm_sm += splitmaps_profile_src->pe_sm_sm;
  splitmaps_profile_dst->pe_sm_rm += splitmaps_profile_src->pe_sm_rm;
  splitmaps_profile_dst->pe_rm_rm += splitmaps_profile_src->pe_rm_rm;
}

/*
 * Calculate Distributions
 */
void gt_stats_get_misms(uint64_t* const misms_array,const uint64_t misms_num,const uint64_t read_length) {
  register const double one_per_cent = read_length/100.0;
  if (misms_num==0) {
    misms_array[MISMS_RANGE_0]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(1,one_per_cent)) {
    misms_array[MISMS_RANGE_1]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(2,one_per_cent)) {
    misms_array[MISMS_RANGE_2]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(3,one_per_cent)) {
    misms_array[MISMS_RANGE_3]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(4,one_per_cent)) {
    misms_array[MISMS_RANGE_4]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(5,one_per_cent)) {
    misms_array[MISMS_RANGE_5]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(6,one_per_cent)) {
    misms_array[MISMS_RANGE_6]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(7,one_per_cent)) {
    misms_array[MISMS_RANGE_7]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(8,one_per_cent)) {
    misms_array[MISMS_RANGE_8]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(9,one_per_cent)) {
    misms_array[MISMS_RANGE_9]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(10,one_per_cent)) {
    misms_array[MISMS_RANGE_10]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(20,one_per_cent)) {
    misms_array[MISMS_RANGE_20]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(50,one_per_cent)) {
    misms_array[MISMS_RANGE_50]++;
  } else {
    misms_array[MISMS_RANGE_BEHOND]++;
  }
}
void gt_stats_get_mmap_distribution(uint64_t* const mmap,const uint64_t num_maps) {
  if (num_maps>0) {
    if (num_maps<=1) {
      ++mmap[MMAP_RANGE_1];
    } else if (num_maps<=5) {
      ++mmap[MMAP_RANGE_5];
    } else if (num_maps<=10) {
      ++mmap[MMAP_RANGE_10];
    } else if (num_maps<=50) {
      ++mmap[MMAP_RANGE_50];
    } else if (num_maps<=100) {
      ++mmap[MMAP_RANGE_100];
    } else if (num_maps<=500) {
      ++mmap[MMAP_RANGE_500];
    } else if (num_maps<=1000) {
      ++mmap[MMAP_RANGE_1000];
    } else {
      ++mmap[MMAP_RANGE_BEHOND];
    }
  }
}
void gt_stats_get_uniq_distribution(uint64_t* const uniq,const uint64_t uniq_degree) {
  if (uniq_degree==GT_NO_STRATA) {
    ++uniq[UNIQ_RANGE_X];
  } else if (uniq_degree==0) {
    ++uniq[UNIQ_RANGE_0];
  } else if (uniq_degree<=1) {
    ++uniq[UNIQ_RANGE_1];
  } else if (uniq_degree<=2) {
    ++uniq[UNIQ_RANGE_2];
  } else if (uniq_degree<=3) {
    ++uniq[UNIQ_RANGE_3];
  } else if (uniq_degree<=10) {
    ++uniq[UNIQ_RANGE_10];
  } else if (uniq_degree<=50) {
    ++uniq[UNIQ_RANGE_50];
  } else if (uniq_degree<=100) {
    ++uniq[UNIQ_RANGE_100];
  } else if (uniq_degree<=500) {
    ++uniq[UNIQ_RANGE_500];
  } else {
    ++uniq[UNIQ_RANGE_BEHOND];
  }
}
void gt_stats_get_inss_distribution(uint64_t* const inss,const int64_t insert_size) {
  if (insert_size<=(-100)) {
    ++inss[INSS_RANGE_NEG];
  } else if (insert_size<=0) {
    ++inss[INSS_RANGE_OVER];
  } else if (insert_size<=100) {
    ++inss[INSS_RANGE_100];
  } else if (insert_size<=200) {
    ++inss[INSS_RANGE_200];
  } else if (insert_size<=300) {
    ++inss[INSS_RANGE_300];
  } else if (insert_size<=400) {
    ++inss[INSS_RANGE_400];
  } else if (insert_size<=500) {
    ++inss[INSS_RANGE_500];
  } else if (insert_size<=600) {
    ++inss[INSS_RANGE_600];
  } else if (insert_size<=700) {
    ++inss[INSS_RANGE_700];
  } else if (insert_size<=800) {
    ++inss[INSS_RANGE_800];
  } else if (insert_size<=900) {
    ++inss[INSS_RANGE_900];
  } else if (insert_size<=1000) {
    ++inss[INSS_RANGE_1000];
  } else if (insert_size<=2000) {
    ++inss[INSS_RANGE_2000];
  } else if (insert_size<=5000) {
    ++inss[INSS_RANGE_5000];
  } else {
    ++inss[INSS_RANGE_BEHOND];
  }
}

#define NUM_JUNCTION_1      0
#define NUM_JUNCTION_2      1
#define NUM_JUNCTION_3      2
#define NUM_JUNCTION_BEHOND 3
#define NUM_JUNCTION_RANGE  4

#define LEN_JUNCTION_100    0
#define LEN_JUNCTION_1000   1
#define LEN_JUNCTION_5000   2
#define LEN_JUNCTION_10000  3
#define LEN_JUNCTION_50000  4
#define LEN_JUNCTION_BEHOND 5
#define LEN_JUNCTION_RANGE  6

void gt_stats_get_juntions_distribution(uint64_t* const junctions,const uint64_t num_junctions) {
  //if (num_junctions>0) {
    if (num_junctions<=1) {
      ++junctions[NUM_JUNCTION_1];
    } else if (num_junctions==2) {
      ++junctions[NUM_JUNCTION_2];
    } else if (num_junctions==3) {
      ++junctions[NUM_JUNCTION_3];
    } else {
      ++junctions[NUM_JUNCTION_BEHOND];
    }
  //}
}
void gt_stats_get_juntions_length_distribution(uint64_t* const junctions_length,const uint64_t length) {
  //if (length>0) {
    if (length<=100) {
      ++junctions_length[LEN_JUNCTION_100];
    } else if (length<=1000) {
      ++junctions_length[LEN_JUNCTION_1000];
    } else if (length<=5000) {
      ++junctions_length[LEN_JUNCTION_5000];
    } else if (length<=10000) {
      ++junctions_length[LEN_JUNCTION_10000];
    } else if (length<=50000) {
      ++junctions_length[LEN_JUNCTION_50000];
    } else {
      ++junctions_length[LEN_JUNCTION_BEHOND];
    }
  //}
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
  stats->total_bases_unaligned=0;
  // Mapped/Maps/MMaps/Uniq...
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  stats->mmap = calloc(MMAP_RANGE,sizeof(uint64_t)); // MMaps
  stats->uniq = calloc(UNIQ_RANGE,sizeof(uint64_t)); // Uniq
  // Maps Error Profile
  stats->maps_error_profile = gt_maps_error_profile_new();
  // Split maps Profile
  stats->splitmaps_profile = gt_splitmaps_profile_new();
  return stats;
}
GT_INLINE void gt_stats_clear(gt_stats *stats) {
  // Length
  stats->min_length=UINT64_MAX;
  stats->max_length=0;
  stats->total_bases=0;
  stats->total_bases_unaligned=0;
  // Mapped/Maps/MMaps/Uniq...
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  memset(stats->mmap,0,MMAP_RANGE*sizeof(uint64_t)); // MMaps
  memset(stats->uniq,0,UNIQ_RANGE*sizeof(uint64_t)); // Uniq
  // Maps Error Profile
  gt_maps_error_profile_clear(stats->maps_error_profile);
  // Split maps Profile
  gt_splitmaps_profile_clear(stats->splitmaps_profile);
}
GT_INLINE void gt_stats_delete(gt_stats *stats) {
  free(stats->mmap);
  free(stats->uniq);
  gt_maps_error_profile_delete(stats->maps_error_profile);
  gt_splitmaps_profile_delete(stats->splitmaps_profile);
  free(stats);
}

/*
 * STATS Merge
 */
void gt_stats_merge(gt_stats** const stats,const uint64_t stats_array_size) {
  register uint64_t i, j;
  for (i=1;i<stats_array_size;++i) {
    // Length
    stats[0]->min_length = GT_MIN(stats[0]->min_length,stats[i]->min_length);
    stats[0]->max_length = GT_MAX(stats[0]->max_length,stats[i]->max_length);
    stats[0]->total_bases += stats[i]->total_bases;
    stats[0]->total_bases_unaligned += stats[i]->total_bases_unaligned;
    // Mapped/Maps
    stats[0]->num_blocks += stats[i]->num_blocks;
    stats[0]->num_alignments += stats[i]->num_alignments;
    stats[0]->num_maps += stats[i]->num_maps;
    stats[0]->num_mapped += stats[i]->num_mapped;
    for (j=0;j<MMAP_RANGE;++j) stats[0]->mmap[j] += stats[i]->mmap[j];
    for (j=0;j<UNIQ_RANGE;++j) stats[0]->uniq[j] += stats[i]->uniq[j];
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
void gt_stats_get_maps_profile(
    gt_maps_error_profile *maps_error_profile,gt_template* const template,
    uint64_t const read_length,gt_map** const mmap) {
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t total_mismatches=0, total_indel_length=0, total_errors_events=0, trimmed=0;
  GT_MULTIMAP_ITERATE_BLOCKS(mmap,num_blocks,map,end_pos) {
    register gt_alignment* const alignment = gt_template_get_block(template,end_pos);
    register gt_string* const read = alignment->read;
    register const bool has_qualities = gt_alignment_has_qualities(alignment);
    register gt_string* const quals = alignment->qualities;
    register char quality_misms = '0';
    GT_MAP_ITERATE(map,map_block) {
      GT_MISMS_ITERATE(map_block,misms) {
        ++total_errors_events;
        // Records position of misms/indel
        if (misms->position < LARGE_READ_POS_RANGE) maps_error_profile->error_position[misms->position]++;
        // Record quality of misms/indel
        if (has_qualities) {
          quality_misms = gt_string_get_string(quals)[misms->position];
          maps_error_profile->qual_score_errors[(uint8_t)quality_misms]++;
        }
        switch (misms->misms_type) {
          case MISMS:
            ++total_mismatches;
            maps_error_profile->misms_transition[
               (gt_cdna_encode[(uint8_t)gt_string_get_string(read)[misms->position]]*MISMS_BASE_RANGE)+
                gt_cdna_encode[(uint8_t)misms->base]]++;
            if (has_qualities) maps_error_profile->qual_score_misms[(uint8_t)quality_misms]++;
            break;
          case DEL:
            // Record trim info
            if (misms->position==0 || misms->position+misms->size==map_block->base_length) {
              trimmed+=misms->size;
            }
          case INS:
            total_indel_length+=misms->size;
            break;
        }
      }
    }
  }
  // Store stats
  maps_error_profile->total_bases_aligned += read_length-trimmed;
  maps_error_profile->total_bases_trimmed += trimmed;
  maps_error_profile->total_mismatches += total_mismatches;
  maps_error_profile->total_indel_length += total_indel_length;
  maps_error_profile->total_errors_events += total_errors_events;
  // Record stats
  gt_stats_get_misms(maps_error_profile->mismatches,total_mismatches,read_length);
  gt_stats_get_misms(maps_error_profile->indel_length,total_indel_length,read_length);
  gt_stats_get_misms(maps_error_profile->errors_events,total_errors_events,read_length);
}

void gt_stats_mmaps_profile(
    gt_stats* const stats,gt_template* const template,const uint64_t total_read_length,const bool best_map) {
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
    // Insert Size Distribution
    if (paired_map) {
      gt_stats_get_inss_distribution(maps_error_profile->inss,gt_template_get_insert_size(mmap));
    }
    // Error Profile
    gt_stats_get_maps_profile(maps_error_profile,template,total_read_length,mmap);
    // SM block stats
    bool has_sm[2] = {true, true};
    GT_MULTIMAP_ITERATE(mmap,map,end_pos) {
      register const uint64_t num_blocks = gt_map_get_num_blocks(map);
      // Calculate total_junctions & total_splitmaps
      if (num_blocks>1) {
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
          if (juntion_position < SHORT_READ_POS_RANGE) splitmaps_profile->junction_position[juntion_position]++;
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
    // Break if we just one the best map
    if (best_map) break;
  } // GT_TEMPLATE_ITERATE_END;
  // Global SM stats about the whole alignment
  if (has_splitsmaps) {
    splitmaps_profile->num_mapped_with_splitmaps++;
    if (only_splitsmaps) splitmaps_profile->num_mapped_only_splitmaps++;
  }
}
GT_INLINE void gt_stats_calculate_template_stats(gt_stats* const stats,gt_template* const template,const bool best_map) {
  /*
   * Blocks/Alignments STATS
   */
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t num_maps = gt_template_get_num_mmaps(template);
  if (best_map && num_maps>0) num_maps = 1;
  ++stats->num_alignments;
  stats->num_blocks += num_blocks;
  /*
   * Length STATS
   */
  register uint64_t total_length = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    register uint64_t read_length =
        (gt_alignment_get_read_length(alignment)==0 && gt_alignment_get_num_maps(alignment)>0) ?
         gt_map_get_base_length(gt_alignment_get_map(alignment,0)) : gt_alignment_get_read_length(alignment);
    total_length += read_length;
    stats->min_length = GT_MIN(stats->min_length,read_length);
    stats->max_length = GT_MAX(stats->max_length,read_length);
  }
  stats->total_bases += total_length;
  /*
   * Maps/MMaps STATS
   */
  if (num_maps > 0 || gt_template_is_mapped(template)){ // Is mapped?
    ++stats->num_mapped; // Mapped/Maps
    stats->num_maps += num_maps;
    gt_stats_get_mmap_distribution(stats->mmap,num_maps); // MMap Distribution
  } else {
    stats->total_bases_unaligned += total_length;
  }
  /*
   * Uniq Distribution
   */
  gt_stats_get_uniq_distribution(stats->uniq,gt_template_get_uniq_degree(template));
  /*
   * MMaps Profile {Insert Size Distribution, Error Profile, SM Profile, ...}
   */
  if (num_maps > 0) gt_stats_mmaps_profile(stats,template,total_length,best_map);
}
