/*
 * PROJECT: GEM-Tools library
 * FILE: gt.stats.c
 * DATE: 02/08/2012
 * DESCRIPTION: Application to retrieve very naive stats from {MAP,SAM,FASTQ} files
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

#define MMAP_RANGE_1 0
#define MMAP_RANGE_5 1
#define MMAP_RANGE_10 2
#define MMAP_RANGE_50 3
#define MMAP_RANGE_100 4
#define MMAP_RANGE_500 5
#define MMAP_RANGE_1000 6
#define MMAP_RANGE_BEHOND 7
#define MMAP_RANGE 8

#define INSS_RANGE_100 0
#define INSS_RANGE_200 1
#define INSS_RANGE_300 2
#define INSS_RANGE_400 3
#define INSS_RANGE_500 4
#define INSS_RANGE_600 5
#define INSS_RANGE_700 6
#define INSS_RANGE_800 7
#define INSS_RANGE_900 8
#define INSS_RANGE_1000 9
#define INSS_RANGE_2000 10
#define INSS_RANGE_5000 11
#define INSS_RANGE_10000 12
#define INSS_RANGE_BEHOND 13
#define INSS_RANGE 14

#define MISMS_RANGE_0 0
#define MISMS_RANGE_1 1
#define MISMS_RANGE_2 2
#define MISMS_RANGE_3 3
#define MISMS_RANGE_4 4
#define MISMS_RANGE_5 5
#define MISMS_RANGE_6 6
#define MISMS_RANGE_7 7
#define MISMS_RANGE_8 8
#define MISMS_RANGE_9 9
#define MISMS_RANGE_10 10
#define MISMS_RANGE_20 11
#define MISMS_RANGE_50 12
#define MISMS_RANGE_BEHOND 13
#define MISMS_RANGE 14

#define UNIQ_RANGE_0 0
#define UNIQ_RANGE_1 1
#define UNIQ_RANGE_2 2
#define UNIQ_RANGE_3 3
#define UNIQ_RANGE_10 4
#define UNIQ_RANGE_50 5
#define UNIQ_RANGE_100 6
#define UNIQ_RANGE_500 7
#define UNIQ_RANGE_BEHOND 8
#define UNIQ_RANGE_X 9
#define UNIQ_RANGE 10

#define ERROR_POS_RANGE 1000

typedef struct {
  char *name_input_file;
  bool mmap_input;
  bool paired_end;
  uint64_t num_reads;
  uint64_t num_threads;
  bool compact;
  bool verbose;
  bool quiet;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .num_reads=0,
    .num_threads=1,
    .compact = false,
    .verbose=false,
    .quiet=false,
};

typedef struct {
  // Mismatch/Indel Distribution
  uint64_t *misms;
  uint64_t *indel;
  uint64_t *errors;
  uint64_t total_hamming;
  uint64_t total_indels;
  uint64_t total_misms_events;
  uint64_t *error_position;
  // Trim stats
  uint64_t total_bases_aligned;
  uint64_t total_bases_trimmed;
  // Insert Size Distribution
  uint64_t *inss;
} gt_maps_profile;

typedef struct {
  // Mapped/Maps
  uint64_t num_blocks;
  uint64_t num_alignments;
  uint64_t num_maps;
  uint64_t num_mapped;
  // MMap Distribution
  uint64_t *mmap;
  // Uniq Distribution
  uint64_t *uniq;
  // Error Profile
  gt_maps_profile *best_match_ep;
  gt_maps_profile *mmaps_match_ep;
  // Length
  uint64_t min_length;
  uint64_t max_length;
  uint64_t total_bases;
} gt_stats;

/*
 * STATS handler functions
 */
void gt_maps_profile_init(gt_maps_profile *maps_profile) {
  maps_profile->error_position = malloc(sizeof(uint64_t)*ERROR_POS_RANGE);
  maps_profile->misms = malloc(sizeof(uint64_t)*MISMS_RANGE);
  maps_profile->indel = malloc(sizeof(uint64_t)*MISMS_RANGE);
  maps_profile->errors = malloc(sizeof(uint64_t)*MISMS_RANGE);
  maps_profile->inss = malloc(sizeof(uint64_t)*INSS_RANGE);
  // Mismatch/Indel Distribution
  register uint64_t j;
  for (j=0;j<MISMS_RANGE;++j) {
    maps_profile->misms[j] = 0;
    maps_profile->indel[j] = 0;
    maps_profile->errors[j] = 0;
  }
  for (j=0;j<ERROR_POS_RANGE;++j) {
    maps_profile->error_position[j] = 0;
  }
  for (j=0;j<INSS_RANGE;++j) maps_profile->inss[j] = 0;
  maps_profile->total_hamming=0;
  maps_profile->total_indels=0;
  maps_profile->total_misms_events=0;
  // Trim stats
  maps_profile->total_bases_aligned=0;
  maps_profile->total_bases_trimmed=0;
}
GT_INLINE void gt_stats_init(gt_stats** _stats) {
  // Allocate handler
  *_stats = malloc(sizeof(gt_stats));
  gt_stats* stats = *_stats;
  gt_cond_fatal_error(!stats,MEM_HANDLER);
  // Mapped/Maps
  stats->num_blocks=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  // MMaps
  register uint64_t j;
  stats->mmap = malloc(sizeof(uint64_t)*MMAP_RANGE);
  for (j=0;j<MMAP_RANGE;++j) stats->mmap[j] = 0;
  // Uniq
  stats->uniq = malloc(sizeof(uint64_t)*UNIQ_RANGE);
  for (j=0;j<UNIQ_RANGE;++j) stats->uniq[j] = 0;
  // Length
  stats->min_length=UINT64_MAX;
  stats->max_length=0;
  stats->total_bases=0;
  // Error Profile
  stats->best_match_ep = malloc(sizeof(gt_maps_profile));
  stats->mmaps_match_ep = malloc(sizeof(gt_maps_profile));
  gt_maps_profile_init(stats->best_match_ep);
  gt_maps_profile_init(stats->mmaps_match_ep);
}
void gt_stats_delete(gt_stats *stats) {
  free(stats->mmap);
  free(stats->uniq);
  free(stats->best_match_ep->error_position);
  free(stats->best_match_ep->misms);
  free(stats->best_match_ep->indel);
  free(stats->best_match_ep->errors);
  free(stats->best_match_ep->inss);
  free(stats->mmaps_match_ep->error_position);
  free(stats->mmaps_match_ep->misms);
  free(stats->mmaps_match_ep->indel);
  free(stats->mmaps_match_ep->errors);
  free(stats->mmaps_match_ep->inss);
  free(stats->best_match_ep);
  free(stats->mmaps_match_ep);
  free(stats);
}
void gt_maps_profile_merge(gt_maps_profile *maps_profile_dst,gt_maps_profile *maps_profile_src) {
  register uint64_t j;
  // Mismatch/Indel Distribution
  for (j=0;j<MISMS_RANGE;++j) {
    maps_profile_dst->misms[j] += maps_profile_src->misms[j];
    maps_profile_dst->indel[j] += maps_profile_src->indel[j];
    maps_profile_dst->errors[j] += maps_profile_src->errors[j];
  }
  maps_profile_dst->total_hamming += maps_profile_src->total_hamming;
  maps_profile_dst->total_indels += maps_profile_src->total_indels;
  maps_profile_dst->total_misms_events += maps_profile_src->total_misms_events;
  // Inss Range
  for (j=0;j<INSS_RANGE;++j) {
    maps_profile_dst->inss[j] += maps_profile_src->inss[j];
  }
  // Error Pos
  for (j=0;j<ERROR_POS_RANGE;++j) {
    maps_profile_dst->error_position[j] += maps_profile_src->error_position[j];
  }
  // Trim stats
  maps_profile_dst->total_bases_aligned += maps_profile_src->total_bases_aligned;
  maps_profile_dst->total_bases_trimmed += maps_profile_src->total_bases_trimmed;
}
void gt_stats_merge(gt_stats** const stats,const uint64_t num_stats) {
  register uint64_t i, j;
  for (i=1;i<num_stats;++i) {
    stats[0]->num_blocks += stats[i]->num_blocks;
    stats[0]->num_alignments += stats[i]->num_alignments;
    stats[0]->num_maps += stats[i]->num_maps;
    stats[0]->num_mapped += stats[i]->num_mapped;
    for (j=0;j<MMAP_RANGE;++j) {
      stats[0]->mmap[j] += stats[i]->mmap[j];
    }
    for (j=0;j<UNIQ_RANGE;++j) {
      stats[0]->uniq[j] += stats[i]->uniq[j];
    }
    // Error Profile
    gt_maps_profile_merge(stats[0]->best_match_ep,stats[i]->best_match_ep);
    gt_maps_profile_merge(stats[0]->mmaps_match_ep,stats[i]->mmaps_match_ep);
    // Length
    stats[0]->min_length = GT_MIN(stats[0]->min_length,stats[i]->min_length);
    stats[0]->max_length = GT_MAX(stats[0]->max_length,stats[i]->max_length);
    stats[0]->total_bases += stats[i]->total_bases;
    // Free
    gt_stats_delete(stats[i]);
  }
}
/*
 * STATS analysis functions
 */
#define GT_STATS_GET_PERCENTAGE_ERROR(percentage) (uint64_t)((double)percentage*one_per_cent)
void gt_stats_get_misms(uint64_t* const misms_array,const uint64_t misms_num,const double one_per_cent) {
  if (misms_num==0) {
    misms_array[MISMS_RANGE_0]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(1)) {
    misms_array[MISMS_RANGE_1]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(2)) {
    misms_array[MISMS_RANGE_2]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(3)) {
    misms_array[MISMS_RANGE_3]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(4)) {
    misms_array[MISMS_RANGE_4]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(5)) {
    misms_array[MISMS_RANGE_5]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(6)) {
    misms_array[MISMS_RANGE_6]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(7)) {
    misms_array[MISMS_RANGE_7]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(8)) {
    misms_array[MISMS_RANGE_8]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(9)) {
    misms_array[MISMS_RANGE_9]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(10)) {
    misms_array[MISMS_RANGE_10]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(20)) {
    misms_array[MISMS_RANGE_20]++;
  } else if (misms_num<=GT_STATS_GET_PERCENTAGE_ERROR(50)) {
    misms_array[MISMS_RANGE_50]++;
  } else {
    misms_array[MISMS_RANGE_BEHOND]++;
  }
}
void gt_stats_get_maps_profile(
    gt_maps_profile *match_ep,uint64_t const read_length,
    gt_map** const mmap,const uint64_t num_blocks) {
  uint64_t num_hamming=0, num_indels=0, num_misms=0, trimmed=0;
  GT_MULTIMAP_ITERATE_BLOCKS(mmap,num_blocks,map,end_pos) {
    GT_MAP_ITERATE(map,map_block) {
      GT_MISMS_ITERATE(map_block,misms) {
        ++num_misms;
        if (misms->position < ERROR_POS_RANGE) match_ep->error_position[misms->position]++;
        switch (misms->misms_type) {
          case MISMS:
            ++num_hamming;
            break;
          case DEL:
            if (misms->position==0 || misms->position+misms->size==map_block->base_length) {
              trimmed+=misms->size;
            }
          case INS:
            num_indels+=misms->size;
            break;
        }
      }
    }
  }
  // Store stats
  match_ep->total_bases_aligned += read_length-trimmed;
  match_ep->total_bases_trimmed += trimmed;
  match_ep->total_indels += num_indels;
  match_ep->total_hamming += num_hamming;
  match_ep->total_indels += num_indels;
  match_ep->total_misms_events += num_misms;
  // Record stats
  const double one_per_cent = (double)read_length/100.0;
  gt_stats_get_misms(match_ep->misms,num_hamming,one_per_cent);
  gt_stats_get_misms(match_ep->indel,num_indels,one_per_cent);
  gt_stats_get_misms(match_ep->errors,num_misms,one_per_cent);
}
void gt_stats_get_mmap_distribution(gt_stats *stats,const uint64_t num_maps) {
  if (num_maps>0) {
    if (num_maps<=1) {
      ++stats->mmap[MMAP_RANGE_1];
    } else if (num_maps<=5) {
      ++stats->mmap[MMAP_RANGE_5];
    } else if (num_maps<=10) {
      ++stats->mmap[MMAP_RANGE_10];
    } else if (num_maps<=50) {
      ++stats->mmap[MMAP_RANGE_50];
    } else if (num_maps<=100) {
      ++stats->mmap[MMAP_RANGE_100];
    } else if (num_maps<=500) {
      ++stats->mmap[MMAP_RANGE_500];
    } else if (num_maps<=1000) {
      ++stats->mmap[MMAP_RANGE_1000];
    } else {
      ++stats->mmap[MMAP_RANGE_BEHOND];
    }
  }
}
void gt_stats_get_uniq_distribution(gt_stats *stats,gt_template* const template) {
  register const uint64_t uniq_degree = gt_template_get_uniq_degree(template);
  if (uniq_degree==GT_NO_STRATA) {
    ++stats->uniq[UNIQ_RANGE_X];
  } else if (uniq_degree==0) {
    ++stats->uniq[UNIQ_RANGE_0];
  } else if (uniq_degree<=1) {
    ++stats->uniq[UNIQ_RANGE_1];
  } else if (uniq_degree<=2) {
    ++stats->uniq[UNIQ_RANGE_2];
  } else if (uniq_degree<=3) {
    ++stats->uniq[UNIQ_RANGE_3];
  } else if (uniq_degree<=10) {
    ++stats->uniq[UNIQ_RANGE_10];
  } else if (uniq_degree<=50) {
    ++stats->uniq[UNIQ_RANGE_50];
  } else if (uniq_degree<=100) {
    ++stats->uniq[UNIQ_RANGE_100];
  } else if (uniq_degree<=500) {
    ++stats->uniq[UNIQ_RANGE_500];
  } else {
    ++stats->uniq[UNIQ_RANGE_BEHOND];
  }
}
void gt_stats_get_inss_distribution(uint64_t* const inss,const uint64_t insert_size) {
  if (insert_size<=100) {
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
GT_INLINE void gt_stats_analize_template(gt_stats* const stats,gt_template* const template,const bool paired_map) {
  // Blocks/Alignments
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t num_maps = (num_blocks>1 || !paired_map) ?
      gt_template_get_num_mmaps(template) : 0;
  ++stats->num_alignments;
  stats->num_blocks += num_blocks;
  // Length
  register uint64_t total_length = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    register uint64_t read_length = (gt_alignment_get_num_maps(alignment)>0) ?
        gt_map_get_base_length(gt_alignment_get_map(alignment,0)) : gt_alignment_get_read_length(alignment);
    total_length += read_length;
    stats->min_length = GT_MIN(stats->min_length,read_length);
    stats->max_length = GT_MAX(stats->max_length,read_length);
  }
  stats->total_bases += total_length;
  // Is mapped?
  if ( (paired_map &&  (num_maps > 0 || gt_template_is_mapped(template))) ||
       (!paired_map && (num_maps > 0 || gt_alignment_is_mapped(gt_template_get_block(template,0)))) ){
    ++stats->num_mapped; // Mapped/Maps
  }
  stats->num_maps += num_maps;
  gt_stats_get_mmap_distribution(stats,num_maps); // MMap Distribution
  // Uniq Distribution
  gt_stats_get_uniq_distribution(stats,template);
  // Insert Size Distribution +  Error Profile
  if (num_maps > 0) {
    register gt_map** const best_map = gt_template_get_mmap(template,0,NULL);
    // Best map error profile
    gt_stats_get_maps_profile(stats->best_match_ep,total_length,best_map,num_blocks);
    if (paired_map) { // Insert Size distribution
      gt_stats_get_inss_distribution(stats->best_match_ep->inss,abs((int64_t)best_map[0]->position-(int64_t)best_map[1]->position));
    }
    // MMaps error profile
    GT_TEMPLATE_ITERATE_(template,mmap) {
      if (paired_map) { // Insert Size distribution
        gt_stats_get_inss_distribution(stats->mmaps_match_ep->inss,abs((int64_t)mmap[0]->position-(int64_t)mmap[1]->position));
      }
      gt_stats_get_maps_profile(stats->mmaps_match_ep,total_length,mmap,num_blocks);
    }
  }
}
void gt_stats_print_mmap_distribution(uint64_t* const mmap,const uint64_t num_alignments,const uint64_t num_mapped) {
#define GT_STATS_PRINT_MMAP_FORMAT "%6lu \t %1.3f%%\n"
#define GT_STATS_PRINT_MMAP(RANGE) mmap[RANGE],100.0*(float)mmap[RANGE]/(float)num_alignments
  fprintf(stderr,"MMap.ranges\n");
  fprintf(stderr,"  -->        [0] \t=> "GT_STATS_PRINT_MMAP_FORMAT,(num_alignments-num_mapped),100.0*(float)(num_alignments-num_mapped)/(float)num_alignments);
  fprintf(stderr,"  -->        [1] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_1));
  fprintf(stderr,"  -->      (1,5] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_5));
  fprintf(stderr,"  -->     (5,10] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_10));
  fprintf(stderr,"  -->    (10,50] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_50));
  fprintf(stderr,"  -->   (50,100] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_100));
  fprintf(stderr,"  -->  (100,500] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_500));
  fprintf(stderr,"  --> (500,1000] \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_1000));
  fprintf(stderr,"  --> (1000,inf) \t=> "GT_STATS_PRINT_MMAP_FORMAT,GT_STATS_PRINT_MMAP(MMAP_RANGE_BEHOND));
}
void gt_stats_print_uniq_distribution(uint64_t* const uniq,const uint64_t num_alignments) {
#define GT_STATS_PRINT_UNIQ_FORMAT "%8lu \t %1.3f%%\n"
#define GT_STATS_PRINT_UNIQ(RANGE) uniq[RANGE],100.0*(float)uniq[RANGE]/(float)num_alignments
  fprintf(stderr,"Uniq.ranges\n");
  fprintf(stderr,"  -->        [X] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_X));
  fprintf(stderr,"  -->        [0] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_0));
  fprintf(stderr,"  -->        [1] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_1));
  fprintf(stderr,"  -->        [2] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_2));
  fprintf(stderr,"  -->        [3] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_3));
  fprintf(stderr,"  -->     (3,10] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_10));
  fprintf(stderr,"  -->    (10,50] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_50));
  register const uint64_t accum = uniq[UNIQ_RANGE_100]+uniq[UNIQ_RANGE_500]+uniq[UNIQ_RANGE_BEHOND];
  fprintf(stderr,"  -->   (50,inf] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,accum,100.0*(float)(accum)/(float)num_alignments);
}
void gt_stats_print_inss_distribution(uint64_t* const inss_best,uint64_t* const inss_all,const uint64_t num_mapped,const uint64_t num_maps) {
#define GT_STATS_PRINT_INSS_FORMAT "%8lu/%8lu \t %1.3f%%/%1.3f%%\n"
#define GT_STATS_PRINT_INSS(RANGE) inss_best[RANGE],inss_all[RANGE],100.0*(float)inss_best[RANGE]/(float)num_mapped,100.0*(float)inss_all[RANGE]/(float)num_maps
  fprintf(stderr,"InsS.ranges\n");
  fprintf(stderr,"  -->       [0,100] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_100));
  fprintf(stderr,"  -->     (100,200] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_200));
  fprintf(stderr,"  -->     (200,300] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_300));
  fprintf(stderr,"  -->     (300,400] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_400));
  fprintf(stderr,"  -->     (400,500] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_500));
  fprintf(stderr,"  -->     (500,600] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_600));
  fprintf(stderr,"  -->     (600,700] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_700));
  fprintf(stderr,"  -->     (700,800] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_800));
  fprintf(stderr,"  -->     (800,900] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_900));
  fprintf(stderr,"  -->    (900,1000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_1000));
  fprintf(stderr,"  -->   (1000,2000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_2000));
  fprintf(stderr,"  -->   (2000,5000] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_5000));
  fprintf(stderr,"  -->    (5000,inf] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_BEHOND));
}
void gt_stats_print_misms_distribution(uint64_t* const misms_best,uint64_t* const misms_all,const uint64_t num_mapped,const uint64_t num_maps,const char* const header) {
#define GT_STATS_PRINT_MISMS_FORMAT "%8lu/%8lu \t %1.3f%%/%1.3f%%\n"
#define GT_STATS_PRINT_MISMS(RANGE) misms_best[RANGE],misms_all[RANGE],100.0*(float)misms_best[RANGE]/(float)num_mapped,100.0*(float)misms_all[RANGE]/(float)num_maps
#define GT_STATS_PRINT_MISMS_VAL(misms_best_val,misms_all_val) misms_best_val,misms_all_val,100.0*(float)misms_best_val/(float)num_mapped,100.0*(float)misms_all_val/(float)num_maps
  fprintf(stderr,"%s\n",header);
  fprintf(stderr,"  -->         [0]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_0));
  fprintf(stderr,"  -->         [1]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_1));
  fprintf(stderr,"  -->         [2]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_2));
  register const uint64_t misms_best_accum = misms_best[MISMS_RANGE_3]+misms_best[MISMS_RANGE_4]+
      misms_best[MISMS_RANGE_5]+misms_best[MISMS_RANGE_6]+misms_best[MISMS_RANGE_7]+
      misms_best[MISMS_RANGE_8]+misms_best[MISMS_RANGE_9]+misms_best[MISMS_RANGE_10];
  register const uint64_t misms_all_accum = misms_all[MISMS_RANGE_3]+misms_all[MISMS_RANGE_4]+
      misms_all[MISMS_RANGE_5]+misms_all[MISMS_RANGE_6]+misms_all[MISMS_RANGE_7]+
      misms_all[MISMS_RANGE_8]+misms_all[MISMS_RANGE_9]+misms_all[MISMS_RANGE_10];
  fprintf(stderr,"  -->      (2,10]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS_VAL(misms_best_accum,misms_all_accum));
  fprintf(stderr,"  -->     (10,20]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_20));
  fprintf(stderr,"  -->     (20,50]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_50));
  fprintf(stderr,"  -->    (50,100]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_BEHOND));
}
uint64_t gt_stats_sum_misms_pos(uint64_t* const pos_error,uint64_t const begin,uint64_t const end) {
  register uint64_t i, accum = 0;
  for (i=begin;i<end;++i) {
    accum += pos_error[i];
  }
  return accum;
}
void gt_stats_print_misms_positions(uint64_t* const pos_error_best,uint64_t* const pos_error_all,
    uint64_t const num_errors_best,uint64_t const num_errors_all,uint64_t const max_length) {
  register const uint64_t max_length_stats = GT_MIN(max_length,ERROR_POS_RANGE);
  register const uint64_t ten_per_cent = max_length_stats/10;
  register uint64_t i, centinel, error_sum_best, error_sum_all;
  fprintf(stderr,"Error.position [0,%lu]nt\n",max_length_stats);
  for (i=0,centinel=0;i<100;i+=10,centinel+=ten_per_cent) {
    register const uint64_t top = (i==90)?max_length_stats:centinel+ten_per_cent;
    error_sum_best = gt_stats_sum_misms_pos(pos_error_best,centinel,top);
    error_sum_all = gt_stats_sum_misms_pos(pos_error_all,centinel,top);
    fprintf(stderr,"  -->   [%3lu-%3lu)nt \t=> %1.3f%% / %1.3f%%\n",
        centinel,top,100.0*(float)error_sum_best/(float)num_errors_best,
        100.0*(float)error_sum_all/(float)num_errors_all);
  }
}
void gt_stats_print_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  fprintf(stderr,"[GENERAL.STATS]\n");
  /*
   * Reads
   */
  register const uint64_t eff_num_reads = paired_end?num_reads>>1:num_reads;
  fprintf(stderr,"Num.reads %lu\n",num_reads);
  fprintf(stderr,"  --> Length.(min,avg,max) (%lu,%lu,%lu)\n",
      stats->min_length,stats->total_bases/stats->num_alignments,stats->max_length);
  fprintf(stderr,"  --> Num.mapped %lu (%2.3f%%)\n",stats->num_mapped,
      100.0*(float)stats->num_mapped/(float)eff_num_reads);
  fprintf(stderr,"Num.alignments %lu\n",stats->num_alignments);
  /*
   * Total bases (aligned/trimmed/unaligned)
   */
  register const uint64_t total_bases_mmaps = stats->mmaps_match_ep->total_bases_aligned+stats->mmaps_match_ep->total_bases_trimmed;
  fprintf(stderr,"  --> Num.bases %lu\n",stats->total_bases);
  fprintf(stderr,"    --> Num.bases.aligned %lu/%lu (%2.3f%%/%2.3f%%)\n",
      stats->best_match_ep->total_bases_aligned,stats->mmaps_match_ep->total_bases_aligned,
      100.0*(float)stats->best_match_ep->total_bases_aligned/(float)stats->total_bases,
      100.0*(float)stats->mmaps_match_ep->total_bases_aligned/(float)total_bases_mmaps);
  fprintf(stderr,"    --> Num.bases.trimmed %lu/%lu (%2.3f%%/%2.3f%%)\n",
      stats->best_match_ep->total_bases_trimmed,stats->mmaps_match_ep->total_bases_trimmed,
      100.0*(float)stats->best_match_ep->total_bases_trimmed/(float)stats->total_bases,
      100.0*(float)stats->mmaps_match_ep->total_bases_trimmed/(float)total_bases_mmaps);
  register const uint64_t total_bases_unaligned_best=stats->total_bases -
      (stats->best_match_ep->total_bases_aligned+stats->best_match_ep->total_bases_trimmed);
  fprintf(stderr,"    --> Num.bases.unaligned %lu (%2.3f%%)\n",
      total_bases_unaligned_best,100.0*(float)total_bases_unaligned_best/(float)stats->total_bases);
  /*
   * Maps
   */
  fprintf(stderr,"  --> Num.maps %lu (%2.3f map/alg)\n",stats->num_maps,(float)stats->num_maps/(float)stats->num_mapped);
  gt_stats_print_mmap_distribution(stats->mmap,eff_num_reads,stats->num_mapped);
  if (paired_end) {
    gt_stats_print_inss_distribution(
        stats->best_match_ep->inss,stats->mmaps_match_ep->inss,stats->num_mapped,stats->num_maps);
  }
  gt_stats_print_uniq_distribution(stats->uniq,eff_num_reads);
  gt_stats_print_misms_distribution(stats->best_match_ep->misms,stats->mmaps_match_ep->misms,
      stats->num_mapped,stats->num_maps,"Mismatches");
  gt_stats_print_misms_distribution(stats->best_match_ep->indel,stats->mmaps_match_ep->indel,
      stats->num_mapped,stats->num_maps,"Indels.length");
  gt_stats_print_misms_distribution(stats->best_match_ep->errors,stats->mmaps_match_ep->errors,
      stats->num_mapped,stats->num_maps,"Error.events");
  gt_stats_print_misms_positions(stats->best_match_ep->error_position,stats->mmaps_match_ep->error_position,
      stats->best_match_ep->total_misms_events,stats->mmaps_match_ep->total_misms_events,stats->max_length);
}
void gt_stats_print_stats_compact(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  /*  #mapped, %mapped, #unmapped, %unmapped, MMap(maps/alg), Bases.aligned(%), Bases.trimmed(%)
   *  mapped(%), MMap(maps/alg), Bases.aligned(%), Bases.trimmed(%)
   */
  // #mapped, %mapped
  register const uint64_t eff_num_reads = paired_end?num_reads>>1:num_reads;
  fprintf(stderr,"%lu,",stats->num_mapped);
  fprintf(stderr,"%2.3f,",100.0*(float)stats->num_mapped/(float)eff_num_reads);
  // #unmapped, %unmapped
  register const uint64_t unmapped = eff_num_reads-stats->num_mapped;
  fprintf(stderr,"%lu,",unmapped);
  fprintf(stderr,"%2.3f,",100.0*(float)unmapped/(float)eff_num_reads);
  // MMap(maps/alg)
  fprintf(stderr,"%2.3f,",(float)stats->num_maps/(float)stats->num_mapped);
  // Bases.aligned(%)
  register const uint64_t total_bases_mmaps = stats->mmaps_match_ep->total_bases_aligned+stats->mmaps_match_ep->total_bases_trimmed;
  fprintf(stderr,"%2.3f,",100.0*(float)stats->mmaps_match_ep->total_bases_aligned/(float)total_bases_mmaps);
  // Bases.trimmed(%)
  fprintf(stderr,"%2.3f,",100.0*(float)stats->mmaps_match_ep->total_bases_trimmed/(float)total_bases_mmaps);
  // Uniq-0
  fprintf(stderr,"%2.3f\n",100.0*(float)stats->uniq[UNIQ_RANGE_0]/(float)eff_num_reads);
}

/*
 * CORE functions
 */
void gt_stats_parallel_generate_stats() {
  // Stats info
  gt_stats** stats = malloc(parameters.num_threads*sizeof(gt_stats*));

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    uint64_t tid = omp_get_thread_num();
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

    gt_status error_code;
    gt_template *template = gt_template_new();
    gt_stats_init(stats+tid);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,parameters.paired_end))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // Extract stats
      gt_stats_analize_template(stats[tid],template,parameters.paired_end);
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }

  // Merge stats
  gt_stats_merge(stats,parameters.num_threads);

  // Print Statistics
  if (!parameters.quiet) {
    if (!parameters.compact) {
      gt_stats_print_stats(stats[0],(parameters.num_reads>0)?
          parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
    } else {
      gt_stats_print_stats_compact(stats[0],(parameters.num_reads>0)?
          parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
    }
  }

  // Clean
  free(stats);
  gt_input_file_close(input_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.stats [ARGS]...\n"
                  "        --input|-i [FILE]\n"
                  "        --mmap-input\n"
                  "        --paired-end|p\n"
                  "        --num-reads|n\n"
                  "        --threads|t\n"
                  "        --compact|c\n"
                  "        --verbose|v\n"
                  "        --quiet|q\n"
                  "        --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    { "num-reads", no_argument, 0, 'n' },
    { "threads", no_argument, 0, 't' },
    { "compact", no_argument, 0, 'c' },
    { "verbose", no_argument, 0, 'v' },
    { "quiet", no_argument, 0, 'q' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:n:t:pchqv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 1:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 'n':
      parameters.num_reads = atol(optarg);
      break;
    case 'c':
      parameters.compact = true;
      break;
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'q':
      parameters.quiet = true;
      break;
    case 'h':
      usage();
      exit(1);
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  // Extract stats
  gt_stats_parallel_generate_stats();

  return 0;
}


