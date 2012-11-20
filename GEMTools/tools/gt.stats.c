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

typedef struct {
  char *name_input_file;
  bool mmap_input;
  bool paired_end;
  uint64_t num_reads;
  uint64_t num_threads;
  bool verbose;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .num_reads=0,
    .num_threads=1,
    .verbose=false,
};

typedef struct {
  uint64_t num_reads;
  uint64_t num_alignments;
  // Mapped/Maps
  uint64_t num_maps;
  uint64_t num_mapped;
  // MMap Distribution
  uint64_t mmap[MMAP_RANGE];
  // Uniq Distribution
  uint64_t uniq[UNIQ_RANGE];
  // Insert Size Distribution
  uint64_t inss[INSS_RANGE];
  // Mismatch/Indel Distribution
  uint64_t misms[MISMS_RANGE];
  uint64_t indel[MISMS_RANGE];
  uint64_t errors[MISMS_RANGE];
  uint64_t error_pos[200];
  uint64_t total_misms;
  uint64_t total_indels;
} gt_stats;

/*
 * STATS functions
 */
void gt_source_stats_init(gt_stats *stats) {
  // Mapped/Maps
  stats->num_reads=0;
  stats->num_alignments=0;
  stats->num_maps=0;
  stats->num_mapped=0;
  // MMaps
  register uint64_t j;
  for (j=0;j<MMAP_RANGE;++j) {
    stats->mmap[j] = 0;
  }
  for (j=0;j<INSS_RANGE;++j) {
    stats->inss[j] = 0;
  }
  // Misms
  for (j=0;j<MISMS_RANGE;++j) {
    stats->misms[j] = 0;
    stats->indel[j] = 0;
    stats->errors[j] = 0;
  }
  stats->total_misms=0;
  stats->total_indels=0;
}
#define GT_STATS_GET_PERCENTAGE_ERROR(percentage) (uint64_t)((double)percentage*one_per_cent)
void gt_source_score_misms(uint64_t* const misms_array,const uint64_t misms_num,const double one_per_cent) {
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
void gt_stats_get_misms_distribution(
    gt_stats* const stats,uint64_t const read_length,gt_map** const map,const uint64_t num_blocks) {
  uint64_t num_misms=0, num_indels=0, total=0, i;
  for (i=0;i<num_blocks;++i) {
    GT_MAP_ITERATE(map[i],map_block) {
      GT_MISMS_ITERATE(map_block,misms) {
        total++;
        switch (misms->misms_type) {
          case MISMS:
            num_misms++;
            break;
          case INS:
          case DEL:
            num_indels+=misms->size;
            break;
        }
      }
    }
  }
  // Record stats
  const double one_per_cent = (double)read_length/100.0;
  gt_source_score_misms(stats->misms,num_misms,one_per_cent);
  gt_source_score_misms(stats->indel,num_indels,one_per_cent);
  gt_source_score_misms(stats->errors,num_misms+num_indels,one_per_cent);
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
void gt_stats_get_inss_distribution(gt_stats *stats,const uint64_t insert_size) {
  if (insert_size<=100) {
    ++stats->inss[INSS_RANGE_100];
  } else if (insert_size<=200) {
    ++stats->inss[INSS_RANGE_200];
  } else if (insert_size<=300) {
    ++stats->inss[INSS_RANGE_300];
  } else if (insert_size<=400) {
    ++stats->inss[INSS_RANGE_400];
  } else if (insert_size<=500) {
    ++stats->inss[INSS_RANGE_500];
  } else if (insert_size<=600) {
    ++stats->inss[INSS_RANGE_600];
  } else if (insert_size<=700) {
    ++stats->inss[INSS_RANGE_700];
  } else if (insert_size<=800) {
    ++stats->inss[INSS_RANGE_800];
  } else if (insert_size<=900) {
    ++stats->inss[INSS_RANGE_900];
  } else if (insert_size<=1000) {
    ++stats->inss[INSS_RANGE_1000];
  } else if (insert_size<=2000) {
    ++stats->inss[INSS_RANGE_2000];
  } else if (insert_size<=5000) {
    ++stats->inss[INSS_RANGE_5000];
  } else {
    ++stats->inss[INSS_RANGE_BEHOND];
  }
}
GT_INLINE void gt_source_analize_template(gt_stats* const stats,gt_template* const template,const bool paired_map) {
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register const uint64_t num_maps = (num_blocks>1 || !paired_map) ?
      gt_template_get_num_mmaps(template) : 0;
  ++stats->num_alignments;
  stats->num_reads += num_blocks;
  if ( (paired_map &&  (num_maps > 0 || gt_template_is_mapped(template))) ||
       (!paired_map && (num_maps > 0 || gt_alignment_is_mapped(gt_template_get_block(template,0)))) ){
    ++stats->num_mapped;
    if (num_maps > 0) {
      stats->num_maps += num_maps;
      gt_stats_get_mmap_distribution(stats,num_maps); // MMaps classification
    }
  }
  gt_stats_get_uniq_distribution(stats,template);
  if (num_maps > 0) {
    gt_stats_get_misms_distribution(stats,gt_template_get_total_length(template),
        gt_template_get_mmap(template,0,NULL),num_blocks);
    if (paired_map) { // Insert Size distribution
      GT_TEMPLATE_ITERATE_(template,mmap) {
        gt_stats_get_inss_distribution(stats,
            abs((int64_t)mmap[0]->position-(int64_t)mmap[1]->position));
      }
    }
  }
}
void gt_stats_merge_stats(gt_stats* const stats,const uint64_t num_stats) {
  register uint64_t i, j;
  for (i=1;i<num_stats;++i) {
    stats->num_reads += stats[i].num_reads;
    stats->num_alignments += stats[i].num_alignments;
    stats->num_maps += stats[i].num_maps;
    stats->num_mapped += stats[i].num_mapped;
    for (j=0;j<MMAP_RANGE;++j) {
      stats->mmap[j] += stats[i].mmap[j];
    }
    for (j=0;j<UNIQ_RANGE;++j) {
      stats->uniq[j] += stats[i].uniq[j];
    }
    for (j=0;j<INSS_RANGE;++j) {
      stats->inss[j] += stats[i].inss[j];
    }
    for (j=0;j<MISMS_RANGE;++j) {
      stats->misms[j] += stats[i].misms[j];
      stats->indel[j] += stats[i].indel[j];
      stats->errors[j] += stats[i].errors[j];
    }
    stats->total_misms += stats[i].total_misms;
    stats->total_indels += stats[i].total_indels;
  }

}
void gt_stats_print_mmap_distribution(gt_stats* const stats,const uint64_t num_reads) {
  fprintf(stderr,"MMap.ranges\n");
  fprintf(stderr,"  -->        [0] \t=> %lu \t %f%%\n",(num_reads-stats->num_mapped),100.0*(float)(num_reads-stats->num_mapped)/(float)num_reads);
  fprintf(stderr,"  -->        [1] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_1],100.0*(float)stats->mmap[MMAP_RANGE_1]/(float)num_reads);
  fprintf(stderr,"  -->      (1,5] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_5],100.0*(float)stats->mmap[MMAP_RANGE_5]/(float)num_reads);
  fprintf(stderr,"  -->     (5,10] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_10],100.0*(float)stats->mmap[MMAP_RANGE_10]/(float)num_reads);
  fprintf(stderr,"  -->    (10,50] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_50],100.0*(float)stats->mmap[MMAP_RANGE_50]/(float)num_reads);
  fprintf(stderr,"  -->   (50,100] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_100],100.0*(float)stats->mmap[MMAP_RANGE_100]/(float)num_reads);
  fprintf(stderr,"  -->  (100,500] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_500],100.0*(float)stats->mmap[MMAP_RANGE_500]/(float)num_reads);
  fprintf(stderr,"  --> (500,1000] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_1000],100.0*(float)stats->mmap[MMAP_RANGE_1000]/(float)num_reads);
  fprintf(stderr,"  --> (1000,inf) \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_BEHOND],100.0*(float)stats->mmap[MMAP_RANGE_BEHOND]/(float)num_reads);
}
void gt_stats_print_uniq_distribution(gt_stats* const stats,const uint64_t num_reads) {
  fprintf(stderr,"Uniq.ranges\n");
  fprintf(stderr,"  -->        [X] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_X],100.0*(float)(stats->uniq[UNIQ_RANGE_X])/(float)num_reads);
  fprintf(stderr,"  -->        [0] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_0],100.0*(float)(stats->uniq[UNIQ_RANGE_0])/(float)num_reads);
  fprintf(stderr,"  -->        [1] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_1],100.0*(float)(stats->uniq[UNIQ_RANGE_1])/(float)num_reads);
  fprintf(stderr,"  -->        [2] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_2],100.0*(float)(stats->uniq[UNIQ_RANGE_2])/(float)num_reads);
  fprintf(stderr,"  -->        [3] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_3],100.0*(float)(stats->uniq[UNIQ_RANGE_3])/(float)num_reads);
  fprintf(stderr,"  -->     (3,10] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_10],100.0*(float)(stats->uniq[UNIQ_RANGE_10])/(float)num_reads);
  fprintf(stderr,"  -->    (10,50] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_50],100.0*(float)(stats->uniq[UNIQ_RANGE_50])/(float)num_reads);
  register const uint64_t accum = stats->uniq[UNIQ_RANGE_100]+stats->uniq[UNIQ_RANGE_500]+stats->uniq[UNIQ_RANGE_BEHOND];
  fprintf(stderr,"  -->   (50,inf] \t=> %lu \t %f%%\n",accum,100.0*(float)(accum)/(float)num_reads);
//  fprintf(stderr,"  -->   (50,100] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_100],100.0*(float)(stats->uniq[UNIQ_RANGE_100])/(float)num_reads);
//  fprintf(stderr,"  -->  (100,500] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_500],100.0*(float)(stats->uniq[UNIQ_RANGE_500])/(float)num_reads);
//  fprintf(stderr,"  -->  (500,inf] \t=> %lu \t %f%%\n",stats->uniq[UNIQ_RANGE_BEHOND],100.0*(float)(stats->uniq[UNIQ_RANGE_BEHOND])/(float)num_reads);
}
void gt_stats_print_inss_distribution(gt_stats* const stats,const uint64_t num_reads) {
  fprintf(stderr,"InsS.ranges\n");
  fprintf(stderr,"  -->       [0,100] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_100],100.0*(float)stats->inss[INSS_RANGE_100]/(float)num_reads);
  fprintf(stderr,"  -->     (100,200] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_200],100.0*(float)stats->inss[INSS_RANGE_200]/(float)num_reads);
  fprintf(stderr,"  -->     (200,300] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_300],100.0*(float)stats->inss[INSS_RANGE_300]/(float)num_reads);
  fprintf(stderr,"  -->     (300,400] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_400],100.0*(float)stats->inss[INSS_RANGE_400]/(float)num_reads);
  fprintf(stderr,"  -->     (400,500] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_500],100.0*(float)stats->inss[INSS_RANGE_500]/(float)num_reads);
  fprintf(stderr,"  -->     (500,600] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_600],100.0*(float)stats->inss[INSS_RANGE_600]/(float)num_reads);
  fprintf(stderr,"  -->     (600,700] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_700],100.0*(float)stats->inss[INSS_RANGE_700]/(float)num_reads);
  fprintf(stderr,"  -->     (700,800] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_800],100.0*(float)stats->inss[INSS_RANGE_800]/(float)num_reads);
  fprintf(stderr,"  -->     (800,900] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_900],100.0*(float)stats->inss[INSS_RANGE_900]/(float)num_reads);
  fprintf(stderr,"  -->    (900,1000] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_1000],100.0*(float)stats->inss[INSS_RANGE_1000]/(float)num_reads);
  fprintf(stderr,"  -->   (1000,2000] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_2000],100.0*(float)stats->inss[INSS_RANGE_2000]/(float)num_reads);
  fprintf(stderr,"  -->   (2000,5000] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_5000],100.0*(float)stats->inss[INSS_RANGE_5000]/(float)num_reads);
  fprintf(stderr,"  -->    (5000,inf] \t=> %lu \t %f%%\n",stats->inss[INSS_RANGE_BEHOND],100.0*(float)stats->inss[INSS_RANGE_BEHOND]/(float)num_reads);
}
void gt_stats_print_misms_distribution(uint64_t* const stats,const uint64_t num_reads,const char* const header) {
  fprintf(stderr,"%s\n",header);
  fprintf(stderr,"  -->         [0]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_0],100.0*(float)stats[MISMS_RANGE_0]/(float)num_reads);
  fprintf(stderr,"  -->         [1]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_1],100.0*(float)stats[MISMS_RANGE_1]/(float)num_reads);
  fprintf(stderr,"  -->         [2]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_2],100.0*(float)stats[MISMS_RANGE_2]/(float)num_reads);
  register const uint64_t accum = stats[MISMS_RANGE_3]+stats[MISMS_RANGE_4]+
      stats[MISMS_RANGE_5]+stats[MISMS_RANGE_6]+stats[MISMS_RANGE_7]+
      stats[MISMS_RANGE_8]+stats[MISMS_RANGE_9]+stats[MISMS_RANGE_10];
  fprintf(stderr,"  -->      (2,10]%% \t=> %lu \t %f%%\n",accum,100.0*(float)accum/(float)num_reads);
//  fprintf(stderr,"  -->         [3]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_3],100.0*(float)stats[MISMS_RANGE_3]/(float)num_reads);
//  fprintf(stderr,"  -->         [4]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_4],100.0*(float)stats[MISMS_RANGE_4]/(float)num_reads);
//  fprintf(stderr,"  -->         [5]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_5],100.0*(float)stats[MISMS_RANGE_5]/(float)num_reads);
//  fprintf(stderr,"  -->         [6]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_6],100.0*(float)stats[MISMS_RANGE_6]/(float)num_reads);
//  fprintf(stderr,"  -->         [7]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_7],100.0*(float)stats[MISMS_RANGE_7]/(float)num_reads);
//  fprintf(stderr,"  -->         [8]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_8],100.0*(float)stats[MISMS_RANGE_8]/(float)num_reads);
//  fprintf(stderr,"  -->         [9]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_9],100.0*(float)stats[MISMS_RANGE_9]/(float)num_reads);
//  fprintf(stderr,"  -->        [10]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_10],100.0*(float)stats[MISMS_RANGE_10]/(float)num_reads);
  fprintf(stderr,"  -->     (10,20]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_20],100.0*(float)stats[MISMS_RANGE_20]/(float)num_reads);
  fprintf(stderr,"  -->     (20,50]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_50],100.0*(float)stats[MISMS_RANGE_50]/(float)num_reads);
  fprintf(stderr,"  -->    (50,100]%% \t=> %lu \t %f%%\n",stats[MISMS_RANGE_BEHOND],100.0*(float)stats[MISMS_RANGE_BEHOND]/(float)num_reads);
}
void gt_stats_print_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  fprintf(stderr,"[GENERAL.STATS]\n");
  fprintf(stderr,"Num.reads %lu\n",num_reads);
  fprintf(stderr,"Num.alignments %lu\n",stats->num_alignments);
  fprintf(stderr,"  --> Num.mapped %lu (%2.3f%% / %2.3f%%)\n",stats->num_mapped,
      100.0*(float)stats->num_mapped/(float)stats->num_alignments,
      100.0*(float)stats->num_mapped/(float)(paired_end?num_reads>>1:num_reads));
  fprintf(stderr,"  --> Num.maps %lu (%2.3f map/alg)\n",stats->num_maps,(float)stats->num_maps/(float)stats->num_mapped);
  gt_stats_print_mmap_distribution(stats,stats->num_alignments);
  if (paired_end) gt_stats_print_inss_distribution(stats,stats->num_maps);
  gt_stats_print_uniq_distribution(stats,stats->num_mapped);
  gt_stats_print_misms_distribution(stats->misms,stats->num_mapped,"Mismatches");
  gt_stats_print_misms_distribution(stats->indel,stats->num_mapped,"Indels");
  gt_stats_print_misms_distribution(stats->errors,stats->num_mapped,"Errors");
}

/*
 * CORE functions
 */
void gt_stats_generate_stats() {
  // Stats info
  gt_stats stats;
  gt_source_stats_init(&stats);

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

  gt_status error_code;
  gt_template *template = gt_template_new();
  while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,parameters.paired_end))) {
    if (error_code!=GT_IMP_OK) {
      gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
    }

    // Extract stats
    gt_source_analize_template(&stats,template,parameters.paired_end);
  }
  // Print Statistics
  gt_stats_print_stats(&stats,(parameters.num_reads>0)?
      parameters.num_reads:stats.num_reads,parameters.paired_end);

  // Clean
  gt_template_delete(template,true,true);
  gt_buffered_input_file_close(buffered_input);
  gt_input_file_close(input_file);
}

void gt_stats_parallel_generate_stats() {
  // Stats info
  gt_stats* stats = malloc(parameters.num_threads*sizeof(gt_stats));

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
    gt_source_stats_init(stats+tid);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,parameters.paired_end))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // Extract stats
      gt_source_analize_template(stats+tid,template,parameters.paired_end);
    }

    // Clean
    gt_template_delete(template,true,true);
    gt_buffered_input_file_close(buffered_input);
  }

  // Merge stats
  gt_stats_merge_stats(stats,parameters.num_threads);

  // Print Statistics
  gt_stats_print_stats(stats,(parameters.num_reads>0)?
      parameters.num_reads:stats->num_alignments,parameters.paired_end);

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
                  "        --verbose|v\n"
                  "        --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    { "num-reads", no_argument, 0, 'n' },
    { "threads", no_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:n:t:phv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 0:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 'n':
      parameters.num_reads = atol(optarg);
      break;
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'v':
      parameters.verbose = true;
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
  //gt_stats_generate_stats();
  gt_stats_parallel_generate_stats();

  return 0;
}


