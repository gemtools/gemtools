/*
 * PROJECT:
 * FILE: gem-source.c
 * DATE:
 * DESCRIPTION:
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE false
#define MAX_FILES 100

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
#define MISMS_RANGE_100 13
#define MISMS_RANGE_BEHOND 14
#define MISMS_RANGE 15

typedef struct {
  char *name_input_file;
  uint64_t num_files;
  uint64_t num_reads;
  char *name_output_file;
  char *name_src_positions_file;
  uint64_t num_threads;
  bool pe;
  bool compact;
  bool output_diffs;
  bool soap_template;
  bool verbose;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .num_files=0,
    .num_reads=0,
    .name_input_file=NULL,
    .name_output_file=NULL,
    .name_src_positions_file=NULL,
    .num_threads=1,
    .pe=false,
    .compact=false,
    .verbose=false,
    .soap_template=false,
    .output_diffs=true
};

typedef struct {
  uint64_t num_alignments;
  // Maps
  uint64_t num_maps;
  // Mapped
  uint64_t num_mapped;
  // Hits (wrt source positions)
  uint64_t num_true_positives;
  uint64_t num_oracles; // (first map)
  uint64_t num_best_oracles;
  uint64_t num_inss_oracles;
  // MMap Distribution
  uint64_t mmap[MMAP_RANGE];
  // Insert Size Distribution
  uint64_t inss[INSS_RANGE];
  // Mismatch/Indel Distribution
  uint64_t misms[MISMS_RANGE];
  uint64_t indel[MISMS_RANGE];
  uint64_t errors[MISMS_RANGE];
  uint64_t error_pos[200];
  uint64_t total_misms;
  uint64_t total_indels;
} gt_source_stats;

typedef struct {
  uint64_t num_exclusive_src;
  uint64_t num_exclusive_cmp;
  uint64_t num_intersection;
  // GOD markers
  uint64_t discrete_god_hits;
  float float_god_hits;
  // Relative GOD markers (wrt MMaps distribution)
  float float_relative_god_hits;
} gt_intersection_stats;

/*
 * CMP functions
 */
int64_t gt_source_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  return gt_map_range_cmp(map_1,map_2,50);
}
int64_t gt_source_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  return gt_mmap_range_cmp(map_1,map_2,num_maps,50);
}
bool gt_mmap_is_proper_pair(gt_map* map_end1,gt_map* map_end2) {
  return (gt_string_eq(map_end1->seq_name,map_end2->seq_name) &&
      abs((int64_t)map_end1->position-(int64_t)map_end2->position) <= 1000);
}
bool gt_template_is_properly_mapped(gt_template * template) {
  if (gt_template_get_num_blocks(template)==1) return false;
  _GT_TEMPLATE_ITERATE(template,map_array) {
    if (gt_mmap_is_proper_pair(map_array[0],map_array[1])) return true;
  }
  return false;
}

/*
 * STATS functions
 */
void gt_source_stats_init(gt_source_stats *stats) {
  // General fields
  stats->num_alignments=0;
  // Maps
  stats->num_maps=0;
  // Mapped
  stats->num_mapped=0;
  // Hits (wrt source positions)
  stats->num_true_positives=0;
  stats->num_oracles=0; // (first map)
  stats->num_best_oracles=0; // (lower levenshtein distance)
  stats->num_inss_oracles=0;
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
void gt_intersection_stats_init(gt_intersection_stats *stats) {
  stats->num_exclusive_src=0;
  stats->num_exclusive_cmp=0;
  stats->num_intersection=0;
  // GOD markers
  stats->discrete_god_hits=0;
  stats->float_god_hits=0.0;
  // Relative GOD markers (wrt MMaps distribution)
  stats->float_relative_god_hits=0.0;
}
void gt_source_score_misms(uint64_t* misms_array,uint64_t misms_num) {
  if (misms_num==0) {
    misms_array[MISMS_RANGE_0]++;
  } else if (misms_num<=1) {
    misms_array[MISMS_RANGE_1]++;
  } else if (misms_num<=2) {
    misms_array[MISMS_RANGE_2]++;
  } else if (misms_num<=3) {
    misms_array[MISMS_RANGE_3]++;
  } else if (misms_num<=4) {
    misms_array[MISMS_RANGE_4]++;
  } else if (misms_num<=5) {
    misms_array[MISMS_RANGE_5]++;
  } else if (misms_num<=6) {
    misms_array[MISMS_RANGE_6]++;
  } else if (misms_num<=7) {
    misms_array[MISMS_RANGE_7]++;
  } else if (misms_num<=8) {
    misms_array[MISMS_RANGE_8]++;
  } else if (misms_num<=9) {
    misms_array[MISMS_RANGE_9]++;
  } else if (misms_num<=10) {
    misms_array[MISMS_RANGE_10]++;
  } else if (misms_num<=20) {
    misms_array[MISMS_RANGE_20]++;
  } else if (misms_num<=50) {
    misms_array[MISMS_RANGE_50]++;
  } else if (misms_num<=100) {
    misms_array[MISMS_RANGE_100]++;
  } else {
    misms_array[MISMS_RANGE_BEHOND]++;
  }
}
void gt_source_get_misms_distribution(gt_source_stats *stats,gt_map *map) {
  uint64_t num_misms=0, num_indels=0, total=0;
  GT_MISMS_ITERATE(map,misms) {
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
  // Record stats
  gt_source_score_misms(stats->misms,num_misms);
  gt_source_score_misms(stats->indel,num_indels);
  gt_source_score_misms(stats->errors,num_misms+num_indels);
}
void gt_source_get_mmap_distribution(gt_source_stats *stats,const uint64_t num_maps) {
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
void gt_source_get_inss_distribution(gt_source_stats *stats,const uint64_t insert_size) {
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
GT_INLINE void gt_source_analize_template(gt_source_stats *stats,gt_template *template,bool paired_map) {
  register const uint64_t num_maps = gt_template_get_num_mmap(template);
  stats->num_maps += num_maps;
  if (num_maps>0 || gt_template_is_mapped(template)) ++stats->num_mapped;
  // MMaps clasification
  gt_source_get_mmap_distribution(stats,num_maps);
  if (gt_template_get_num_mmap(template) > 0) {
    gt_source_get_misms_distribution(stats,
        gt_alignment_get_map(gt_template_get_block(template,0),0));
  }
  // Insert Size distribution
  if (paired_map) {
    _GT_TEMPLATE_ITERATE(template,mmap) {
      gt_source_get_inss_distribution(stats,
          abs((int64_t)mmap[0]->position-(int64_t)mmap[1]->position));
    }
  }
}
void gt_source_print_mmap_distribution(gt_source_stats *stats,uint64_t num_reads) {
  fprintf(stderr,"MMap.ranges\n");
  fprintf(stderr,"  -->        [1] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_1],100.0*(float)stats->mmap[MMAP_RANGE_1]/(float)num_reads);
  fprintf(stderr,"  -->      (1,5] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_5],100.0*(float)stats->mmap[MMAP_RANGE_5]/(float)num_reads);
  fprintf(stderr,"  -->     (5,10] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_10],100.0*(float)stats->mmap[MMAP_RANGE_10]/(float)num_reads);
  fprintf(stderr,"  -->    (10,50] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_50],100.0*(float)stats->mmap[MMAP_RANGE_50]/(float)num_reads);
  fprintf(stderr,"  -->   (50,100] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_100],100.0*(float)stats->mmap[MMAP_RANGE_100]/(float)num_reads);
  fprintf(stderr,"  -->  (100,500] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_500],100.0*(float)stats->mmap[MMAP_RANGE_500]/(float)num_reads);
  fprintf(stderr,"  --> (500,1000] \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_1000],100.0*(float)stats->mmap[MMAP_RANGE_1000]/(float)num_reads);
  fprintf(stderr,"  --> (1000,inf) \t=> %lu \t %f%%\n",stats->mmap[MMAP_RANGE_BEHOND],100.0*(float)stats->mmap[MMAP_RANGE_BEHOND]/(float)num_reads);
}
void gt_source_print_inss_distribution(gt_source_stats *stats,uint64_t num_reads) {
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
void gt_source_print_misms_distribution(uint64_t* stats,uint64_t num_reads,char *header) {
  fprintf(stderr,"%s\n",header);
  fprintf(stderr,"  -->         [0] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_0],100.0*(float)stats[MISMS_RANGE_0]/(float)num_reads);
  fprintf(stderr,"  -->         [1] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_1],100.0*(float)stats[MISMS_RANGE_1]/(float)num_reads);
  fprintf(stderr,"  -->         [2] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_2],100.0*(float)stats[MISMS_RANGE_2]/(float)num_reads);
  fprintf(stderr,"  -->         [3] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_3],100.0*(float)stats[MISMS_RANGE_3]/(float)num_reads);
  fprintf(stderr,"  -->         [4] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_4],100.0*(float)stats[MISMS_RANGE_4]/(float)num_reads);
  fprintf(stderr,"  -->         [5] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_5],100.0*(float)stats[MISMS_RANGE_5]/(float)num_reads);
  fprintf(stderr,"  -->         [6] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_6],100.0*(float)stats[MISMS_RANGE_6]/(float)num_reads);
  fprintf(stderr,"  -->         [7] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_7],100.0*(float)stats[MISMS_RANGE_7]/(float)num_reads);
  fprintf(stderr,"  -->         [8] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_8],100.0*(float)stats[MISMS_RANGE_8]/(float)num_reads);
  fprintf(stderr,"  -->         [9] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_9],100.0*(float)stats[MISMS_RANGE_9]/(float)num_reads);
  fprintf(stderr,"  -->        [10] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_10],100.0*(float)stats[MISMS_RANGE_10]/(float)num_reads);
  fprintf(stderr,"  -->     (10,20] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_20],100.0*(float)stats[MISMS_RANGE_20]/(float)num_reads);
  fprintf(stderr,"  -->     (20,50] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_50],100.0*(float)stats[MISMS_RANGE_50]/(float)num_reads);
  fprintf(stderr,"  -->    (50,100] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_100],100.0*(float)stats[MISMS_RANGE_100]/(float)num_reads);
  fprintf(stderr,"  -->   (100,inf] \t=> %lu \t %f%%\n",stats[MISMS_RANGE_BEHOND],100.0*(float)stats[MISMS_RANGE_BEHOND]/(float)num_reads);
}
void gt_source_print_general_stats(gt_source_stats *stats,uint64_t num_reads,bool paired_end) {
  if (parameters.compact) {
    // CMP_MAPPED($1)    PER_CMP_MAPPED($2)    NUM_MAPS($7)
    fprintf(stdout,"\t%lu\t%2.5f\t%lu\n",
        stats->num_mapped,(float)stats->num_mapped/(float)num_reads,stats->num_maps);
  } else {
    fprintf(stderr,"[GENERAL.STATS]\n");
    fprintf(stderr,"Num.reads %lu\n",num_reads);
    fprintf(stderr,"Num.alignments %lu\n",stats->num_alignments);
    fprintf(stderr,"  --> Num.mapped %lu (%2.3f%%)\n",stats->num_mapped,100.0*(float)stats->num_mapped/(float)num_reads);
    fprintf(stderr,"  --> Num.maps %lu (%2.3f%%)\n",stats->num_maps,100.0*(float)stats->num_maps/(float)stats->num_alignments);
    gt_source_print_mmap_distribution(stats,num_reads);
    if (paired_end) {
      gt_source_print_inss_distribution(stats,num_reads);
    }
    gt_source_print_misms_distribution(stats->misms,num_reads,"Mismatches");
    gt_source_print_misms_distribution(stats->indel,num_reads,"Indels");
    gt_source_print_misms_distribution(stats->errors,num_reads,"Errors");
  }
}
void gt_source_test_print_stats(gt_source_stats *stats_src,gt_source_stats *stats_cmp,uint64_t num_reads,bool paired_end) {
  if (parameters.compact) {
    // CMP_MAPPED($1)    PER_CMP_MAPPED($2)   CMP_HITS($3)    PER_CMP_HITS($4)   ORACLES($5)
    //    PER_ORACLES($6)    BEST_ORACLES($5)   PER_BEST_ORACLES($6)   NUM_MAPS($7)
    fprintf(stdout,"\t%lu\t%2.5f\t%lu\t%2.5f\t%lu\t%2.5f\t%lu\t%2.5f\t%lu\n",
        stats_cmp->num_mapped,(float)stats_cmp->num_mapped/(float)num_reads,
        stats_cmp->num_true_positives,(float)stats_cmp->num_true_positives/(float)num_reads,
        stats_cmp->num_oracles,(float)stats_cmp->num_oracles/(float)num_reads,
        stats_cmp->num_best_oracles,(float)stats_cmp->num_best_oracles/(float)num_reads,
        stats_cmp->num_maps);
  } else {
    fprintf(stderr,"[SOURCE.STATS]\n");
    fprintf(stderr,"Num.reads %lu\n",num_reads);
    fprintf(stderr,"Num.SRC.mapped %lu (%2.3f%%)\n",stats_src->num_mapped,100.0*(float)stats_src->num_mapped/(float)num_reads);
    fprintf(stderr,"  --> Num.SRC.maps %lu (%2.3f%%)\n",stats_src->num_maps,100.0*(float)stats_src->num_maps/(float)num_reads);
    fprintf(stderr,"Num.CMP.mapped %lu (%2.3f%%)\n",stats_cmp->num_mapped,100.0*(float)stats_cmp->num_mapped/(float)num_reads);
    fprintf(stderr,"  --> Num.CMP.maps %lu (%2.3f%%)\n",stats_cmp->num_maps,100.0*(float)stats_cmp->num_maps/(float)num_reads);
    fprintf(stderr,"    --> Num.CMP.TP %lu (%2.3f%%)\n",stats_cmp->num_true_positives,100.0*(float)stats_cmp->num_true_positives/(float)num_reads);
    fprintf(stderr,"    --> Num.CMP.ORACLE %lu (%2.3f%%)\n",stats_cmp->num_oracles,100.0*(float)stats_cmp->num_oracles/(float)num_reads);
    fprintf(stderr,"    --> Num.CMP.BEST_ORACLE %lu (%2.3f%%)\n",stats_cmp->num_best_oracles,100.0*(float)stats_cmp->num_best_oracles/(float)num_reads);
    fprintf(stderr,"    --> Num.CMP.INS_ORACLE %lu (%2.3f%%)\n",stats_cmp->num_inss_oracles,100.0*(float)stats_cmp->num_inss_oracles/(float)num_reads);
    fprintf(stderr,"[SRC.DISTRIBUTIONS]\n");
    gt_source_print_mmap_distribution(stats_src,num_reads);
    if (paired_end) gt_source_print_inss_distribution(stats_src,num_reads);
    fprintf(stderr,"[CMP.DISTRIBUTIONS]\n");
    gt_source_print_mmap_distribution(stats_cmp,num_reads);
    if (paired_end) gt_source_print_inss_distribution(stats_cmp,num_reads);
  }
}

void gt_inclusion_test_print_stats(gt_source_stats *stats_src,gt_source_stats *stats_cmp,gt_intersection_stats *intersection_stats) {

}

/*
 * Helper functions
 */
GT_INLINE gt_status gt_source_input_generic_parser_get_template(
    gt_buffered_input_file* buffered_input,gt_template *template,bool paired_end,bool soap) {
  gt_status error_code;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      if ((error_code=gt_input_map_parser_get_template(buffered_input,template))==GT_BMI_FAIL) return GT_BMI_FAIL;
      if (gt_template_get_num_blocks(template)==1 && paired_end) {
        if ((error_code=gt_input_map_parser_get_alignment(buffered_input,gt_template_dyn_get_block(template,1)))==GT_BMI_FAIL) return GT_BMI_FAIL;
      }
      break;
    case SAM:
    default: /* SAM */
      if (paired_end) {
        if (soap) {
          if ((error_code=gt_input_sam_parser_get_soap_template(buffered_input,template))==GT_BMI_FAIL) return GT_BMI_FAIL;
        } else {
          if ((error_code=gt_input_sam_parser_get_template(buffered_input,template))==GT_BMI_FAIL) return GT_BMI_FAIL;
        }
      } else {
        if ((error_code=gt_input_sam_parser_get_alignment(buffered_input,gt_template_dyn_get_block(template,0)))==GT_BMI_FAIL) return GT_BMI_FAIL;
      }
      break;
  }
  return error_code;
}
#define GT_SOURCE_READ_INPUT(buffered_input,template) \
    if ((error_code=gt_source_input_generic_parser_get_template( \
        buffered_input,template,parameters.pe,parameters.soap_template))==GT_BMI_FAIL) { gt_error_msg("Bad template"); exit(0); }
gt_status gt_source_read__synch_files(
    gt_buffered_input_file* buffered_input_src,gt_buffered_input_file* buffered_input_cmp,
    gt_template *template_src,gt_template *template_cmp,
    gt_source_stats *stats_src,gt_source_stats *stats_cmp) {
  gt_status error_code;
  // Read SRC
  GT_SOURCE_READ_INPUT(buffered_input_src,template_src);
  if (error_code==0) return 0;
  if (stats_src) ++stats_src->num_alignments;
  // Read CMP
  GT_SOURCE_READ_INPUT(buffered_input_cmp,template_cmp);
  if (error_code==0) return 0;
  if (stats_cmp) ++stats_cmp->num_alignments;
  while (!gt_string_eq(gt_template_get_tag(template_cmp),gt_template_get_tag(template_cmp))) {
    GT_SOURCE_READ_INPUT(buffered_input_src,template_src);
    if (error_code==0) return 0;
    if (stats_src) ++stats_src->num_alignments;
  }
  return error_code;
}

/*
 * CORE functions
 */
void gt_source_translate() {
  // Stats info
  gt_source_stats stats;
  gt_source_stats_init(&stats);

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

  gt_template *template = gt_template_new();
  gt_status error_code;
  register bool keep_reading = true;
  while (keep_reading) {
    // Read template
    GT_SOURCE_READ_INPUT(buffered_input,template);
    if (error_code==0) break;
    // Print template
    if (parameters.pe) {
      gt_output_map_fprint_template(stdout,template,GT_ALL,false);
    } else {
      gt_output_map_fprint_alignment(stdout,gt_template_get_block(template,0),GT_ALL,false);
    }
  }
  // Close buffered readers & files
  gt_buffered_input_file_close(buffered_input);
  gt_input_file_close(input_file);
}
void gt_source_generate_stats() {
  // Stats info
  gt_source_stats stats;
  gt_source_stats_init(&stats);

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

  gt_template *template = gt_template_new();
  gt_status error_code;
  register bool keep_reading = true;
  while (keep_reading) {
    // Read template
    GT_SOURCE_READ_INPUT(buffered_input,template);
    if (error_code==0) break;
    ++stats.num_alignments;
    // Extract stats
    gt_source_analize_template(&stats,template,parameters.pe);
  }
  // Print Statistics
  gt_source_print_general_stats(
      &stats,(parameters.num_reads>0)?parameters.num_reads:stats.num_alignments,parameters.pe);

  // Close buffered readers & files
  gt_buffered_input_file_close(buffered_input);
  gt_input_file_close(input_file);
}
gt_map** gt_source_get_best_mmap(gt_template* template) {
  gt_map** mmap = gt_template_get_mmap(template,0,NULL);
  _GT_TEMPLATE_ITERATE(template,mmap_it) {
    register const uint64_t lev_cmp = gt_map_get_levenshtein_distance(mmap[0])+gt_map_get_levenshtein_distance(mmap[1]);
    register const uint64_t lev_it = gt_map_get_levenshtein_distance(mmap_it[0])+gt_map_get_levenshtein_distance(mmap_it[1]);
    if (lev_it<lev_cmp) {
      mmap = mmap_it;
    }
  }
  return mmap;
}
gt_map** gt_source_get_inss_mmap(gt_template* template) {
  gt_map** mmap = gt_template_get_mmap(template,0,NULL);
  _GT_TEMPLATE_ITERATE(template,mmap_it) {
    register const uint64_t inss_cmp = abs((int64_t)mmap[0]->position-(int64_t)mmap[1]->position);
    register const uint64_t inss_it = abs((int64_t)mmap_it[0]->position-(int64_t)mmap_it[1]->position);
    if (inss_it<inss_cmp) {
      mmap = mmap_it;
    }
  }
  return mmap;
}
void gt_source_test() {
  // Stats info
  gt_source_stats stats_src;
  gt_source_stats stats_cmp;
  gt_source_stats_init(&stats_src);
  gt_source_stats_init(&stats_cmp);

  // Open files & buffers
  gt_input_file* input_file_src = (parameters.name_src_positions_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_src_positions_file,false);
  gt_input_file* input_file_cmp = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_buffered_input_file* buffered_input_src = gt_buffered_input_file_new(input_file_src);
  gt_buffered_input_file* buffered_input_cmp = gt_buffered_input_file_new(input_file_cmp);

  // Template
  gt_template *template_src = gt_template_new();
  gt_template *template_cmp = gt_template_new();
  gt_status error_code;
  register bool keep_reading = true;
  while (keep_reading) {
    /*
     * Read & Synch files
     */
    error_code=gt_source_read__synch_files(
        buffered_input_src,buffered_input_cmp,template_src,template_cmp,&stats_src,&stats_cmp);
    if (error_code==0) break;
    /*
     * Source-test
     */
    gt_source_analize_template(&stats_cmp,template_cmp,parameters.pe);
    gt_source_analize_template(&stats_src,template_src,parameters.pe);
    if (parameters.pe) {
      register bool oracle_hit=false, tp=false;
      register uint64_t oracle_pos;
      gt_map** mmap_src = gt_template_get_mmap(template_src,0,NULL);
      // Oracle Hits
      if (gt_template_get_num_mmap(template_cmp)>0) {
        if (gt_source_mmap_cmp(mmap_src,gt_template_get_mmap(template_cmp,0,NULL),2)==0) { ++stats_cmp.num_oracles; oracle_hit = true;}
        if (gt_source_mmap_cmp(mmap_src,gt_source_get_best_mmap(template_cmp),2)==0) {++stats_cmp.num_best_oracles;}
        if (gt_source_mmap_cmp(mmap_src,gt_source_get_inss_mmap(template_cmp),2)==0) {++stats_cmp.num_inss_oracles;}
      }
      // TP Hits
      if (gt_template_is_mmap_contained_fx(gt_source_mmap_cmp,template_cmp,mmap_src)) {
        ++stats_cmp.num_true_positives;
        tp = true; oracle_pos=1;
        _GT_TEMPLATE_ITERATE(template_cmp,mmap) {
          if (gt_source_mmap_cmp(gt_template_get_mmap(template_src,0,NULL),mmap,2)==0) break;
          ++oracle_pos;
        }
      }
      // Print Diffs
      if (parameters.output_diffs && !oracle_hit) {
        fprintf(stdout,"%s\t",oracle_hit?"YES":"NO");
        if (tp) {
          fprintf(stdout,"YES{%lu}\t",oracle_pos);
        } else {
          fprintf(stdout,"NO\t");
        }
        fprintf(stdout,"%s\t",template_cmp->tag);
        gt_output_map_fprint_template_maps(stdout,template_src,GT_ALL,false);
        fprintf(stdout,"\t");
        gt_output_map_fprint_counters(stdout,template_cmp->counters,template_cmp->max_complete_strata,false);
        fprintf(stdout,"\t");
        gt_output_map_fprint_template_maps(stdout,template_cmp,GT_ALL,false);
        fprintf(stdout,"\n");
      }
    } else {
      register gt_alignment* alignment_src = gt_template_get_block(template_src,0);
      register gt_alignment* alignment_cmp = gt_template_get_block(template_cmp,0);
      if (gt_alignment_get_num_maps(alignment_cmp)>0 &&
          gt_source_map_cmp(gt_alignment_get_map(alignment_src,0),
                            gt_alignment_get_map(alignment_cmp,0))==0) {
        ++stats_cmp.num_oracles;
      } else if (parameters.output_diffs) {
        fprintf(stdout,"%s\t",alignment_cmp->tag);
        gt_output_map_fprint_alignment_maps(stdout,alignment_src,GT_ALL,false);
        fprintf(stdout,"\t");
        gt_output_map_fprint_counters(stdout,alignment_cmp->counters,alignment_cmp->max_complete_strata,false);
        fprintf(stdout,"\t");
        gt_output_map_fprint_alignment_maps(stdout,alignment_cmp,GT_ALL,false);
        fprintf(stdout,"\n");
      }
      if (gt_alignment_is_map_contained_fx(gt_source_map_cmp,alignment_cmp,gt_alignment_get_map(alignment_src,0))) {
        ++stats_cmp.num_true_positives;
      }
    }
  }
  // Read the rest
  while (keep_reading) {
    // Read template
    GT_SOURCE_READ_INPUT(buffered_input_src,template_src);
    if (error_code==0) break;
    ++stats_src.num_alignments;
    // Extract stats
    gt_source_analize_template(&stats_src,template_src,parameters.pe);
  }


  // Close files & buffered readers
  gt_buffered_input_file_close(buffered_input_src);
  gt_buffered_input_file_close(buffered_input_cmp);
  gt_input_file_close(input_file_src);
  gt_input_file_close(input_file_cmp);

  // Print stats
  gt_source_test_print_stats(&stats_src,&stats_cmp,stats_src.num_alignments,parameters.pe);
}

void gt_source_merge() {
  // Open files & buffered readers
  gt_input_file *input_file_base, *input_file_src;
  gt_buffered_input_file *buffered_input_base, *buffered_input_src;
  input_file_base = gt_input_stream_open(stdin);
  input_file_src  = gt_input_file_open(parameters.name_input_file,false);
  buffered_input_base = gt_buffered_input_file_new(input_file_base);
  buffered_input_src = gt_buffered_input_file_new(input_file_src);

  // Declare templates
  gt_template *template_base = gt_template_new();
  gt_template *template_src = gt_template_new();
  gt_status error_code;
  register bool keep_reading = true;
  while (keep_reading) {
    /*
     * Read & Synch files
     */
    // Read SRC
    GT_SOURCE_READ_INPUT(buffered_input_base,template_base);
    if (error_code==0) break;
    // Read CMP
    GT_SOURCE_READ_INPUT(buffered_input_src,template_src);
    if (error_code==0) break;
    while (!gt_string_eq(gt_template_get_tag(template_base),gt_template_get_tag(template_src))) {
      GT_SOURCE_READ_INPUT(buffered_input_base,template_base);
      if (error_code==0){ gt_error_msg("Bad Master Synch"); exit(0);}
    }
    if (error_code==0) break;
    /*
     * Merge Templates
     */
    if (parameters.pe) {
      // Merge templates
      gt_template* merged_template = gt_template_copy(template_base);
      gt_template_merge_template_mmaps_fx(gt_source_mmap_cmp,merged_template,template_src);
      // Output merged
      gt_output_map_fprint_template(stdout,merged_template,GT_ALL,false);
      gt_template_delete_handler(merged_template);
    } else {
      register gt_alignment* alignment_base = gt_template_get_block(template_base,0);
      register gt_alignment* alignment_src = gt_template_get_block(template_src,0);
      // Merge alignments
      gt_alignment* merged_alignment = gt_alignment_copy(alignment_base);
      gt_alignment_merge_alignment_maps_fx(gt_source_map_cmp,merged_alignment,alignment_src);
      // Output merged
      gt_output_map_fprint_alignment(stdout,merged_alignment,GT_ALL,false);
      gt_alignment_delete_handler(merged_alignment);
    }
  }
  // Read the remaining reads
  while (keep_reading) {
    GT_SOURCE_READ_INPUT(buffered_input_base,template_base);
    if (error_code==0) break;
    if (parameters.pe) {
      gt_output_map_fprint_template(stdout,template_base,GT_ALL,false);
    } else {
      register gt_alignment* alignment_base = gt_template_get_block(template_base,0);
      gt_output_map_fprint_alignment(stdout,alignment_base,GT_ALL,false);
    }
  }
  // Fire to the hole!!
}

void gt_source_inclusion() {
  // Stats info
  gt_source_stats stats_src;
  gt_source_stats stats_cmp;
  gt_intersection_stats intersection_stats;
  gt_source_stats_init(&stats_src);
  gt_source_stats_init(&stats_cmp);

  // Open files & buffers
  gt_input_file* input_file_src = (parameters.name_src_positions_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_src_positions_file,false);
  gt_input_file* input_file_cmp = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_buffered_input_file* buffered_input_src = gt_buffered_input_file_new(input_file_src);
  gt_buffered_input_file* buffered_input_cmp = gt_buffered_input_file_new(input_file_cmp);

  // Auxiliary variables
  gt_alignment *excl_src_alignment = gt_alignment_new();
  gt_alignment *excl_cmp_alignment = gt_alignment_new();
  gt_alignment *intersect_alignment = gt_alignment_new();
  gt_template *excl_src_template = gt_template_new();
  gt_template *excl_cmp_template = gt_template_new();
  gt_template *intersect_template = gt_template_new();

  // Template
  gt_template *template_src = gt_template_new();
  gt_template *template_cmp = gt_template_new();
  gt_status error_code;
  register bool keep_reading = true;
  while (keep_reading) {
    /*
     * Read & Synch files
     */
    error_code=gt_source_read__synch_files(
        buffered_input_src,buffered_input_cmp,template_src,template_cmp,&stats_cmp,&stats_src);
    /*
     * Exclusion/Inclusion
     */
    gt_source_analize_template(&stats_src,template_src,parameters.pe);
    gt_source_analize_template(&stats_cmp,template_cmp,parameters.pe);
    if (parameters.pe) {
      /*
       * Exclusion/Inclusion
       */
      gt_template_handler_dup(intersect_template,template_src);
      gt_template_handler_dup(excl_src_template,template_src);
      gt_template_handler_dup(excl_cmp_template,template_src);
      // Intersection with SRC
      gt_template_intersect_template_mmaps_fx(gt_source_mmap_cmp,intersect_template,template_src,template_cmp);
      // Exclusion sets
      gt_template_subtract_template_mmaps_fx(gt_source_mmap_cmp,excl_src_template,template_src,template_cmp);
      gt_template_subtract_template_mmaps_fx(gt_source_mmap_cmp,excl_cmp_template,template_cmp,template_src);
      /*
       * Alignment stats
       */
      intersection_stats.num_intersection += gt_template_get_num_mmap(intersect_template);
      intersection_stats.num_exclusive_src += gt_template_get_num_mmap(excl_src_template);
      intersection_stats.num_exclusive_cmp += gt_template_get_num_mmap(excl_cmp_template);
      register const uint64_t num_maps_src = gt_template_get_num_mmap(template_src);
      register const uint64_t num_maps_cmp = gt_template_get_num_mmap(template_cmp);
      if (num_maps_src>0) {
        register const uint64_t num_maps_intersection = gt_template_get_num_mmap(intersect_template);
        if (num_maps_intersection==num_maps_src) ++intersection_stats.discrete_god_hits;
        float mmaps_proportion = (float)num_maps_intersection/(float)num_maps_src;
        intersection_stats.float_god_hits += mmaps_proportion;
        intersection_stats.float_relative_god_hits += num_maps_src*mmaps_proportion;
      }
//      if (parameters.output_diffs && gt_template_get_num_mmap(excl_src_template)>0) {
      if (parameters.output_diffs && num_maps_src>0 && num_maps_cmp==0) {
        fprintf(stderr,"%s\t",template_src->tag);
        gt_output_map_fprint_counters(stderr,template_src->counters,0,false);
        fprintf(stderr,"\t");
        gt_output_map_fprint_template_maps(stderr,excl_src_template,GT_ALL,false);
        fprintf(stderr,"\t");
        gt_output_map_fprint_template_maps(stderr,template_cmp,GT_ALL,false);
        fprintf(stderr,"\n");
      }
      // Clean aux alingments
      gt_template_clear_handler(intersect_template);
      gt_template_clear_handler(excl_src_template);
      gt_template_clear_handler(excl_cmp_template);
    } else {
      register gt_alignment* alignment_src = gt_template_get_block(template_src,0);
      register gt_alignment* alignment_cmp = gt_template_get_block(template_cmp,0);
      // Intersection with SRC
      gt_alignment_intersect_alignment_maps_fx(gt_source_map_cmp,intersect_alignment,alignment_src,alignment_cmp);
      // Exclusion sets
      gt_alignment_subtract_alignment_maps_fx(gt_source_map_cmp,excl_src_alignment,alignment_src,alignment_cmp);
      gt_alignment_subtract_alignment_maps_fx(gt_source_map_cmp,excl_cmp_alignment,alignment_cmp,alignment_src);
      /*
       * Alignment stats
       */
      intersection_stats.num_intersection += gt_alignment_get_num_maps(intersect_alignment);
      intersection_stats.num_exclusive_src += gt_alignment_get_num_maps(excl_src_alignment);
      intersection_stats.num_exclusive_cmp += gt_alignment_get_num_maps(excl_cmp_alignment);
      register const uint64_t num_maps_src = gt_alignment_get_num_maps(alignment_src);
      register const uint64_t num_maps_cmp = gt_alignment_get_num_maps(alignment_cmp);
      if (num_maps_src>0) {
        register const uint64_t num_maps_intersection = gt_alignment_get_num_maps(intersect_alignment);
        if (num_maps_intersection==num_maps_src) ++intersection_stats.discrete_god_hits;
        float mmaps_proportion = (float)num_maps_intersection/(float)num_maps_src;
        intersection_stats.float_god_hits += mmaps_proportion;
        intersection_stats.float_relative_god_hits += num_maps_src*mmaps_proportion;
      }
//      if (parameters.output_diffs && gt_alignment_get_num_maps(excl_src_alignment)>0) {
      if (parameters.output_diffs && num_maps_src>0 && alignment_cmp==0) {
        fprintf(stderr,"%s\t",alignment_src->tag);
        gt_output_map_fprint_counters(stderr,alignment_src->counters,0,false);
        fprintf(stderr,"\t");
        gt_output_map_fprint_alignment_maps(stderr,excl_src_alignment,GT_ALL,false);
        fprintf(stderr,"\t");
        gt_output_map_fprint_alignment_maps(stderr,alignment_cmp,GT_ALL,false);
        fprintf(stderr,"\n");
      }
      // Clean aux alingments
      gt_alignment_clear_handler(intersect_alignment);
      gt_alignment_clear_handler(excl_src_alignment);
      gt_alignment_clear_handler(excl_cmp_alignment);
    }
    // Progress
    if (parameters.verbose && stats_src.num_alignments%100000==0) {
      fprintf(stderr,".. done %lu\n",stats_src.num_alignments);
    }
  }

  // Close files & buffered readers
  gt_buffered_input_file_close(buffered_input_src);
  gt_buffered_input_file_close(buffered_input_cmp);
  gt_input_file_close(input_file_src);
  gt_input_file_close(input_file_cmp);

  // Print stats
  gt_inclusion_test_print_stats(&stats_src,&stats_cmp,&intersection_stats);
}


void usage() {
  fprintf(stderr, "USE: ./gem-source [COMMAND] [ARGS]...\n"
                  "      Commands::\n"
                  "        tr -i [FILE]\n"
                  "        source-test --src [FILE] -i [FILE]\n"
                  "        stats -i [FILE] \n"
                  "        merge --src [FILE] -i [FILE] \n"
                  "        inclusion -i [FILE] \n"
                  "      Options::\n"
                  "        --input|-i [FILE]\n"
                  "        --src|-s   [FILE]\n"
                  "        --paired-end|-p \n"
                  "        --soap \n"
                  "        --no-diffs \n"
                  "        --compact|c \n"
                  "        --verbose|v \n"
                  "        --help|h    \n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "src", required_argument, 0, 's' },
    { "no-diffs", no_argument, 0, 1 },
    { "soap", no_argument, 0, 2 },
    { "num-reads", no_argument, 0, 'r' },
    { "paired-end", no_argument, 0, 'p' },
    { "compact", no_argument, 0, 'c' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:s:r:T:hcpv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 1: /* --no-diffs */
      parameters.output_diffs = false;
      break;
    case 2: /* --soap */
      parameters.soap_template = true;
      break;
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 's':
      parameters.name_src_positions_file = optarg;
      break;
    case 'p':
      parameters.pe = true;
      break;
    case 'T':
      parameters.num_threads = atol(optarg);
      break;
    case 'r':
      parameters.num_reads = atol(optarg);
      break;
    case 'c':
      parameters.compact = true;
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  register char* const example_name = argv[1];
  if (argc==1) { usage(); exit(0); }
  parse_arguments(argc-1, argv+1);

  // X -> MAP
  if (strcmp(example_name,"tr")==0) { // stats
    gt_source_translate();
    return 0;
  }

  // General Stats
  if (strcmp(example_name,"stats")==0) { // stats
    gt_source_generate_stats();
    return 0;
  }

  // Stats wrt source
  if (strcmp(example_name,"source-test")==0) {
    gt_source_test();
    return 0;
  }

  // Merge alignments
  if (strcmp(example_name,"merge")==0) {
    gt_source_merge();
    return 0;
  }

  // Intersections Stats
  if (strcmp(example_name,"inclusion")==0) {
    gt_source_inclusion();
    return 0;
  }

  gt_error_msg("Incorrect test name provided");
  usage();
  return -1;
}


