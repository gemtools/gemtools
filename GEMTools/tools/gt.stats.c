/*
 * PROJECT: GEM-Tools library
 * FILE: gt.stats.c
 * DATE: 02/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to retrieve very naive stats from {MAP,SAM,FASTQ} files
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef struct {
  /* [Input] */
  char *name_input_file;
  bool mmap_input;
  bool paired_end;
  uint64_t num_reads;
  bool best_map;
  /* [Tests] */
  bool error_profile;
  bool mismatch_transitions;
  bool mismatch_quality;
  bool splitmaps_profile;
  /* [Output] */
  bool compact;
  bool verbose;
  bool quiet;
  /* [Misc] */
  uint64_t num_threads;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .num_reads=0,
    .best_map=false,
    .error_profile = false,
    .mismatch_transitions = false,
    .mismatch_quality = false,
    .splitmaps_profile = false,
    .compact = false,
    .verbose=false,
    .quiet=false,
    .num_threads=1,
};

/*
 * STATS Print results
 */
void gt_stats_print_mmap_distribution(uint64_t* const mmap,const uint64_t num_alignments,const uint64_t num_mapped) {
#define GT_STATS_PRINT_MMAP_FORMAT "%6" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_MMAP(RANGE) mmap[RANGE],100.0*(float)mmap[RANGE]/(float)num_alignments
  if(!num_alignments) return;
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
#define GT_STATS_PRINT_UNIQ_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_UNIQ(RANGE) uniq[RANGE],100.0*(float)uniq[RANGE]/(float)num_alignments
  if(!num_alignments) return;
  fprintf(stderr,"Uniq.ranges\n");
  fprintf(stderr,"  -->        [X] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_X));
  fprintf(stderr,"  -->        [0] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_0));
  fprintf(stderr,"  -->        [1] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_1));
  fprintf(stderr,"  -->        [2] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_2));
  fprintf(stderr,"  -->        [3] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_3));
  fprintf(stderr,"  -->     (3,10] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_10));
  fprintf(stderr,"  -->    (10,50] \t=> "GT_STATS_PRINT_UNIQ_FORMAT,GT_STATS_PRINT_UNIQ(UNIQ_RANGE_50));
  register const uint64_t accum = uniq[UNIQ_RANGE_100]+uniq[UNIQ_RANGE_500]+uniq[UNIQ_RANGE_BEHOND];
  fprintf(stderr,"  -->   (50,inf) \t=> "GT_STATS_PRINT_UNIQ_FORMAT,accum,100.0*(float)(accum)/(float)num_alignments);
}
void gt_stats_print_inss_distribution(uint64_t* const inss,const uint64_t num_maps) {
#define GT_STATS_PRINT_INSS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_INSS(RANGE) inss[RANGE],100.0*(float)inss[RANGE]/(float)num_maps
  if(!num_maps) return;
  fprintf(stderr,"InsS.ranges\n");
  fprintf(stderr,"  -->   (-inf,-100] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_NEG));
  fprintf(stderr,"  -->      (-100,0] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_OVER));
  fprintf(stderr,"  -->       (0,100] \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_100));
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
  fprintf(stderr,"  -->    (5000,inf) \t=> "GT_STATS_PRINT_INSS_FORMAT,GT_STATS_PRINT_INSS(INSS_RANGE_BEHOND));
}
void gt_stats_print_misms_distribution(uint64_t* const misms,const uint64_t num_maps,const char* const header) {
#define GT_STATS_PRINT_MISMS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_MISMS(RANGE) misms[RANGE],100.0*(float)misms[RANGE]/(float)num_maps
#define GT_STATS_PRINT_MISMS_VAL(misms_val) misms_val,100.0*(float)misms_val/(float)num_maps
  if(!num_maps) return;
  fprintf(stderr,"%s\n",header);
  fprintf(stderr,"  -->         [0]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_0));
  fprintf(stderr,"  -->         [1]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_1));
  fprintf(stderr,"  -->         [2]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS(MISMS_RANGE_2));
  register const uint64_t misms_accum = misms[MISMS_RANGE_3]+misms[MISMS_RANGE_4]+
      misms[MISMS_RANGE_5]+misms[MISMS_RANGE_6]+misms[MISMS_RANGE_7]+
      misms[MISMS_RANGE_8]+misms[MISMS_RANGE_9]+misms[MISMS_RANGE_10];
  fprintf(stderr,"  -->      (2,10]%% \t=> "GT_STATS_PRINT_MISMS_FORMAT,GT_STATS_PRINT_MISMS_VAL(misms_accum));
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
void gt_stats_print_read_event_positions(
    uint64_t* const pos_error,uint64_t const num_errors,uint64_t const max_length) {
  register const uint64_t max_length_stats = GT_MIN(max_length,LARGE_READ_POS_RANGE);
  register const uint64_t ten_per_cent = max_length_stats/10;
  register uint64_t i, centinel, error_sum;
  if(!num_errors) return;
  fprintf(stderr,"Error.position [0,%" PRIu64 "]nt\n",max_length_stats);
  for (i=0,centinel=0;i<100;i+=10,centinel+=ten_per_cent) {
    register const uint64_t top = (i==90)?max_length_stats:centinel+ten_per_cent;
    error_sum = gt_stats_sum_misms_pos(pos_error,centinel,top);
    fprintf(stderr,"  -->   [%3" PRIu64 " - %3" PRIu64 ")nt \t=> %1.3f%%\n",
        centinel,top,100.0*(float)error_sum/(float)num_errors);
  }
}
void gt_stats_print_num_junctions_distribution(uint64_t* const num_junctions,uint64_t const total) {
#define GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_NUM_JUNCTIONS(RANGE) num_junctions[RANGE],100.0*(float)num_junctions[RANGE]/(float)total
  fprintf(stderr,"  -->           [1]%% \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(NUM_JUNCTION_1));
  fprintf(stderr,"  -->           [2]%% \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(NUM_JUNCTION_2));
  fprintf(stderr,"  -->           [3]%% \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(NUM_JUNCTION_3));
  fprintf(stderr,"  -->       (3,inf)%% \t=> "GT_STATS_PRINT_NUM_JUNCTIONS_FORMAT,GT_STATS_PRINT_NUM_JUNCTIONS(NUM_JUNCTION_BEHOND));
}
void gt_stats_print_length_junctions_distribution(uint64_t* const length_junctions,uint64_t const total_junctions) {
#define GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT "%8" PRIu64 " \t %1.3f%%\n"
#define GT_STATS_PRINT_LENGTH_JUNCTION(RANGE) length_junctions[RANGE],100.0*(float)length_junctions[RANGE]/(float)total_junctions
  fprintf(stderr,"SM.Length.Junctions\n");
  fprintf(stderr,"  -->       [0,100]%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_100));
  fprintf(stderr,"  -->    (100,1000]%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_1000));
  fprintf(stderr,"  -->   (1000,5000]%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_5000));
  fprintf(stderr,"  -->  (5000,10000]%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_10000));
  fprintf(stderr,"  --> (10000,50000]%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_50000));
  fprintf(stderr,"  -->   (50000,inf)%% \t=> "GT_STATS_PRINT_LENGTH_JUNCTION_FORMAT,GT_STATS_PRINT_LENGTH_JUNCTION(LEN_JUNCTION_BEHOND));
}
void gt_stats_print_junction_position_distribution(uint64_t* const junction_position,uint64_t const total_junctions,uint64_t const max_length) {
  register const uint64_t max_length_stats = GT_MIN(max_length,SHORT_READ_POS_RANGE);
  register const uint64_t ten_per_cent = max_length_stats/10;
  register uint64_t i, centinel;
  if(!total_junctions) return;
  fprintf(stderr,"Juntion.position [0,%" PRIu64 "]nt\n",max_length_stats);
  for (i=0,centinel=0;i<100;i+=10,centinel+=ten_per_cent) {
    register const uint64_t top = (i==90)?max_length_stats:centinel+ten_per_cent;
    register const uint64_t num_junctions = gt_stats_sum_misms_pos(junction_position,centinel,top);
    fprintf(stderr,"  -->   [%3" PRIu64 " - %3" PRIu64")nt \t=> %1.3f%%\n",centinel,top,100.0*(float)num_junctions/(float)total_junctions);
  }
}

void gt_stats_print_qualities_error_distribution(uint64_t* const qualities_error,uint64_t const total_error,const char* const header) {
  register uint64_t i;
  if(!total_error) return;
  fprintf(stderr,"%s\n",header);
//  // Print Header
//  fprintf(stderr,"    ");
//  for (j=0;j<16;++j) fprintf(stderr,"[%4lu]",j);
//  fprintf(stderr,"\n");
//  // Print Values
//  for (i=0;i<16;++i) {
//    fprintf(stderr,"%3lu ",qual);
//    for (j=0;j<16;++j,++qual) {
//      fprintf(stderr,"[%4.1f]",100.0*((double)qualities_error[qual]/(double)total_error));
//    }
//    fprintf(stderr,"\n");
//  }
  // CSV, Compact format
  for (i=0;i<256;++i) {
    fprintf(stderr,"\"%04.01f\",",100.0*((double)qualities_error[i]/(double)total_error));
  }
  fprintf(stderr,"\n");
}
void gt_stats_print_misms_transition_table(uint64_t* const misms_trans,uint64_t const total_misms,const char* const header) {
  register uint64_t i, pos=0;
  if(!total_misms) return;
  fprintf(stderr,"%s\n",header);
  // Print Header
  fprintf(stderr,"    [  A  ][  C  ][  G  ][  T  ][  N  ]\n");
  fprintf(stderr,"[A] ");
  for (i=0;i<MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stderr,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"[C] ");
  for (i=0;i<MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stderr,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"[G] ");
  for (i=0;i<MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stderr,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"[T] ");
  for (i=0;i<MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stderr,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"[N] ");
  for (i=0;i<MISMS_BASE_RANGE;++i,++pos) {
    fprintf(stderr,"[%5.2f]",100.0*(double)misms_trans[pos]/(double)total_misms);
  }
  fprintf(stderr,"\n");
}
void gt_stats_print_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  fprintf(stderr,"[GENERAL.STATS]\n");
  /*
   * Reads
   */
  register const uint64_t eff_num_reads = paired_end?num_reads>>1:num_reads;
  if(stats->min_length>stats->max_length) stats->min_length=0; // For the case of zero input lines
  fprintf(stderr,"Num.reads %" PRIu64 "\n",num_reads);
  fprintf(stderr,"  --> Length.(min,avg,max) (%" PRIu64 ",%" PRIu64 ",%" PRIu64 ")\n",
	  stats->min_length,stats->num_blocks?stats->total_bases/stats->num_blocks:(uint64_t)0,stats->max_length);
  fprintf(stderr,"  --> Num.mapped %" PRIu64 " (%2.3f%%)\n",stats->num_mapped,
	  eff_num_reads?100.0*(float)stats->num_mapped/(float)eff_num_reads:0.0);
  fprintf(stderr,"Num.alignments %" PRIu64 "\n",stats->num_alignments);
  /*
   * Total bases (aligned/trimmed/unaligned)
   */
  register const uint64_t total_bases_mmaps =
      stats->maps_error_profile->total_bases_aligned+stats->maps_error_profile->total_bases_trimmed;
  fprintf(stderr,"  --> Num.bases %" PRIu64 "\n",stats->total_bases);
  fprintf(stderr,"    --> Num.bases.unaligned %" PRIu64 " (%2.3f%%)\n",
	  stats->total_bases_unaligned,stats->total_bases?100.0*(float)stats->total_bases_unaligned/(float)stats->total_bases:0.0);
  fprintf(stderr,"    --> Num.bases.allMaps.aligned %" PRIu64 " (%2.3f%%)\n",
      stats->maps_error_profile->total_bases_aligned,
	  total_bases_mmaps?100.0*(float)stats->maps_error_profile->total_bases_aligned/(float)total_bases_mmaps:0.0);
  fprintf(stderr,"    --> Num.bases.allMaps.trimmed %" PRIu64 " (%2.3f%%)\n",
      stats->maps_error_profile->total_bases_trimmed,
	  total_bases_mmaps?100.0*(float)stats->maps_error_profile->total_bases_trimmed/(float)total_bases_mmaps:0.0);
  /*
   * Maps
   */
  fprintf(stderr,"  --> Num.maps %" PRIu64 " (%2.3f map/alg)\n",stats->num_maps,stats->num_mapped?(float)stats->num_maps/(float)stats->num_mapped:0.0);
  gt_stats_print_mmap_distribution(stats->mmap,eff_num_reads,stats->num_mapped);
  if (paired_end) {
    gt_stats_print_inss_distribution(stats->maps_error_profile->inss,stats->num_maps);
  }
  gt_stats_print_uniq_distribution(stats->uniq,eff_num_reads);
  if (parameters.error_profile) {
    fprintf(stderr,"[ERROR.PROFILE]\n");
    gt_stats_print_misms_distribution(stats->maps_error_profile->mismatches,stats->num_maps,"Mismatches");
    gt_stats_print_misms_distribution(stats->maps_error_profile->indel_length,stats->num_maps,"Indels.length");
    gt_stats_print_misms_distribution(stats->maps_error_profile->errors_events,stats->num_maps,"Error.events");
    gt_stats_print_read_event_positions(
        stats->maps_error_profile->error_position,stats->maps_error_profile->total_errors_events,stats->max_length);
  }
  /*
   * Print Quality Scores vs Errors/Misms
   */
  if (parameters.mismatch_quality) {
    fprintf(stderr,"[MISMATCH.QUALITY]\n");
    gt_stats_print_qualities_error_distribution(
        stats->maps_error_profile->qual_score_misms,stats->maps_error_profile->total_mismatches,"Qualities.Misms");
    fprintf(stderr,"[ERRORS.QUALITY]\n");
    gt_stats_print_qualities_error_distribution(
        stats->maps_error_profile->qual_score_errors,stats->maps_error_profile->total_errors_events,"Qualities.Error");
  }
  /*
   * Print Mismatch transition table
   */
  if (parameters.mismatch_transitions) {
    fprintf(stderr,"[MISMATCH.TRANSITIONS]\n");
    gt_stats_print_misms_transition_table(
        stats->maps_error_profile->misms_transition,stats->maps_error_profile->total_mismatches,"MismsTransitions");
  }
  /*
   * Print Splitmaps profile
   */
  if (parameters.splitmaps_profile) {
    gt_splitmaps_profile* const splitmap_stats = stats->splitmaps_profile;
    fprintf(stderr,"[SPLITMAPS.PROFILE]\n");
    if (splitmap_stats->total_splitmaps==0) {
      fprintf(stderr,"SM.Total \t 0\n");
    } else {
      fprintf(stderr,"SM.Num.mapped.withSM \t %" PRIu64 " (%2.3f%%)\n",splitmap_stats->num_mapped_with_splitmaps,
	      stats->num_mapped?100.0*(float)splitmap_stats->num_mapped_with_splitmaps/(float)stats->num_mapped:0.0);
      fprintf(stderr,"SM.Num.mapped.onlyBySM \t %" PRIu64 " (%2.3f%%)\n",splitmap_stats->num_mapped_only_splitmaps,
	      stats->num_mapped?100.0*(float)splitmap_stats->num_mapped_only_splitmaps/(float)stats->num_mapped:0.0);

      register const uint64_t num_single_maps = stats->num_maps*(paired_end?2:1);
      fprintf(stderr,"SM.Num.Splitted.Segments \t %" PRIu64 " (%2.3f%% SplitMaps/Maps)\n",splitmap_stats->total_splitmaps,
	      num_single_maps?100.0*(float)splitmap_stats->total_splitmaps/(float)num_single_maps:0.0);
      fprintf(stderr,"SM.Num.Junctions \t %" PRIu64 " (%2.3f Juntions/SplitMaps)\n",splitmap_stats->total_junctions,
	      splitmap_stats->total_splitmaps?(float)splitmap_stats->total_junctions/(float)splitmap_stats->total_splitmaps:0.0);

      gt_stats_print_num_junctions_distribution(splitmap_stats->num_junctions,splitmap_stats->total_splitmaps);
      gt_stats_print_length_junctions_distribution(splitmap_stats->length_junctions,splitmap_stats->total_junctions);
      gt_stats_print_junction_position_distribution(splitmap_stats->junction_position,splitmap_stats->total_junctions,stats->max_length);
      if (parameters.paired_end) {
        fprintf(stderr,"SM.PE.Combinations\n");
        fprintf(stderr,"  --> SM+SM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_sm_sm,100.0*(float)splitmap_stats->pe_sm_sm/(float)stats->num_maps);
        fprintf(stderr,"  --> SM+RM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_sm_rm,100.0*(float)splitmap_stats->pe_sm_rm/(float)stats->num_maps);
        fprintf(stderr,"  --> RM+RM %" PRIu64 " (%2.3f)\n",splitmap_stats->pe_rm_rm,100.0*(float)splitmap_stats->pe_rm_rm/(float)stats->num_maps);
      }
    }
  }
}
void gt_stats_print_stats_compact(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  /*
   *  #mapped, %mapped, #unmapped, %unmapped, MMap(maps/alg), Bases.aligned(%), Bases.trimmed(%)
   *  mapped(%), MMap(maps/alg), Bases.aligned(%), Bases.trimmed(%), #Uniq-0, %Uniq-0
   */
  // #mapped, %mapped
  register const uint64_t eff_num_reads = paired_end?num_reads>>1:num_reads;
  fprintf(stderr,"%" PRIu64 ",",stats->num_mapped);
  fprintf(stderr,"%2.3f,",eff_num_reads?100.0*(float)stats->num_mapped/(float)eff_num_reads:0.0);
  // #unmapped, %unmapped
  register const uint64_t unmapped = eff_num_reads-stats->num_mapped;
  fprintf(stderr,"%" PRIu64 ",",unmapped);
  fprintf(stderr,"%2.3f,",eff_num_reads?100.0*(float)unmapped/(float)eff_num_reads:0.0);
  // MMap(maps/alg)
  fprintf(stderr,"%2.3f,",stats->num_mapped?(float)stats->num_maps/(float)stats->num_mapped:0.0);
  // Bases.aligned(%)
  register const uint64_t total_bases_mmaps =
      stats->maps_error_profile->total_bases_aligned+stats->maps_error_profile->total_bases_trimmed;
  fprintf(stderr,"%2.3f,",total_bases_mmaps?100.0*(float)stats->maps_error_profile->total_bases_aligned/(float)total_bases_mmaps:0.0);
  // Bases.trimmed(%)
  fprintf(stderr,"%2.3f,",total_bases_mmaps?100.0*(float)stats->maps_error_profile->total_bases_trimmed/(float)total_bases_mmaps:0.0);
  // #Uniq-0, %Uniq-0
  fprintf(stderr,"%" PRIu64 ",",stats->uniq[UNIQ_RANGE_0]);
  fprintf(stderr,"%2.3f\n",eff_num_reads?100.0*(float)stats->uniq[UNIQ_RANGE_0]/(float)eff_num_reads:0.0);
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
    stats[tid] = gt_stats_new();
    gt_generic_parser_attr generic_parser_attr = GENERIC_PARSER_ATTR_DEFAULT(parameters.paired_end);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,&generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // Extract stats
      gt_stats_calculate_template_stats(stats[tid],template,parameters.best_map);
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
  gt_stats_delete(stats[0]); free(stats);
  gt_input_file_close(input_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.stats [ARGS]...\n"
                  "       [Input]\n"
                  "        --input|-i [FILE]\n"
                  "        --mmap-input\n"
                  "        --paired-end|p\n"
                  "        --num-reads|n\n"
                  "        --best-map\n"
                  "       [Tests]\n"
                  "        --error-profile|E\n"
                  "        --mismatch-transitions|T\n"
                  "        --mismatch-quality|Q\n"
                  "        --splitmaps-profile|S\n"
                  "       [Output]\n"
                  "        --compact|c\n"
                  "        --verbose|v\n"
                  "        --quiet|q\n"
                  "       [Misc]\n"
                  "        --threads|t\n"
                  "        --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* [Input] */
    { "input", required_argument, 0, 'i' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    { "num-reads", no_argument, 0, 'n' },
    { "best-map", no_argument, 0, 2 },
    /* [Tests] */
    { "error-profile", no_argument, 0, 'E' },
    { "mismatch-transitions", no_argument, 0, 'T' },
    { "mismatch-quality", no_argument, 0, 'Q' },
    { "splitmaps-profile", no_argument, 0, 'S' },
    /* [Output] */
    { "compact", no_argument, 0, 'c' },
    { "verbose", no_argument, 0, 'v' },
    { "quiet", no_argument, 0, 'q' },
    /* [Misc] */
    { "threads", no_argument, 0, 't' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:pn:ETQScvqt:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /* [Input] */
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
    case 2:
      parameters.best_map = true;
      break;
    /* [Tests] */
    case 'E': // --error-profile
      parameters.error_profile = true;
      break;
    case 'T': // --mismatch-transitions
      parameters.mismatch_transitions = true;
      break;
    case 'Q': // --mismatch-quality
      parameters.mismatch_quality = true;
      break;
    case 'S': // --splitmaps-profile
      parameters.splitmaps_profile = true;
      break;
    /* [Output] */
    case 'c':
      parameters.compact = true;
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'q':
      parameters.quiet = true;
      break;
    /* [Misc] */
    case 't':
      parameters.num_threads = atol(optarg);
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


