/*
 * PROJECT: GEM-Tools library
 * FILE: gt.stats.c
 * DATE: 02/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to retrieve very naive stats from {MAP,SAM,FASTQ} files
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

#define GT_STATS_OUT_FILE stdout

typedef struct {
  /* [Input] */
  char *name_input_file;
  char *name_reference_file;
  char *name_output_file;
  FILE* output_file;
  FILE* output_file_json;
  bool mmap_input;
  bool paired_end;
  uint64_t num_reads;
  /* [Tests] */
  bool first_map;
  bool maps_profile;
  bool mismatch_transitions;
  bool mismatch_quality;
  bool splitmaps_profile;
  bool indel_profile;
  bool population_profile;
  /* [MAP Specific] */
  bool use_only_decoded_maps;
  /* [Output] */
  bool verbose;
  bool compact; // FIXME Deleteme
  bool print_json;
  bool print_both;
  /* [Misc] */
  uint64_t num_threads;
} gt_stats_args;

gt_stats_args parameters = {
    /* [Input] */
    .name_input_file=NULL,
    .name_reference_file=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .num_reads=0,
    .output_file=NULL,
    .output_file_json=NULL,
    /* [Tests] */
    .first_map=false,
    .maps_profile = false,
    .mismatch_transitions = false,
    .mismatch_quality = false,
    .splitmaps_profile = false,
    .indel_profile = false,
    .population_profile = false,
    /* [MAP Specific] */
    .use_only_decoded_maps = false,
    /* [Output] */
    .verbose=false,
    .compact = false,
    /* [Misc] */
    .num_threads=1,
    .print_json=false,
    .print_both=false,
};
/*
 * STATS Print results
 */
JsonNode* gt_stats_print_json_general_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end){
  /**
   * General stats
   */
  const uint64_t num_templates = paired_end ? num_reads/2 : num_reads;
  // For the case of zero input lines
  if(stats->min_length>stats->max_length) stats->min_length=0;
  if(stats->mapped_min_length>stats->mapped_max_length) stats->mapped_min_length=0;

  JsonNode* general_stats = json_mkobject();
  json_append_member(general_stats, "num_reads", json_mknumber(num_reads));
  json_append_member(general_stats, "num_templates", json_mknumber(num_templates));
  json_append_member(general_stats, "read_lenght_min", json_mknumber(stats->min_length));
  json_append_member(general_stats, "read_lenght_avg", json_mknumber(GT_DIV(stats->total_bases,stats->num_blocks)));
  json_append_member(general_stats, "read_lenght_max", json_mknumber(stats->max_length));
  json_append_member(general_stats, "templates_mapped", json_mknumber(stats->num_mapped));
  json_append_member(general_stats, "reads_mapped", json_mknumber(stats->num_mapped_reads));
  json_append_member(general_stats, "reads_mapped_length_min", json_mknumber(stats->mapped_min_length));
  json_append_member(general_stats, "reads_mapped_length_avg", json_mknumber(GT_DIV(stats->total_bases_aligned,(paired_end)?stats->num_mapped*2:stats->num_mapped)));
  json_append_member(general_stats, "reads_mapped_length_max", json_mknumber(stats->mapped_max_length));
  json_append_member(general_stats, "num_bases", json_mknumber(stats->total_bases));
  json_append_member(general_stats, "num_bases_aligned", json_mknumber(stats->total_bases_aligned));
  json_append_member(general_stats, "bases_prop", gt_json_int_named_tuple(5,
      "A", stats->nt_counting[0],
      "C", stats->nt_counting[1],
      "G", stats->nt_counting[2],
      "T", stats->nt_counting[3],
      "N", stats->nt_counting[4]
  ));
  json_append_member(general_stats, "num_alignments", json_mknumber(stats->num_alignments));
  json_append_member(general_stats, "read_length_ranges", gt_json_int_named_tuple(
      GT_STATS_LENGTH_RANGE,
      "[0,5]", stats->length[GT_STATS_LENGTH_RANGE_5],
      "(5,40]", stats->length[GT_STATS_LENGTH_RANGE_40],
      "(40,80]", stats->length[GT_STATS_LENGTH_RANGE_80],
      "(80,100]", stats->length[GT_STATS_LENGTH_RANGE_100],
      "(100,150]", stats->length[GT_STATS_LENGTH_RANGE_150],
      "(150,300]", stats->length[GT_STATS_LENGTH_RANGE_300],
      "(300,800]", stats->length[GT_STATS_LENGTH_RANGE_800],
      "(800,1000]", stats->length[GT_STATS_LENGTH_RANGE_1000],
      "(1000,2000]", stats->length[GT_STATS_LENGTH_RANGE_2000],
      "(2000,5000]", stats->length[GT_STATS_LENGTH_RANGE_5000],
      "(5000,inf)", stats->length[GT_STATS_LENGTH_RANGE_BEHOND]
  ));
  json_append_member(general_stats, "read_length_ranges_mapped", gt_json_int_named_tuple(
      GT_STATS_LENGTH_RANGE,
      "[0,5]", stats->length_mapped[GT_STATS_LENGTH_RANGE_5],
      "(5,40]", stats->length_mapped[GT_STATS_LENGTH_RANGE_40],
      "(40,80]", stats->length_mapped[GT_STATS_LENGTH_RANGE_80],
      "(80,100]", stats->length_mapped[GT_STATS_LENGTH_RANGE_100],
      "(100,150]", stats->length_mapped[GT_STATS_LENGTH_RANGE_150],
      "(150,300]", stats->length_mapped[GT_STATS_LENGTH_RANGE_300],
      "(300,800]", stats->length_mapped[GT_STATS_LENGTH_RANGE_800],
      "(800,1000]", stats->length_mapped[GT_STATS_LENGTH_RANGE_1000],
      "(1000,2000]", stats->length_mapped[GT_STATS_LENGTH_RANGE_2000],
      "(2000,5000]", stats->length_mapped[GT_STATS_LENGTH_RANGE_5000],
      "(5000,inf)", stats->length_mapped[GT_STATS_LENGTH_RANGE_BEHOND]
  ));

  json_append_member(general_stats, "read_qualities_avg", gt_json_int_array(32, 192, stats->avg_quality));

  JsonNode* read_length_quals = json_mkobject();
  json_append_member(read_length_quals, "[0,5]", gt_json_int_array(GT_STATS_LENGTH_RANGE_5*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(5,40]", gt_json_int_array(GT_STATS_LENGTH_RANGE_40*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(40,80]", gt_json_int_array(GT_STATS_LENGTH_RANGE_80*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(80,100]", gt_json_int_array(GT_STATS_LENGTH_RANGE_100*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(100,150]", gt_json_int_array(GT_STATS_LENGTH_RANGE_150*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(300,800]", gt_json_int_array(GT_STATS_LENGTH_RANGE_300*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(800,1000]", gt_json_int_array(GT_STATS_LENGTH_RANGE_800*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(1000,2000]", gt_json_int_array(GT_STATS_LENGTH_RANGE_1000*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(2000,5000]", gt_json_int_array(GT_STATS_LENGTH_RANGE_2000*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(read_length_quals, "(5000,inf)", gt_json_int_array(GT_STATS_LENGTH_RANGE_BEHOND*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->length__quality));
  json_append_member(general_stats, "read_qualities_per_length", read_length_quals);
  return general_stats;
}

JsonNode* gt_stats_create_error_distribution(uint64_t* data){
  return gt_json_int_named_tuple(
      GT_STATS_MISMS_RANGE,
      "[0]", data[GT_STATS_MISMS_RANGE_0],
      "[1]", data[GT_STATS_MISMS_RANGE_1],
      "[2]", data[GT_STATS_MISMS_RANGE_2],
      "[3]", data[GT_STATS_MISMS_RANGE_3],
      "[4]", data[GT_STATS_MISMS_RANGE_4],
      "[5]", data[GT_STATS_MISMS_RANGE_5],
      "[6]", data[GT_STATS_MISMS_RANGE_6],
      "[7]", data[GT_STATS_MISMS_RANGE_7],
      "[8]", data[GT_STATS_MISMS_RANGE_8],
      "[9]", data[GT_STATS_MISMS_RANGE_9],
      "[10]", data[GT_STATS_MISMS_RANGE_10],
      "(10,20]", data[GT_STATS_MISMS_RANGE_20],
      "(20,50]", data[GT_STATS_MISMS_RANGE_50],
      "(50,100]", data[GT_STATS_MISMS_RANGE_BEHOND]
  );
}
#define GT_STATS_GET_IXD_TRANSITION_1_CTX(a,b,c,i) ((((a*GT_STATS_MISMS_BASE_RANGE+b)*GT_STATS_MISMS_BASE_RANGE)+c)*GT_STATS_MISMS_BASE_RANGE+i)
JsonNode* gt_stats_print_json_maps_profile(gt_stats* const stats,const uint64_t num_reads,const bool paired_end) {
  const gt_maps_profile* const maps_profile = stats->maps_profile;
  JsonNode* profile = json_mkobject();
  json_append_member(profile, "num_bases", json_mknumber(maps_profile->total_bases));
  json_append_member(profile, "num_matching_bases", json_mknumber(maps_profile->total_bases_matching));
  json_append_member(profile, "num_trimmed_bases", json_mknumber(maps_profile->total_bases_trimmed));
  json_append_member(profile, "multi_map_ranges", gt_json_int_named_tuple(
      GT_STATS_MMAP_RANGE,
      "[0]", stats->mmap[GT_STATS_MMAP_RANGE_0],
      "[1]", stats->mmap[GT_STATS_MMAP_RANGE_1],
      "(1,5]", stats->mmap[GT_STATS_MMAP_RANGE_5],
      "(5,10]", stats->mmap[GT_STATS_MMAP_RANGE_10],
      "(10,50]", stats->mmap[GT_STATS_MMAP_RANGE_50],
      "(50,100]", stats->mmap[GT_STATS_MMAP_RANGE_100],
      "(100,500]", stats->mmap[GT_STATS_MMAP_RANGE_500],
      "(500,1000]", stats->mmap[GT_STATS_MMAP_RANGE_1000],
      "(1000, inf)", stats->mmap[GT_STATS_MMAP_RANGE_BEHOND]
  ));
  JsonNode* mmap_read_length = json_mkobject();
  json_append_member(mmap_read_length, "[0]", gt_json_int_array(GT_STATS_MMAP_RANGE_0*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "[1]", gt_json_int_array(GT_STATS_MMAP_RANGE_1*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(1,5]", gt_json_int_array(GT_STATS_MMAP_RANGE_5*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(5,10]", gt_json_int_array(GT_STATS_MMAP_RANGE_10*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(10,50]", gt_json_int_array(GT_STATS_MMAP_RANGE_50*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(50,100]", gt_json_int_array(GT_STATS_MMAP_RANGE_100*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(100,500]", gt_json_int_array(GT_STATS_MMAP_RANGE_500*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(500,1000]", gt_json_int_array(GT_STATS_MMAP_RANGE_1000*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(mmap_read_length, "(1000,inf)", gt_json_int_array(GT_STATS_MMAP_RANGE_BEHOND*GT_STATS_MMAP_RANGE, GT_STATS_MMAP_RANGE, stats->length__mmap));
  json_append_member(profile, "read_length_multi_map_ranges", mmap_read_length);

  JsonNode* mmap_read_quality = json_mkobject();
  json_append_member(mmap_read_quality, "[0]", gt_json_int_array(GT_STATS_MMAP_RANGE_0*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "[1]", gt_json_int_array(GT_STATS_MMAP_RANGE_1*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(1,5]", gt_json_int_array(GT_STATS_MMAP_RANGE_5*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(5,10]", gt_json_int_array(GT_STATS_MMAP_RANGE_10*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(10,50]", gt_json_int_array(GT_STATS_MMAP_RANGE_50*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(50,100]", gt_json_int_array(GT_STATS_MMAP_RANGE_100*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(100,500]", gt_json_int_array(GT_STATS_MMAP_RANGE_500*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(500,1000]", gt_json_int_array(GT_STATS_MMAP_RANGE_1000*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(mmap_read_quality, "(1000,inf)", gt_json_int_array(GT_STATS_MMAP_RANGE_BEHOND*GT_STATS_QUAL_SCORE_RANGE, GT_STATS_QUAL_SCORE_RANGE, stats->mmap__avg_quality));
  json_append_member(profile, "read_quality_multi_map_ranges", mmap_read_quality);

  json_append_member(profile, "unique_ranges", gt_json_int_named_tuple(
      GT_STATS_UNIQ_RANGE,
      "[X]", stats->uniq[GT_STATS_UNIQ_RANGE_X],
      "[0]", stats->uniq[GT_STATS_UNIQ_RANGE_0],
      "[1]", stats->uniq[GT_STATS_UNIQ_RANGE_1],
      "[2]", stats->uniq[GT_STATS_UNIQ_RANGE_2],
      "[3]", stats->uniq[GT_STATS_UNIQ_RANGE_3],
      "(3,10]", stats->uniq[GT_STATS_UNIQ_RANGE_10],
      "(10,50]", stats->uniq[GT_STATS_UNIQ_RANGE_50],
      "(50,100]", stats->uniq[GT_STATS_UNIQ_RANGE_100],
      "(100, 500]", stats->uniq[GT_STATS_UNIQ_RANGE_500],
      "(500, inf)", stats->uniq[GT_STATS_UNIQ_RANGE_BEHOND]
  ));
  if(paired_end){
    json_append_member(profile, "strands", gt_json_int_named_tuple(
        4,
        "F+R", maps_profile->pair_strand_fr,
        "R+F", maps_profile->pair_strand_rf,
        "F+F", maps_profile->pair_strand_ff,
        "R+R", maps_profile->pair_strand_rr
    ));
  }else{
    json_append_member(profile, "strands", gt_json_int_named_tuple(
        2,
        "F", maps_profile->single_strand_f,
        "R", maps_profile->single_strand_r
    ));
  }
  JsonNode* error_profile = json_mkobject();
  json_append_member(error_profile, "total_mismatches", json_mknumber(maps_profile->total_mismatches));
  json_append_member(error_profile, "total_errors", json_mknumber(maps_profile->total_errors_events));
  json_append_member(error_profile, "total_indel_length", json_mknumber(maps_profile->total_indel_length));
  json_append_member(error_profile, "total_levenshtein", json_mknumber(maps_profile->total_levenshtein));
  json_append_member(error_profile, "mismatch_distribution", gt_stats_create_error_distribution(maps_profile->mismatches));
  json_append_member(error_profile, "levenshtein_distribution", gt_stats_create_error_distribution(maps_profile->levenshtein));
  json_append_member(error_profile, "error_distribution", gt_stats_create_error_distribution(maps_profile->errors_events));
  json_append_member(error_profile, "inserts_distribution", gt_stats_create_error_distribution(maps_profile->insertion_length));
  json_append_member(error_profile, "deletions_distribution", gt_stats_create_error_distribution(maps_profile->deletion_length));
  json_append_member(error_profile, "error_positions", gt_json_int_array(
      0,
      GT_MIN(stats->max_length,GT_STATS_LARGE_READ_POS_RANGE),
      maps_profile->error_position
  ));
  json_append_member(profile, "error_profile", error_profile);

  JsonNode* quality_profile = json_mkobject();
  json_append_member(quality_profile, "mismatch_qualities_avg", gt_json_int_array(32, 192, maps_profile->qual_score_misms));
  json_append_member(quality_profile, "error_qualities_avg", gt_json_int_array(32, 192, maps_profile->qual_score_errors));
  json_append_member(profile, "quality_profile", quality_profile);

  JsonNode* transition_profile = json_mkobject();
  json_append_member(transition_profile, "transitions", gt_json_int_array(0,GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE, maps_profile->misms_transition));
  JsonNode* transition_context_order = json_mkarray();
  JsonNode* transition_context = json_mkarray();
  // transition 1context
  char bases[] = {'A','C','G','T','N'};
  char* name = gt_malloc_(4, sizeof(char), false, false);
  uint64_t a,b,c;
  for (b=0;b<4;++b) {
    for (a=0;a<4;++a) {
      for (c=0;c<4;++c) {
        sprintf(name, "%c%c%c", bases[a],bases[b],bases[c]);
        json_append_element(transition_context_order, json_mkstring(name));
        uint64_t i;
        for (i=0;i<GT_STATS_MISMS_BASE_RANGE;++i) {
          json_append_element(transition_context, json_mknumber(maps_profile->misms_1context[GT_STATS_GET_IXD_TRANSITION_1_CTX(a,b,c,i)]));
        }
      }
    }
  }
  gt_free(name);

  json_append_member(transition_profile, "transition_context", transition_context);
  json_append_member(transition_profile, "transition_context_order", transition_context_order);
  json_append_member(profile, "transition_profile", transition_profile);
  return profile;
}
JsonNode* gt_stats_print_json_splits_profile(gt_stats* const stats,const uint64_t num_reads,const bool paired_end) {
  gt_splitmaps_profile* const splitmap_stats = stats->splitmaps_profile;
  JsonNode* profile = json_mkobject();
  json_append_member(profile, "total_reads_with_splitmap", json_mknumber(splitmap_stats->total_splitmaps));
  json_append_member(profile, "total_junctions", json_mknumber(splitmap_stats->total_junctions));
  json_append_member(profile, "mapped_with_sm", json_mknumber(splitmap_stats->num_mapped_with_splitmaps));
  json_append_member(profile, "mapped_only_with_sm", json_mknumber(splitmap_stats->num_mapped_only_splitmaps));
  json_append_member(profile, "num_junctions", gt_json_int_named_tuple(
      GT_STATS_NUM_JUNCTION_RANGE,
      "[1]", splitmap_stats->num_junctions[GT_STATS_NUM_JUNCTION_1],
      "[2]", splitmap_stats->num_junctions[GT_STATS_NUM_JUNCTION_2],
      "[3]", splitmap_stats->num_junctions[GT_STATS_NUM_JUNCTION_3],
      "(3,inf]", splitmap_stats->num_junctions[GT_STATS_NUM_JUNCTION_BEHOND]
  ));
  json_append_member(profile, "junction_length", gt_json_int_named_tuple(
      GT_STATS_LEN_JUNCTION_RANGE,
      "[0,100]", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_100],
      "(100,1000]", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_1000],
      "(1000,5000]", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_5000],
      "(5000,10000]", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_10000],
      "(10000,50000]", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_50000],
      "(50000,inf)", splitmap_stats->length_junctions[GT_STATS_LEN_JUNCTION_BEHOND]
  ));
  json_append_member(profile, "junctions_positions", gt_json_int_array(
       0,
       GT_MIN(stats->max_length,GT_STATS_SHORT_READ_POS_RANGE),
       splitmap_stats->junction_position
   ));

  if (paired_end) {
    json_append_member(profile, "splits_sm_sm", json_mknumber(splitmap_stats->pe_sm_sm));
    json_append_member(profile, "splits_sm_rm", json_mknumber(splitmap_stats->pe_sm_rm));
    json_append_member(profile, "splits_rm_rm", json_mknumber(splitmap_stats->pe_rm_rm));
  }
  return profile;
}

void gt_stats_print_json_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  JsonNode* root = json_mkobject();
  json_append_member(root, "general", gt_stats_print_json_general_stats(stats, num_reads, paired_end));
  json_append_member(root, "maps_profile", gt_stats_print_json_maps_profile(stats, num_reads, paired_end));
  json_append_member(root, "splits_profile", gt_stats_print_json_splits_profile(stats, num_reads, paired_end));

  fprintf(parameters.output_file_json, "%s\n", json_stringify(root, "  "));
  json_delete(root);
}

void gt_stats_print_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  /*
   * General.Stats (Reads,Alignments,...)
   */
  fprintf(parameters.output_file,"[GENERAL.STATS]\n");
  gt_stats_print_general_stats(parameters.output_file,stats,num_reads,paired_end);
  /*
   * Maps
   */
  if (parameters.maps_profile) {
    fprintf(parameters.output_file,"[MAPS.PROFILE]\n");
    gt_stats_print_maps_stats(parameters.output_file,stats,num_reads,paired_end);
  }
  /*
   * Print Quality Scores vs Errors/Misms
   */
  if (parameters.mismatch_quality && num_reads>0) {
    const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(parameters.output_file,"[MISMATCH.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          parameters.output_file,maps_profile->qual_score_misms,maps_profile->total_mismatches);
    }
    if (maps_profile->total_errors_events > 0) {
      fprintf(parameters.output_file,"[ERRORS.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          parameters.output_file,maps_profile->qual_score_errors,maps_profile->total_errors_events);
    }
  }
  /*
   * Print Mismatch transition table
   */
  if (parameters.mismatch_transitions && num_reads>0) {
    const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(parameters.output_file,"[MISMATCH.TRANSITIONS]\n");
      fprintf(parameters.output_file,"MismsTransitions\n");
      gt_stats_print_misms_transition_table(
          parameters.output_file,maps_profile->misms_transition,maps_profile->total_mismatches);
      fprintf(parameters.output_file,"MismsTransitions.1-Nucleotide.Context\n");
      gt_stats_print_misms_transition_table_1context(
          parameters.output_file,maps_profile->misms_1context,maps_profile->total_mismatches);
    }
  }
  /*
   * Print Splitmaps profile
   */
  if (parameters.splitmaps_profile) {
    fprintf(parameters.output_file,"[SPLITMAPS.PROFILE]\n");
    gt_stats_print_split_maps_stats(parameters.output_file,stats,parameters.paired_end);
  }
  /*
   * Print Population profile
   */
  if (parameters.population_profile) {
    fprintf(parameters.output_file,"[POPULATION.PROFILE]\n");
    gt_stats_print_population_stats(parameters.output_file,stats,num_reads,parameters.paired_end);
  }
}
void gt_stats_print_stats_compact(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  // #mapped, %mapped
  const uint64_t num_templates = paired_end ? num_reads>>1 : num_reads; // SE => 1 template. PE => 1 template
  fprintf(parameters.output_file,"%" PRIu64 ",",stats->num_mapped);
  fprintf(parameters.output_file,"%2.3f,",num_templates?100.0*(float)stats->num_mapped/(float)num_templates:0.0);
  // #unmapped, %unmapped
  const uint64_t unmapped = num_templates-stats->num_mapped;
  fprintf(parameters.output_file,"%" PRIu64 ",",unmapped);
  fprintf(parameters.output_file,"%2.3f,",num_templates?100.0*(float)unmapped/(float)num_templates:0.0);
  // MMap(maps/alg)
  fprintf(parameters.output_file,"%2.3f,",stats->num_mapped?(float)stats->num_maps/(float)stats->num_mapped:0.0);
  // Bases.aligned(%)
  fprintf(parameters.output_file,"%2.3f,",GT_GET_PERCENTAGE(stats->maps_profile->total_bases_matching,stats->maps_profile->total_bases));
  // Bases.trimmed(%)
  fprintf(parameters.output_file,"%2.3f,",GT_GET_PERCENTAGE(stats->maps_profile->total_bases_trimmed,stats->maps_profile->total_bases));
  // #Uniq-0, %Uniq-0
  const uint64_t all_uniq = stats->uniq[GT_STATS_UNIQ_RANGE_0]+
      stats->uniq[GT_STATS_UNIQ_RANGE_1]+stats->uniq[GT_STATS_UNIQ_RANGE_2]+
      stats->uniq[GT_STATS_UNIQ_RANGE_3]+stats->uniq[GT_STATS_UNIQ_RANGE_10]+
      stats->uniq[GT_STATS_UNIQ_RANGE_50]+stats->uniq[GT_STATS_UNIQ_RANGE_100]+
      stats->uniq[GT_STATS_UNIQ_RANGE_500]+stats->uniq[GT_STATS_UNIQ_RANGE_BEHOND];
  fprintf(parameters.output_file,"%" PRIu64 ",",all_uniq);
  fprintf(parameters.output_file,"%2.3f\n",num_templates?100.0*(float)all_uniq/(float)num_templates:0.0);
}

/*
 * CORE functions
 */
void gt_stats_parallel_generate_stats() {
  // Stats info
  gt_stats_analysis stats_analysis = GT_STATS_ANALYSIS_DEFAULT();
  gt_stats** stats = gt_calloc(parameters.num_threads,gt_stats*,false);

  // Select analysis
  stats_analysis.first_map = parameters.first_map;
  stats_analysis.maps_profile = parameters.maps_profile|parameters.mismatch_quality|parameters.mismatch_transitions;
  stats_analysis.nucleotide_stats = true;
  stats_analysis.splitmap_profile = parameters.splitmaps_profile;
  stats_analysis.indel_profile = parameters.indel_profile;
  stats_analysis.population_profile = parameters.population_profile;
  stats_analysis.use_map_counters = !parameters.use_only_decoded_maps;

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);

  gt_sequence_archive* sequence_archive = NULL;
  if (stats_analysis.indel_profile) {
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      fprintf(stderr,"\n");
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
#ifdef HAVE_OPENMP
    uint64_t tid = omp_get_thread_num();
#else
    uint64_t tid = 0;
#endif

    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

    gt_status error_code;
    gt_template *template = gt_template_new();
    stats[tid] = gt_stats_new();
    gt_generic_parser_attributes* generic_parser_attribute = gt_input_generic_parser_attributes_new(parameters.paired_end);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attribute))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // Extract stats
      gt_stats_calculate_template_stats(stats[tid],template,sequence_archive,&stats_analysis);
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }

  // Merge stats
  gt_stats_merge(stats,parameters.num_threads);

  /*
   * Print Statistics
   *   Use stats[0]->num_blocks as the number of blocks in a MAP/SAM/FASTA/FASTQ file
   *   is the number of reads in a FASTA/FASTQ
   */
  if(parameters.print_json){
    gt_stats_print_json_stats(stats[0],(parameters.num_reads>0)?
        parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
  }
  if(!parameters.print_json || parameters.print_both){
    if (!parameters.compact) {
      gt_stats_print_stats(stats[0],(parameters.num_reads>0)?
          parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
    } else {
      gt_stats_print_stats_compact(stats[0],(parameters.num_reads>0)?
          parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
    }
  }

  // Clean
  gt_stats_delete(stats[0]); gt_free(stats);
  gt_input_file_close(input_file);
}

void parse_arguments(int argc,char** argv) {
  struct option* gt_stats_getopt = gt_options_adaptor_getopt(gt_stats_options);
  gt_string* const gt_stats_short_getopt = gt_options_adaptor_getopt_short(gt_stats_options);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_stats_short_getopt),gt_stats_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i': // input
      parameters.name_input_file = optarg;
      break;
    case 200: // mmap-input
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'r': // reference
      parameters.name_reference_file = optarg;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'I': // gem-index
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'p': // paired-end
      parameters.paired_end = true;
      break;
    case 'n': // num-reads
      parameters.num_reads = atol(optarg);
      break;
    case 'o': // output
      parameters.name_output_file = optarg;
      break;
    case 'f':
      if(strcmp("report", optarg) == 0){
        parameters.print_json = false;
      }else if(strcmp("json", optarg) == 0){
        parameters.print_json = true;
      }else if(strcmp("both", optarg) == 0){
        parameters.print_json = true;
        parameters.print_both = true;
      }else{
        gt_fatal_error_msg("Unknown format %s, supported formats are 'report' or 'json' or 'both'", optarg);
      }
      break;
    /* Analysis */
    case 300: // first-map
      parameters.first_map = true;
      break;
    case 'a': // all-tests
      parameters.maps_profile = true;
      parameters.mismatch_transitions = true;
      parameters.mismatch_quality = true;
      parameters.splitmaps_profile = true;
      parameters.population_profile = true;
      break;
    case 'M': // maps-profile
      parameters.maps_profile = true;
      break;
    case 'T': // mismatch-transitions
      parameters.mismatch_transitions = true;
      break;
    case 'Q': // mismatch-quality
      parameters.mismatch_quality = true;
      break;
    case 'R': // rna-profile // FIXME name
      parameters.splitmaps_profile = true;
      break;
    case 'P': // population-profile
      parameters.population_profile = true;
      break;
    case 'D': // indel-profile
      gt_fatal_error(NOT_IMPLEMENTED);
      parameters.indel_profile = true;
      break;
    /* MAP Specific */
    case 400:
      parameters.use_only_decoded_maps = true;
      break;
    /* Misc */
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      fprintf(stderr, "USE: ./gt.stats [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_stats_options,gt_stats_groups,false,false);
      exit(1);
    case 'H':
      fprintf(stderr, "USE: ./gt.stats [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_stats_options,gt_stats_groups,false,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_stats_options,gt_stats_groups,false,true);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Checks
   */
  if (parameters.indel_profile && parameters.name_reference_file==NULL) {
    gt_error_msg("To generate the indel-profile, a reference file(.fa/.fasta) or GEMindex(.gem) is required");
  }
  // Free
  gt_string_delete(gt_stats_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  parameters.output_file = stdout;
  parameters.output_file_json = stderr;

  // init output paramters
  if(parameters.name_output_file != NULL){
    parameters.output_file = fopen(parameters.name_output_file, "w");
    gt_cond_fatal_error(parameters.output_file==NULL,FILE_OPEN,parameters.name_output_file);
  }
  if(parameters.print_json && !parameters.print_both){
    parameters.output_file_json = parameters.output_file;
  }
  // Extract stats
  gt_stats_parallel_generate_stats();
  // close output
  if(parameters.name_output_file != NULL){
    fclose(parameters.output_file);
  }

  return 0;
}


