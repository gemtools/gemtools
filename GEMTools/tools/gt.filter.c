/*
 * PROJECT: GEM-Tools library
 * FILE: gt.filter.c
 * DATE: 02/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Application to filter {MAP,SAM,FASTQ} files and output the filtered result
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef struct {
  /* I/O */
  char* name_input_file;
  char* name_output_file;
  char* name_reference_file;
  char* name_gem_index_file;
  bool mmap_input;
  bool paired_end;
  bool no_output;
  /* Filter Template/Alignments */
  bool mapped;
  bool unmapped;
  /* Filter Maps */
  bool no_split_maps;
  bool only_split_maps;
  bool best_map;
  uint64_t max_matches;
  bool make_counters;
  bool only_unmapped;
  bool only_mapped;
  uint64_t min_event_distance;
  uint64_t max_event_distance;
  uint64_t min_levenshtein_distance;
  uint64_t max_levenshtein_distance;
  gt_vector* filter_map_ids;
  /* Filter-pairs */
  int64_t max_inss;
  int64_t min_inss;
  bool filter_by_strand;
  bool allow_strand_rf;
  bool allow_strand_fr;
  bool allow_strand_ff;
  bool allow_strand_rr;
  /* Filter-Realign */
  bool mismatch_recovery;
  bool realign_hamming;
  bool realign_levenshtein;
  /* Checking/Report */
  bool check;
  /* Hidden */
  bool special_functionality;
  bool error_plot;
  bool insert_size_plot;
  bool show_sequence_list;
  bool display_pretty;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool perform_map_filter;
  bool load_index;
} gt_stats_args;

gt_stats_args parameters = {
    /* I/O */
    .name_input_file=NULL,
    .name_output_file=NULL,
    .name_reference_file=NULL,
    .name_gem_index_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .no_output=false,
    /* Filter Template/Alignments */
    .mapped=false,
    .unmapped=false,
    /* Filter Maps */
    .no_split_maps=false,
    .only_split_maps=false,
    .best_map=false,
    .max_matches=GT_ALL,
    .make_counters=false,
    .only_unmapped=false,
    .only_mapped=false,
    .min_event_distance=0,
    .max_event_distance=UINT64_MAX,
    .min_levenshtein_distance=0,
    .max_levenshtein_distance=UINT64_MAX,
    .filter_map_ids=NULL,
    /* Filter-pairs */
    .max_inss=INT64_MAX,
    .min_inss=INT64_MIN,
    .filter_by_strand=false,
    .allow_strand_rf=false,
    .allow_strand_fr=false,
    .allow_strand_ff=false,
    .allow_strand_rr=false,
    /* Filter-Realign */
    .mismatch_recovery=false,
    .realign_hamming=false,
    .realign_levenshtein=false,
    /* Checking/Report */
    .check = false,
    /* Hidden */
    .special_functionality = false,
    .error_plot = false,
    .insert_size_plot = false,
    .show_sequence_list = false,
    .display_pretty = false,
    /* Misc */
    .num_threads=1,
    .verbose=false,
    /* Control flags */
    .perform_map_filter=false,
    .load_index=false
};

void gt_filter_delete_map_ids(gt_vector* filter_map_ids) {
  // Free vector
  if (filter_map_ids!=NULL) {
    GT_VECTOR_ITERATE(filter_map_ids,map_id,pos,gt_string*) {
      gt_string_delete(*map_id);
    }
    gt_vector_delete(filter_map_ids);
  }
}

GT_INLINE bool gt_filter_is_sequence_name_allowed(gt_string* const seq_name) {
  GT_VECTOR_ITERATE(parameters.filter_map_ids,map_id,pos,gt_string*) {
    if (gt_string_equals(seq_name,*map_id)) return true;
  }
  return false;
}

void gt_template_filter(gt_template* template_dst,gt_template* template_src) {
  /*
   * SE
   */
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_REDUCTION(template_dst,alignment_dst);
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      // Check sequence name
      if (parameters.filter_map_ids!=NULL) {
        if (!gt_filter_is_sequence_name_allowed(map->seq_name)) continue;
      }
      // Check SM contained
      register const uint64_t num_blocks = gt_map_get_num_blocks(map);
      if (parameters.no_split_maps && num_blocks>1) continue;
      if (parameters.only_split_maps && num_blocks==1) continue;
      // Check strata
      if (parameters.min_event_distance > 0 || parameters.max_event_distance < UINT64_MAX) {
        register const int64_t total_distance = gt_map_get_global_distance(map);
        if (parameters.min_event_distance > total_distance || total_distance > parameters.max_event_distance) continue;
      }
      // Check levenshtein distance
      if (parameters.min_levenshtein_distance > 0 || parameters.max_levenshtein_distance < UINT64_MAX) {
        register const int64_t total_distance = gt_map_get_global_levenshtein_distance(map);
        if (parameters.min_levenshtein_distance > total_distance || total_distance > parameters.max_levenshtein_distance) continue;
      }
      // Insert the map
      gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
      // Skip the rest if best
      if (parameters.best_map) return;
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  /*
   * PE
   */
  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  GT_TEMPLATE__ATTR_ITERATE(template_src,mmap,mmap_attr) {
    // Check sequence name
    if (parameters.filter_map_ids!=NULL) {
      if (!gt_filter_is_sequence_name_allowed(mmap[0]->seq_name)) continue;
      if (!gt_filter_is_sequence_name_allowed(mmap[1]->seq_name)) continue;
    }
    // Check SM contained
    register uint64_t has_sm = false;
    if (parameters.no_split_maps || parameters.only_split_maps) {
      GT_MULTIMAP_ITERATE(mmap,map,end_p) {
        if (gt_map_get_num_blocks(map)>1) {
          has_sm = true; break;
        }
      }
    }
    if (parameters.no_split_maps && has_sm) continue;
    if (parameters.only_split_maps && !has_sm) continue;
    // Check strata
    if (parameters.min_event_distance > 0 || parameters.max_event_distance < UINT64_MAX) {
      register const int64_t total_distance = gt_map_get_global_distance(mmap[0])+gt_map_get_global_distance(mmap[1]);
      if (parameters.min_event_distance > total_distance || total_distance > parameters.max_event_distance) continue;
    }
    // Check levenshtein distance
    if (parameters.min_levenshtein_distance > 0 || parameters.max_levenshtein_distance < UINT64_MAX) {
      register const int64_t total_distance = gt_map_get_global_levenshtein_distance(mmap[0])+gt_map_get_global_levenshtein_distance(mmap[1]);
      if (parameters.min_levenshtein_distance > total_distance || total_distance > parameters.max_levenshtein_distance) continue;
    }
    // Check inss
    if (parameters.min_inss > INT64_MIN || parameters.max_inss < INT64_MAX) {
      uint64_t gt_err;
      register const int64_t inss = gt_template_get_insert_size(mmap,&gt_err);
      if (parameters.min_inss > inss || inss > parameters.max_inss) continue;
    }
    // Check strandness
    if (parameters.filter_by_strand) {
      if (mmap[0]->strand==FORWARD && mmap[1]->strand==REVERSE && !parameters.allow_strand_fr) continue;
      if (mmap[0]->strand==REVERSE && mmap[1]->strand==FORWARD && !parameters.allow_strand_rf) continue;
      if (mmap[0]->strand==FORWARD && mmap[1]->strand==FORWARD && !parameters.allow_strand_ff) continue;
      if (mmap[0]->strand==REVERSE && mmap[1]->strand==REVERSE && !parameters.allow_strand_rr) continue;
    }
    // Add the mmap
    register gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
    gt_template_insert_mmap(template_dst,mmap_copy,mmap_attr);
    free(mmap_copy);
    // Skip the rest if best
    if (parameters.best_map) return;
  }
}

void gt_filter_hidden_options(FILE* stream,gt_template* const template,gt_sequence_archive* const sequence_archive) {
  if (parameters.error_plot) {
    /*
     * Print levenshtein distance of the maps
     */
    if (parameters.best_map)  {
      uint64_t best_distance = UINT64_MAX;
      GT_TEMPLATE_ITERATE_(template,mmap) {
        register const uint64_t dist = gt_map_get_global_levenshtein_distance(*mmap);
        if (dist < best_distance) best_distance = dist;
      }
      if (best_distance < UINT64_MAX) fprintf(stream,"%lu\n",best_distance);
    } else {
      GT_TEMPLATE_ITERATE_(template,mmap) {
        fprintf(stream,"%lu\n",gt_map_get_global_levenshtein_distance(*mmap));
      }
    }
  } else if (parameters.insert_size_plot && gt_template_get_num_blocks(template)>1) {
    /*
     * Print insert size
     */
    uint64_t error_code;
    GT_TEMPLATE_ITERATE_(template,mmap) {
      fprintf(stream,"%lu\n",gt_template_get_insert_size(mmap,&error_code));
      if (parameters.best_map) break;
    }
  }
}

GT_INLINE void gt_filter_mismatch_recovery_maps(
    char* const name_input_file,const uint64_t current_line_num,
    gt_template* const template,gt_sequence_archive* const sequence_archive) {
  // Unfolded as to report errors in the recovery
  register gt_status error_code;
  register uint64_t alignment_pos = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    register uint64_t map_pos = 0;
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if ((error_code=gt_map_recover_mismatches_sa(map,alignment->read,sequence_archive))) {
        gt_error_msg("Unrecoverable Alignment '%s':%"PRIu64"\n\tREAD::'"PRIgts"':%"PRIu64":%"PRIu64" ",
            name_input_file,current_line_num,PRIgts_content(template->tag),alignment_pos,map_pos);
        gt_output_map_fprint_map_pretty_sa(stderr,map,alignment->read,sequence_archive);
      }
      ++map_pos;
    }
    gt_alignment_recalculate_counters(alignment);
    ++alignment_pos;
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}

GT_INLINE void gt_filter_check_maps(
    char* const name_input_file,const uint64_t current_line_num,
    gt_template* const template,gt_sequence_archive* const sequence_archive,
    uint64_t* const total_algs_checked,uint64_t* const total_algs_correct,
    uint64_t* const total_maps_checked,uint64_t* const total_maps_correct) {
  register bool alignment_correct=true;
  register gt_status error_code;
  register uint64_t alignment_pos = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    register uint64_t map_pos = 0;
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if ((error_code=gt_map_check_alignment_sa(map,alignment->read,sequence_archive))) {
        gt_error_msg("Wrong Alignment '%s':%"PRIu64"\n\tREAD::'"PRIgts"':%"PRIu64":%"PRIu64" ",
            name_input_file,current_line_num,PRIgts_content(template->tag),alignment_pos,map_pos);
        gt_output_map_fprint_map_pretty_sa(stderr,map,alignment->read,sequence_archive);
        alignment_correct = false;
      } else {
        ++(*total_maps_correct);
      }
      ++(*total_maps_checked);
      ++map_pos;
    }
    ++alignment_pos;
  }
  ++(*total_algs_checked);
  if (alignment_correct) ++(*total_algs_correct);
}

gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  gt_log("Loading reference file ...");
  if (parameters.name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(parameters.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  gt_log("Done.");
  return sequence_archive;
}

void gt_filter_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* output_file;

  // Open out file
  if (!parameters.no_output) {
    output_file = (parameters.name_output_file==NULL) ?
          gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
  }

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) {
    sequence_archive = gt_filter_open_sequence_archive(true);
  }

  // Parallel reading+process
  uint64_t total_algs_checked=0, total_algs_correct=0, total_maps_checked=0, total_maps_correct=0;
  #pragma omp parallel num_threads(parameters.num_threads) \
                       reduction(+:total_algs_checked,total_algs_correct,total_maps_checked,total_maps_correct)
  {
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = NULL;
    if (!parameters.no_output) {
      buffered_output = gt_buffered_output_file_new(output_file);
      gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
    }
    // I/O attributes
    gt_generic_parser_attr generic_parser_attr = GENERIC_PARSER_ATTR_DEFAULT(parameters.paired_end);
    gt_output_map_attributes output_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
    generic_parser_attr.map_parser_attr.max_parsed_maps = parameters.max_matches; // Limit max-matches

    gt_template* template = gt_template_new();
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,&generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
      }

      // Consider mapped/unmapped
      register const bool is_mapped = gt_template_is_mapped(template);
      if (parameters.mapped && !is_mapped) continue;
      if (parameters.unmapped && is_mapped) continue;

      // Hidden options (aborts the rest)
      if (parameters.special_functionality) {
        gt_filter_hidden_options(stdout,template,sequence_archive);
      } else {
        // First realign
        if (parameters.realign_levenshtein) {
          gt_template_realign_levenshtein(template,sequence_archive);
        } else if (parameters.realign_hamming) {
          gt_template_realign_hamming(template,sequence_archive);
        } else if (parameters.mismatch_recovery) {
          gt_filter_mismatch_recovery_maps(parameters.name_input_file,buffered_input->current_line_num-1,template,sequence_archive);
        }

        // Map level filtering
        if (parameters.perform_map_filter) {
          gt_template *template_filtered = gt_template_copy(template,false,false);
          gt_template_filter(template_filtered,template);
          gt_template_delete(template);
          template = template_filtered;
        }

        // Make counters
        if (parameters.make_counters || parameters.perform_map_filter) {
          gt_template_recalculate_counters(template);
        }

        // Check
        if (parameters.check) {
          gt_filter_check_maps(parameters.name_input_file,buffered_input->current_line_num-1,
              template,sequence_archive,&total_algs_checked,&total_algs_correct,&total_maps_checked,&total_maps_correct);
        }

        // Print template
        if (!parameters.no_output && gt_output_map_bofprint_template(buffered_output,template,&output_attributes)) {
          gt_error_msg("Fatal error outputting read '"PRIgts"'(InputLine:%"PRIu64")\n",
              PRIgts_content(gt_template_get_string_tag(template)),buffered_input->current_line_num-1);
        }
      }
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
    if (!parameters.no_output) gt_buffered_output_file_close(buffered_output);
  }

  // Print check report
   if (parameters.check) {
     gt_log("Checked %lu alignments. Total.Correct %lu (%2.3f %%). Total.Maps.Correct %lu (%2.3f %%)",
         total_algs_checked,total_algs_correct,GT_GET_PERCENTAGE(total_algs_correct,total_algs_checked),
         total_maps_correct,GT_GET_PERCENTAGE(total_maps_correct,total_maps_checked));
   }

  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_filter_delete_map_ids(parameters.filter_map_ids);
  gt_input_file_close(input_file);
  if (!parameters.no_output) gt_output_file_close(output_file);
}

void usage(const bool print_hidden) {
  fprintf(stderr, "USE: ./gt.filter [ARGS]...\n"
                  "         [I/O]\n"
                  "           --input|-i [FILE]\n"
                  "           --output|-o [FILE]\n"
                  "           --reference|-r [FILE] (MultiFASTA/FASTA)\n"
                  "           --gem-index|-I [FILE] (GEM2-Index)\n"
                  "           --mmap-input\n"
                  "           --paired-end|p\n"
                  "           --no_output\n"
                  "         [Filter]\n"
                  "           --unmapped|--mapped\n"
                  "           --no-split-maps|--only-split-maps\n"
                  "           --best-map\n"
                  "           --max-matches <number>\n"
                  "           --make-counters <number>\n"
                  "           --min-strata <number>\n"
                  "           --max-strata <number>\n"
                  "           --min-levenshtein-error <number>\n"
                  "           --max-levenshtein-error <number>\n"
                  "           --map-id [SequenceId],... (Eg 'Chr1','Chr2')\n"
                  "         [Filter-pairs]\n"
                  "           --pair-strandness [COMB],...\n"
                  "               [COMB] := 'FR'|'RF'|'FF'|'RR'\n"
                  "           --min-inss <number>\n"
                  "           --max-inss <number>\n"
                  "         [Filter-Realign/Check]\n"
                  "           --check\n"
                  "           --mismatch-recovery\n"
                  "           --hamming-realign\n"
                  "           --levenshtein-realign\n"
                  "         [Misc]\n"
                  "           --threads|t\n"
                  "           --verbose|v\n"
                  "           --help|h\n");
  if (print_hidden) {
  fprintf(stderr, "         [Hidden]\n"
                  "           --sequence-list\n"
                  "           --display-pretty\n"
                  "           --error-plot\n"
                  "           --insert-size-plot\n");
  }
}

void gt_filter_get_argument_pair_strandness(char* const strandness_opt) {
  char *opt;
  opt = strtok(strandness_opt,",");
  while (opt!=NULL) {
    if (gt_streq(opt,"FR")) {
      parameters.allow_strand_fr = true;
    } else if (gt_streq(opt,"RF")) {
      parameters.allow_strand_rf = true;
    } else if (gt_streq(opt,"FF")) {
      parameters.allow_strand_ff = true;
    } else if (gt_streq(opt,"RR")) {
      parameters.allow_strand_rr = true;
    } else {
      gt_fatal_error_msg("Strandness option not recognized '%s'\n",opt);
    }
    opt = strtok(NULL,","); // Reload
  }
  parameters.filter_by_strand = true;
}

void gt_filter_get_argument_map_id(char* const strandness_opt) {
  // Allocate vector
  parameters.filter_map_ids = gt_vector_new(20,sizeof(gt_string*));
  // Add all the valid map Ids (sequence names)
  char *opt;
  opt = strtok(strandness_opt,",");
  while (opt!=NULL) {
    // Get id
    gt_string* map_id = gt_string_new(0);
    gt_string_set_string(map_id,opt);
    // Add to the vector
    gt_vector_insert(parameters.filter_map_ids,map_id,gt_string*);
    // Next
    opt = strtok(NULL,","); // Reload
  }
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* I/O */
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "reference", required_argument, 0, 'r' },
    { "gem-index", required_argument, 0, 'I' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    { "no_output", no_argument, 0, 14 },
    /* Filter Template/Alignments */
    { "mapped", no_argument, 0, 2 },
    { "unmapped", no_argument, 0, 3 },
    /* Filter Maps */
    { "no-split-maps", no_argument, 0, 4 },
    { "only-split-maps", no_argument, 0, 5 },
    { "best-map", no_argument, 0, 6 },
    { "max-matches", required_argument, 0, 7 },
    { "make-counters", no_argument, 0, 8 },
    { "min-strata", required_argument, 0, 9 },
    { "max-strata", required_argument, 0, 10 },
    { "min-levenshtein-error", required_argument, 0, 11 },
    { "max-levenshtein-error", required_argument, 0, 12 },
    { "map-id", required_argument, 0, 13 },
    /* Filter-pairs */
    { "pair-strandness", required_argument, 0, 30 },
    { "min-inss", required_argument, 0, 31 },
    { "max-inss", required_argument, 0, 32 },
    /* Filter-Realign */
    { "mismatch-recovery", no_argument, 0, 40 },
    { "hamming-realign", no_argument, 0, 41 },
    { "levenshtein-realign", no_argument, 0, 42 },
    /* Checking/Report */
    { "check", no_argument, 0, 50 },
    /* Hidden */
    { "error-plot", no_argument, 0, 55 },
    { "insert-size-plot", no_argument, 0, 56 },
    { "sequence-list", no_argument, 0, 57 },
    { "display-pretty", no_argument, 0, 58 },
    /* Misc */
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:r:I:pt:hHv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      if (gt_streq(optarg,"null")) parameters.no_output = true;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      break;
    case 1:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 14:
      parameters.no_output = true;
      break;
    /* Filter Template/Alignments */
    case 2: // mapped
      parameters.mapped = true;
      break;
    case 3:
      parameters.unmapped = true;
      break;
    /* Filter Maps */
    case 4:
      parameters.perform_map_filter = true;
      parameters.no_split_maps = true;
      break;
    case 5:
      parameters.perform_map_filter = true;
      parameters.only_split_maps = true;
      break;
    case 6:
      parameters.perform_map_filter = true;
      parameters.best_map = true;
      break;
    case 7:
      parameters.max_matches = atoll(optarg);
      break;
    case 8:
      parameters.make_counters = true;
      break;
    case 9:
      parameters.perform_map_filter = true;
      parameters.min_event_distance = atoll(optarg);
      break;
    case 10:
      parameters.perform_map_filter = true;
      parameters.max_event_distance = atoll(optarg);
      break;
    case 11:
      parameters.perform_map_filter = true;
      parameters.min_levenshtein_distance = atoll(optarg);
      break;
    case 12:
      parameters.perform_map_filter = true;
      parameters.max_levenshtein_distance = atoll(optarg);
      break;
    case 13: // map-id
      parameters.perform_map_filter = true;
      gt_filter_get_argument_map_id(optarg);
      break;
    /* Filter-pairs */
    case 30: // pair-strandness
      parameters.perform_map_filter = true;
      gt_filter_get_argument_pair_strandness(optarg);
      break;
    case 31:
      parameters.perform_map_filter = true;
      parameters.min_inss = atoll(optarg);
      break;
    case 32:
      parameters.perform_map_filter = true;
      parameters.max_inss = atoll(optarg);
      break;
    /* Filter-Realign */
    case 40:
      parameters.load_index = true;
      parameters.mismatch_recovery = true;
      break;
    case 41:
      parameters.load_index = true;
      parameters.realign_hamming = true;
      break;
    case 42:
      parameters.load_index = true;
      parameters.realign_levenshtein = true;
      break;
    /* Checking/Report */
    case 50:
      parameters.load_index = true;
      parameters.check = true;
      break;
    /* Hidden */
    case 55:
      parameters.special_functionality = true;
      parameters.error_plot = true;
      break;
    case 56:
      parameters.special_functionality = true;
      parameters.insert_size_plot = true;
      break;
    case 57:
      parameters.special_functionality = true;
      parameters.load_index = true;
      parameters.show_sequence_list = true;
      break;
    case 58:
      parameters.special_functionality = true;
      parameters.load_index = true;
      parameters.display_pretty = true;
      break;
    /* Misc */
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      usage(false);
      exit(1);
    case 'H':
      usage(true);
      exit(1);
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters check
   */
  if (parameters.load_index && parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  /*
   * Select functionality
   */
  if (parameters.show_sequence_list) {
    /*
     * Show sequence archive summary
     */
    gt_sequence_archive* sequence_archive = gt_filter_open_sequence_archive(false);
    gt_sequence_archive_iterator sequence_archive_it;
    gt_sequence_archive_new_iterator(sequence_archive,&sequence_archive_it);
    register gt_segmented_sequence* seq;
    while ((seq=gt_sequence_archive_iterator_next(&sequence_archive_it))) {
      fprintf(stderr,"%s\t%lu\n",seq->seq_name->buffer,seq->sequence_total_length);
    }
  } else {
    /*
     * Filter !!
     */
    gt_filter_read__write();
  }
  return 0;
}

