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
  bool mmap_input;
  bool paired_end;
  /* Filter */
  bool no_split_maps;
  bool best_map;
  uint64_t max_matches;
  bool make_counters;
  bool realign_hamming;
  bool realign_levenshtein;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .name_reference_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    /* Filter */
    .no_split_maps=false,
    .best_map=false,
    .max_matches=GT_ALL,
    .make_counters=false,
    .realign_hamming=false,
    .realign_levenshtein=false,
    /* Misc */
    .num_threads=1,
    .verbose=false,
};

void gt_template_filter(gt_template *template_dst,gt_template *template_src) {
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_REDUCTION(template_dst,alignment_dst);
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      if (!parameters.no_split_maps || gt_map_get_num_blocks(map)==1) {
        gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
        if (parameters.best_map) return;
      }
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  GT_TEMPLATE__ATTR_ITERATE(template_src,mmap,mmap_attr) {
    // Check SM contained
    if (parameters.no_split_maps) {
      register bool has_sm;
      GT_MULTIMAP_ITERATE(mmap,map,end_p) {
        has_sm = gt_map_get_num_blocks(map)>1;
        if (has_sm) break;
      }
      if (has_sm) continue;
    }
    // Add the mmap
    register gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
    gt_template_add_mmap(template_dst,mmap_copy,mmap_attr);
    free(mmap_copy);
    // Skip the rest if best
    if (parameters.best_map) return;
  }
}

void gt_filter_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Open reference file
  register gt_sequence_archive* sequence_archive;
  if (parameters.name_reference_file!=NULL && (parameters.realign_hamming || parameters.realign_levenshtein)) {
    sequence_archive = gt_sequence_archive_new();
    register gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    fprintf(stderr,"Loading reference file ...");
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      fprintf(stderr,"\n");
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
    fprintf(stderr," done! \n");
  }

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
    gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
    gt_status error_code;
    gt_generic_parser_attr generic_parser_attr = GENERIC_PARSER_ATTR_DEFAULT(parameters.paired_end);

    // Limit max-matches
    generic_parser_attr.map_parser_attr.max_parsed_maps = parameters.max_matches;

    gt_template* template = gt_template_new();
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,&generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
      }

      // First realign
      if (parameters.realign_levenshtein) {
        gt_template_realign_levenshtein(template,sequence_archive);
      } else if (parameters.realign_hamming) {
        gt_template_realign_hamming(template,sequence_archive);
      }

      // Pick up best-map || erase splitmaps
      if (parameters.best_map || parameters.no_split_maps) {
        gt_template *template_best = gt_template_copy(template,false,false);
        gt_template_filter(template_best,template);
      }

      // Make counters
      if (parameters.make_counters) {
        gt_template_recalculate_counters(template);
      }

      // Print template
      gt_output_map_bofprint_template(buffered_output,template,output_attributes);
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.filter [ARGS]...\n"
                  "         [I/O]\n"
                  "           --input|-i [FILE]\n"
                  "           --output|-o [FILE]\n"
                  "           --reference|-r [FILE]\n"
                  "           --mmap-input\n"
                  "           --paired-end|p\n"
                  "         [Filter]\n"
                  "           --no-split-maps\n"
                  "           --best-map\n"
                  "           --max-matches <number>\n"
                  "           --make-counters <number>\n"
                  "           --hamming-realign\n"
                  "           --levenshtein-realign\n"
                  "         [Misc]\n"
                  "           --threads|t\n"
                  "           --verbose|v\n"
                  "           --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* I/O */
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "reference", required_argument, 0, 'r' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    /* Filter */
    { "no-split-maps", no_argument, 0, 2 },
    { "best-map", no_argument, 0, 3 },
    { "max-matches", required_argument, 0, 4 },
    { "make-counters", no_argument, 0, 5 },
    { "hamming-realign", no_argument, 0, 6 },
    { "levenshtein-realign", no_argument, 0, 7 },
    /* Misc */
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:r:pt:hv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      break;
    case 1:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    /* Filter */
    case 2:
      parameters.no_split_maps = true;
      break;
    case 3:
      parameters.best_map = true;
      break;
    case 4:
      parameters.max_matches = atoll(optarg);
      break;
    case 5:
      parameters.make_counters = true;
      break;
    case 6:
      parameters.realign_hamming = true;
      break;
    case 7:
      parameters.realign_levenshtein = true;
      break;
    /* Misc */
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
  /*
   * Parameters check
   */
  if (parameters.realign_hamming || parameters.realign_levenshtein) {
    if (parameters.name_reference_file==NULL) gt_error_msg("Reference file required to realign");
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  // Filter !
  gt_filter_read__write();

  return 0;
}

