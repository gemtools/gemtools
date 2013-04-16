/*
 * PROJECT: GEM-Tools library
 * FILE: gt.map.2.sam.c
 * DATE: 02/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Converter from MAP to SAM
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
  /* Headers */
  /* XXX */
  bool compact_format;
  /* Optional Fields */
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
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
  /* Headers */
  /* XXX */
  .compact_format=false,
  /* Optional Fields */
  /* Misc */
  .num_threads=1,
  .verbose=false,
  /* Control flags */
  .load_index=false
};

gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  register gt_sequence_archive* sequence_archive = gt_sequence_archive_new();
  gt_log("Loading reference file ...");
  if (parameters.name_gem_index_file!=NULL) { // Load GEM-IDX
    gt_gemIdx_load_archive(parameters.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    register gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  gt_log("Done.");
  return sequence_archive;
}

void gt_map2sam_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) sequence_archive = gt_filter_open_sequence_archive();

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    // I/O attributes
    gt_map_parser_attr input_attributes = GT_MAP_PARSER_ATTR_DEFAULT(parameters.paired_end);
    gt_output_sam_attributes output_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();

    gt_template* template = gt_template_new();
    while ((error_code=gt_input_map_parser_get_template_g(buffered_input,template,&generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
      }

      // Print template
      if (gt_output_map_bofprint_template(buffered_output,template,&output_attributes)) {
        gt_error_msg("Fatal error outputting read '"PRIgts"'(InputLine:%"PRIu64")\n",
            PRIgts_content(gt_template_get_string_tag(template)),buffered_input->current_line_num-1);
      }
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.filter [ARGS]...\n"
                  "         [I/O]\n"
                  "           --input|-i [FILE]\n"
                  "           --output|-o [FILE]\n"
                  "           --reference|-r [FILE] (MultiFASTA/FASTA)\n"
                  "           --gem-index|-I [FILE] (GEM2-Index)\n"
               //   "           --mmap-input\n"
                  "           --paired-end|p\n"
                  "         [Headers]\n"
                  "           --sorting 'unknown'|'unsorted'|'queryname'|'coordinate'\n"
                  "           --read-group  \n"
                  "           --program  \n"
                  "           --comment <string>"
                  "         [Alignments]\n"
                  "           --compact|c\n"
                  "           --score-alignments|s"
                  "         [Optional Fields]\n"
                  "           -- "
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
    { "gem-index", required_argument, 0, 'I' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    /* Headers */

    /* XXX */
    { "compact", no_argument, 0, 'c' },

    /* Optional Fields */

    /* Misc */
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:r:I:pct:hHv",long_options,&option_index);
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
    /* Headers */

    /* XXX */
    case 'c':
      parameters.compact_format = true;
      break;
    /* Optional Fields */

    /* Misc */
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
    case 'H':
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
  if (parameters.load_index && parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  // map2sam !!
  gt_map2sam_read__write();

  return 0;
}

