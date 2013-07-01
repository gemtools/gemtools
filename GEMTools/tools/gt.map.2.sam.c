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
  gt_qualities_offset_t quality_format;
  /* Headers */

  /* SAM format */
  bool compact_format;
  /* Optional Fields */
  bool optional_field_NH;
  bool optional_field_NM;
  bool optional_field_XT;
  bool optional_field_XS;
  bool optional_field_md;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool load_index;
  bool load_index_sequences;
} gt_stats_args;

gt_stats_args parameters = {
  /* I/O */
  .name_input_file=NULL,
  .name_output_file=NULL,
  .name_reference_file=NULL,
  .name_gem_index_file=NULL,
  .mmap_input=false,
  .paired_end=false,
  .quality_format=GT_QUALS_OFFSET_33,
  /* Headers */
  /* SAM format */
  .compact_format=false,
  /* Optional Fields */
  .optional_field_NH=false,
  .optional_field_NM=false,
  .optional_field_XT=false,
  .optional_field_XS=false,
  .optional_field_md=false,
  /* Misc */
  .num_threads=1,
  .verbose=false,
  /* Control flags */
  .load_index=false,
  .load_index_sequences=false
};

gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
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
  return sequence_archive;
}

void gt_map2sam_read__write() {
  // Open file IN/OUT
  gt_input_file* const input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* const output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
  gt_sam_headers* const sam_headers = gt_sam_header_new(); // SAM headers

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) {
    sequence_archive = gt_filter_open_sequence_archive(parameters.load_index_sequences);
    gt_sam_header_set_sequence_archive(sam_headers,sequence_archive);
  }

  // Print SAM headers
  gt_output_sam_ofprint_headers_sh(output_file,sam_headers);

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    // I/O attributes
    gt_map_parser_attributes* const input_map_attributes = gt_input_map_parser_attributes_new(parameters.paired_end);
    gt_output_sam_attributes* const output_sam_attributes = gt_output_sam_attributes_new();
    // Set out attributes
    gt_output_sam_attributes_set_compact_format(output_sam_attributes,parameters.compact_format);
    gt_output_sam_attributes_set_qualities_offset(output_sam_attributes,parameters.quality_format);
    if (parameters.optional_field_NH) gt_sam_attributes_add_tag_NH(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_NM) gt_sam_attributes_add_tag_NM(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XT) gt_sam_attributes_add_tag_XT(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_md) gt_sam_attributes_add_tag_md(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XS) gt_sam_attributes_add_tag_XS(output_sam_attributes->sam_attributes);
    gt_template* template = gt_template_new();
    while ((error_code=gt_input_map_parser_get_template(buffered_input,template,input_map_attributes))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
        continue;
      }

      // Print SAM template
      gt_output_sam_bofprint_template(buffered_output,template,output_sam_attributes);
    }

    // Clean
    gt_template_delete(template);
    gt_input_map_parser_attributes_delete(input_map_attributes);
    gt_output_sam_attributes_delete(output_sam_attributes);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_sam_header_delete(sam_headers);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.map.2.sam [ARGS]...\n"
                  "         [I/O]\n"
                  "           --input|-i [FILE]\n"
                  "           --output|-o [FILE]\n"
                  "           --reference|-r [FILE] (MultiFASTA/FASTA)\n"
                  "           --gem-index|-I [FILE] (GEM2-Index)\n"
               // "           --mmap-input\n"
                  "           --paired-end|-p\n"
                  "           --quality-format|-q 'offset-33'|'offset-64'\n"
                  "         [Headers]\n"
               // "           --comment <string>"
                  "         [SAM format]\n"
                  "           --compact|-c\n"
                  "         [Optional Fields]\n"
                  "           --NH\n"
                  "           --NM\n"
                  "           --XT\n"
                  "           --XS\n"
                  "           --md\n"
                  "         [Misc]\n"
                  "           --threads|-t\n"
                  "           --verbose|-v\n"
                  "           --help|-h\n");
  /*
   * Pending ...
   *        --score-alignments|s"
   *
   *
   */
  //                  "           --RG \n" // TODO: Bufff RG:Z:0 NH:i:16 XT:A:U

  // "           --headers [FILE] (Only {@RG,@PG,@CO} lines)\n"
  // "           --sorting 'unknown'|'unsorted'|'queryname'|'coordinate'\n"
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
    { "quality-format", no_argument, 0, 'q' },
    /* Headers */

    /* SAM format */
    { "compact", no_argument, 0, 'c' },

    /* Optional Fields */
    { "NH", no_argument, 0, 50 },
    { "NM", no_argument, 0, 51 },
    { "XT", no_argument, 0, 52 },
    { "md", no_argument, 0, 53 },
    { "XS", no_argument, 0, 54 },
    /* Misc */
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:r:I:pq:ct:hHv",long_options,&option_index);
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
      parameters.load_index = true;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      parameters.load_index = true;
      break;
    case 1:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 'q':
      if (gt_streq(optarg,"offset-64")) {
        parameters.quality_format=GT_QUALS_OFFSET_64;
      } else if (gt_streq(optarg,"offset-33")) {
        parameters.quality_format=GT_QUALS_OFFSET_33;
      } else {
        gt_fatal_error_msg("Quality format not recognized: '%s'",optarg);
      }
      break;
    /* Headers */

    /* SAM format */
    case 'c':
      parameters.compact_format = true;
      break;
    /* Optional Fields */
    case 50: // NH
      parameters.optional_field_NH = true;
      break;
    case 51: // NM
      parameters.optional_field_NM = true;
      break;
    case 52: // XT
      parameters.optional_field_XT = true;
      break;
    case 53: // md
      parameters.optional_field_md = true;
      break;
    case 54: // XS
      parameters.optional_field_XS = true;
      break;
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
  if(!parameters.load_index && parameters.optional_field_XS){
    gt_fatal_error_msg("Reference file required to compute XS field in SAM");
  }
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // map2sam !!
  gt_map2sam_read__write();

  return 0;
}

