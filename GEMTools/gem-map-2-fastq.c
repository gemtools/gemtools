/*
 * PROJECT: GEM-Tools library
 * FILE: gem-map-2-fastq.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef struct {
  char *name_input_file;
  char *name_output_file;
  uint64_t num_threads;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .num_threads=1
};

void usage() {
  fprintf(stderr, "USE: ./gem-map-filter -i input -o output \n"
                  "        --input     <File>\n"
                  "        --output    <File>\n"
                  "        --help|h    \n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:T:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'T':
      parameters.num_threads = atol(optarg);
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
  parse_arguments(argc, argv);

  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);

  // Parallel working threads
  gt_buffered_map_input* map_input = gt_buffered_map_input_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  while ((error_code=gt_buffered_map_input_get_template(map_input,template))) {
    if (error_code==GT_BMI_FAIL) continue;
    register const uint64_t num_blocks = gt_template_get_num_blocks(template);
    register uint64_t i;
    for (i=0;i<num_blocks;++i) {
      register gt_alignment* alignment = gt_template_get_block(template,0);
      printf("@%s\n%s\n+\n%s\n",gt_template_get_tag(template),
          gt_alignment_get_read(alignment),gt_alignment_get_qualities(alignment));
    }
  }
  gt_buffered_map_input_close(map_input);

  // Close files
  gt_input_file_close(input_file);

  return 0;
}


