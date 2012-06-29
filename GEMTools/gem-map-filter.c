/*
 * PROJECT: GEM-Tools library
 * FILE: gem-map-filter.c
 * DATE: 01/06/2012
 * DESCRIPTION: Basic tool to perform simple map file conversions, filters and treatment in general
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
    // Iterate over all elements
    GT_ALIGNMENT_ITERATE(template,alignment) {
      printf(">> %s\n",gt_template_get_tag(template));
      // coutners();
      GT_MAPS_ITERATE(alignment,map) {
        GT_MAP_BLOCKS_ITERATE(map,map_block) {
          printf("\t%s\t",gt_map_get_seq_name(map));
          printf("init_pos=%lu\t",gt_map_get_position(map));
          printf("end_pos=%lu\t",gt_map_get_position(map)+gt_map_get_length(map));
          printf("length=%lu\t",gt_map_get_length(map));
          printf("strand=%c\t",gt_map_get_direction(map)==FORWARD?'F':'R');
          printf("distance=%lu\n",gt_map_get_distance(map));
          printf("\t\t{ ");
          GT_MISMS_ITERATE(map_block,misms) {
            if (gt_misms_get_type(misms)==MISMS) {
              printf("(M,%lu,%c) ",gt_misms_get_position(misms),gt_misms_get_base(misms));
            } else {
              printf("(%c,%lu,%lu) ",gt_misms_get_type(misms)==INS?'I':'D',
                  gt_misms_get_position(misms),gt_misms_get_size(misms));
            }
          }
          printf("}\n");
        }
      }
    }
  }
  gt_buffered_map_input_close(map_input);

  // Close files
  gt_input_file_close(input_file);

  return 0;
}


