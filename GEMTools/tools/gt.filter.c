/*
 * PROJECT: GEM-Tools library
 * FILE: gt.filter.c
 * DATE: 02/08/2012
 * DESCRIPTION: Application to filter {MAP,SAM} files and output the filtered result
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE false

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
  char* name_input_file;
  char* name_output_file;
  bool mmap_input;
  uint64_t num_threads;
  bool verbose;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .num_threads=1,
    .verbose=false,
};

void gt_filter_read__write() {
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
    while ((error_code=gt_input_map_parser_get_template(buffered_input,template))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // DO STH // TODO
    }

    // Clean
    gt_template_delete(template,true);
    gt_buffered_input_file_close(buffered_input);
  }

  // Clean
  gt_input_file_close(input_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.stats [ARGS]...\n"
                  "        --input|-i [FILE]\n"
                  "        --output|-o [FILE]\n"
                  "        --mmap-input\n"
                  "        --threads|t\n"
                  "        --verbose|v\n"
                  "        --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "mmap-input", no_argument, 0, 1 },
    { "threads", no_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:t:phv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 0:
      parameters.mmap_input = true;
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

  // Filter !
  gt_filter_read__write();

  return 0;
}


