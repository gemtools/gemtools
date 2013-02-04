/*
 * PROJECT: GEM-Tools library
 * FILE: gt.mapset.c
 * DATE: 08/11/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to perform set operations {UNION,INTERSECTION,DIFFERENCE} over alignment files {MAP,SAM}
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef enum { GT_MAP_SET_INTERSECTION, GT_MAP_SET_UNION, GT_MAP_SET_DIFFERENCE } gt_operation;

typedef struct {
  /* [I/O] */
  char* name_input_file_1;
  char* name_input_file_2;
  char* name_output_file;
  bool mmap_input;
  bool paired_end;
  bool files_contain_same_reads;
  /* [Misc] */
  uint64_t num_threads;
  bool verbose;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file_1=NULL,
    .name_input_file_2=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .files_contain_same_reads=false,
    .num_threads=1,
    .verbose=false,
};

void gt_merge_map_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Mutex
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_merge_map_files(&input_mutex,input_file_1,input_file_2,
        parameters.paired_end,parameters.files_contain_same_reads,output_file);
  }

  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.merge.map [ARGS]...\n"
                  "       [ARGS]\n"
                  "         --i1 [FILE]\n"
                  "         --i2 [FILE]\n"
                  "         --output|-o [FILE]\n"
                  "         --paired-end|-p\n"
                  "         --files-same-reads|-s\n"
                  "       [Misc]\n"
                  "         --threads|t\n"
                  "         --verbose|v\n"
                  "         --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  // Parse arguments
  struct option long_options[] = {
    { "i1", required_argument, 0, 1 },
    { "i2", required_argument, 0, 2 },
    { "mmap-input", no_argument, 0, 3 },
    { "output", required_argument, 0, 'o' },
    { "paired-end", no_argument, 0, 'p' },
    { "files-same-reads", no_argument, 0, 's' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { "threads", required_argument, 0, 't' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:t:pshv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 1:
      parameters.name_input_file_1 = optarg;
      break;
    case 2:
      parameters.name_input_file_2 = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 3:
      parameters.mmap_input = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 's':
      parameters.files_contain_same_reads = true;
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
  // Check parameters
  if (!parameters.name_input_file_1) {
    gt_fatal_error_msg("Input file 1 required (--i1)\n");
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  // Filter !
  gt_merge_map_read__write();

  return 0;
}


