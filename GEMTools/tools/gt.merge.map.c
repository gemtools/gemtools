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
    .num_threads=1,
    .verbose=false,
};

pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;


GT_INLINE gt_status gt_merge_map_read_template_sync(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_buffered_output_file* buffered_output,gt_template* const template_master,gt_template* const template_slave) {
  register gt_status error_code_master, error_code_slave;
  do {
    // Read Synch blocks
    error_code_master=gt_input_map_parser_synch_blocks_by_subset(&input_mutex,buffered_input_master,buffered_input_slave,parameters.paired_end);
    if (error_code_master==GT_IMP_EOF) return GT_IMP_EOF;
    if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error synchronizing files");
    // Read master (always guaranteed)
    if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))==GT_IMP_FAIL) {
      gt_fatal_error_msg("Fatal error parsing file <<Master>>");
    }
    // Check slave
    if (gt_buffered_input_file_eob(buffered_input_slave)) { // Slave exhausted. Dump master & return EOF
      gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
      while ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))) {
        if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
      }
    } else {
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_slave,template_slave))==GT_IMP_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
      }
      // Synch loop
      while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
        // Print non correlative master's template
        gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
        // Fetch next master's template
        if (gt_buffered_input_file_eob(buffered_input_master)) gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
        if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))!=GT_IMP_OK) {
          gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        }
      }
      return GT_IMP_OK;
    }
  } while (true);
}

void gt_merge_map_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    // Parallel reading+process
    gt_buffered_input_file* buffered_input_master = gt_buffered_input_file_new(input_file_1);
    gt_buffered_input_file* buffered_input_slave  = gt_buffered_input_file_new(input_file_2);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input_master,buffered_output);

    gt_template *template_1 = gt_template_new();
    gt_template *template_2 = gt_template_new();
    while (gt_merge_map_read_template_sync(buffered_input_master,buffered_input_slave,buffered_output,template_1,template_2)) {
      // Merge maps
      register gt_template *ptemplate = gt_template_union_template_mmaps(template_1,template_2);
      // Print template
      gt_output_map_bofprint_gem_template(buffered_output,ptemplate,GT_ALL,true);
      // Delete template
      gt_template_delete(ptemplate);
    }

    // Clean
    gt_template_delete(template_1);
    gt_template_delete(template_2);
    gt_buffered_input_file_close(buffered_input_master);
    gt_buffered_input_file_close(buffered_input_slave);
    gt_buffered_output_file_close(buffered_output);
  }

  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.merge.map [OP] [ARGS]...\n"
                  "       [ARGS]\n"
                  "         --i1 [FILE]\n"
                  "         --i2 [FILE]\n"
                  "         --output|-o [FILE]\n"
                  "         --paired-end|p\n"
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
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { "threads", required_argument, 0, 't' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:t:phv",long_options,&option_index);
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


