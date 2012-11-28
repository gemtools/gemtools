/*
 * PROJECT: GEM-Tools library
 * FILE: gt.mapset.c
 * DATE: 08/11/2012
 * DESCRIPTION: Application to perform set operations {UNION,INTERSECTION,DIFFERENCE} over alignment files {MAP,SAM}
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef enum { GT_MAP_SET_INTERSECTION, GT_MAP_SET_UNION, GT_MAP_SET_DIFFERENCE } gt_operation;

typedef struct {
  char* name_input_file_1;
  char* name_input_file_2;
  char* name_output_file;
  gt_operation operation;
  bool mmap_input;
  bool paired_end;
  uint64_t eq_threshold;
  bool verbose;
} gt_stats_args;

gt_stats_args parameters = {
    .name_input_file_1=NULL,
    .name_input_file_2=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .eq_threshold=50,
    .verbose=false,
};

int64_t gt_mapset_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  return gt_map_range_cmp(map_1,map_2,parameters.eq_threshold);
}
int64_t gt_mapset_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  return gt_mmap_range_cmp(map_1,map_2,num_maps,parameters.eq_threshold);
}

GT_INLINE gt_status gt_mapset_read_template_sync(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_template* const template_master,gt_template* const template_slave,
    gt_buffered_output_file* const buffered_output) {
  register bool synch = false, slave_read = false;
  // Read master
  register gt_status error_code_master, error_code_slave;
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,parameters.paired_end))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,parameters.paired_end))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  if (error_code_master==GT_IMP_EOF) {
    if (error_code_slave!=GT_IMP_EOF) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
    return GT_IMP_EOF;
  }
  while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
    // Print non correlative master's template
    gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
    // Fetch next master's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_master,template_master,parameters.paired_end))!=GT_IMP_OK) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
  }
  return GT_IMP_OK;
}

void gt_mapset_read__write() {
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Parallel reading+process
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  gt_status error_code;
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  while (gt_mapset_read_template_sync(buffered_input_1,buffered_input_2,template_1,template_2,buffered_output)) {
    // Apply operation
    register gt_template *ptemplate;
    switch (parameters.operation) {
      case GT_MAP_SET_UNION:
        ptemplate=gt_template_union_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        break;
      case GT_MAP_SET_INTERSECTION:
        ptemplate=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        break;
      case GT_MAP_SET_DIFFERENCE:
        ptemplate=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
    // Print template
    gt_output_map_bofprint_gem_template(buffered_output,ptemplate,GT_ALL,true);
    // Delete template
    gt_template_delete(ptemplate);
  }

  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);

  // Clean
  gt_input_file_close(input_file_1);
  gt_output_file_close(output_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.mapset [OP] [ARGS]...\n"
                  "       [OP]\n"
                  "         union\n"
                  "         intersection\n"
                  "         difference\n"
                  "       [ARGS]\n"
                  "         --i1 [FILE]\n"
                  "         --i2 [FILE]\n"
                  "         --mmap-input\n"
                  "         --output|-o [FILE]\n"
                  "         --paired-end|p\n"
                  "         --verbose|v\n"
                  "         --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  // Parse operation
  if (argc<=1) gt_fatal_error_msg("Please specify operation {union,intersection}");
  if (gt_streq(argv[1],"INTERSECCTION") || gt_streq(argv[1],"Intersection") || gt_streq(argv[1],"intersection")) {
    parameters.operation = GT_MAP_SET_INTERSECTION;
  } else if (gt_streq(argv[1],"UNION") || gt_streq(argv[1],"Union") || gt_streq(argv[1],"union")) {
    parameters.operation = GT_MAP_SET_UNION;
  } else if (gt_streq(argv[1],"DIFFERENCE") || gt_streq(argv[1],"Difference") || gt_streq(argv[1],"difference")) {
    parameters.operation = GT_MAP_SET_DIFFERENCE;
  } else {
    if (argv[1][0]=='I' || argv[1][0]=='i') {
      fprintf(stderr,"\tAssuming 'Intersection' ...");
      parameters.operation = GT_MAP_SET_INTERSECTION;
    } else if (argv[1][0]=='I' || argv[1][0]=='i') {
      fprintf(stderr,"\tAssuming 'Union' ...");
      parameters.operation = GT_MAP_SET_UNION;
    } else if (argv[1][0]=='D' || argv[1][0]=='d') {
      fprintf(stderr,"\tAssuming 'Difference' ...");
      parameters.operation = GT_MAP_SET_DIFFERENCE;
    } else {
      gt_fatal_error_msg("Unknown operation '%s'\n",argv[1]);
    }
  }
  argc--; argv++;
  // Parse arguments
  struct option long_options[] = {
    { "i2", required_argument, 0, 1 },
    { "i1", required_argument, 0, 2 },
    { "mmap-input", no_argument, 0, 3 },
    { "output", required_argument, 0, 'o' },
    { "paired-end", no_argument, 0, 'p' },
    { "eq-th", no_argument, 0, 4 },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:phv",long_options,&option_index);
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
    case 4:
      parameters.eq_threshold = atol(optarg);
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
  gt_mapset_read__write();

  return 0;
}


