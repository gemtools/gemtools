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

typedef enum { GT_MAP_SET_UNKNOWN,
               GT_MAP_SET_INTERSECTION, GT_MAP_SET_UNION, GT_MAP_SET_DIFFERENCE,
               GT_MAP_SET_JOIN, GT_MAP_SET_COMPARE,
               GT_MERGE_MAP, GT_DISPLAY_COMPACT_MAP} gt_operation;

typedef struct {
  gt_operation operation;
  char* name_input_file_1;
  char* name_input_file_2;
  char* name_output_file;
  bool mmap_input;
  bool paired_end;
  bool files_contain_same_reads;
  double eq_threshold;
  bool strict;
  bool verbose;
  uint64_t num_threads;
} gt_stats_args;

gt_stats_args parameters = {
    .operation=GT_MAP_SET_UNKNOWN,
    .name_input_file_1=NULL,
    .name_input_file_2=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .files_contain_same_reads=true,
    .eq_threshold=0.5,
    .strict=false,
    .verbose=false,
    .num_threads=1
};
uint64_t current_read_length;

int64_t gt_mapset_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  return parameters.strict ? gt_map_cmp(map_1,map_2) : gt_map_range_cmp(map_1,map_2,eq_threshold);
}
int64_t gt_mapset_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  return parameters.strict ? gt_mmap_cmp(map_1,map_2,num_maps) : gt_mmap_range_cmp(map_1,map_2,num_maps,eq_threshold);
}

GT_INLINE gt_status gt_mapset_read_template_sync(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_buffered_output_file* const buffered_output,gt_template* const template_master,gt_template* const template_slave,
    const gt_operation operation) {
  // Read master
  gt_status error_code_master, error_code_slave;
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired_end);
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  // Check EOF conditions
  if (error_code_master==GT_IMP_EOF) {
    if (error_code_slave!=GT_IMP_EOF) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
    return GT_IMP_EOF;
  } else if (error_code_slave==GT_IMP_EOF) { // Slave exhausted. Dump master & return EOF
    do {
      if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
      if (operation==GT_MAP_SET_UNION || operation==GT_MAP_SET_DIFFERENCE) {
        gt_output_map_bofprint_template(buffered_output,template_master,output_attributes);
      }
    } while ((error_code_master=gt_input_generic_parser_get_template(
                buffered_input_master,template_master,generic_parser_attr)));
    return GT_IMP_EOF;
  }
  // Synch loop
  while (gt_string_cmp(gt_template_get_string_tag(template_master),
      gt_template_get_string_tag(template_slave))) {
    // Print non correlative master's template
    if (operation==GT_MAP_SET_UNION || operation==GT_MAP_SET_DIFFERENCE) {
      gt_output_map_bofprint_template(buffered_output,template_master,output_attributes);
    }
    // Fetch next master's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_master,template_master,generic_parser_attr))!=GT_IMP_OK) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
  }
  return GT_IMP_OK;
}

GT_INLINE gt_status gt_mapset_read_template_get_commom_map(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_template* const template_master,gt_template* const template_slave) {
  gt_status error_code_master, error_code_slave;
  gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired_end);
  // Read master
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  if (error_code_master==GT_IMP_EOF) return GT_IMP_EOF;
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  if (error_code_slave==GT_IMP_EOF) { // Check EOF conditions
    gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
        PRIgts_content(gt_template_get_string_tag(template_master)));
  }
  // Synch loop
  while (gt_string_cmp(gt_template_get_string_tag(template_master),gt_template_get_string_tag(template_slave))) {
    // Fetch next slave's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_slave,template_slave,generic_parser_attr))!=GT_IMP_OK) {
      gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
          PRIgts_content(gt_template_get_string_tag(template_master)));
    }
  }
  return GT_IMP_OK;
}

void gt_mapset_perform_set_operations() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  // Template I/O (synch)
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_sync(buffered_input_1,buffered_input_2,
      buffered_output,template_1,template_2,parameters.operation)) {
    // Record current read length
    current_read_length = gt_template_get_total_length(template_1);
    // Apply operation
    gt_template *ptemplate;
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
    gt_output_map_bofprint_template(buffered_output,ptemplate,output_attributes);
    // Delete template
    gt_template_delete(ptemplate);
  }

  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}

void gt_mapset_perform_cmp_operations() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  // Template I/O (synch)
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  gt_output_map_attributes* output_map_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_get_commom_map(buffered_input_1,buffered_input_2,template_1,template_2)) {
    // Record current read length
    current_read_length = gt_template_get_total_length(template_1);
    // Apply operation
    switch (parameters.operation) {
      case GT_MAP_SET_JOIN:
        // Print Master's TAG+Counters+Maps
        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_1),
            template_1->attributes,output_map_attributes); // Master's Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_2),
            template_1->attributes,output_map_attributes); // Slave's Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_1,output_map_attributes); // Master's Maps
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_2,output_map_attributes); // Slave's Maps
        gt_bofprintf(buffered_output,"\n");
        break;
      case GT_MAP_SET_COMPARE: {
        // Perform simple cmp operations
        gt_template *template_master_minus_slave=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        gt_template *template_slave_minus_master=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_2,template_1);
        gt_template *template_intersection=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        /*
         * Print results :: (TAG (Master-Slave){COUNTER MAPS} (Slave-Master){COUNTER MAPS} (Intersection){COUNTER MAPS})
         */
        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
        // Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_master_minus_slave),
            template_master_minus_slave->attributes,output_map_attributes); // (Master-Slave){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_slave_minus_master),
            template_slave_minus_master->attributes,output_map_attributes); // (Slave-Master){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_intersection),
            template_intersection->attributes,output_map_attributes); // (Intersection){COUNTER}
        // Maps
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_master_minus_slave,output_map_attributes); // (Master-Slave){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_slave_minus_master,output_map_attributes); // (Slave-Master){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_intersection,output_map_attributes); // (Intersection){COUNTER}
        gt_bofprintf(buffered_output,"\n");
        // Delete templates
        gt_template_delete(template_master_minus_slave);
        gt_template_delete(template_slave_minus_master);
        gt_template_delete(template_intersection);
        }
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
  }

  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}

void gt_mapset_perform_merge_map() {
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
    if (parameters.files_contain_same_reads) {
      gt_merge_synch_map_files(&input_mutex,parameters.paired_end,output_file,input_file_1,input_file_2);
    } else {
      gt_merge_unsynch_map_files(&input_mutex,input_file_1,input_file_2,parameters.paired_end,output_file);
    }
  }

  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
void gt_mapset_display_compact_map() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file_1==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_output_map_attributes* const output_map_attributes = gt_output_map_attributes_new();
    output_map_attributes->compact = true;
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        // Print compact summary
        gt_bofprintf(buffered_output,"End1::"PRIgts"[%"PRIu64"]\t",PRIgts_content(alignment->tag),gt_string_get_length(alignment->read));
        gt_output_map_bofprint_counters(buffered_output,alignment->counters,alignment->attributes,output_map_attributes);
        gt_bofprintf(buffered_output,"\t");
        uint64_t printed = 0;
        GT_ALIGNMENT_ITERATE(alignment,map) {
          if (printed>0) {
            gt_bofprintf(buffered_output,","PRIgts,PRIgts_content(map->seq_name));
          } else {
            gt_bofprintf(buffered_output,PRIgts,PRIgts_content(map->seq_name));
          }
          ++printed;
        }
        gt_bofprintf(buffered_output,"\n");
      }
    } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
    // Clean
    gt_output_map_attributes_delete(output_map_attributes);
  }

  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

#define GT_MAPSET_OPERATIONS "union,intersection,difference,compare,join,merge-map,display-compact"

void gt_filter_parse_operation(char* const string_operation) {
  if (gt_streq(string_operation,"INTERSECCTION") || gt_streq(string_operation,"Intersection") || gt_streq(string_operation,"intersection")) {
    parameters.operation = GT_MAP_SET_INTERSECTION;
  } else if (gt_streq(string_operation,"UNION") || gt_streq(string_operation,"Union") || gt_streq(string_operation,"union")) {
    parameters.operation = GT_MAP_SET_UNION;
  } else if (gt_streq(string_operation,"DIFFERENCE") || gt_streq(string_operation,"Difference") || gt_streq(string_operation,"difference")) {
    parameters.operation = GT_MAP_SET_DIFFERENCE;
  } else if (gt_streq(string_operation,"COMPARE") || gt_streq(string_operation,"Compare") || gt_streq(string_operation,"compare")) {
    parameters.operation = GT_MAP_SET_COMPARE;
  } else if (gt_streq(string_operation,"JOIN") || gt_streq(string_operation,"Join") || gt_streq(string_operation,"join")) {
    parameters.operation = GT_MAP_SET_JOIN;
  } else if (gt_streq(string_operation,"MERGE-MAP") || gt_streq(string_operation,"Merge-map") || gt_streq(string_operation,"merge-map")) {
    parameters.operation = GT_MERGE_MAP;
  } else if (gt_streq(string_operation,"DISPLAY-COMPACT") || gt_streq(string_operation,"Display-compact") || gt_streq(string_operation,"display-compact")) {
    parameters.operation = GT_DISPLAY_COMPACT_MAP;
  } else {
    if (string_operation[0]=='I' || string_operation[0]=='i') {
      fprintf(stderr,"\tAssuming 'Intersection' ...\n");
      parameters.operation = GT_MAP_SET_INTERSECTION;
    } else if (string_operation[0]=='U' || string_operation[0]=='u') {
      fprintf(stderr,"\tAssuming 'Union' ...\n");
      parameters.operation = GT_MAP_SET_UNION;
    } else if (string_operation[0]=='D' || string_operation[0]=='d') {
      fprintf(stderr,"\tAssuming 'Difference' ...\n");
      parameters.operation = GT_MAP_SET_DIFFERENCE;
    } else if (string_operation[0]=='C' || string_operation[0]=='c') {
      fprintf(stderr,"\tAssuming 'Compare' ...\n");
      parameters.operation = GT_MAP_SET_COMPARE;
    } else if (string_operation[0]=='P' || string_operation[0]=='p') {
      fprintf(stderr,"\tAssuming 'Join' ...\n");
      parameters.operation = GT_MAP_SET_JOIN;
    } else if (string_operation[0]=='M' || string_operation[0]=='m') {
      fprintf(stderr,"\tAssuming 'Merge-map' ...\n");
      parameters.operation = GT_MERGE_MAP;
    } else {
      gt_fatal_error_msg("Unknown operation '%s' in {"GT_MAPSET_OPERATIONS"}",string_operation);
    }
  }
}

void parse_arguments(int argc,char** argv) {
  struct option* gt_mapset_getopt = gt_options_adaptor_getopt(gt_mapset_options);
  gt_string* const gt_mapset_short_getopt = gt_options_adaptor_getopt_short(gt_mapset_options);
  int option, option_index;
  while (true) {
    // Get option & Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_mapset_short_getopt),gt_mapset_getopt,&option_index))==-1) break;
    // c=getopt_long(argc,argv,"i:o:psht:v",long_options,&option_index);
    switch (option) {
    /* Operations */
    case 'C':
      gt_filter_parse_operation(optarg);
      break;
    /* I/O */
    case 300:
      parameters.name_input_file_1 = optarg;
      break;
    case 301:
      parameters.name_input_file_2 = optarg;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 302:
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    /* Compare Function */
    case 's': // files-with-same-reads
      parameters.files_contain_same_reads = true;
      break;
    case 400: // eq-th
      parameters.eq_threshold = atof(optarg);
      break;
    case 401: // strict
      parameters.strict = true;
      break;
    /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'h':
      fprintf(stderr, "USE: ./gt.mapset [OPERATION] [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_mapset_options,gt_mapset_groups,false,false);
      exit(1);
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  // Check parameters
  if (parameters.operation==GT_MAP_SET_UNKNOWN) {
    gt_fatal_error_msg("Please specify operation {"GT_MAPSET_OPERATIONS"}");
  }
  if (parameters.operation!=GT_DISPLAY_COMPACT_MAP && !parameters.name_input_file_1) {
    gt_fatal_error_msg("Input file 1 required (--i1)\n");
  }
  // Free
  gt_string_delete(gt_mapset_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // Do it !
  if (parameters.operation==GT_MERGE_MAP) {
    gt_mapset_perform_merge_map();
  } else if (parameters.operation==GT_DISPLAY_COMPACT_MAP) {
    gt_mapset_display_compact_map();
  } else if (parameters.operation==GT_MAP_SET_INTERSECTION ||
      parameters.operation==GT_MAP_SET_UNION ||
      parameters.operation==GT_MAP_SET_DIFFERENCE) {
    gt_mapset_perform_set_operations();
  } else {
    gt_mapset_perform_cmp_operations();
  }

  return 0;
}


