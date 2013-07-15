/*
 * PROJECT: GEM-Tools library
 * FILE: gt.gtfcount.c
 * DATE: 10/07/2013
 * AUTHOR(S): Thasso Griebel <thasso.griebel@gmail.com>
 * DESCRIPTION: Annotation map a file against a reference annotation
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

typedef struct {
  char *input_file;
  char *output_file;
  char *annotation;
  char *gene_id;
  bool paired;
  uint64_t num_threads;
} gt_region_args;


gt_region_args parameters = {
    .input_file=NULL,
    .output_file=NULL,
    .annotation=NULL,
    .gene_id=NULL,
    .paired=false,
    .num_threads=1
};


GT_INLINE void gt_region_read(gt_gtf* const gtf) {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.input_file,false);
  gt_output_file* output_file = (parameters.output_file==NULL) ?
            gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.output_file,SORTED_FILE);

  // Parallel I/O
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_vector* hits = gt_vector_new(16, sizeof(gt_gtf_entry*));
    gt_output_map_attributes* const output_map_attributes = gt_output_map_attributes_new();
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired,buffered_output,template){
      gt_vector_clear(hits);
      gt_gtf_search_template(gtf, hits, template);
      GT_VECTOR_ITERATE(hits, v, c, gt_gtf_entry*){
        gt_gtf_entry* e = *v;
        if(parameters.gene_id != NULL && e->gene_id != NULL){
          if(strcmp(e->gene_id->buffer, parameters.gene_id) == 0){
            //gt_output_map_bofprint_gem_template(buffered_output, template, output_map_attributes);
            gt_output_map_fprint_gem_template(stdout, template, output_map_attributes);
            break;
          }
        }
      }
    }GT_END_READING_WRITING_LOOP(input_file,output_file,template);
    // cleanup per threads
    gt_vector_delete(hits);
    gt_output_map_attributes_delete(output_map_attributes);
  }
  // Clean global
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}


void parse_arguments(int argc,char** argv) {
  struct option* gt_region_getopt = gt_options_adaptor_getopt(gt_region_options);
  gt_string* const gt_region_short_getopt = gt_options_adaptor_getopt_short(gt_region_options);

  int option, option_index;
  while (true) {
    // Get option & Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_region_short_getopt),gt_region_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i':
      parameters.input_file = optarg;
      break;
    case 'o':
      parameters.output_file = optarg;
      break;
    case 'a':
      parameters.annotation = optarg;
      break;
    case 'p':
      parameters.paired = true;
      break;
    /* Misc */
    case 'g':
      parameters.gene_id = optarg;
      break;
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'h':
      fprintf(stderr, "USE: gt.gtfcount [OPERATION] [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_region_options,gt_region_groups,false,false);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_region_options,gt_region_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  // Check parameters
  if (parameters.annotation==NULL) {
    gt_fatal_error_msg("Please specify a reference annotation");
  }
  // Free
  gt_string_delete(gt_region_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();
  parse_arguments(argc,argv);

  // read gtf file
  gt_gtf* const gtf = gt_gtf_read_from_file(parameters.annotation, parameters.num_threads);
  gt_region_read(gtf);
  return 0;
}


