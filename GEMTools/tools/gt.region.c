/*
 * PROJECT: GEM-Tools library
 * FILE: gt.gtfcount.c
 * DATE: 10/07/2013
 * AUTHOR(S): Thasso Griebel <thasso.griebel@gmail.com>
 * DESCRIPTION: Annotation map a file against a reference annotation
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

typedef struct {
  char *input_file;
  char *output_file;
  char *annotation;
  gt_vector *gene_ids;
  gt_vector *ref_ids;
  bool paired;
  uint64_t num_threads;
} gt_region_args;


gt_region_args parameters = {
    .input_file=NULL,
    .output_file=NULL,
    .annotation=NULL,
    .gene_ids=NULL,
    .ref_ids=NULL,
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
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_vector* hits = gt_vector_new(16, sizeof(gt_gtf_entry*));
    gt_output_map_attributes* const output_map_attributes = gt_output_map_attributes_new();
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired,buffered_output,template){
      bool printed = false;
      if (gtf != NULL){
        gt_vector_clear(hits);
        gt_gtf_search_template(gtf, hits, template);
        GT_VECTOR_ITERATE(hits, v, c, gt_gtf_entry*){
          gt_gtf_entry* e = *v;
          if(parameters.gene_ids != NULL && e->gene_id != NULL){
            GT_VECTOR_ITERATE(parameters.gene_ids, gene_id, pos, gt_string*) {
              if(gt_string_equals(e->gene_id, *gene_id)){
                //gt_output_map_bofprint_gem_template(buffered_output, template, output_map_attributes);
                gt_output_map_fprint_gem_template(stdout, template, output_map_attributes);
                printed = true;
                break;
              }
            }
          }
          if (printed){
            break; 
          }
        }
      }
      if (printed){
        continue;
      }
      if(parameters.ref_ids != NULL){
        GT_TEMPLATE_ITERATE_ALIGNMENT(template, ali){
          GT_ALIGNMENT_ITERATE(ali, map){
            GT_VECTOR_ITERATE(parameters.ref_ids, ref_id, pos, gt_string*) {
              if(gt_string_equals(gt_map_get_string_seq_name(map), *ref_id)){
                //gt_output_map_bofprint_gem_template(buffered_output, template, output_map_attributes);
                gt_output_map_fprint_gem_template(stdout, template, output_map_attributes);
                printed = true;
                break;
              }
            }
            if (printed){break;}
          }
          if (printed){break;}
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

void gt_region_get_argument_gene_id(char* const gene_ids) {
  // Allocate vector
  parameters.gene_ids = gt_vector_new(20,sizeof(gt_string*));
  // Add all the valid map Ids (sequence names)
  char *opt;
  opt = strtok(gene_ids, ",");
  while (opt!=NULL) {
    // Get id
    gt_string* gene_id = gt_string_new(0);
    gt_string_set_string(gene_id,opt);
    // Add to the vector
    gt_vector_insert(parameters.gene_ids,gene_id,gt_string*);
    // Next
    opt = strtok(NULL,","); // Reload
  }
}

void gt_region_get_argument_ref_id(char* const ref_ids) {
  // Allocate vector
  parameters.ref_ids = gt_vector_new(20,sizeof(gt_string*));
  // Add all the valid map Ids (sequence names)
  char *opt;
  opt = strtok(ref_ids, ",");
  while (opt!=NULL) {
    // Get id
    gt_string* ref_id = gt_string_new(0);
    gt_string_set_string(ref_id,opt);
    // Add to the vector
    gt_vector_insert(parameters.ref_ids,ref_id,gt_string*);
    // Next
    opt = strtok(NULL,","); // Reload
  }
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
      gt_region_get_argument_gene_id(optarg);
      break;
    case 'r':
      gt_region_get_argument_ref_id(optarg);
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
  if (parameters.annotation==NULL && parameters.gene_ids != NULL) {
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
  gt_gtf* gtf = NULL; 
  if(parameters.annotation != NULL){
    gtf = gt_gtf_read_from_file(parameters.annotation, parameters.num_threads);
  }
  gt_region_read(gtf);
  return 0;
}


