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
  bool paired;
  bool verbose;
  uint64_t num_threads;
} gt_gtfcount_args;

gt_gtfcount_args parameters = {
    .input_file=NULL,
    .output_file=NULL,
    .annotation=NULL,
    .paired=false,
    .verbose=false,
    .num_threads=1
};

GT_INLINE void gt_gtfcount_merge_counts_(gt_shash* const source, gt_shash* const target){
  GT_SHASH_BEGIN_ITERATE(source, key, value, uint64_t){
    if(!gt_shash_is_contained(target, key)){
      uint64_t* v = gt_malloc_uint64();
      *v = *value;
      gt_shash_insert(target, key, v, uint64_t);
    }else{
      uint64_t* v = gt_shash_get(target,key,uint64_t);
      *v += (*value);
    }

  }GT_SHASH_END_ITERATE;
}



uint64_t gt_gtfcount_read(gt_gtf* const gtf, gt_shash* const gene_counts, gt_shash* const type_counts) {
  // Open file IN/OUT
  gt_input_file* input_file = NULL;
  uint64_t i = 0;
  if(parameters.input_file != NULL){
    input_file = gt_input_file_open(parameters.input_file,false);
  }else{
    input_file = gt_input_stream_open(stdin);
  }
  // create maps for the threads
  gt_shash** gene_counts_list = gt_calloc(parameters.num_threads, gt_shash*, true);
  gt_shash** type_counts_list = gt_calloc(parameters.num_threads, gt_shash*, true);
  for(i=0; i<parameters.num_threads; i++){
    gene_counts_list[i] = gt_shash_new();
    type_counts_list[i] = gt_shash_new();
  }
  uint64_t* read_counts = malloc(parameters.num_threads*sizeof(uint64_t));
  for(i=0; i<parameters.num_threads; i++){
    read_counts[i] = 0;
  }
  // Parallel reading+process
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_status error_code;
    gt_template* template = gt_template_new();
    gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired);

    // local maps
    uint64_t tid = omp_get_thread_num();
    gt_shash* l_gene_counts = gene_counts_list[tid];
    gt_shash* l_type_counts = type_counts_list[tid];

    while ((error_code = gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attr))) {
      if (error_code != GT_IMP_OK) {
        gt_fatal_error_msg("Fatal error parsing file \n");
      }

      if (gt_template_get_num_blocks(template)==1){
        GT_TEMPLATE_REDUCTION(template,alignment);
        if(gt_alignment_is_mapped(alignment) && gt_alignment_get_num_maps(alignment) == 1) read_counts[tid]++;
        gt_gtf_count_alignment(gtf, alignment, l_type_counts, l_gene_counts);
      } else {
        if (!gt_template_is_mapped(template)) {
          GT_TEMPLATE_REDUCE_BOTH_ENDS(template,alignment_end1,alignment_end2);

          if(gt_alignment_is_mapped(alignment_end1)&& gt_alignment_get_num_maps(alignment_end1) == 1) read_counts[tid]++;
          gt_gtf_count_alignment(gtf, alignment_end1, l_type_counts, l_gene_counts);
          if(gt_alignment_is_mapped(alignment_end2)&& gt_alignment_get_num_maps(alignment_end2) == 1) read_counts[tid]++;
          gt_gtf_count_alignment(gtf, alignment_end2, l_type_counts, l_gene_counts);
        } else {
          if(gt_template_is_mapped(template) && gt_template_get_num_mmaps(template) == 1) read_counts[tid] += 2;
          gt_gtf_count_template(gtf, template, l_type_counts, l_gene_counts);
        }
      }
    }
    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }
  uint64_t total_counts =0;
  // merge the count tables and delete them
  for(i=0; i<parameters.num_threads; i++){
    gt_gtfcount_merge_counts_(gene_counts_list[i], gene_counts);
    gt_gtfcount_merge_counts_(type_counts_list[i], type_counts);
    gt_shash_delete(gene_counts_list[i], true);
    gt_shash_delete(type_counts_list[i], true);
    total_counts += read_counts[i];
  }
  // Clean
  gt_input_file_close(input_file);
  return total_counts;
}


void parse_arguments(int argc,char** argv) {
  struct option* gt_gtfcount_getopt = gt_options_adaptor_getopt(gt_gtfcount_options);
  gt_string* const gt_gtfcount_short_getopt = gt_options_adaptor_getopt_short(gt_gtfcount_options);

  int option, option_index;
  while (true) {
    // Get option & Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_gtfcount_short_getopt),gt_gtfcount_getopt,&option_index))==-1) break;
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
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
      parameters.num_threads = atol(optarg);
      break;
    case 'h':
      fprintf(stderr, "USE: gt.gtfcount [OPERATION] [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_gtfcount_options,gt_gtfcount_groups,false,false);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_gtfcount_options,gt_gtfcount_groups,true,false);
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
  gt_string_delete(gt_gtfcount_short_getopt);
}

GT_INLINE void gt_gtfcount_warn(const char* const msg){
  if(parameters.verbose){
    fprintf(stderr, "%s", msg);
  }
}

GT_INLINE uint64_t gt_gtfcount_get_count_(gt_shash* const table, char* const element){
  if(!gt_shash_is_contained(table, element)){
    return 0;
  }
  uint64_t* v = gt_shash_get(table,element,uint64_t);
  return *v;
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();
  parse_arguments(argc,argv);

  // read gtf file
  gt_gtfcount_warn("Reading GTF...");
  gt_gtf* const gtf = gt_gtf_read_from_file(parameters.annotation, parameters.num_threads);
  gt_gtfcount_warn("Done\n");

  // local maps
  gt_shash* gene_counts = gt_shash_new();
  gt_shash* type_counts = gt_shash_new();
  gt_gtfcount_warn("Counting...");
  uint64_t total_counts = gt_gtfcount_read(gtf, gene_counts, type_counts);
  gt_gtfcount_warn("Done\n");

  printf("Annotation type counts\n");
  printf("----------------------\n");
  GT_SHASH_BEGIN_ITERATE(type_counts, key, e, uint64_t){
    printf("  %15s: %"PRIu64" (%.2f%%)\n", key, *e, (((float)*e/(float)total_counts) * 100.0));
  }GT_SHASH_END_ITERATE
  printf("----------------------\n");
  printf("Total counts: %"PRIu64"\n", total_counts);

  FILE* output = stdout;
  GT_SHASH_BEGIN_ITERATE(gene_counts, key, e, uint64_t){
    fprintf(output, "%s\t%"PRIu64"\n", key, *e );
  }GT_SHASH_END_ITERATE


  gt_shash_delete(gene_counts, true);
  gt_shash_delete(type_counts, true);
  return 0;
}


