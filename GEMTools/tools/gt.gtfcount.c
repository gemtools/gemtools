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
  char *gene_counts_file;
  char *annotation;
  bool shell;
  bool paired;
  bool weighted_counts;
  bool unique_only;
  bool single_hit_only;
  bool verbose;
  uint64_t num_threads;
} gt_gtfcount_args;

typedef struct {
  uint64_t single_genes;
  uint64_t multi_genes;
  uint64_t no_genes;
  uint64_t single_reads;
} gt_gtfcount_pair_counts;


gt_gtfcount_args parameters = {
    .input_file=NULL,
    .output_file=NULL,
    .gene_counts_file=NULL,
    .annotation=NULL,
    .paired=false,
    .unique_only=true,
    .weighted_counts=false,
    .single_hit_only=true,
    .shell=false,
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

GT_INLINE void gt_gtfcount_merge_counts_weighted_(gt_shash* const source, gt_shash* const target){
  GT_SHASH_BEGIN_ITERATE(source, key, value, double){
    if(!gt_shash_is_contained(target, key)){
      double* v = malloc(sizeof(double*));
      *v = *value;
      gt_shash_insert(target, key, v, double);
    }else{
      double* v = gt_shash_get(target,key,double);
      *v += (*value);
    }

  }GT_SHASH_END_ITERATE;
}
GT_INLINE void gt_gtfcount_count_alignment(gt_gtf* gtf, gt_alignment* alignment, uint64_t* single_ends, uint64_t* read_counts,
                                           gt_shash* l_type_counts, gt_shash* l_gene_counts, gt_shash* pattern_counts,
                                           gt_shash* private_gene_counts){
  // increase read counts for unique hits
  uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  // count the alignment
  // in case we use weighted counts, the weight is set to 1 for single end reads and to 0.5 for paired end reads
  // in case no weighting should be applied, the count is set to < 0
  double weight = parameters.weighted_counts ? (parameters.paired ? 1.0 : 1.0) : -1.0;

  // clear the private counter hash
  gt_shash_clear(private_gene_counts, true);
  uint64_t hits = gt_gtf_count_alignment(gtf, alignment, num_maps == 1 ? pattern_counts : NULL, private_gene_counts);
  // now we have the full counts for this alignment and we have to add weighted counts to the l_gene_counts
  if(hits > 0 && ((!parameters.unique_only || num_maps == 1) && (!parameters.single_hit_only || hits == 1))){
    *single_ends += 1;
    *read_counts += 1;
    GT_SHASH_BEGIN_ITERATE(private_gene_counts, key, e, double){
      if(weight < 0.0){
        // unweighted counts
        gt_gtf_count_(l_gene_counts, key);
      }else{
        // weighted count
        double v = ((*e)/(double)num_maps) * weight;
        gt_gtf_count_weight_(l_gene_counts, key, v);
      }
    }GT_SHASH_END_ITERATE;
  }
  if(num_maps == 1 && hits == 1){
    // add type counts for unique reads with single gene hits
    GT_SHASH_BEGIN_KEY_ITERATE(private_gene_counts, key){
      gt_gtf_entry* gene = gt_gtf_get_gene_by_id(gtf, key);
      if(gene != NULL && gene->gene_type != NULL){
        gt_gtf_count_(l_type_counts, gene->gene_type->buffer);
      }
    }GT_SHASH_END_ITERATE;
  }
}

/**
 * This call does the stats count and the gene counts and returns the total number of reads that were taken into account
 * for the stats counts (NOT for the read counts, that depdends on the weighting scheme)
 */
uint64_t gt_gtfcount_read(gt_gtf* const gtf,
                          gt_shash* const gene_counts,
                          gt_shash* const type_counts,
                          gt_shash* const single_patterns_counts,
                          gt_shash* const pair_patterns_counts,
                          gt_gtfcount_pair_counts* pair_counts) {
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
  gt_shash** single_patterns_list = gt_calloc(parameters.num_threads, gt_shash*, true);
  gt_shash** pair_patterns_list = gt_calloc(parameters.num_threads, gt_shash*, true);
  uint64_t* singel_gene_pairs = malloc(parameters.num_threads * sizeof(uint64_t));
  uint64_t* multi_gene_pairs = malloc(parameters.num_threads * sizeof(uint64_t));
  uint64_t* no_gene_pairs = malloc(parameters.num_threads * sizeof(uint64_t));
  uint64_t* single_ends = malloc(parameters.num_threads * sizeof(uint64_t));

  for(i=0; i<parameters.num_threads; i++){
    gene_counts_list[i] = gt_shash_new();
    type_counts_list[i] = gt_shash_new();
    single_patterns_list[i] = gt_shash_new();
    pair_patterns_list[i] = gt_shash_new();
    singel_gene_pairs[i] = 0;
    multi_gene_pairs[i] = 0;
    no_gene_pairs[i] = 0;
    single_ends[i] = 0;
  }
  // general read counts array per thread
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
    gt_shash* l_single_patterns = single_patterns_list[tid];
    gt_shash* l_pair_patterns = pair_patterns_list[tid];

    gt_shash* private_gene_counts = gt_shash_new();

    while ((error_code = gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attr))) {
      if (error_code != GT_IMP_OK) {
        gt_fatal_error_msg("Fatal error parsing file \n");
      }
      // clear the private counter hash
      gt_shash_clear(private_gene_counts, true);

      if (gt_template_get_num_blocks(template)==1){
        // single end alignments
        GT_TEMPLATE_REDUCTION(template,alignment);
        gt_gtfcount_count_alignment(gtf, alignment, &single_ends[tid], &read_counts[tid],
                                    l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts);
      } else {
        if (!gt_template_is_mapped(template)) {
          // paired-end reads with unpaired alignments
          GT_TEMPLATE_REDUCE_BOTH_ENDS(template,alignment_end1,alignment_end2);
          gt_gtfcount_count_alignment(gtf, alignment_end1, &single_ends[tid], &read_counts[tid],
                                      l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts);
          gt_gtfcount_count_alignment(gtf, alignment_end2, &single_ends[tid], &read_counts[tid],
                                      l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts);
        } else {
          // paired-end alignments
          // increase read counts for unique hits
          uint64_t num_maps = gt_template_get_num_mmaps(template);
          // count the alignment
          // in case we use weighted counts, the weight is set to 1 for single end reads and to 0.5 for paired end reads
          // in case no weighting should be applied, the count is set to < 0
          double weight = parameters.weighted_counts ? 1.0 : -1.0;
          uint64_t hits = gt_gtf_count_template(gtf, template, num_maps == 1 ? l_pair_patterns : NULL, private_gene_counts);
          // now we have the full counts for this alignment and we have to add weighted counts to the l_gene_counts
          if(hits > 0 && ((!parameters.unique_only || num_maps == 1) && (!parameters.single_hit_only || hits == 1))){
            read_counts[tid] += 1;
            if(hits == 1){
              singel_gene_pairs[tid]++;
            }else if(hits > 1){
              multi_gene_pairs[tid]++;
            }
            GT_SHASH_BEGIN_ITERATE(private_gene_counts, key, e, double){
              if(weight < 0.0){
                // unweighted counts
                gt_gtf_count_(l_gene_counts, key);
              }else{
                // weighted count
                double v = ((*e)/(double)num_maps) * weight;
                gt_gtf_count_weight_(l_gene_counts, key, v);
              }
            }GT_SHASH_END_ITERATE;
          }else{
            no_gene_pairs[tid]++;
          }
          if(num_maps == 1 && hits == 1){
            // add type counts for unique reads with single gene hits
            GT_SHASH_BEGIN_KEY_ITERATE(private_gene_counts, key){
              gt_gtf_entry* gene = gt_gtf_get_gene_by_id(gtf, key);
              if(gene != NULL && gene->gene_type != NULL){
                gt_gtf_count_(l_type_counts, gene->gene_type->buffer);
              }
            }GT_SHASH_END_ITERATE;
          }
        }
      }
    }
    // Clean
    gt_shash_delete(private_gene_counts, true);
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }
  uint64_t total_counts =0;
  // merge the count tables and delete them
  for(i=0; i<parameters.num_threads; i++){
    if(parameters.weighted_counts){
      gt_gtfcount_merge_counts_weighted_(gene_counts_list[i], gene_counts);
    }else{
      gt_gtfcount_merge_counts_(gene_counts_list[i], gene_counts);
    }
    gt_gtfcount_merge_counts_(type_counts_list[i], type_counts);
    gt_gtfcount_merge_counts_(single_patterns_list[i], single_patterns_counts);
    gt_gtfcount_merge_counts_(pair_patterns_list[i], pair_patterns_counts);
    gt_shash_delete(gene_counts_list[i], true);
    gt_shash_delete(type_counts_list[i], true);
    gt_shash_delete(pair_patterns_list[i], true);
    gt_shash_delete(single_patterns_list[i], true);
    total_counts += read_counts[i];

    pair_counts->multi_genes += multi_gene_pairs[i];
    pair_counts->single_genes += singel_gene_pairs[i];
    pair_counts->no_genes += no_gene_pairs[i];
    pair_counts->single_reads += single_ends[i];
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
    case 'g':
      parameters.gene_counts_file = optarg;
      break;
    case 'a':
      parameters.annotation = optarg;
      break;
    case 'p':
      parameters.paired = true;
      break;
    /* Counts */
    case 'w':
      parameters.weighted_counts = true;
      break;
    case 'm':
      parameters.unique_only = false;
      break;
    case 's':
      parameters.single_hit_only = false;
      break;
    /* Misc */
    case 500:
      parameters.shell = true;
      break;
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

/***
 * SHELL and shell parser utilities
 */
GT_INLINE char* gt_gtfcount_shell_parse_ref(char** line){
  char* ref = *line;
  GT_READ_UNTIL(line, **line==':');
  if(GT_IS_EOL(line))return NULL;
  **line = EOS;
  GT_NEXT_CHAR(line);
  return ref;
}
GT_INLINE uint64_t gt_gtfcount_shell_parse_start(char** line){
  char* ref = *line;
  uint64_t n = 0;
  GT_READ_UNTIL(line, **line=='-' || **line=='\n');
  **line = EOS;
  n = atol(ref);
  GT_NEXT_CHAR(line);
  return n;
}
GT_INLINE uint64_t gt_gtfcount_shell_parse_end(char** line){
  if(**line == '\n') return 0;
  char* ref = *line;
  uint64_t n = 0;
  GT_READ_UNTIL(line, **line==' ' || **line=='\n');
  **line = EOS;
  n = atol(ref);
  GT_NEXT_CHAR(line);
  return n;
}
GT_INLINE char* gt_gtfcount_shell_parse_type(char** line){
  if(**line == '\n') return 0;
  while(**line == ' ') GT_NEXT_CHAR(line);
  char* ref = *line;
  GT_READ_UNTIL(line, **line=='\n');
  **line = EOS;
  GT_NEXT_CHAR(line);
  if(strlen(ref)==0) return NULL;
  return ref;
}
GT_INLINE void gt_gtfcount_run_shell(gt_gtf* const gtf){
  fprintf(stdout, "Search the annotation with queries like : <ref>:<start>[-<end>] [type]\n");
  fprintf(stdout, ">");
  gt_vector* hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  size_t buf_size = 1024;
  ssize_t read;
  char* line = malloc(buf_size * sizeof(char));
  uint64_t start, end=0 ;
  char* type = NULL;
  while((read = getline(&line, &buf_size, stdin)) != -1){
    // parse the line
    char* ref = gt_gtfcount_shell_parse_ref(&line);
    if(ref==NULL){
      fprintf(stdout, "Unable to parse reference name.\n");
      fprintf(stdout, ">");
      continue;
    }
    start = gt_gtfcount_shell_parse_start(&line);
    end = gt_gtfcount_shell_parse_end(&line);
    if(end != 0 && start > end){
      fprintf(stdout, "start > end not allowed!\n");
      fprintf(stdout, ">");
      continue;
    }
    type = gt_gtfcount_shell_parse_type(&line);


    if(end == 0) end = start;
    uint64_t num_results = gt_gtf_search(gtf,hits, ref, start, end, true);
    if(num_results == 0){
      fprintf(stdout, "Nothing found :(\n");
    }else{
      GT_VECTOR_ITERATE(hits, v, c, gt_gtf_entry*){
        gt_gtf_entry* e = *v;
        if(type != NULL && (e->type == NULL || strcmp(type, e->type->buffer) != 0)) continue;
        gt_gtf_print_entry_(e, NULL);
      }
    }
    fprintf(stdout, ">");
  }

}
int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();
  parse_arguments(argc,argv);

  // read gtf file
  gt_gtfcount_warn("Reading GTF...");
  gt_gtf* const gtf = gt_gtf_read_from_file(parameters.annotation, parameters.num_threads);
  gt_gtfcount_warn("Done\n");

  // run the shell
  if(parameters.shell){
    gt_gtfcount_run_shell(gtf);
    exit(0);
  }

  // local counting maps for types and genes counts
  gt_shash* gene_counts = gt_shash_new(); // gene counts
  gt_shash* type_counts = gt_shash_new(); // type (exon/intron/utr/rRNA/...)
  gt_shash* pair_pattern_counts = gt_shash_new(); // patterns for pairs (exon/intron/unknown/na)
  gt_shash* single_pattern_counts = gt_shash_new(); // patterns for single (exon/intron/unknown/na)

  // general counting stats struct
  gt_gtfcount_pair_counts pair_counts = {
      .single_genes = 0,
      .multi_genes = 0,
      .no_genes= 0,
      .single_reads =0
  };

  /// MAIN CALL TO COUNTING
  gt_gtfcount_warn("Counting...");
  uint64_t total_counts = gt_gtfcount_read(gtf, gene_counts, type_counts, single_pattern_counts, pair_pattern_counts, &pair_counts);
  gt_gtfcount_warn("Done\n");

  // write output for stats
  FILE* output = stdout;
  if(parameters.output_file != NULL){
    output = fopen(parameters.output_file, "w");
    if(output == NULL){
      gt_perror();
      exit(1);
    }
  }


  /*
   * Print type counts
   */
  uint64_t total_types = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(type_counts, e, uint64_t){
      total_types += *e;
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");
  fprintf(output, "Type counts (%"PRIu64" (%.2f%%))\n", total_types, (((float)total_types/(float)total_counts) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");
  GT_SHASH_BEGIN_ITERATE(type_counts, key, e, uint64_t){
    fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", key, *e, (((float)*e/(float)total_types) * 100.0));
  }GT_SHASH_END_ITERATE


  uint64_t total_type_counts = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(pair_pattern_counts, e, uint64_t){
      total_type_counts += *e;
  }GT_SHASH_END_ITERATE;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(single_pattern_counts, e, uint64_t){
      total_type_counts += *e;
  }GT_SHASH_END_ITERATE;
  /*
   * Print type counts for hits to single genes and hits to multiple genes
   */
  uint64_t type_counts_single_total = 0;
  uint64_t type_counts_multi_total = 0;
  GT_SHASH_BEGIN_ITERATE(pair_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") == NULL){
      type_counts_single_total += *e;
    }else{
      type_counts_multi_total += *e;
    }
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");

  fprintf(output, "PE Annotation type counts for single gene hits (Single: %"PRIu64" (%.2f%%))\n", type_counts_single_total, (((float)type_counts_single_total/total_type_counts) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");
  GT_SHASH_BEGIN_ITERATE(pair_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") == NULL){
      fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", key, *e, (((float)*e/(float)type_counts_single_total) * 100.0));
    }
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");
  fprintf(output, "PE Annotation type counts for multi gene hits (Multi: %"PRIu64" (%.2f%%))\n", type_counts_multi_total, (((float)type_counts_multi_total/total_type_counts) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");
  GT_SHASH_BEGIN_ITERATE(pair_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") != NULL){
      fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", key, *e, (((float)*e/(float)type_counts_multi_total) * 100.0));
    }
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");


  /*
   * Print type counts for hits to single genes and hits to multiple genes
   */
  type_counts_single_total = 0;
  type_counts_multi_total = 0;
  GT_SHASH_BEGIN_ITERATE(single_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") == NULL){
      type_counts_single_total += *e;
    }else{
      type_counts_multi_total += *e;
    }
  }GT_SHASH_END_ITERATE

  fprintf(output, "SE Annotation type counts for single gene hits (Single: %"PRIu64" (%.2f%%))\n", type_counts_single_total, (((float)type_counts_single_total/total_type_counts) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");
  GT_SHASH_BEGIN_ITERATE(single_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") == NULL){
      fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", key, *e, (((float)*e/(float)type_counts_single_total) * 100.0));
    }
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");
  fprintf(output, "SE Annotation type counts for multi gene hits (Multi: %"PRIu64" (%.2f%%))\n", type_counts_multi_total, (((float)type_counts_multi_total/total_type_counts) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");
  GT_SHASH_BEGIN_ITERATE(single_pattern_counts, key, e, uint64_t){
    if(strstr(key, "_mg") != NULL){
      fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", key, *e, (((float)*e/(float)type_counts_multi_total) * 100.0));
    }
  }GT_SHASH_END_ITERATE
  fprintf(output, "-----------------------------------------------------------------------\n");

  /*
   * Print Pair patterns
   */
  uint64_t paired_total = pair_counts.multi_genes + pair_counts.single_genes + pair_counts.no_genes + pair_counts.single_reads;
  fprintf(output, "Unique PE-reads Gene-Matches (pairs: %"PRIu64" singles: %"PRIu64" total: %"PRIu64")\n", paired_total-pair_counts.single_reads, pair_counts.single_reads, paired_total);
  fprintf(output, "-----------------------------------------------------------------------\n");
  fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", "Single end reads", pair_counts.single_reads, (((float)pair_counts.single_reads/(float)paired_total) * 100.0));
  fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", "Pair not mapped to gene", pair_counts.no_genes, (((float)pair_counts.no_genes/(float)paired_total) * 100.0));
  fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", "Pair mapped to single gene", pair_counts.single_genes, (((float)pair_counts.single_genes/(float)paired_total) * 100.0));
  fprintf(output, "  %40s: %"PRIu64" (%.5f%%)\n", "Pair mapped to multiple genes", pair_counts.multi_genes, (((float)pair_counts.multi_genes/(float)paired_total) * 100.0));
  fprintf(output, "-----------------------------------------------------------------------\n");


  if(parameters.output_file != NULL){
    fclose(output);
  }

  // print gene counts
  if(parameters.gene_counts_file != NULL){
    output = fopen(parameters.gene_counts_file, "w");
    if(output == NULL){
      gt_perror();
      exit(1);
    }
    // print gene table in case we use weighted counts, print as floats
    // otherwise as ints
    if(!parameters.weighted_counts){
      GT_SHASH_BEGIN_ITERATE(gene_counts, key, e, uint64_t){
          fprintf(output, "%s\t%"PRIu64"\n", key, *e );
      }GT_SHASH_END_ITERATE
    }else{
      GT_SHASH_BEGIN_ITERATE(gene_counts, key, e, double){
          fprintf(output, "%s\t%.4f\n", key, *e );
      }GT_SHASH_END_ITERATE
    }
    fclose(output);
  }


  gt_shash_delete(gene_counts, true);
  gt_shash_delete(type_counts, true);
  gt_shash_delete(pair_pattern_counts, true);
  gt_shash_delete(single_pattern_counts, true);
  return 0;
}


