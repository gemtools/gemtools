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

#define GT_GTFCOUNT_HEAD(output,title) fprintf(output, "---%*s%*s---\n",(int)(30+strlen(title)/2) ,title, (int)(30-strlen(title)/2), "");\
  uint64_t __title_length=0;for(__title_length=(30+strlen(title)/2)+(30-strlen(title)/2)+6; __title_length>0; __title_length--)fprintf(output, "-");fprintf(output, "\n")
#define GT_GTFCOUNT_PRINT(output,space,name,num) fprintf(output, "  %"space"s: %10"PRIu64"\n", name, num)
#define GT_GTFCOUNT_PRINT_P(output,space,name,num,tot) fprintf(output, "  %"space"s: %10"PRIu64" (%7.3f%%)\n", name, num, tot?(((double)num/(double)tot )*100.0):0.0)
#define GT_GTFCOUNT_PRINT_R(output,space,name,num) fprintf(output, "  %"space"s: %10.3f\n", name, num)

typedef struct {
  char *input_file;
  char *name_output_file;
  char *gene_counts_file;
  char *annotation;
  FILE *output_file;
  FILE *output_file_json;
  bool shell;
  bool paired;
  bool coverage_profiles;
  bool weighted_counts;
  bool unique_only;
  bool single_hit_only;
  bool single_end_counts;
  bool verbose;
  bool print_json;
  bool print_both;
  bool count_bases;
  float exon_overlap;
  uint64_t num_threads;
} gt_gtfcount_args;

typedef struct {
  // general stats from the dataset
  uint64_t num_templates; // total number of templates
  uint64_t num_pe_reads; // number properly paired paired-end reads (counts 2 per template)
  uint64_t num_se_reads; // number of single end reads + number of unpaired PE reds
  uint64_t num_pe_unmapped; // number of unmapped paired-end reads (counts 2 per template)
  uint64_t num_se_unmapped; // number of unmapped se-templates + number of unpaired unmapped PE reads
  uint64_t num_pe_mappings; // total number of mappings PE-template mappings (counts 1 per template)
  uint64_t num_se_mappings; // total number of mappings from SE reads or unpaired PE-reads (counts per unmapped read)
  uint64_t num_unique_pe_maps; // # uniquely mapped pe-templates (counts 2 per template)
  uint64_t num_unique_se_maps; // # uniquely mapped se reads (counts 1 per read)
  uint64_t num_junctions; // # total number of junctions for single mapping reads
  uint64_t num_annotated_junctions; // # total number of junctions for single mapping reads that are found in the annotation

  uint64_t counted_reads; // total number of reads counted for gene counts (raw reads, counts 2 for a pe-template)
  uint64_t considered_pe_mappings; // total number of pe-mappings taken into account for counting (count 2 per template)
  uint64_t considered_se_mappings; // total number of se-mappings taken into account for counting
  uint64_t counted_pe_single_gene; // # PE-templates matched to single gene (2 per template)
  uint64_t counted_pe_multi_gene; // # PE-templates matched to multiple genes (2 per template)
  uint64_t counted_se_single_gene; // # SE-templates or unpaired PE-reads matched to single gene
  uint64_t counted_se_multi_gene; // # SE-templates or unpaired PE-reads matched to multiple genes
  uint64_t* single_transcript_coverage; // # coverage profiles for exons
  uint64_t* gene_body_coverage; // # coverage profiles for gene body
} gt_gtfcount_count_stats;


gt_gtfcount_args parameters = {
    .input_file=NULL,
    .name_output_file=NULL,
    .output_file=NULL,
    .output_file_json=NULL,
    .gene_counts_file=NULL,
    .coverage_profiles=false,
    .annotation=NULL,
    .paired=false,
    .unique_only=true,
    .weighted_counts=false,
    .single_hit_only=true,
    .single_end_counts=false,
    .shell=false,
    .exon_overlap=0.0,
    .print_both=false,
    .print_json=false,
    .count_bases=false,
    .verbose=false,
    .num_threads=1
};

GT_INLINE gt_gtfcount_count_stats* gt_gtfcount_count_stats_new(void){
  gt_gtfcount_count_stats* s = malloc(sizeof(gt_gtfcount_count_stats));
  s->num_templates = 0;
  s->num_pe_reads = 0;
  s->num_se_reads = 0;
  s->num_pe_unmapped = 0;
  s->num_se_unmapped = 0;
  s->num_pe_mappings = 0;
  s->num_se_mappings = 0;
  s->num_unique_pe_maps = 0;
  s->num_unique_se_maps = 0;
  s->num_junctions = 0;
  s->num_annotated_junctions = 0;
  s->counted_reads = 0;
  s->counted_pe_single_gene = 0;
  s->counted_pe_multi_gene = 0;
  s->counted_se_single_gene = 0;
  s->counted_se_multi_gene = 0;
  s->considered_pe_mappings = 0;
  s->considered_se_mappings = 0;
  s->single_transcript_coverage = NULL;
  s->gene_body_coverage = NULL;
  return s;
}

GT_INLINE void gt_gtfcount_count_stats_delete(gt_gtfcount_count_stats* stats){
  if(stats->single_transcript_coverage != NULL){
    free(stats->single_transcript_coverage);
  }
  if(stats->gene_body_coverage != NULL){
    free(stats->gene_body_coverage);
  }
  free(stats);
}

GT_INLINE void gt_gtfcount_count_stats_merge(gt_gtfcount_count_stats* target, gt_gtfcount_count_stats* source){
  target->num_templates += source->num_templates;
  target->num_pe_reads += source->num_pe_reads;
  target->num_se_reads += source->num_se_reads;
  target->num_pe_unmapped += source->num_pe_unmapped;
  target->num_se_unmapped += source->num_se_unmapped;
  target->num_pe_mappings += source->num_pe_mappings;
  target->num_se_mappings += source->num_se_mappings;
  target->num_unique_pe_maps += source->num_unique_pe_maps;
  target->num_unique_se_maps += source->num_unique_se_maps;
  target->num_junctions += source->num_junctions;
  target->num_annotated_junctions += source->num_annotated_junctions;
  target->counted_reads += source->counted_reads;
  target->counted_pe_single_gene += source->counted_pe_single_gene;
  target->counted_pe_multi_gene += source->counted_pe_multi_gene;
  target->counted_se_single_gene += source->counted_se_single_gene;
  target->counted_se_multi_gene += source->counted_se_multi_gene;
  target->considered_pe_mappings += source->considered_pe_mappings;
  target->considered_se_mappings += source->considered_se_mappings;
}


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
GT_INLINE void gt_gtfcount_count_alignment(gt_gtf* gtf, gt_alignment* alignment, gt_gtfcount_count_stats* stats,
                                           gt_shash* l_type_counts, gt_shash* l_gene_counts, gt_shash* pattern_counts,
                                           gt_shash* private_gene_counts, gt_gtf_count_parms* params){
  // increase read counts for unique hits
  uint64_t num_maps = gt_alignment_get_num_maps(alignment);

  stats->num_se_mappings += num_maps;
  stats->num_se_reads++;
  switch(num_maps){
    case 0: stats->num_se_unmapped++; break;
    case 1: stats->num_unique_se_maps++; break;
  }

  // count the alignment
  // in case we use weighted counts, the weight is set to 1 for single end reads and to 0.5 for paired end reads
  // in case no weighting should be applied, the count is set to < 0
  double weight = parameters.weighted_counts ? (parameters.paired ? 1.0 : 1.0) : -1.0;

  // clear the private counter hash
  gt_shash_clear(private_gene_counts, true);
  uint64_t hits = gt_gtf_count_alignment(gtf, alignment, num_maps == 1 ? pattern_counts : NULL, private_gene_counts, params);


  // now we have the full counts for this alignment and we have to add weighted counts to the l_gene_counts
  if(hits > 0 && ((!parameters.unique_only || hits == 1)) && (!parameters.paired || parameters.single_end_counts)){
    stats->counted_reads++;
    stats->considered_se_mappings += num_maps;
    switch(hits){
      case 1: stats->counted_se_single_gene++; break;
      default: stats->counted_se_multi_gene++; break;
    }
    GT_SHASH_BEGIN_KEY_ITERATE(private_gene_counts, key){
      if(weight < 0.0){
        // unweighted counts        
        if(params->count_bases){
          gt_gtf_count_custom_(l_gene_counts, key, gt_alignment_get_read_length(alignment));
        }else{
          gt_gtf_count_(l_gene_counts, key);
        }
      }else{
        //double v = ((*e)/(double)hits) * weight;
        if(params->count_bases){
          uint64_t l = gt_alignment_get_read_length(alignment);
          double v = ((double)l/(double)(hits * l)) * weight;
          gt_gtf_count_weight_(l_gene_counts, key, v);
        }else{
          double v = (1.0/(double)hits) * weight;
          gt_gtf_count_weight_(l_gene_counts, key, v);
        }
      }
    }GT_SHASH_END_ITERATE;
  }

  if(num_maps == 1 && hits == 1){
    // add type counts for unique reads with single gene hits
    GT_SHASH_BEGIN_KEY_ITERATE(private_gene_counts, key){
      gt_gtf_entry* gene = gt_gtf_get_gene_by_id(gtf, key);
      if(gene != NULL && gene->gene_type != NULL){
        if(params->count_bases){
          gt_gtf_count_custom_(l_type_counts, gene->gene_type->buffer, gt_alignment_get_read_length(alignment));
        }else{
          gt_gtf_count_(l_type_counts, gene->gene_type->buffer);
        }
      }
    }GT_SHASH_END_ITERATE;
  }
}

/**
 * This call does the stats count and the gene counts and returns the total number of reads that were taken into account
 * for the stats counts (NOT for the read counts, that depdends on the weighting scheme)
 */
GT_INLINE void gt_gtfcount_read(gt_gtf* const gtf,
                                gt_shash* const gene_counts,
                                gt_shash* const type_counts,
                                gt_shash* const single_patterns_counts,
                                gt_shash* const pair_patterns_counts,
                                gt_gtfcount_count_stats* pair_counts) {
  // Open file IN/OUT
  gt_input_file* input_file = NULL;
  uint64_t i = 0;
  uint64_t j = 0;
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
  gt_gtfcount_count_stats** stats_list =  gt_calloc(parameters.num_threads, gt_gtfcount_count_stats*, true);
  gt_gtf_count_parms** thread_params = gt_calloc(parameters.num_threads, gt_gtf_count_parms*, true);

  for(i=0; i<parameters.num_threads; i++){
    gene_counts_list[i] = gt_shash_new();
    type_counts_list[i] = gt_shash_new();
    single_patterns_list[i] = gt_shash_new();
    pair_patterns_list[i] = gt_shash_new();
    stats_list[i] = gt_gtfcount_count_stats_new();
  }

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_status error_code;
    gt_template* template = gt_template_new();
    gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired);

    // local maps
#ifdef HAVE_OPENMP
    uint64_t tid = omp_get_thread_num();
#else
	 uint64_t tid = 0;
#endif
    gt_shash* l_gene_counts = gene_counts_list[tid];
    gt_shash* l_type_counts = type_counts_list[tid];
    gt_shash* l_single_patterns = single_patterns_list[tid];
    gt_shash* l_pair_patterns = pair_patterns_list[tid];

    gt_gtf_count_parms* params = gt_gtf_count_params_new(parameters.coverage_profiles);
    params->count_bases = parameters.count_bases;
    thread_params[tid] = params;
    params->exon_overlap = parameters.exon_overlap;
    gt_shash* private_gene_counts = gt_shash_new();

    while ((error_code = gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attr))) {
      if (error_code != GT_IMP_OK) {
        gt_fatal_error_msg("Fatal error parsing file \n");
      }
      // clear the private counter hash
      gt_shash_clear(private_gene_counts, true);
      // count general stats
      stats_list[tid]->num_templates++;

      if (gt_template_get_num_blocks(template)==1){
        // single end alignments
        GT_TEMPLATE_REDUCTION(template,alignment);
        gt_gtfcount_count_alignment(gtf, alignment, stats_list[tid],
                                    l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts, params);
      } else {
        if (!gt_template_is_mapped(template) &&
            (gt_alignment_get_num_maps(gt_template_get_block(template, 0)) > 0
             || gt_alignment_get_num_maps(gt_template_get_block(template, 1)) > 0)) {
          // paired-end reads with unpaired alignments
          GT_TEMPLATE_REDUCE_BOTH_ENDS(template,alignment_end1,alignment_end2);
          gt_gtfcount_count_alignment(gtf, alignment_end1, stats_list[tid],
                                      l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts, params);
          gt_gtfcount_count_alignment(gtf, alignment_end2, stats_list[tid],
                                      l_type_counts, l_gene_counts, l_single_patterns, private_gene_counts, params);
        } else {
          // paired-end alignments
          gt_gtfcount_count_stats* stats = stats_list[tid];
          uint64_t num_maps = gt_template_get_num_mmaps(template);
          stats->num_pe_mappings += num_maps;
          stats->num_pe_reads += 2;
          switch(num_maps){
            case 0: stats->num_pe_unmapped += 2; break;
            case 1: stats->num_unique_pe_maps += 2; break;
          }

          // count the alignment
          // in case we use weighted counts, the weight is set to 1 for single end reads and to 0.5 for paired end reads
          // in case no weighting should be applied, the count is set to < 0
          double weight = parameters.weighted_counts ? 1.0 : -1.0;
          uint64_t hits = gt_gtf_count_template(gtf, template, num_maps == 1 ? l_pair_patterns : NULL, private_gene_counts, params);
          // now we have the full counts for this alignment and we have to add weighted counts to the l_gene_counts
          if(hits > 0 && ((!parameters.unique_only || hits == 1) )){
            stats->counted_reads +=2;
            stats->considered_pe_mappings += (num_maps * 2);
            switch(hits){
              case 1: stats->counted_pe_single_gene += 2; break;
              default: stats->counted_pe_multi_gene += 2; break;
            }
            GT_SHASH_BEGIN_ITERATE(private_gene_counts, key, e, double){
              if(weight < 0.0){
                uint64_t v = parameters.single_end_counts ? *e : 2;
                if(params->count_bases){
                  if(v == 2){
                    v = gt_template_get_total_length(template);
                  }else{
                    v = gt_template_get_total_length(template) / 2;
                  }
                }
                // unweighted counts
                gt_gtf_count_sum_(l_gene_counts, key, v);
              }else{
                uint64_t v = parameters.single_end_counts ? *e : 2;
                double diff = (double) hits;
                double vv = (double)v/diff;
                if(params->count_bases){
                  if(v == 2){
                    v = gt_template_get_total_length(template);
                  }else{
                    v = gt_template_get_total_length(template) / 2;
                  }
                  diff = gt_template_get_total_length(template) * hits;
                  vv = (double) v * ((double)v/diff);
                }
                // weighted count
                gt_gtf_count_weight_(l_gene_counts, key, vv);
              }
            }GT_SHASH_END_ITERATE;
          }
          if(num_maps == 1 && hits == 1){
            // add type counts for unique reads with single gene hits
            GT_SHASH_BEGIN_KEY_ITERATE(private_gene_counts, key){
              gt_gtf_entry* gene = gt_gtf_get_gene_by_id(gtf, key);
              if(gene != NULL && gene->gene_type != NULL){
                gt_gtf_count_sum_(l_type_counts, gene->gene_type->buffer, 2);
              }
            }GT_SHASH_END_ITERATE;
          }
        }
      }
    }

    stats_list[tid]->num_junctions += params->num_junctions;
    stats_list[tid]->num_annotated_junctions += params->num_annotated_junctions;
    // Clean

    gt_shash_delete(private_gene_counts, true);
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }


  // merge the count tables and delete them
  // and merge and delete coverage
  uint64_t* single_transcript_coverage = NULL;
  uint64_t* gene_body = NULL;
  if(parameters.coverage_profiles){
    single_transcript_coverage = GT_GTF_INIT_COVERAGE();
    gene_body = GT_GTF_INIT_COVERAGE();
  }
  for(i=0; i<parameters.num_threads; i++){
    if(parameters.weighted_counts){
      gt_gtfcount_merge_counts_weighted_(gene_counts_list[i], gene_counts);
    }else{
      gt_gtfcount_merge_counts_(gene_counts_list[i], gene_counts);
    }
    gt_gtfcount_merge_counts_(type_counts_list[i], type_counts);
    gt_gtfcount_merge_counts_(single_patterns_list[i], single_patterns_counts);
    gt_gtfcount_merge_counts_(pair_patterns_list[i], pair_patterns_counts);
    gt_gtfcount_count_stats_merge(pair_counts, stats_list[i]);

    gt_shash_delete(gene_counts_list[i], true);
    gt_shash_delete(type_counts_list[i], true);
    gt_shash_delete(pair_patterns_list[i], true);
    gt_shash_delete(single_patterns_list[i], true);
    if(parameters.coverage_profiles){
      for(j=0; j<GT_GTF_COVERAGE_LENGTH;j++){
        single_transcript_coverage[j] += thread_params[i]->single_transcript_coverage[j];
      }
      for(j=0; j<GT_GTF_COVERAGE_LENGTH;j++){
        gene_body[j] += thread_params[i]->gene_body_coverage[j];
      }
      gt_gtf_count_params_delete(thread_params[i]);
    }

    gt_gtfcount_count_stats_delete(stats_list[i]);
  }
  if(parameters.coverage_profiles){
    pair_counts->single_transcript_coverage = single_transcript_coverage;
    pair_counts->gene_body_coverage = gene_body;
  }

  // Clean
  gt_input_file_close(input_file);
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
      parameters.name_output_file = optarg;
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
    case 'e':
      parameters.exon_overlap = atof(optarg);
      break;
    case 'm':
      parameters.unique_only = false;
      break;
    case 's':
      parameters.single_end_counts = true;
      break;
    case 'f':
      if(strcmp("report", optarg) == 0){
        parameters.print_json = false;
      }else if(strcmp("json", optarg) == 0){
        parameters.print_json = true;
      }else if(strcmp("both", optarg) == 0){
        parameters.print_json = true;
        parameters.print_both = true;
      }else{
        gt_fatal_error_msg("Unknown format %s, supported formats are 'report' or 'json' or 'both'", optarg);
      }
      break;
    case 400:
      parameters.count_bases = true;
      break;
    /* Misc */
    case 500:
      parameters.shell = true;
      break;
    case 'c':
      parameters.coverage_profiles = true;
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
        gt_gtf_print_entry_(stdout,e, NULL);
      }
    }
    fprintf(stdout, ">");
  }

}

GT_INLINE void gt_gtfcount_print_general(FILE* output, gt_gtfcount_count_stats* stats){

  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  uint64_t total_unmapped_reads = stats->num_pe_unmapped + stats->num_se_unmapped;
  uint64_t total_mapped_reads = total_reads-total_unmapped_reads;
  GT_GTFCOUNT_HEAD(output, "General stats");
  fprintf(output, "General:\n");
  GT_GTFCOUNT_PRINT(output, "40", "Templates", stats->num_templates);
  GT_GTFCOUNT_PRINT(output, "40", "Reads", total_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Paired reads", stats->num_pe_reads, total_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Single reads", stats->num_se_reads, total_reads);
  fprintf(output, "General mapping:\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Mapped reads",(total_mapped_reads), total_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Unmapped reads",(total_unmapped_reads), total_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Mapped Paired reads",(stats->num_pe_reads-stats->num_pe_unmapped), stats->num_pe_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Unmapped Paired reads",(stats->num_pe_unmapped), stats->num_pe_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Mapped Single reads",(stats->num_se_reads-stats->num_se_unmapped), stats->num_se_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Unmapped Single reads",(stats->num_se_unmapped), stats->num_se_reads);
  fprintf(output, "Multimaps:\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Uniquely mapped reads",(stats->num_unique_pe_maps+stats->num_unique_se_maps), total_mapped_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Ambiguously mapped reads",(total_mapped_reads-(stats->num_unique_pe_maps+stats->num_unique_se_maps)), total_mapped_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Uniquely mapped Paired reads",(stats->num_unique_pe_maps), total_mapped_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Ambiguously mapped Paired reads",(stats->num_pe_reads-stats->num_unique_pe_maps-stats->num_pe_unmapped), total_mapped_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Uniquely mapped Single reads",(stats->num_unique_se_maps), total_mapped_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Ambiguously mapped Single reads",(stats->num_se_reads-stats->num_unique_se_maps-stats->num_se_unmapped), total_mapped_reads);
}

GT_INLINE JsonNode* gt_gtfcount_print_general_json(gt_gtfcount_count_stats* stats){
  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  uint64_t total_unmapped_reads = stats->num_pe_unmapped + stats->num_se_unmapped;
  uint64_t total_mapped_reads = total_reads-total_unmapped_reads;
  JsonNode* node = json_mkobject();
  json_append_member(node, "templates", json_mknumber(stats->num_templates));
  json_append_member(node, "reads", json_mknumber(total_reads));
  json_append_member(node, "paired_reads", json_mknumber(stats->num_pe_reads));
  json_append_member(node, "single_reads", json_mknumber(stats->num_se_reads));
  json_append_member(node, "mapped_reads", json_mknumber(total_mapped_reads));
  json_append_member(node, "unmapped_reads", json_mknumber(total_unmapped_reads));

  json_append_member(node, "mapped_paired_reads", json_mknumber(stats->num_pe_reads-stats->num_pe_unmapped));
  json_append_member(node, "unmapped_paired_reads", json_mknumber(stats->num_pe_unmapped));
  json_append_member(node, "mapped_single_reads", json_mknumber(stats->num_se_reads-stats->num_se_unmapped));
  json_append_member(node, "unmapped_single_reads", json_mknumber(stats->num_se_unmapped));

  json_append_member(node, "mapped_unique_reads", json_mknumber(stats->num_unique_pe_maps+stats->num_unique_se_maps));
  json_append_member(node, "mapped_ambigous_reads", json_mknumber(total_mapped_reads-(stats->num_unique_pe_maps+stats->num_unique_se_maps)));

  json_append_member(node, "mapped_unique_paired_reads", json_mknumber(stats->num_unique_pe_maps));
  json_append_member(node, "mapped_ambigous_paired_reads", json_mknumber(stats->num_pe_reads - stats->num_unique_pe_maps - stats->num_pe_unmapped));

  json_append_member(node, "mapped_unique_single_reads", json_mknumber(stats->num_unique_se_maps));
  json_append_member(node, "mapped_ambigous_single_reads", json_mknumber(stats->num_se_reads - stats->num_unique_se_maps - stats->num_se_unmapped));
  return node;
}

GT_INLINE JsonNode* gt_gtf_count_coverage_array_(uint64_t* data, uint64_t range){
  JsonNode* node = json_mkarray();
  uint64_t i = 0;
  for(i=0; i<GT_GTF_COVERAGE_BUCKETS;i++){
    json_append_element(node, json_mknumber(data[GT_GTF_COVERGAGE_GET_BUCKET(range, i)]));
  }
  return node;
}

GT_INLINE JsonNode* gt_gtf_count_print_coverage_json(gt_gtfcount_count_stats* stats){
  JsonNode* node = json_mkobject();
  JsonNode* single_transcripts = json_mkobject();
  JsonNode* gene_body = json_mkobject();
  // gene body
  json_append_member(gene_body, "all", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_ALL));
  json_append_member(gene_body, "150", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_150));
  json_append_member(gene_body, "250", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_250));
  json_append_member(gene_body, "500", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_500));
  json_append_member(gene_body, "1000", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_1000));
  json_append_member(gene_body, "2500", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_2500));
  json_append_member(gene_body, "5000", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_5000));
  json_append_member(gene_body, "7500", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_7500));
  json_append_member(gene_body, "10000", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_10000));
  json_append_member(gene_body, "15000", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_15000));
  json_append_member(gene_body, "20000", gt_gtf_count_coverage_array_(stats->gene_body_coverage, GT_GTF_COVERAGE_LENGTH_20000));
  // single transcript
  json_append_member(single_transcripts, "all", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_ALL));
  json_append_member(single_transcripts, "150", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_150));
  json_append_member(single_transcripts, "250", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_250));
  json_append_member(single_transcripts, "500", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_500));
  json_append_member(single_transcripts, "1000", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_1000));
  json_append_member(single_transcripts, "2500", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_2500));
  json_append_member(single_transcripts, "5000", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_5000));
  json_append_member(single_transcripts, "7500", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_7500));
  json_append_member(single_transcripts, "10000", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_10000));
  json_append_member(single_transcripts, "15000", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_15000));
  json_append_member(single_transcripts, "20000", gt_gtf_count_coverage_array_(stats->single_transcript_coverage, GT_GTF_COVERAGE_LENGTH_20000));

  json_append_member(node, "single_transcript", single_transcripts);
  json_append_member(node, "gene_body", gene_body);
  return node;
}


GT_INLINE void gt_gtfcount_print_count_stats(FILE* output, gt_gtfcount_count_stats* stats){
  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  GT_GTFCOUNT_HEAD(output, "Gene count stats");
  GT_GTFCOUNT_PRINT_P(output, "40", "Total counted reads", stats->counted_reads, total_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Total considered single mappings", stats->considered_se_mappings, ((stats->num_pe_mappings*2)+stats->num_se_mappings));
  GT_GTFCOUNT_PRINT_P(output, "40", "Total considered paired mappings", stats->considered_pe_mappings, ((stats->num_pe_mappings*2)+stats->num_se_mappings));
  GT_GTFCOUNT_PRINT_R(output, "40", "Mappings/Reads ratio", (((double)(stats->considered_se_mappings + stats->considered_pe_mappings))/(double)stats->counted_reads));
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Reads fall in single gene", (stats->counted_pe_single_gene+stats->counted_se_single_gene), stats->counted_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Reads fall in multiple genes", (stats->counted_pe_multi_gene+stats->counted_se_multi_gene), stats->counted_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Paired Reads fall in single gene", (stats->counted_pe_single_gene), stats->counted_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Single Reads fall in single gene", (stats->counted_se_single_gene), stats->counted_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Paired Reads fall in multiple gene", (stats->counted_pe_multi_gene), stats->counted_reads);
  GT_GTFCOUNT_PRINT_P(output, "40", "Single Reads fall in multiple gene", (stats->counted_se_multi_gene), stats->counted_reads);
  fprintf(output, "\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Annotated junction hits", (stats->num_annotated_junctions), stats->num_junctions);
  GT_GTFCOUNT_PRINT_P(output, "40", "Denovo junction hits", (stats->num_junctions-stats->num_annotated_junctions), stats->num_junctions);
}
GT_INLINE JsonNode* gt_gtfcount_print_count_stats_json(gt_gtfcount_count_stats* stats){
  JsonNode* node = json_mkobject();
  json_append_member(node, "total_counted_reads", json_mknumber(stats->counted_reads));
  json_append_member(node, "total_considered_single_reads", json_mknumber(stats->considered_se_mappings));
  json_append_member(node, "total_considered_paired_reads", json_mknumber(stats->considered_pe_mappings));
  json_append_member(node, "reads_single_gene", json_mknumber(stats->counted_pe_single_gene+stats->counted_se_single_gene));
  json_append_member(node, "reads_multiple_genes", json_mknumber(stats->counted_pe_multi_gene+stats->counted_se_multi_gene));
  json_append_member(node, "reads_single_gene_se", json_mknumber(stats->counted_se_single_gene));
  json_append_member(node, "reads_multiple_genes_se", json_mknumber(stats->counted_se_multi_gene));
  json_append_member(node, "reads_single_gene_pe", json_mknumber(stats->counted_pe_single_gene));
  json_append_member(node, "reads_multiple_genes_pe", json_mknumber(stats->counted_pe_multi_gene));
  json_append_member(node, "annotated_junction_hits", json_mknumber(stats->num_annotated_junctions));
  json_append_member(node, "denovo_junction_hits", json_mknumber(stats->num_junctions - stats->num_annotated_junctions));
  return node;
}

GT_INLINE void gt_gtfcount_print_type_counts(FILE* output, gt_gtfcount_count_stats* stats, gt_shash* type_counts){
  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  uint64_t total_types = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(type_counts, e, uint64_t){
      total_types += *e;
  }GT_SHASH_END_ITERATE

  GT_GTFCOUNT_HEAD(output, "Type Counts");
  fprintf(output, "The type counts are generated only from uniquely\n"
                  "mapping reads which hit a single gene.\n\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Total Reads used for type counts", total_types, total_reads);
  fprintf(output, "\nTypes:\n");
  GT_SHASH_BEGIN_ITERATE(type_counts, key, e, uint64_t){
    GT_GTFCOUNT_PRINT_P(output, "40", key, *e, total_types);
  }GT_SHASH_END_ITERATE
}
GT_INLINE void gt_gtfcount_print_pair_patterns(FILE* output, gt_gtfcount_count_stats* stats, gt_shash* patterns){
  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  uint64_t total_patterns = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(patterns, e, uint64_t){
      total_patterns += *e;
  }GT_SHASH_END_ITERATE

  GT_GTFCOUNT_HEAD(output, "Paired reads patterns");
  fprintf(output, "The paired reads patterns are generated only from the uniquely mapping\n"
                  "paired reads, counting 1 for a pair.\n"
                  "Junctions are indicated by '^' and pairs are split by '|'.\n\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Total Reads used for paired counts", (total_patterns*2), total_reads);
  fprintf(output, "\nPatterns:\n");
  GT_SHASH_BEGIN_ITERATE(patterns, key, e, uint64_t){
    GT_GTFCOUNT_PRINT_P(output, "40", key, *e, total_patterns);
  }GT_SHASH_END_ITERATE
}
GT_INLINE void gt_gtfcount_print_single_patterns(FILE* output, gt_gtfcount_count_stats* stats, gt_shash* patterns){
  uint64_t total_reads = stats->num_pe_reads + stats->num_se_reads;
  uint64_t total_patterns = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(patterns, e, uint64_t){
      total_patterns += *e;
  }GT_SHASH_END_ITERATE

  GT_GTFCOUNT_HEAD(output, "Single reads patterns");
  fprintf(output, "The single reads patterns are generated only from the uniquely mapping\n"
                  "single reads, counting 1 for each read.\n"
                  "Junctions are indicated by '^'.\n\n");
  GT_GTFCOUNT_PRINT_P(output, "40", "Total Reads used for single counts", total_patterns, total_reads);
  fprintf(output, "\nPatterns:\n");
  GT_SHASH_BEGIN_ITERATE(patterns, key, e, uint64_t){
    GT_GTFCOUNT_PRINT_P(output, "40", key, *e, total_patterns);
  }GT_SHASH_END_ITERATE
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
  gt_gtfcount_count_stats* counting_stats = gt_gtfcount_count_stats_new();

  /// MAIN CALL TO COUNTING
  gt_gtfcount_read(gtf, gene_counts, type_counts, single_pattern_counts, pair_pattern_counts, counting_stats);

  // init output paramters
  parameters.output_file = stdout;
  parameters.output_file_json = stderr;
  if(parameters.name_output_file != NULL){
    parameters.output_file = fopen(parameters.name_output_file, "w");
    gt_cond_fatal_error(parameters.output_file==NULL,FILE_OPEN,parameters.name_output_file);
  }
  if(parameters.print_json && !parameters.print_both){
    parameters.output_file_json = parameters.output_file;
  }

  if(parameters.print_both || !parameters.print_json){
    /**
     * Print human readable stats
     */
    fprintf(parameters.output_file, "\n");
    gt_gtfcount_print_general(parameters.output_file, counting_stats);
    fprintf(parameters.output_file, "\n");
    gt_gtfcount_print_count_stats(parameters.output_file, counting_stats);
    fprintf(parameters.output_file, "\n");
    gt_gtfcount_print_type_counts(parameters.output_file, counting_stats, type_counts);
    if(gt_shash_get_num_elements(pair_pattern_counts) > 0){
      fprintf(parameters.output_file, "\n");
      gt_gtfcount_print_pair_patterns(parameters.output_file, counting_stats, pair_pattern_counts);
    }
    if(gt_shash_get_num_elements(single_pattern_counts) > 0){
      fprintf(parameters.output_file, "\n");
      gt_gtfcount_print_single_patterns(parameters.output_file, counting_stats, single_pattern_counts);
    }
    fprintf(parameters.output_file, "\n");
  }
  if(parameters.print_both || parameters.print_json){
    // print JSON stats
    JsonNode* root = json_mkobject();
    json_append_member(root, "general", gt_gtfcount_print_general_json(counting_stats));
    json_append_member(root, "counts", gt_gtfcount_print_count_stats_json(counting_stats));
    json_append_member(root, "type_counts", gt_json_int_hash(type_counts));
    json_append_member(root, "pair_patterns", gt_json_int_hash(pair_pattern_counts));
    json_append_member(root, "single_patterns", gt_json_int_hash(single_pattern_counts));
    if(parameters.coverage_profiles){
      json_append_member(root, "coverage", gt_gtf_count_print_coverage_json(counting_stats));
    }

    fprintf(parameters.output_file_json, "%s\n", json_stringify(root, "  "));
    json_delete(root);
  }

  if(parameters.name_output_file != NULL){
    fclose(parameters.output_file);
  }



  // print gene counts
  if(parameters.gene_counts_file != NULL){
    FILE* output = fopen(parameters.gene_counts_file, "w");
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

  gt_gtfcount_count_stats_delete(counting_stats);
  gt_shash_delete(gene_counts, true);
  gt_shash_delete(type_counts, true);
  gt_shash_delete(pair_pattern_counts, true);
  gt_shash_delete(single_pattern_counts, true);
  return 0;
}


