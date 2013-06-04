#include <getopt.h>
#include <omp.h>
#include <math.h>

#include "gem_tools.h"

#define is_good_mapq(qual) (qual >= 78 && (qual <= 90 || qual >= 119))
#define get_mapq(score) ((int)floor((sqrt(score)/256.0)*255))

typedef struct {
  /* I/O */
  char* name_input_file;
  char* name_output_file;
  char* name_reference_file;
  bool mmap_input;
  bool paired_end;
  /* Output */
  bool count_alignments;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
} gt_quality_args;

gt_quality_args params = {
  /* I/O */
  .name_input_file=NULL,
  .name_output_file=NULL,
  .mmap_input=false,
  .paired_end=false,
  /* Output */
  .count_alignments=false,
  /* Misc */
  .num_threads=1,
  .verbose=false,
};

/*
 * CORE functions
 */

void gt_template_filter(gt_template* template_dst,gt_template* template_src) {
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src, alignment_src) {
    GT_TEMPLATE_REDUCTION(template_dst, alignment_dst);
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      if (is_good_mapq(get_mapq(map->gt_score))) {
        gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
      }
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;

  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template_src,mmap,mmap_attr) {
    if (is_good_mapq(get_mapq(mmap_attr->gt_score))) {
      register gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
      gt_template_add_mmap_array(template_dst,mmap_copy,mmap_attr);
      free(mmap_copy);
    }
  }
}

void gt_get_by_score() {

  // Open file
  gt_input_file* input_file = (params.name_input_file == NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(params.name_input_file,params.mmap_input);
  gt_output_file* out_file = (params.name_output_file==NULL) ?
      gt_output_stream_new(stdout, SORTED_FILE) : gt_output_file_new(params.name_output_file,SORTED_FILE);

  gt_output_map_attributes attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();

  // Parallel reading+process
  #pragma omp parallel num_threads(params.num_threads)
  {

    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(out_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    gt_status error_code;
    gt_template *template = gt_template_new();
    gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(params.paired_end);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",params.name_input_file);
      }
      register const bool is_mapped = gt_template_is_mapped(template);

      if (is_mapped) {
        //Filter template
        gt_template *template_filtered = gt_template_dup(template,false,false);
        gt_template_filter(template_filtered, template);
        gt_template_delete(template);
        template = template_filtered;

        //Recompute counters
        gt_template_recalculate_counters(template);
      }

      // Print
      if (!params.count_alignments) {
        if (gt_output_map_bofprint_template(buffered_output,template,&attributes)) {
          gt_error_msg("Fatal error outputting read '"PRIgts"'(InputLine:%"PRIu64")\n",
              PRIgts_content(gt_template_get_string_tag(template)),buffered_input->current_line_num-1);
        }
      }

    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }


  gt_input_file_close(input_file);
  gt_output_file_close(out_file);
}

void usage() {
  fprintf(stderr, "USE: ./gt.quality [ARGS]...\n"
                  "         [I/O]\n"
                  "           --input|-i [FILE]\n"
                  "           --output|-o [FILE]\n"
                  "           --mmap-input\n"
                  "           --paired-end|p\n"
                 // "         [Output]\n"
                 // "           --count-alignments|-c\n"
                  "         [Misc]\n"
                  "           --threads|t\n"
                  "           --verbose|v\n"
                  "           --help|h\n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* I/O */
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "mmap-input", no_argument, 0, 1 },
    { "paired-end", no_argument, 0, 'p' },
    /* Output */
    { "count-alignments", no_argument, 0, 'c' },
    /* Misc */
    { "threads", required_argument, 0, 't' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:r:pct:hv",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /* I/O */
    case 'i':
      params.name_input_file = optarg;
      break;
    case 'o':
      params.name_output_file = optarg;
      break;
    case 1:
      params.mmap_input = true;
      break;
    case 'p':
      params.paired_end = true;
      break;
    /* Output */
    case 'c':
      params.count_alignments = true;
      break;
    /* Misc */
    case 't':
      params.num_threads = atol(optarg);
      break;
    case 'v':
      params.verbose = true;
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

  parse_arguments(argc, argv);

  gt_get_by_score();

  return 0;
}
