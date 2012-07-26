/*
 * PROJECT: GEM-Tools library
 * FILE: gem-tools-examples.c
 * DATE: 01/06/2012
 * DESCRIPTION: Basic tool to perform simple map file conversions, filters and treatment in general
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE true

typedef struct {
  char *name_input_file;
  char *name_output_file;
  uint64_t num_threads;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .num_threads=1
};

void usage() {
  fprintf(stderr, "USE: ./gem-tools-examples <command> -i input -o output \n"
                  "      Commands::\n"
                  "        map-2-fastq\n"
                  "        parse-map\n"
                  "        map-2-map\n"
                  "        parse-sam\n"
                  "      Options::\n"
                  "        --input     <File>\n"
                  "        --output    <File>\n"
                  "        --help|h    \n");
}

void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:T:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'T':
      parameters.num_threads = atol(optarg);
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
}


void gt_example_map_2_fastq() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_template(buffered_input,template))) {
    if (error_code==GT_BMI_FAIL) continue;
    register const uint64_t num_blocks = gt_template_get_num_blocks(template);
    register uint64_t i;
    for (i=0;i<num_blocks;++i) {
      register gt_alignment* alignment = gt_template_get_block(template,i);
      printf("@%s\n%s\n+\n%s\n",gt_template_get_tag(template),
          gt_alignment_get_read(alignment),gt_alignment_get_qualities(alignment));
    }
  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
}

/*
 * Displays Template's content
 */
void gt_example_display_template(gt_template* template) {
  // Iterate over template's maps
  //   Eg. Single end {map1,map2,map3,...}
  //   Eg. Paired Alignment {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
  GT_TEMPLATE_ITERATE(template,map_array) {
    register const uint64_t template_num_blocks = gt_template_get_num_blocks(template);
    printf(">> %s  => Is %s (contains %"PRIu64" blocks)\n", gt_template_get_tag(template),
      (template_num_blocks==1)?"Single-end":(template_num_blocks==2?"Paired-end":"Weird"),
      template_num_blocks);
    GT_MAP_ARRAY_ITERATE(map_array,map,end_position) {
      gt_alignment* alignment = gt_template_get_block(template,end_position);
      // As maps can contain more than one block (due to split maps) we iterate over all of them
      printf("\t BEGIN_MAPS_BLOCK [TotalDistance=%"PRIu64"] { ", gt_map_get_global_distance(map));
      GT_MAP_BLOCKS_ITERATE(map,map_block) {
        printf("\n\t\t%s\t",gt_map_get_seq_name(map_block));
        /// IMPORTANT NOTE: Positions are base-1 (Genomic coordinates)
        printf("InitPos=%"PRIu64"\t",gt_map_get_position(map_block));
        printf("EndPos=%"PRIu64"\t",gt_map_get_position(map_block)+gt_map_get_length(map_block));
        printf("Len=%"PRIu64"\t",gt_map_get_length(map_block));
        printf("Strand=%c\t",gt_map_get_direction(map_block)==FORWARD?'F':'R');
        printf("Dist=%"PRIu64"\t",gt_map_get_distance(map_block));
        printf("LevDist=%"PRIu64"\t",gt_map_get_levenshtein_distance(map_block));

        printf("Misms{ ");
        GT_MISMS_ITERATE(map_block,misms) {
          /// IMPORTANT NOTE: Mismatch' Positions are base-0 (like strings in C)
          if (gt_misms_get_type(misms)==MISMS) {
            printf("(M,%"PRIu64",%c)[%c] ",gt_misms_get_position(misms),gt_misms_get_base(misms),
                gt_alignment_get_qualities(alignment)[gt_misms_get_position(misms)]);
          } else {
            printf("(%c,%"PRIu64",%"PRIu64")[%c] ",gt_misms_get_type(misms)==INS?'I':'D',
                gt_misms_get_position(misms),gt_misms_get_size(misms),
                gt_alignment_get_qualities(alignment)[gt_misms_get_position(misms)]);
          }
        }
        printf("}");

      }
      printf("\n\t } END_MAP_BLOCKS\n");
    }
  }
}

/*
 * Single thread MAP file parsing
 */
void gt_example_map_parsing() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_template(buffered_input,template))) {
    if (error_code==GT_BMI_FAIL) continue;

    // Print template's content (separately...)
    gt_example_display_template(template);

  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
}

/*
 * Single thread MAP file parsing and dump to a MAP file
 *   POST: input.map == output.map
 */
void gt_example_map_dump_map() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Open output file
  gt_buffered_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_buffered_output_stream_new(stdout,SORTED_FILE) :
      gt_buffered_output_file_new(parameters.name_output_file,SORTED_FILE);

  /*** BEGIN PPAR REGION */

  // Create I/O buffers
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_output_buffer* output_buffer = gt_buffered_output_file_request_buffer(output_file);
  // Read all templates
  gt_template* template = gt_template_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_template__sync_output(buffered_input,template,&output_buffer))) {
    if (error_code==GT_BMI_FAIL) continue;

    // Print MAP template
    gt_output_map_bprint_template(output_buffer,template,GT_ALL,true);

    // NOTE:: In case the parsing is not synch
    //    if (gt_buffered_input_file_eob(map_input_buffer)) {
    //      output_buffer = gt_buffered_output_file_dump_buffer(output_file,output_buffer);
    //    }
  }
  // Close I/O buffers
  gt_buffered_input_file_close(buffered_input);
  gt_buffered_output_file_release_buffer(output_file,output_buffer);

  /*** END PPAR REGION */

  // Close files
  gt_buffered_output_file_close(output_file);
  gt_input_file_close(input_file);
}

void gt_example_map_dump_map_parallel() {
//  // Parallel working threads
//  #pragma omp parallel num_threads(parameters.num_threads)
//  {
//
//  }
}

void gt_example_sam_parsing() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
//  while ((error_code=gt_buffered_sam_input_get_template(input_file,template))) {
//    if (error_code==GT_BMI_FAIL) continue;
//
//    // Print template's content (separately...)
//    gt_example_display_template(template);
//
//  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
}

int main(int argc,char** argv) {
  // Parsing command-line options
  register char* const example_name = argv[1];
  if (argc==1) { usage(); exit(0); }
  parse_arguments(argc-1, argv+1);

  //
  // Examples
  //

  // MAP to FASTQ conversion
  if (strcmp(example_name,"map-2-fastq")==0) {
    gt_example_map_2_fastq();
    return 0;
  }

  // Single thread MAP file parsing
  if (strcmp(example_name,"parse-map")==0) {
    gt_example_map_parsing();
    return 0;
  }

  // Single thread MAP file parsing and dump to a MAP file
  if (strcmp(example_name,"map-2-map")==0) {
    gt_example_map_dump_map();
    return 0;
  }

  // Single thread SAM file parsing
  if (strcmp(example_name,"parse-sam")==0) {
    gt_example_sam_parsing();
    return 0;
  }

  gt_error_msg("Incorrect test name provided");
  usage();
  return -1;
}


