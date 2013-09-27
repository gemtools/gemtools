/*
 * PROJECT: GEM-Tools library
 * FILE: gt.construct.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Constructor. Used to load examples/debug/testing/deploy/develop/...
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE false

typedef struct {
  char *name_input_file;
  char *name_output_file;
  uint64_t option;
  uint64_t number;
  uint64_t param1;
  uint64_t param2;
  uint64_t param3;
  uint64_t param4;
  uint64_t param5;
  uint64_t param6;
  uint64_t param7;
  uint64_t num_threads;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .option=0,
    .number=0,
    .num_threads=1
};

void gt_map_2_fastq() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_template(buffered_input,template,NULL))) {
    if (error_code==GT_BMI_FAIL) continue;
    const uint64_t num_blocks = gt_template_get_num_blocks(template);
    uint64_t i;
    for (i=0;i<num_blocks;++i) {
      gt_alignment* alignment = gt_template_get_block(template,i);
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
  const uint64_t template_num_blocks = gt_template_get_num_blocks(template);
  printf(">> %s  => Is %s (contains %"PRIu64" blocks)\n",
      gt_template_get_tag(template),
      (template_num_blocks==1)?"Single-end":(template_num_blocks==2?"Paired-end":"Weird"),
      template_num_blocks);
  GT_TEMPLATE_ITERATE(template,map_array) {
    GT_MMAP_ITERATE(map_array,map,end_position) {
      gt_alignment* alignment = gt_template_get_block(template,end_position);
      // As maps can contain more than one block (due to split maps) we iterate over all of them
      printf("\t BEGIN_MAPS_BLOCK [TotalDistance=%"PRIu64"] { ", gt_map_get_global_distance(map));
      GT_MAP_ITERATE(map,map_block) {
        printf("\n\t\t%s\t",gt_map_get_seq_name(map_block));
        /// IMPORTANT NOTE: Positions are base-1 (Genomic coordinates)
        printf("InitPos=%"PRIu64"\t",gt_map_get_position(map_block));
        printf("EndPos=%"PRIu64"\t",gt_map_get_position(map_block)+gt_map_get_length(map_block));
        printf("Len=%"PRIu64"\t",gt_map_get_length(map_block));
        printf("Strand=%c\t",gt_map_get_strand(map_block)==FORWARD?'F':'R');
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
void gt_example_map_parsing() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_map_parser_attributes map_parser_attributes = GT_MAP_PARSER_ATTR_DEFAULT(false);
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  gt_input_map_parser_attributes_set_skip_model(&map_parser_attributes,true);
  while ((error_code=gt_input_map_parser_get_template(buffered_input,template,&map_parser_attributes))) {
    if (error_code!=GT_IMP_OK) {
      gt_error_msg("Parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
      continue;
    }

    gt_example_display_template(template);
    // Output the same line (parsed + printed)
    gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
    gt_output_map_fprint_template(stdout,template,&output_map_attributes);
  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
}
void gt_remove_maps_with_n_or_more_mismatches() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_input_file_open("bug",false);
  gt_output_file* output_file = gt_output_stream_new(stdout, SORTED_FILE);
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(4)
#endif
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    gt_status error_code;
    gt_alignment *alignment_src = gt_alignment_new();
    gt_alignment *alignment_dst;
    gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(false);
    while ((error_code = gt_input_generic_parser_get_alignment(buffered_input,alignment_src,generic_parser_attr))) {
      if (error_code != GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file \n");
      }

      // Filter by number of mismatches
      alignment_dst = gt_alignment_copy(alignment_src,false);
      GT_ALIGNMENT_ITERATE(alignment_src,map) {
        // Count number of mismatches
        int num_misms = 0;
        GT_MISMS_ITERATE(map,misms) {
          if (misms->misms_type == MISMS)
            num_misms++;
        }
        // Add the map
        if (num_misms <= 3) {
          gt_alignment_insert_map(alignment_dst,gt_map_copy(map), true);
        }
      }

      // Print template
      gt_output_map_bofprint_alignment(buffered_output,alignment_dst,output_attributes);
      gt_alignment_delete(alignment_dst);
    }

    // Clean
    gt_alignment_delete(alignment_src);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
void gt_load_reference__dump_it() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Load reference
  gt_sequence_archive* sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
  if (gt_input_multifasta_parser_get_archive(input_file,sequence_archive)!=GT_IFP_OK) {
    gt_error_msg("Fatal error parsing reference\n");
  }

  // Dump summary of reference file
  if (gt_input_fasta_is_multifasta(input_file)) {
    fprintf(stderr,"File '%s' is MULTIFASTA\n",input_file->file_name);
  }
  gt_sequence_archive_iterator seq_arch_it;
  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
  gt_segmented_sequence* seq;
  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
    fprintf(stderr,"SEQUENCE '%s' [length=%"PRIu64"]\n",seq->seq_name->buffer,seq->sequence_total_length);
  }
  fprintf(stderr,"\n");

  // TEST 1

//  // Dump the content of the reference file
//  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
//  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
//    fprintf(stderr,">%s\n",seq->seq_name->buffer);
//    uint64_t i;
//    for (i=0; i<seq->sequence_total_length; ++i) {
//      fprintf(stderr,"%c",gt_segmented_sequence_get_char_at(seq,i));
//    }
//    if (seq->sequence_total_length) fprintf(stderr,"\n");
//  }

  // TEST 2

//  // Dump the content of the reference file (using Iterators)
//  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
//  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
//    fprintf(stderr,">%s\n",seq->seq_name->buffer);
//    gt_segmented_sequence_iterator sequence_iterator;
//    gt_segmented_sequence_new_iterator(seq,0,GT_ST_FORWARD,&sequence_iterator);
//    if (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
//      while (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
//        fprintf(stderr,"%c",gt_segmented_sequence_iterator_next(&sequence_iterator));
//      }
//      fprintf(stderr,"\n");
//    }
//  }

  // Freedom !!
  gt_sequence_archive_delete(sequence_archive);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
void gt_filter_fastq() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
    gt_output_fasta_attributes* output_attributes = gt_output_fasta_attributes_new();

    gt_status error_code;
    gt_dna_read* dna_read = gt_dna_read_new();
    //gt_output_fasta_attributes_set_format(output_attributes, gt_input_fasta_get_format(input_file));

    while ((error_code=gt_input_fasta_parser_get_read(buffered_input,dna_read))) {
      if (error_code!=GT_IFP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
      }
      gt_output_fasta_bofprint_dna_read(buffered_output,dna_read, output_attributes);
    }

    // Clean
    gt_dna_read_delete(dna_read);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
    gt_output_fasta_attributes_delete(output_attributes);
  }

  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
void gt_debug_mem_leak_parsing_map() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_alignment* alignment = gt_alignment_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_alignment(buffered_input,alignment,NULL))) {
    if (error_code==GT_BMI_FAIL) continue;
  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
}

void gt_mm_write_seq_mem(void* mem,uint64_t size) {
  uint64_t num=1, i;
  for (i=0;i<size;i+=8) {
    (*((uint64_t*)mem)) = num;
    num*=3;
    mem+=8;
  }
}
uint64_t gt_mm_read_seq_mem(void* mem,uint64_t size) {
  uint64_t sum = 0, i;
  for (i=0;i<size;i+=8) {
    sum = sum ^ (*((uint64_t*)mem));
    sum += i;
    mem+=8;
  }
  return sum;
}
void gt_mm_write_rand_mem(void* mem,uint64_t size,uint64_t num_rand_its) {
  const uint64_t num_words = (size/8)-1;
  uint64_t i, num=1;
  gt_srand();
  for (i=0;i<num_rand_its;++i) {
    const uint64_t pos = gt_rand(0,num_words);
    *(((uint64_t*)mem)+pos) = num;
    num*=3;
  }
}
uint64_t gt_mm_read_rand_mem(void* mem,uint64_t size,uint64_t num_rand_its) {
  const uint64_t num_words = (size/8)-1;
  uint64_t i, sum = 0;
  gt_srand();
  for (i=0;i<num_rand_its;++i) {
    const uint64_t pos = gt_rand(0,num_words);
    sum = sum ^ (*(((uint64_t*)mem)+pos));
    sum += i;
  }
  return sum;
}

void gt_mm_performance_test() {
  gt_handle_error_signals();
  gt_profiler_init();
  gt_mm_set_tmp_folder("/scratch_tmp/");
  /*
   * parameters.number = numBytes to allocate
   */
  const uint64_t num_rand_its = 100000000;
  uint64_t seq_check = 0, rand_check = 0;
  gt_mm* mm;
  switch (parameters.option) {
    case 0: /* gt_mm_bulk_malloc */
    case 1: /* gt_mm_bulk_mmalloc_temp */
    case 2: /* gt_mm_bulk_mmalloc */
      if (parameters.option==0) {
        fprintf(stderr,"[GTMalloc.profile]\n");
        // Allocate
        GT_START_TIMER(0);
        mm = gt_mm_bulk_malloc(parameters.number,false);
        GT_STOP_TIMER(0);
      } else if (parameters.option==1) {
        fprintf(stderr,"[GTMapMallocTMP.profile]\n");
        // Allocate
        GT_START_TIMER(0);
        mm = gt_mm_bulk_mmalloc_temp(parameters.number);
        GT_STOP_TIMER(0);
      } else { // parameters.option==2
        fprintf(stderr,"[GTMapMalloc.profile]\n");
        // Allocate
        GT_START_TIMER(0);
        mm = gt_mm_bulk_mmalloc(parameters.number,false);
        GT_STOP_TIMER(0);
      }
      // Read sequential
      GT_START_TIMER(1);
      gt_mm_write_seq_mem(gt_mm_get_base_mem(mm),parameters.number);
      GT_STOP_TIMER(1);
      // Read sequential
      GT_START_TIMER(2);
      seq_check = gt_mm_read_seq_mem(gt_mm_get_base_mem(mm),parameters.number);
      GT_STOP_TIMER(2);
      // Write random
      GT_START_TIMER(3);
      gt_mm_write_rand_mem(gt_mm_get_base_mem(mm),parameters.number,num_rand_its);
      GT_STOP_TIMER(3);
      // Read random
      GT_START_TIMER(4);
      rand_check = gt_mm_read_rand_mem(gt_mm_get_base_mem(mm),parameters.number,num_rand_its);
      GT_STOP_TIMER(4);
      // Report
      fprintf(stderr,"--> Time.alloc           %2.3f       \n",GT_GET_TIMER(0));
      fprintf(stderr,"--> Time.writeSeq.all    %2.3f       \n",GT_GET_TIMER(1));
      fprintf(stderr,"--> Time.readSeq.all     %2.3f  (%"PRIu64"X)\n",GT_GET_TIMER(2),seq_check);
      fprintf(stderr,"--> Time.writeRand.100M  %2.3f       \n",GT_GET_TIMER(3));
      fprintf(stderr,"--> Time.readRand.100M   %2.3f  (%"PRIu64"X)\n",GT_GET_TIMER(4),rand_check);
      // Free
      gt_mm_free(mm);
      break;
    case 3: /* gt_mm_bulk_mmap_file */
    case 4: /* gt_mm_bulk_load_file */
    case 5: /* gt_mm_bulk_mload_file */
    case 6: /* gt_mm_bulk_mmap_file_populate */
      if (parameters.option==3) {
        fprintf(stderr,"[GTMapMallocFILE.profile]\n");
        GT_START_TIMER(0);
        mm = gt_mm_bulk_mmap_file(parameters.name_input_file,GT_MM_READ_ONLY,false);
        GT_STOP_TIMER(0);
      } else if (parameters.option==4) {
        fprintf(stderr,"[GTLoadFILE.profile]\n");
        GT_START_TIMER(0);
        mm = gt_mm_bulk_load_file(parameters.name_input_file,parameters.param1);
        GT_STOP_TIMER(0);
      } else if (parameters.option==5) {
        fprintf(stderr,"[GTMapLoadFILE.profile]\n");
        GT_START_TIMER(0);
        mm = gt_mm_bulk_mload_file(parameters.name_input_file,parameters.param1);
        GT_STOP_TIMER(0);
      } else {
        fprintf(stderr,"[GTMapMallocFILE.POPULATE.profile]\n");
        GT_START_TIMER(0);
        mm = gt_mm_bulk_mmap_file(parameters.name_input_file,GT_MM_READ_ONLY,true);
        GT_STOP_TIMER(0);
      }
      // Read sequential
      GT_START_TIMER(2);
      seq_check = gt_mm_read_seq_mem(gt_mm_get_base_mem(mm),mm->allocated);
      GT_STOP_TIMER(2);
      // Read random
      GT_START_TIMER(4);
      rand_check = gt_mm_read_rand_mem(gt_mm_get_base_mem(mm),mm->allocated,num_rand_its);
      GT_STOP_TIMER(4);
      // Report
      fprintf(stderr,"--> Time.alloc           %2.3f       \n",GT_GET_TIMER(0));
      fprintf(stderr,"--> Time.readSeq.all     %2.3f  (%"PRIu64"X)\n",GT_GET_TIMER(2),seq_check);
      fprintf(stderr,"--> Time.readRand.10M    %2.3f  (%"PRIu64"X)\n",GT_GET_TIMER(4),rand_check);
      // Free
      gt_mm_free(mm);
      break;
    default:
      break;
  }
  gt_profiler_exit();
}

void usage() {
  fprintf(stderr, "USE: ./gem-tools-examples -i input -o output \n"
                  "      Options::\n"
                  "        --input|i     <File>\n"
                  "        --output|o    <File>\n"
                  "        --select|s    <Number>\n"
                  "        --number|n    <Number>\n"
                  "        --help|h\n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "select", required_argument, 0, 's' },
    { "number", required_argument, 0, 'n' },
    // parameters
    { "param1", required_argument, 0, '1' },
    { "param2", required_argument, 0, '2' },
    { "param3", required_argument, 0, '3' },
    { "param4", required_argument, 0, '4' },
    { "param5", required_argument, 0, '5' },
    { "param6", required_argument, 0, '6' },
    { "param7", required_argument, 0, '7' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:s:n:1:2:3:4:5:6:7:T:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 's':
      parameters.option = atol(optarg);
      break;
    case 'n':
     parameters.number = atol(optarg);
     break;
    /* Numbered Params */
    case '1':
     parameters.param1 = atol(optarg);
     break;
    case '2':
     parameters.param2 = atol(optarg);
     break;
    case '3':
     parameters.param3 = atol(optarg);
     break;
    case '4':
     parameters.param4 = atol(optarg);
     break;
    case '5':
     parameters.param5 = atol(optarg);
     break;
    case '6':
     parameters.param6 = atol(optarg);
     break;
    case '7':
     parameters.param7 = atol(optarg);
     break;
    case 'T':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
}

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  /*
   * Useful constructs
   */
  //gt_load_reference__dump_it();

  /*
   * Load it!
   */
  gt_mm_performance_test();

  return 0;
}


