/*
 * PROJECT: GEM-Tools library
 * FILE: gt.construct.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Constructor. Used to load examples/debug/testing/deploy/develop/...
 */

#include <getopt.h>
#include <omp.h>

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE false

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

/*
 * MAP to FASTQ conversion
 */
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
  register const uint64_t template_num_blocks = gt_template_get_num_blocks(template);
  printf(">> %s  => Is %s (contains %"PRIu64" blocks)\n",
      gt_template_get_tag(template),
      (template_num_blocks==1)?"Single-end":(template_num_blocks==2?"Paired-end":"Weird"),
      template_num_blocks);
  GT_TEMPLATE_ITERATE(template,map_array) {
    GT_MULTIMAP_ITERATE(map_array,map,end_position) {
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

/*
 * Single thread MAP file parsing
 */
void gt_example_map_parsing() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_map_parser_attr map_parser_attributes = GT_MAP_PARSER_ATTR_DEFAULT(false);
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_template* template = gt_template_new();
  gt_status error_code;
  gt_input_map_parser_attributes_set_skip_model(&map_parser_attributes,true);
  while ((error_code=gt_input_map_parser_get_template_g(buffered_input,template,&map_parser_attributes))) {
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

void gt_dummy_example() {
  register gt_status error_code = 0;
  register gt_alignment* alignment = gt_alignment_new();
  error_code+=gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20):::65472,"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20):::12224,"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20):::1984,"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20):::1984,"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20):::1984",alignment);

  fprintf(stdout, "[%d] the alignment has %" PRIu64 " maps\n", error_code, gt_alignment_get_num_maps(alignment));
}


void gt_example_map_string_parsing() {
  register gt_status error_code = 0;

  /*
   * Parsing single maps
   */
  register gt_map* map = gt_map_new();
  // Parse old map
  error_code+=gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map);
  // Parse new SE-map
  error_code+=gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map);
  // Parse new PE-map
  error_code+=gt_input_map_parse_map("chr15:-:102516634:66G9::chr15:+:102516358:66>1+10:::7936",map);

  /*
   * Parsing list of maps
   */
  register gt_vector* map_list = gt_vector_new(10,sizeof(gt_map*));
  // Parse multiple old split-map
  error_code+=gt_input_map_parse_map_list("[31;35]=chr16:R[2503415;2503411]~chr16:R2503271",map_list);
  error_code+=gt_input_map_parse_map_list("[30;34]=chr10:F74776624~chr10:F[74790025;74790029]",map_list);
  error_code+=gt_input_map_parse_map_list("[23-50]=chr1:F[188862944-188868041]~chr19:F53208292",map_list);
  error_code+=gt_input_map_parse_map_list("[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list);
  error_code+=gt_input_map_parse_map_list("[26]=chr7:R1203797~chr7:R1203108",map_list);
  error_code+=gt_input_map_parse_map_list(
      "chrM:F6598<+1>16@0/0,[23]=chr6:R31322884~chr6:R31237276,"
      "[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list);
  // Parse multiple old SE-map
  error_code+=gt_input_map_parse_map_list("chr1:F8926499@0/0,chr12:R7027116G39A42@77/2",map_list);
  // Parse multiple new SE-map
  error_code+=gt_input_map_parse_map_list(
      "chrX:-:155255234:1T36A37,chrY:-:59358240:1T36A37:200,"
      "chr15:-:102516664:1>1-28>5+8A37,chr16:+:64108:3>1-30>1+1>4+3A37,"
      "chr9:+:14540:3>1-34A33A2>1-",map_list);
  // Parse multiple new PE-map
  error_code+=gt_input_map_parse_map_list(
      "chr15:-:102516742:(3)3GCA67::chr15:+:102516611:76:::7936,"
      "chr16:+:64114:1>1-26>1+1>4+47::chr16:-:64196:68A6T:::12224,"
      "chr1:+:16731:(5)35>92*16(20),chrY:-:59355959:(5)6G4G24>1-3A1CA1AAA1>1-1(20),"
      "chrX:-:155252953:(5)6G4G24>1-3A1CA1AAA1>1-1(20)",map_list);

  /*
   * Parsing counters
   */

  /*
   * Parsing Alignments
   */
  register gt_alignment* alignment = gt_alignment_new();
  error_code+=gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20):65472,"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20):12224,"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20):1984,"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20):1984,"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20):1984",alignment);

  error_code+=gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "-",alignment);

  /*
   * Parsing Templates
   */
  register gt_template* template = gt_template_new();
  error_code+=gt_input_map_parse_template(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20),"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20),"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20),"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20),"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20)",template);
  printf("MCD = %" PRIu64 "\n",gt_template_get_mcs(template));

  printf("Test %s\n",error_code==0?"Passed":"Failed");
}


void gt_example_sam_parsing() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,false);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
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

void gt_constructor_merge_template() {
  gt_status error_code;
  gt_template* templateA = gt_template_new();
  gt_template* templateB = gt_template_new();

  /*
   * Alignment
   */
  gt_input_map_parse_template(
      "ID\tACGT\t####\t1\tchr1:-:20:4",templateA);
  gt_input_map_parse_template(
      "ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1",templateB);
  // Merge into source
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  gt_output_map_fprint_template(stdout,templateA,output_attributes);
  gt_output_map_fprint_template(stdout,templateB,output_attributes);
  gt_template_merge_template_mmaps(templateA,templateB);
  gt_output_map_fprint_template(stdout,templateA,output_attributes);


  /*
   * Template
   */
  error_code = gt_input_map_parse_template(
      "TAG_A\t"
      "ACTCCAGTCA AA\t"
      "########## ##\t"
      "1:0:1:0\t"
      "chr5:-:2001:10::chr5:+:3500:2,"
      "chr14:-:740:(5)4A::chr2:+:9087:2",templateA);
  //      "chr14:-:74094161:(5)4A::chr2:+:9087:2",templateA);
  if (error_code) {
    gt_fatal_error_msg("Error parsing A");
  }

  error_code = gt_input_map_parse_template(
      "TAG_B\t"
      "ACTCCAGTCA AA\t"
      "########## ##\t"
      "0:0:1:1\t"
      "chr14:-:740:(5)4A::chr2:+:9087:2,"
      "chr8:-:1009:(5)2C1A::chr1:+:9087:2"
      ,templateB);
  if (error_code) {
    gt_fatal_error_msg("Error parsing B");
  }

  gt_output_map_fprint_template(stdout,templateA,output_attributes);
  gt_output_map_fprint_template(stdout,templateB,output_attributes);

//  gt_template_merge_template_mmaps(templateA,templateB);
//  gt_output_map_fprint_template(stdout,templateA,GT_ALL,true);

  gt_template* tunion = gt_template_union_template_mmaps(templateA,templateB);
  gt_output_map_fprint_template(stdout,tunion,output_attributes);

  gt_template_delete(templateA);
  gt_template_delete(templateB);
  gt_template_delete(tunion);
}

void usage() {
  fprintf(stderr, "USE: ./gem-tools-examples -i input -o output \n"
                  "      Options::\n"
                  "        --input     <File>\n"
                  "        --output    <File>\n"
                  "        --help|h    \n");
}

void gt_constructor_copy_template() {
  gt_template* source = gt_template_new();
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  // test no maps no mmaps
  gt_input_map_parse_template(
      "ID\tACGT\t####\t1\tchr1:-:20:4",source);
  gt_template* copy = gt_template_copy(source, false, false);
  gt_string* string = gt_string_new(1024);
  gt_output_map_sprint_template(string, copy, output_attributes);
  // convert to string for simple check
  char * line = gt_string_get_string(string);
  gt_streq(line,"ID\tACGT\t####\t0\t-\n");
  gt_template_delete(copy);

  // test maps no mmaps
  gt_string_clear(string);
  copy = gt_template_copy(source, true, false);
  string = gt_string_new(1024);
  gt_output_map_sprint_template(string, copy, output_attributes);
  // convert to string for simple check
  line = gt_string_get_string(string);
  gt_streq(line,"ID\tACGT\t####\t1\tchr1:-:20:4\n");
  gt_template_delete(copy);

  // test maps no mmaps with multi map string
  gt_input_map_parse_template(
    "ID\tACGT\t####\t1\tchr1:-:20:4,chr2:-:40:4",source);

  gt_string_clear(string);
  copy = gt_template_copy(source, true, false);
  string = gt_string_new(1024);
  gt_output_map_sprint_template(string, copy, output_attributes);
  // convert to string for simple check
  line = gt_string_get_string(string);
  gt_streq(line,"ID\tACGT\t####\t1\tchr1:-:20:4,chr2:-:40:4\n");
  gt_template_delete(copy);

  // test maps and mmaps
  gt_string_clear(string);
  copy = gt_template_copy(source, true, true);
  string = gt_string_new(1024);
  gt_output_map_sprint_template(string, copy, output_attributes);
  // convert to string for simple check
  line = gt_string_get_string(string);
  gt_streq(line,"ID\tACGT\t####\t1\tchr1:-:20:4,chr2:-:40:4\n");
  gt_template_delete(copy);
}

void gt_constructor_parse_again() {
  gt_map* map = gt_map_new();
  gt_vector* map_list = gt_vector_new(10,sizeof(gt_map*));
  gt_alignment* alignment = gt_alignment_new();
  gt_template* template = gt_template_new();

  /*
   * Parse old map
   */
  gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map);
  /*
   * Parse new map
   */
  gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map);

  /*
   * Parsing single maps
   */
  // Parse old map
  gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map);
  // Parse new SE-map
  gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map);
  // Parse new PE-map
  gt_input_map_parse_map("chr15:-:102516634:66G9::chr15:+:102516358:66>1+10:::7936",map);

  /*
   * Parsing list of maps
   */
  // Parse multiple old split-map
  gt_input_map_parse_map("[26]=chr7:R1203797~chr7:R1203108",map);
  gt_input_map_parse_map_list("[31;35]=chr16:R[2503415;2503411]~chr16:R2503271",map_list);
  gt_input_map_parse_map_list("[30;34]=chr10:F74776624~chr10:F[74790025;74790029]",map_list);
  gt_input_map_parse_map_list("[23-50]=chr1:F[188862944-188868041]~chr19:F53208292",map_list);
  gt_input_map_parse_map_list("[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list);
  gt_input_map_parse_map_list("[26]=chr7:R1203797~chr7:R1203108",map_list);
  gt_input_map_parse_map_list(
      "chrM:F6598<+1>16@0/0,[23]=chr6:R31322884~chr6:R31237276,"
      "[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list);
  // Parse multiple old SE-map
  gt_input_map_parse_map_list("chr1:F8926499@0/0,chr12:R7027116G39A42@77/2",map_list);
  // Parse multiple new SE-map
  gt_input_map_parse_map_list(
      "chrX:-:155255234:1T36A37,chrY:-:59358240:1T36A37:200,"
      "chr15:-:102516664:1>1-28>5+8A37,chr16:+:64108:3>1-30>1+1>4+3A37,"
      "chr9:+:14540:3>1-34A33A2>1-",map_list);
  // Parse multiple new PE-map
  gt_input_map_parse_map_list(
      "chr15:-:102516742:(3)3GCA67::chr15:+:102516611:76:::7936,"
      "chr16:+:64114:1>1-26>1+1>4+47::chr16:-:64196:68A6T:::12224,"
      "chr1:+:16731:(5)35>92*16(20),chrY:-:59355959:(5)6G4G24>1-3A1CA1AAA1>1-1(20),"
      "chrX:-:155252953:(5)6G4G24>1-3A1CA1AAA1>1-1(20)",map_list);

  /*
   * Parsing counters
   */

  /*
   * Parsing Alignments
   */
  gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20):65472,"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20):12224,"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20):1984,"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20):1984,"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20):1984",alignment);
  gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "-",alignment);

  /*
   * Parsing Templates
   */
  gt_input_map_parse_template(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20),"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20),"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20),"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20),"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20)",template);
  /*
   * Check mcs delegation for single end
   */
  gt_template_get_mcs(template);
}

void gt_constructor_template_merge_hung() {
  gt_template* source = gt_template_new();
  gt_template* target = gt_template_new();

  gt_input_map_parse_template(
      "ID/1\t"
      "GAGAGAACAGGCCTCTGAGCCCAAGCCAAGCCATCGCATCCCCTGTGACTTGCCCGTATATATGCCCAGATGGCCTGAAGTAACTGAAGAATCACAAAAGA\t"
      "BPYcceeegegggiiiiiiiiiiiiiiighhhfgghhiiihhhifhfhhhfhhdghhiiihfbggedddeabbcdcccb`ZaW[^^abbcGW[`^`R]`BB\t"
      "0:0:0:0:1:4+27:46:45:34:12:16:3:0:1\t"
      "chr7:+:24435865:1T>5-46A8C38,chr2:+:64479547:>1-1T>1-1T>1-94,chr4:+:16999371:1T>5-35>2-9A8C38,chr9:+:90025610:1T>5-28A17A8C38,chr11:-:130623749:2>1-1T1T46A2>2-43",source);
  gt_input_map_parse_template(
      "ID/1\t"
      "GAGAGAACAGGCCTCTGAGCCCAAGCCAAGCCATCGCATCCCCTGTGACTTGCCCGTATATATGCCCAGATGGCCTGAAGTAACTGAAGAATCACAAAAGA\t"
      "BPYcceeegegggiiiiiiiiiiiiiiighhhfgghhiiihhhifhfhhhfhhdghhiiihfbggedddeabbcdcccb`ZaW[^^abbcGW[`^`R]`BB\t"
      "0:1:0:0:0:0:1:3\t"
      "chr2:+:64479108:10>436*91,chr2:+:64479544:TTC1TGT94,chr13:+:91496959:1GA1TGT27T18A47,chr2:+:216792393:1GA1TGT46A8C38,chr2:+:64478959:AG1T1G1TT1>585*91",target);
  // merge into source
  gt_template_merge_template_mmaps(source,target);
  gt_string* string = gt_string_new(1024);
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  gt_output_map_fprint_counters(stderr,gt_template_get_block(source,0)->counters);
  gt_output_map_sprint_template(string, source, output_attributes);
  // convert to string
  char * line = gt_string_get_string(string);
  printf("%s\n", line);
//  fail_unless(gt_streq(line,
//      "HISEQ8_0071:3:1101:19107:2010#TGACCA/1\t"
//      "NGTCATGAGTGCAAAATGCAAATGCAAGTTTGGCCAGAAGTCCGGTCACCATCCAGGGGAGACTCCACCTCTCATCACCCCAGGCTCAGCCCAAAGCTGAT\t"
//      "BPYcceeegegggiiiiiiiiiiiiiiighhhfgghhiiihhhifhfhhhfhhdghhiiihfbggedddeabbcdcccb`ZaW[^^abbcGW[`^`R]`BB\t"
//      "0:0:0:0:2+0:1\t"
//      "chr19:+:35613736:C7>78*37>316*38>132*18,chr19:+:35613741:(5)3>78*37>316*36(20),chr19:+:35613819:(5)CAG37>316*36(20)\n"));
}

void gt_remove_maps_with_n_or_more_mismatches() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_input_file_open("bug",false);
  gt_output_file* output_file = gt_output_stream_new(stdout, SORTED_FILE);
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  // Parallel reading+process
  #pragma omp parallel num_threads(4)
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    gt_status error_code;
    gt_alignment *alignment_src = gt_alignment_new();
    gt_alignment *alignment_dst;
    gt_generic_parser_attr generic_parser_attr = GENERIC_PARSER_ATTR_DEFAULT(false);
    while ((error_code = gt_input_generic_parser_get_alignment(buffered_input,alignment_src,&generic_parser_attr))) {
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
          gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
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
  register gt_sequence_archive* sequence_archive = gt_sequence_archive_new();
  if (gt_input_multifasta_parser_get_archive(input_file,sequence_archive)!=GT_IFP_OK) {
    gt_error_msg("Fatal error parsing reference\n");
  }

  // Dump summary of reference file
  if (gt_input_fasta_is_multifasta(input_file)) {
    fprintf(stderr,"File '%s' is MULTIFASTA\n",input_file->file_name);
  }
  gt_sequence_archive_iterator seq_arch_it;
  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
  register gt_segmented_sequence* seq;
  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
    fprintf(stderr,"SEQUENCE '%s' [length=%lu]\n",seq->seq_name->buffer,seq->sequence_total_length);
  }
  fprintf(stderr,"\n");

  // TEST 1

//  // Dump the content of the reference file
//  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
//  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
//    fprintf(stderr,">%s\n",seq->seq_name->buffer);
//    register uint64_t i;
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
  #pragma omp parallel num_threads(parameters.num_threads)
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

void gt_debug_mem_leak() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);

  // Buffered reading of the file
  gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
  gt_alignment* alignment = gt_alignment_new();
  gt_status error_code;
  while ((error_code=gt_input_map_parser_get_alignment(buffered_input,alignment))) {
    if (error_code==GT_BMI_FAIL) continue;
  }
  gt_buffered_input_file_close(buffered_input);

  // Close files
  gt_input_file_close(input_file);
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

int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  //
  // Load it!
  //

  // Leo's
  //gt_remove_maps_with_n_or_more_mismatches();

  //gt_constructor_template_merge_hung();

  //gt_constructor_parse_again();
  //gt_constructor_copy_template();

  //gt_constructor_merge_template();
  //gt_example_map_string_parsing();
  //gt_dummy_example();

  //gt_filter_fastq();
  //gt_load_reference__dump_it();

  gt_example_map_parsing();

  return 0;
}


