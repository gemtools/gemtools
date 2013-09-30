/*
 * sel_reads.h
 *
 *  Created on: 8 Jul 2013
 *      Author: heath
 */

#ifndef SEL_READS_H_
#define SEL_READS_H_

#define DEFAULT_INS_CUTOFF 0.01 /* Insert sizes in the upper or lower cutoff percentiles will not be used */
#define MAP_THRESHOLD 3
#define AP_BUF_SIZE 16384
#define QUAL_FASTQ 33
#define QUAL_SOLEXA 64
#define SOLEXA_BAD_QUAL 2
#define DEFAULT_QUAL (QUAL_FASTQ)
#define MISSING_QUAL 40 // Value to use in alignment score if qualities not available
#define INDEL_QUAL 40 // Value to use in alignment score for indels
#define MAX_GT_SCORE 0xFFFF
#define MAX_QUAL 42
#define ALIGN_NORM 0
#define ALIGN_BS_POS 1
#define ALIGN_BS_NEG 2
#define ALIGN_TYPES 3
#define AL_FORWARD 4
#define AL_REVERSE 8
#define AL_DIRECTIONS ((AL_FORWARD)|(AL_REVERSE))
#define AL_USED 128

typedef struct {
  char *input_files[2];
  char *output_file;
  char *dist_file;
  double ins_cutoff;
  bool mmap_input;
  gt_output_file_compression compress;
  gt_generic_parser_attributes *parser_attr;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output[2];
  int64_t min_insert;
  int64_t max_insert;
  double *ins_dist;
  uint8_t *ins_phred;
  int num_threads;
  int indel_quality;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina)
} sr_param;

#endif /* SEL_READS_H_ */
