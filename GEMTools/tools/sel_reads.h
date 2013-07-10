/*
 * sel_reads.h
 *
 *  Created on: 8 Jul 2013
 *      Author: heath
 */

#ifndef SEL_READS_H_
#define SEL_READS_H_

#define DEFAULT_INS_CUTOFF 0.05 /* Insert sizes in the upper or lower cutoff percentiles will not be used */
#define MAP_THRESHOLD 3
#define AP_BUF_SIZE 16384
#define QUAL_FASTQ 33
#define QUAL_SOLEXA 64
#define SOLEXA_BAD_QUAL 2
#define DEFAULT_QUAL (QUAL_FASTQ)
#define DEFAULT_LTRIM 0
#define DEFAULT_RTRIM 0

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
  bool compress;
  bool combined;
  gt_generic_parser_attributes *parser_attr;
  int num_threads;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina)
  int bad_qual; // Normally 2 for Illumina
  int left,right; // Set ends of reads to 'N'
} sr_param;

#endif /* SEL_READS_H_ */
