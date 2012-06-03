/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_FILE_H_
#define GT_INPUT_FILE_H_

#include "gt_commons.h"

typedef enum { STREAM, REGULAR_FILE, MAPPED_FILE } gt_file_type;
typedef enum { FASTA, FASTQ, MAP, MAPQ, MMAP, MMAPQ, UNKNOWN } gt_file_format;
typedef struct {
  /* Input file */
  char* file_name;
  gt_file_type file_type;
  FILE* file;
  int fildes;
  uint8_t* mapped_file;
  bool eof;
  uint64_t file_size;
  gt_file_format file_format;
  pthread_mutex_t input_mutex;
  /* Auxiliary Buffer (for synch purposes) */
  uint8_t* synch_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
} gt_input_file;

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream);
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file);
void gt_input_file_close(gt_input_file* const input_file);
gt_file_format gt_input_file_get_file_format(gt_input_file* const input_file);

/*
 * Advanced I/O
 *   // TODO: Nested conditions would need ...
 */
gt_input_file* gt_input_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t segment_number,const uint64_t total_segments);
gt_input_file* gt_input_reads_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t num_init_line,const uint64_t num_end_line);

/*
 * Mutex functions
 */
inline gt_status gt_input_get_lock(gt_input_file* const input_file);
inline gt_status gt_input_release_lock(gt_input_file* const input_file);

/*
 * Reading from input (NO thread safe, must call mutex functions before)
 */
inline size_t gt_input_get_lines(gt_input_file* const input_file,gt_vector* buffer_dst);

#endif /* GT_INPUT_FILE_H_ */
