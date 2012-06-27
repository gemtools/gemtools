/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_FILE_H_
#define GT_INPUT_FILE_H_

#include "gt_commons.h"

// Codes gt_status
#define GT_INPUT_FILE_OK 0
#define GT_INPUT_FILE_CLOSE_ERR 1
#define GT_INPUT_FILE_EOF 0
#define GT_INPUT_FILE_LINE_READ 1

/*
 * File formats
 */
typedef enum { FASTA, FASTQ, MAP, UNKNOWN } gt_file_format;
// MAP specific info
typedef struct {
  bool contains_qualities;
  char separator;
  // uint64_t num_blocks_template; /* As we mixed files, this can vary */
  // gt_map_version format_version; /* We might even tolerate mixtures */
} gt_map_file_format;
// FASTQ/FASTA specific info
typedef struct {/*TODO*/} gt_fast_file_format;
/* */

// GT Input file
typedef enum { STREAM, REGULAR_FILE, MAPPED_FILE } gt_file_type;
typedef struct {
  /* Input file */
  char* file_name;
  gt_file_type file_type;
  FILE* file;
  int fildes;
  bool eof;
  uint64_t file_size;
  /* File format */
  gt_file_format file_format;
  union {
    gt_map_file_format map_type;
    gt_fast_file_format fast_type;
  };
  pthread_mutex_t input_mutex;
  /* Auxiliary Buffer (for synch purposes) */
  uint8_t* file_buffer;
  uint64_t buffer_size;
  uint64_t buffer_begin;
  uint64_t buffer_pos;
  uint64_t global_pos;
  uint64_t processed_lines;
  /* ID generator */
  uint64_t processed_id;
} gt_input_file;

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream);
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file);
gt_status gt_input_file_close(gt_input_file* const input_file);

/*
 * Advanced I/O
 *   // TODO: Nested conditions would need ...
 */
gt_input_file* gt_input_file_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t segment_number,const uint64_t total_segments);
gt_input_file* gt_input_file_reads_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t num_init_line,const uint64_t num_end_line);

/*
 * Mutex functions
 */
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file);
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file);
GT_INLINE uint64_t gt_input_file_next_id(gt_input_file* const input_file);

/*
 * Format detection
 */
gt_file_format gt_input_file_detect_file_format(gt_input_file* const input_file);

/*
 * Reading from input (NO thread safe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,
    gt_vector* buffer_dst,const uint64_t num_lines);

#endif /* GT_INPUT_FILE_H_ */
