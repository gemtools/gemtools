/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_FILE_H_
#define GT_INPUT_FILE_H_

#include "gt_commons.h"

#include <zlib.h>
#include <bzlib.h>

// Codes gt_status
#define GT_INPUT_FILE_OK 0
#define GT_INPUT_FILE_CLOSE_ERR 1
#define GT_INPUT_FILE_EOF 0
#define GT_INPUT_FILE_LINE_READ 1

/*
 * Checkers
 */
#define GT_INPUT_FILE_CHECK(input_file) \
  GT_NULL_CHECK(input_file); \
  GT_NULL_CHECK(input_file->file_name)

/*
 * File specifics (formats, attributes, ...)
 */
typedef enum { FASTA, MAP, SAM, FILE_FORMAT_UNKNOWN } gt_file_format;
// MAP specific info
typedef struct {
  bool contains_qualities;
} gt_map_file_format;
// FASTQ/FASTA specific info
typedef enum { F_FASTA, F_FASTQ, F_MULTI_FASTA } gt_file_fasta_format;
typedef struct {
  gt_file_fasta_format fasta_format;
} gt_fasta_file_format;
// SAM specific info (headers)
typedef struct {
  // gt_reference_sequences reference_sequences;
  char* program_name;
  char* program_version;
  /* ... */
} gt_sam_headers;
/* */


/*
 * GT Input file
 */
typedef enum { STREAM, REGULAR_FILE, MAPPED_FILE, GZIPPED_FILE, BZIPPED_FILE } gt_file_type;
typedef struct {
  /* Input file */
  char* file_name;
  gt_file_type file_type;
  void* file;
  int fildes;
  bool eof;
  uint64_t file_size;
  /* File format */
  gt_file_format file_format;
  union {
    gt_map_file_format map_type;
    gt_fasta_file_format fasta_type;
    gt_sam_headers sam_headers;
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
 *   // TODO
 */
//gt_input_file* gt_input_file_segmented_file_open(
//    char* const file_name,const bool mmap_file,
//    const uint64_t segment_number,const uint64_t total_segments);
//gt_input_file* gt_input_file_reads_segmented_file_open(
//    char* const file_name,const bool mmap_file,
//    const uint64_t num_init_line,const uint64_t num_end_line);
/* Format detection */
gt_file_format gt_input_file_detect_file_format(gt_input_file* const input_file);

/*
 * Accessors (Mutex,ID,...) functions
 */
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file);
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file);
GT_INLINE uint64_t gt_input_file_next_id(gt_input_file* const input_file);

/*
 * Basic line functions
 */
GT_INLINE size_t gt_input_file_dump_to_buffer(gt_input_file* const input_file,gt_vector* const buffer_dst);
GT_INLINE size_t gt_input_file_fill_buffer(gt_input_file* const input_file);
GT_INLINE size_t gt_input_file_next_line(gt_input_file* const input_file,gt_vector* const buffer_dst);

GT_INLINE size_t gt_input_file_next_map_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,uint64_t* const num_blocks);
GT_INLINE size_t gt_input_file_next_sam_record(
    gt_input_file* const input_file,gt_vector* const buffer_dst,gt_string* const first_field);
GT_INLINE bool gt_input_file_cmp_next_sam_record(gt_input_file* const input_file,gt_string* const reference_tag);

/*
 * Line Readers (thread-unsafe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines);

/*
 * Processing Macros (direct parsing from input file)
 */
#define GT_INPUT_FILE_CHECK_BUFFER(input_file) \
  if (gt_expect_false(input_file->buffer_pos >= input_file->buffer_size)) { \
    gt_input_file_fill_buffer(input_file); \
  }
#define GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst) \
  if (gt_expect_false(input_file->buffer_pos >= input_file->buffer_size)) { \
    if (gt_expect_true(buffer_dst!=NULL)) gt_input_file_dump_to_buffer(input_file,buffer_dst); \
    gt_input_file_fill_buffer(input_file); \
  }

#define GT_INPUT_FILE_NEXT_CHAR(input_file) \
  ++input_file->buffer_pos; \
  GT_INPUT_FILE_CHECK_BUFFER(input_file)
#define GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst) \
  ++input_file->buffer_pos; \
  GT_INPUT_FILE_CHECK_BUFFER__DUMP(input_file,buffer_dst)

#define GT_INPUT_FILE_SKIP_EOL(input_file) \
  if (!input_file->eof) { \
    GT_INPUT_FILE_NEXT_CHAR(input_file); /* Skip EOF/DOS_EOL */  \
    if (gt_expect_false(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==EOL)) { \
      GT_INPUT_FILE_NEXT_CHAR(input_file); /* Skip DOS_EOL */  \
    } \
  }
#define GT_INPUT_FILE_HANDLE_EOL(input_file,buffer_dst) \
  if (!input_file->eof) { \
    GT_INPUT_FILE_NEXT_CHAR__DUMP(input_file,buffer_dst); /* Skip EOF/DOS_EOL */  \
    if (gt_expect_false(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==EOL)) { \
      ++input_file->buffer_pos; \
      if (gt_expect_true(buffer_dst!=NULL)) { \
        gt_input_file_dump_to_buffer(input_file,buffer_dst); \
        gt_vector_dec_used(buffer_dst); \
        *gt_vector_get_last_elm(buffer_dst,char)=EOL; \
      } \
      gt_input_file_fill_buffer(input_file); \
    } \
  }

#define GT_INPUT_FILE_CURRENT_CHAR(input_file) input_file->file_buffer[input_file->buffer_pos]


#endif /* GT_INPUT_FILE_H_ */
