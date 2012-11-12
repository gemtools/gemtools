/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_file.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_FILE_H_
#define GT_OUTPUT_FILE_H_

#include "gt_commons.h"
#include "gt_output_buffer.h"

#define GT_MAX_OUTPUT_BUFFERS 10

typedef enum { SORTED_FILE, UNSORTED_FILE } gt_output_file_type;
typedef struct {
  /* Output file */
  char* file_name;
  FILE* file;
  gt_output_file_type file_type;
  /* Output Buffers */
  gt_output_buffer* buffer[GT_MAX_OUTPUT_BUFFERS];
  uint64_t buffer_busy;
  uint64_t buffer_write_pending;
  /* Block ID (for synchronization purposes) */
  uint32_t mayor_block_id;
  uint32_t minor_block_id;
  /* Mutexes */
  pthread_mutex_t out_buffer_mutex;
  pthread_cond_t  out_buffer_cond;
  pthread_mutex_t out_file_mutex;
} gt_output_file;

// Codes gt_status
#define GT_OUTPUT_FILE_OK 0
#define GT_OUTPUT_FILE_FAIL -1

/*
 * Output File Setup
 */
gt_output_file* gt_output_stream_new(FILE* const file,const gt_output_file_type output_file_type);
gt_output_file* gt_output_file_new(char* const file_name,const gt_output_file_type output_file_type);
gt_status gt_output_file_close(gt_output_file* const output_file);

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* gt_output_file_request_buffer(gt_output_file* const output_file);
GT_INLINE void gt_output_file_release_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer);
GT_INLINE gt_output_buffer* gt_output_file_dump_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer);

#endif /* GT_OUTPUT_FILE_H_ */
