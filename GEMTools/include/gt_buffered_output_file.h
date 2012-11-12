/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_output_file.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_BUFFERED_OUTPUT_FILE_H_
#define GT_BUFFERED_OUTPUT_FILE_H_

#include "gt_commons.h"
#include "gt_output_buffer.h"


typedef struct {
  /* Output file */
  gt_output_file* output_file;
  /* Output Buffer */
  gt_output_buffer* buffer;
} gt_buffered_output_file;

// Codes gt_status
#define GT_BUFFERED_OUTPUT_FILE_OK 0
#define GT_BUFFERED_OUTPUT_FILE_FAIL -1

/*
 * Buffered Output File Setup
 */
gt_buffered_output_file* gt_buffered_output_file_new(gt_output_file* const output_file);
gt_status gt_buffered_output_file_close(gt_buffered_output_file* const buffered_output_file);

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* gt_buffered_output_file_request_buffer(
    gt_buffered_output_file* const buffered_output_file);
GT_INLINE void gt_buffered_output_file_release_buffer(
    gt_buffered_output_file* const buffered_output_file,gt_output_buffer* const output_buffer);
GT_INLINE gt_output_buffer* gt_buffered_output_file_dump_buffer(
    gt_buffered_output_file* const buffered_output_file,gt_output_buffer* const output_buffer);

#endif /* GT_BUFFERED_OUTPUT_FILE_H_ */
