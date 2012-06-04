/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_buffer.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_BUFFER_H_
#define GT_OUTPUT_BUFFER_H_

#include "gt_commons.h"
#include "gt_buffered_output_file.h"

typedef enum { GT_OUTPUT_BUFFER_FREE, GT_OUTPUT_BUFFER_BUSY, GT_OUTPUT_BUFFER_WRITE_PENDING } gt_output_buffer_state;
typedef struct {
  gt_buffered_output_file* output_file;
  /* Buffer */
  uint8_t buffer_id;
  gt_vector* buffer;
  gt_output_buffer_state buffer_state;
} gt_output_buffer;

#endif /* GT_OUTPUT_BUFFER_H_ */
