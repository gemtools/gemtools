/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_buffer.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_BUFFER_H_
#define GT_OUTPUT_BUFFER_H_

#include "gt_commons.h"

typedef enum { GT_OUTPUT_BUFFER_FREE, GT_OUTPUT_BUFFER_BUSY, GT_OUTPUT_BUFFER_WRITE_PENDING } gt_output_buffer_state;
typedef struct {
  uint8_t buffer_id;
  uint32_t block_id;
  gt_vector* buffer;
  gt_output_buffer_state buffer_state;
} gt_output_buffer;

/*
 * Setup
 */
GT_INLINE gt_output_buffer* gem_output_buffer_new(void);
GT_INLINE void gem_output_buffer_clear(gt_output_buffer* output_buffer);
GT_INLINE void gem_output_buffer_delete(gt_output_buffer* output_buffer);

#endif /* GT_OUTPUT_BUFFER_H_ */
