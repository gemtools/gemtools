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

struct _gt_buffered_output_file;
typedef struct {
  /* Block ID (for synchronization purposes) */
  uint32_t mayor_block_id;
  uint32_t minor_block_id;
  bool is_final_block;
  /* Buffer */
//  uint8_t buffer_id; // FIXME
  gt_vector* buffer;
  gt_output_buffer_state buffer_state;
  /* Buffered output file*/
  struct _gt_buffered_output_file* buffered_output_file;
} gt_output_buffer;

// Codes gt_status
#define GT_OUT_BUFFER_OK 0
#define GT_OUT_BUFFER_FAIL -1

/*
 * Checkers
 */
#define GT_OUTPUT_BUFFER_CHECK(output_buffer) gt_fatal_check( \
    (output_buffer)==NULL||(output_buffer)->buffer==NULL,NULL_HANDLER)

/*
 * Setup
 */
GT_INLINE gt_output_buffer* gt_output_buffer_new(void);
GT_INLINE void gt_output_buffer_clear(gt_output_buffer* const output_buffer);
GT_INLINE void gt_output_buffer_initiallize(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state);
GT_INLINE void gt_output_buffer_delete(gt_output_buffer* const output_buffer);

/*
 * Accessors
 */
GT_INLINE void gt_output_buffer_set_state(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state);
GT_INLINE gt_output_buffer_state gt_output_buffer_get_state(gt_output_buffer* const output_buffer);
GT_INLINE void gt_output_buffer_set_mayor_block_id(gt_output_buffer* const output_buffer,const uint32_t mayor_block_id);
GT_INLINE uint32_t gt_output_buffer_get_mayor_block_id(gt_output_buffer* const output_buffer);
GT_INLINE void gt_output_buffer_set_minor_block_id(gt_output_buffer* const output_buffer,const uint32_t minor_block_id);
GT_INLINE void gt_output_buffer_inc_minor_block_id(gt_output_buffer* const output_buffer);
GT_INLINE uint32_t gt_output_buffer_get_minor_block_id(gt_output_buffer* const output_buffer);
GT_INLINE void gt_output_buffer_set_partial_block(gt_output_buffer* const output_buffer);

/*
 * Adaptors
 */
GT_INLINE char* gt_output_buffer_to_string(gt_output_buffer* const output_buffer);
GT_INLINE gt_vector* gt_output_buffer_to_vchar(gt_output_buffer* const output_buffer);

/*
 * Buffer printer
 */
GT_INLINE gt_status gt_vbprintf(gt_output_buffer** const output_buffer,const char *template,va_list v_args);
GT_INLINE gt_status gt_bprintf(gt_output_buffer** const output_buffer,const char *template,...);
// If you happen to know how much memory are you going to use
GT_INLINE gt_status gt_vbprintf_mem(
    gt_output_buffer** const output_buffer,const uint64_t expected_mem_usage,
    const char *template,va_list v_args);
GT_INLINE gt_status gt_bprintf_mem(
    gt_output_buffer** const output_buffer,const uint64_t expected_mem_usage,const char *template,...);

/*
 * Memory usage helper functions
 */
GT_INLINE uint64_t gt_calculate_memory_required_v(const char *template,va_list v_args);
GT_INLINE uint64_t gt_calculate_memory_required_va(const char *template,...);

#endif /* GT_OUTPUT_BUFFER_H_ */
