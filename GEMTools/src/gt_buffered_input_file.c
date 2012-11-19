/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_input_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_buffered_input_file.h"

#define GT_BMI_BUFFER_SIZE GT_BUFFER_SIZE_4M
#define GT_BMI_NUM_LINES GT_NUM_LINES_5K

/*
 * Buffered map file handlers
 */
gt_buffered_input_file* gt_buffered_input_file_new(gt_input_file* const input_file) {
  GT_NULL_CHECK(input_file);
  gt_buffered_input_file* buffered_map_file = malloc(sizeof(gt_buffered_input_file));
  gt_cond_fatal_error(!buffered_map_file,MEM_HANDLER);
  /* Input file */
  buffered_map_file->input_file = input_file;
  /* Block buffer and cursors */
  buffered_map_file->block_id = UINT32_MAX;
  buffered_map_file->block_buffer = gt_vector_new(GT_BMI_BUFFER_SIZE,sizeof(uint8_t));
  buffered_map_file->cursor = (char*) gt_vector_get_mem(buffered_map_file->block_buffer,uint8_t);
  buffered_map_file->current_line_num = UINT64_MAX;
  /* Attached output buffer */
  buffered_map_file->buffered_output_file = NULL;
  return buffered_map_file;
}
gt_status gt_buffered_input_file_close(gt_buffered_input_file* const buffered_input_file) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  gt_vector_delete(buffered_input_file->block_buffer);
  return GT_BMI_OK;
}
GT_INLINE uint64_t gt_buffered_input_file_get_cursor_pos(gt_buffered_input_file* const buffered_input_file) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  GT_NULL_CHECK(buffered_input_file->cursor);
  return buffered_input_file->cursor-gt_vector_get_mem(buffered_input_file->block_buffer,char);
}
GT_INLINE bool gt_buffered_input_file_eob(gt_buffered_input_file* const buffered_input_file) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  return gt_buffered_input_file_get_cursor_pos(buffered_input_file) >= gt_vector_get_used(buffered_input_file->block_buffer);
}
GT_INLINE gt_status gt_buffered_input_file_get_block(
    gt_buffered_input_file* const buffered_input_file,const uint64_t num_lines,const bool use_mutex) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  register gt_input_file* const input_file = buffered_input_file->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  if (gt_expect_true(use_mutex)) gt_input_file_lock(input_file);
  if (input_file->eof) {
    if (gt_expect_true(use_mutex)) gt_input_file_unlock(input_file);
    return GT_BMI_EOF;
  }
  buffered_input_file->block_id = gt_input_file_next_id(input_file) % UINT32_MAX;
  buffered_input_file->current_line_num = input_file->processed_lines+1;
  buffered_input_file->lines_in_buffer = // FIXME:
      gt_input_file_get_lines(input_file,buffered_input_file->block_buffer,
          gt_expect_true(num_lines)?num_lines:GT_BMI_NUM_LINES);
  if (gt_expect_true(use_mutex)) gt_input_file_unlock(input_file);
  // Setup the block
  buffered_input_file->cursor = gt_vector_get_mem(buffered_input_file->block_buffer,char);
  return buffered_input_file->lines_in_buffer;
}
GT_INLINE gt_status gt_buffered_input_file_add_lines_to_block(
    gt_buffered_input_file* const buffered_input_file,const uint64_t num_lines) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  register gt_input_file* const input_file = buffered_input_file->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  register const uint64_t current_position =
      buffered_input_file->cursor - gt_vector_get_mem(buffered_input_file->block_buffer,char);
  register const uint64_t lines_added =
      gt_input_file_get_lines(input_file,buffered_input_file->block_buffer,
          gt_expect_true(num_lines)?num_lines:GT_BMI_NUM_LINES);
  buffered_input_file->lines_in_buffer += lines_added;
  buffered_input_file->cursor = gt_vector_get_elm(buffered_input_file->block_buffer,current_position,char);
  return lines_added;
}

/*
 * Block Synchronization with Output
 */
GT_INLINE void gt_buffered_input_file_attach_buffered_output(
    gt_buffered_input_file* const buffered_input_file,gt_buffered_output_file* const buffered_output_file) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file);
  buffered_input_file->buffered_output_file = buffered_output_file;
}
