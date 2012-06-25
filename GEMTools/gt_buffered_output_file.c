/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_output_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_buffered_output_file.h"

/*
 * Setup
 */
gt_buffered_output_file* gt_buffered_output_stream_new(FILE* file,const gt_output_file_type output_file_type) {
  // TODO
}
gt_buffered_output_file* gt_buffered_output_file_new(char* const file_name,const gt_output_file_type output_file_type) {
  // TODO
}
gt_status gt_buffered_output_file_close(gt_buffered_output_file* const buffered_output_file) {
  // TODO
}

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* gt_buffered_output_file_get_buffer(
    gt_buffered_output_file* const buffered_output_file) {
  // TODO
}
GT_INLINE gt_status gt_buffered_output_file_dump_buffer(
    gt_buffered_output_file* const buffered_output_file,
    gt_output_buffer* const output_buffer) {
  // TODO
}
