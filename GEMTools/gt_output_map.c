/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_output_map.h"

/*
 * Buffered Map Output File
 */
typedef struct {
  uint64_t block_id;
  gt_output_buffer* buffer;
  gt_buffered_output_file* buffered_output_file;
} gt_buffered_map_output;

gt_buffered_map_output* gt_buffered_map_output_new(gt_buffered_output_file* const buffered_output_file) {

}
gt_status gt_buffered_map_output_close(gt_buffered_map_output* const buffered_map_output) {

}
GT_INLINE gt_status gt_buffered_map_output_print_template(
    gt_buffered_map_output* const buffered_map_output,gt_template* const template) {

}
GT_INLINE gt_status gt_buffered_map_output_print_alignment(
    gt_buffered_map_output* const buffered_map_output,gt_alignment* const alignment) {

}

/*
 * MAP Printers
 */
GT_INLINE gt_status gt_output_map_print_template(gt_output_buffer *output_buffer) {

}
GT_INLINE gt_status gt_output_map_print_alignment(gt_output_buffer *output_buffer) {

}
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file) {

}
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file) {

}
