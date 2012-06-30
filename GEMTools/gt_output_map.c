/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_commons.h"
#include "gt_output_map.h"

/*
 * Setup
 */
gt_buffered_map_output* gt_buffered_map_output_new(gt_buffered_output_file* const buffered_output_file) {
  //TODO
}
gt_status gt_buffered_map_output_close(gt_buffered_map_output* const buffered_map_output) {
  //TODO
}

/*
 * MAP building block printers
 */
GT_INLINE gt_status gt_output_map_fprint_counters(
    FILE* file,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact) {
  //TODO
}
GT_INLINE gt_status gt_output_map_fprint_map(FILE* file,gt_map* const map) {
  //TODO
}

/*
 * High-level MAP Printers
 */
GT_INLINE gt_status gt_buffered_map_output_print_template(gt_buffered_map_output* const buffered_map_output,gt_template* const template) {
  //TODO
}
GT_INLINE gt_status gt_buffered_map_output_print_alignment(gt_buffered_map_output* const buffered_map_output,gt_alignment* const alignment) {
  //TODO
}
/* */
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer *output_buffer,gt_template* const template) {
  //TODO
}
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer *output_buffer,gt_alignment* const alignment) {
  //TODO
}
GT_INLINE gt_status gt_output_map_sprint_template(char **line_ptr,gt_template* const template) {
  //TODO
}
GT_INLINE gt_status gt_output_map_sprint_alignment(char **line_ptr,gt_alignment* const alignment) {
  //TODO
}
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template) {
  //TODO
}
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment) {
  //TODO
}
