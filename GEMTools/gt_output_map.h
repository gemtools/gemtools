/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_MAP_H_
#define GT_OUTPUT_MAP_H_

#include "gt_commons.h"
#include "gt_template.h"
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"

/*
 * Buffered Map Output File
 */
typedef struct {
  uint64_t block_id;
  gt_output_buffer* buffer;
  gt_buffered_output_file* buffered_output_file;
} gt_buffered_map_output;

/*
 * Setup
 */
gt_buffered_map_output* gt_buffered_map_output_new(gt_buffered_output_file* const buffered_output_file);
gt_status gt_buffered_map_output_close(gt_buffered_map_output* const buffered_map_output);

/*
 * MAP building block printers
 */
GT_INLINE gt_status gt_output_map_fprint_counters(
    FILE* file,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_fprint_map(FILE* file,gt_map* const map);
// TODO

/*
 * High-level MAP Printers
 */
GT_INLINE gt_status gt_buffered_map_output_print_template(gt_buffered_map_output* const buffered_map_output,gt_template* const template);
GT_INLINE gt_status gt_buffered_map_output_print_alignment(gt_buffered_map_output* const buffered_map_output,gt_alignment* const alignment);
/* */
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer *output_buffer,gt_template* const template);
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer *output_buffer,gt_alignment* const alignment);
GT_INLINE gt_status gt_output_map_sprint_template(char **line_ptr,gt_template* const template);
GT_INLINE gt_status gt_output_map_sprint_alignment(char **line_ptr,gt_alignment* const alignment);
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template);
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment);

#endif /* GT_OUTPUT_MAP_H_ */
