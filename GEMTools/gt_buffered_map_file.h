/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_map_file.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_BUFFERED_MAP_FILE_H_
#define GT_BUFFERED_MAP_FILE_H_

#include "gt_commons.h"
#include "gt_input_file.h"

typedef struct {
  /* Input file */
  gt_input_file* input_file;
  /* Block ID */
  uint64_t block_id;
  /* Block buffer and cursors */
  gt_vector* block_buffer;
  uint8_t* cursor;
  uint64_t cursor_pos;
  uint64_t current_line_num;
} gt_buffered_map_file;

// Codes gt_status
#define GT_BMF_OK 1
#define GT_BMF_FAIL 0

/*
 * Buffered map file handlers
 */
gt_buffered_map_file* gt_buffered_map_file_new(gt_input_file* const input_file);
gt_status gt_buffered_map_file_delete(gt_buffered_map_file* const buffered_map_input);
bool gt_buffered_map_file_eof(gt_buffered_map_file* const buffered_map_input);

/*
 * MAP/MAPQ/MMAP/MMAPQ Simple Parsers
 *   Simple parsing procedures to extract one template/alignment from
 *   the buffered file (reads one line)
 *   Simple extraction involves syntax checking
 */
gt_status gt_buffered_map_file_get_template(gt_buffered_map_file* const buffered_map_input,gt_template* template);
/*
 * If file type is MAP/MAPQ then the alignment can be extracted
 * directly without wrapping it in a template.
 * PRE: FileType is MAP/MAPQ (otherwise an error will be returned)
 */
gt_status gt_buffered_map_file_get_alignment(gt_buffered_map_file* const buffered_map_input,gt_alignment* alignment);

/*
 * MAP/MAPQ/MMAP/MMAPQ Lazy Parsers
 *   Lazy parsing works in 3 steps
 *     (1) Parses TAG(s),READ(s),QUALITIES(s),COUNTERS
 *     (2) Parses MAPS' SEQUENCE,STRAND,POSITION
 *     (3) Parses MAPS' CIGAR string
 */
gt_status gt_buffered_map_file_lazy_get_template(gt_buffered_map_file* const buffered_map_input,gt_template* template);
gt_status gt_buffered_map_file_lazy_get_alignment(gt_buffered_map_file* const buffered_map_input,gt_alignment* alignment);

#endif /* GT_BUFFERED_MAP_FILE_H_ */
