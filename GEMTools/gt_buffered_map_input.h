/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_map_input.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_BUFFERED_MAP_FILE_H_
#define GT_BUFFERED_MAP_FILE_H_

#include "gt_commons.h"
#include "gt_input_file.h"
#include "gt_template.h"

// Codes gt_status
#define GT_BMI_OK 1
#define GT_BMI_FAIL -1
#define GT_BMI_EOF 0
// PE (Parsing Errors)
// TODO

// Constants
#define GT_MAP_MCS '+'
#define GT_MAP_COUNTS_SEP ':'
#define GT_MAP_COUNTS_TIMES 'x'

typedef struct {
  /* Input file */
  gt_input_file* input_file;
  /* Block buffer and cursors */
  uint64_t block_id;
  gt_vector* block_buffer;
  char* cursor;
  uint64_t lines_in_buffer;
  uint64_t current_line_num;
} gt_buffered_map_input;

/*
 * Lazy parsing:
 *   PARSE_READ         Parses TAG,READ,QUALITY,COUNTERS
 *   PARSE_READ__MAPS   Parses TAG,READ,QUALITY,COUNTERS,MAPS
 *   PARSE_ALL          Parses TAG,READ,QUALITY,COUNTERS,MAPS,CIGAR
 */
typedef enum {PARSE_READ, PARSE_READ__MAPS, PARSE_ALL } gt_lazy_parse_mode;

/*
 * Checkers
 */
#define GT_BMI_CHECK(buffered_map_input) gt_fatal_check( \
  buffered_map_input==NULL||buffered_map_input->input_file|| \
  buffered_map_input->block_buffer==NULL||buffered_map_input->cursor==NULL,NULL_HANDLER)

/*
 * Buffered map file handlers
 */
gt_buffered_map_input* gt_buffered_map_input_new(gt_input_file* const input_file);
gt_status gt_buffered_map_input_close(gt_buffered_map_input* const buffered_map_input);
GT_INLINE uint64_t gt_buffered_map_input_get_cursor_pos(gt_buffered_map_input* const buffered_map_input);
GT_INLINE bool gt_buffered_map_input_eob(gt_buffered_map_input* const buffered_map_input);
GT_INLINE gt_status gt_buffered_map_input_get_block(gt_buffered_map_input* const buffered_map_input);

/*
 * MAP File Format test
 */
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors);
GT_INLINE bool gt_buffered_map_input_test_map(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    gt_map_file_format* const map_file_format,const bool show_errors);

/*
 * MAP/MAPQ/MMAP/MMAPQ Lazy Parsers
 *   Lazy parsing works in 3 steps
 *     (1) Parses TAG(s),READ(s),QUALITIES(s),COUNTERS -> gt_template/gt_alignment
 *     (2) Parses MAPS' SEQUENCE,STRAND,POSITION
 *     (3) Parses MAPS' CIGAR string
 */
GT_INLINE void gt_buffered_map_input_parse_error(
    gt_buffered_map_input* const buffered_map_input,
    const uint64_t line_num,const gt_status error_code);
GT_INLINE void gt_buffered_map_input_parse_next_record(gt_buffered_map_input* const buffered_map_input);
// Parse Template/Alignment
GT_INLINE gt_status gt_buffered_map_input_parse_template(
    gt_buffered_map_input* const buffered_map_input,gt_template* const template,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps);
GT_INLINE gt_status gt_buffered_map_input_parse_alignment(
    gt_buffered_map_input* const buffered_map_input,gt_alignment* alignment,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps);
// Parse Maps
GT_INLINE gt_status gt_buffered_map_input_parse_template_maps(gt_template* template,uint64_t num_maps);
GT_INLINE gt_status gt_buffered_map_input_parse_alignment_maps(gt_alignment* alignment,uint64_t num_maps);
// Parse Mismatches
GT_INLINE gt_status gt_buffered_map_input_parse_template_mismatch_string(
    gt_template* template,const gt_map_version map_file_format);
GT_INLINE gt_status gt_buffered_map_input_parse_alignment_mismatch_string(
    gt_alignment* alignment,const gt_map_version map_file_format);

/*
 * MAP/MAPQ/MMAP/MMAPQ High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_buffered_map_input_get_template(
    gt_buffered_map_input* const buffered_map_input,gt_template* const template);
GT_INLINE gt_status gt_buffered_map_input_get_alignment(
    gt_buffered_map_input* const buffered_map_input,gt_alignment* const alignment);

#endif /* GT_BUFFERED_MAP_FILE_H_ */
