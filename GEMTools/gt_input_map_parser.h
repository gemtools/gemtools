/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_parser.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_MAP_PARSER_H_
#define GT_INPUT_MAP_PARSER_H_

#include "gt_commons.h"
#include "gt_template.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_buffered_output_file.h"

#include "gt_input_parser.h"

// Codes gt_status
#define GT_IMP_OK 1
#define GT_IMP_FAIL -1
#define GT_IMP_EOF 0

/*
 * Parsing error/state codes // TODO Move common errors to input_parser.h
 */
#define GT_IMP_PE_WRONG_FILE_FORMAT 10
#define GT_IMP_PE_NOT_IMPLEMENTED 11

#define GT_IMP_PE_PREMATURE_EOL 20
#define GT_IMP_PE_PENDING_BLOCKS 21
#define GT_IMP_PE_EOB 22
#define GT_IMP_PE_BAD_SEPARATOR 23
#define GT_IMP_PE_BAD_NUMBER_OF_BLOCKS 24
#define GT_IMP_PE_BAD_CHARACTER 25

#define GT_IMP_PE_READ_BAD_CHARACTER 30

#define GT_IMP_PE_QUAL_BAD_PREMATURE_EOB 40
#define GT_IMP_PE_QUAL_BAD_CHARACTER 41
#define GT_IMP_PE_QUAL_BAD_SEPARATOR 42

#define GT_IMP_PE_COUNTERS_BAD_CHARACTER 50

#define GT_IMP_PE_MAP_PENDING_MAPS 60
#define GT_IMP_PE_MAP_ALREADY_PARSED 61
#define GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS 62
#define GT_IMP_PE_MAP_BAD_CHARACTER 63

#define GT_IMP_PE_MISMS_ALREADY_PARSED 70
#define GT_IMP_PE_MISMS_BAD_CHARACTER 71
#define GT_IMP_PE_MISMS_BAD_MISMS_POS 72

/*
 * MAP file format constants
 */
#define GT_MAP_MCS '+'
#define GT_MAP_COUNTS_SEP ':'
#define GT_MAP_COUNTS_TIMES 'x'
#define GT_MAP_SEP ':'
#define GT_MAP_NONE '-'
#define GT_MAP_NEXT ','

#define GT_MAP_SEP_S ":"
#define GT_MAP_NEXT_S ","
#define GT_MAP_TEMPLATE_SEP "::"
#define GT_MAP_TEMPLATE_SCORE ":::"
#define GT_MAP_SMCS "+"
#define GT_MAP_COUNTS_STIMES "x"
#define GT_MAP_COUNTS_NOT_UNIQUE '!'
#define GT_MAP_COUNTS_NOT_UNIQUE_S "!"

#define GT_MAP_STRAND_FORWARD_SYMBOL '+'
#define GT_MAP_STRAND_FORWARD_LETTER 'F'
#define GT_MAP_STRAND_REVERSE_SYMBOL '-'
#define GT_MAP_STRAND_REVERSE_LETTER 'R'

#define GT_MAP_INDEL_INSERTION '+'
#define GT_MAP_INDEL_DELETION '-'
#define GT_MAP_INDEL_SPLICE '*'

#define GT_MAP_SKIP_POSITIVE '+'
#define GT_MAP_SKIP_NEGATIVE '-'
#define GT_MAP_SKIP_SPLICE '*'

/*
 * Lazy parsing:
 *   PARSE_READ         Parses TAG,READ,QUALITY,COUNTERS
 *   PARSE_READ__MAPS   Parses TAG,READ,QUALITY,COUNTERS,MAPS
 *   PARSE_ALL          Parses TAG,READ,QUALITY,COUNTERS,MAPS,CIGAR
 */
typedef enum {PARSE_READ, PARSE_READ__MAPS, PARSE_ALL} gt_lazy_parse_mode;

/*
 * MAP File Format test
 */
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors);

/*
 * MAP Lazy Parsers
 *   Lazy parsing works in 3 steps
 *     (1) Parses TAG(s),READ(s),QUALITIES(s),COUNTERS -> gt_template/gt_alignment
 *     (2) Parses MAPS' SEQUENCE,STRAND,POSITION
 *     (3) Parses MAPS' CIGAR string
 */
GT_INLINE void gt_input_map_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_map_parser_next_record(gt_buffered_input_file* const buffered_map_input);
// Parse Mismatches (From lazy parsing)
GT_INLINE gt_status gt_input_map_parser_parse_template_mismatch_string(gt_template* template);
GT_INLINE gt_status gt_input_map_parser_parse_alignment_mismatch_string(gt_alignment* alignment);
// Parse Maps (From lazy parsing)
GT_INLINE gt_status gt_input_map_parser_parse_template_maps(gt_template* template,uint64_t num_maps);
GT_INLINE gt_status gt_input_map_parser_parse_alignment_maps(gt_alignment* alignment,uint64_t num_maps);

/*
 * MAP High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_input_map_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template);
GT_INLINE gt_status gt_input_map_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment);
// Synchronized BMI. Attached to an output file & buffer,
// dumps the output whenever the end-of-input-block is reached
GT_INLINE gt_status gt_input_map_parser_get_template__sync_output(
    gt_buffered_input_file* const buffered_map_input,gt_template* template,gt_output_buffer** const output_buffer);
GT_INLINE gt_status gt_input_map_parser_get_alignment__sync_output(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* alignment,gt_output_buffer** const output_buffer);

#endif /* GT_INPUT_MAP_PARSER_H_ */
