/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_fastq_parser.h
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */


#ifndef GT_INPUT_FASTQ_PARSER_H_
#define GT_INPUT_FASTQ_PARSER_H_

#include "gt_commons.h"
#include "gt_alignment_handling.h"
#include "gt_template_handling.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"

// Codes gt_status
#define GT_IFP_OK 1
#define GT_IFP_FAIL -1
#define GT_IFP_EOF 0

/*
 * Parsing error/state codes
 */
// #define GT_IFP_PE_WRONG_FILE_FORMAT 10

/*
 * FASTQ File Format test
 */
GT_INLINE bool gt_input_file_test_fastq(
    gt_input_file* const input_file,,const bool show_errors);

/*
 * Parsing Building Blocks
 */
GT_INLINE void gt_input_fastq_parser_next_record(gt_buffered_input_file* const buffered_map_input);
GT_INLINE void gt_input_fastq_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_fastq_parser_get_sequence(
    gt_buffered_input_file* const buffered_map_input,char* const sequence,char* const qualities);
GT_INLINE gt_status gt_input_fastq_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment);
GT_INLINE gt_status gt_input_fastq_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,const uint64_t num_blocks);

#endif /* GT_INPUT_FASTQ_PARSER_H_ */
