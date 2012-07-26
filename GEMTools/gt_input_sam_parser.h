/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.h
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */


#ifndef GT_INPUT_SAM_PARSER_H_
#define GT_INPUT_SAM_PARSER_H_

#include "gt_commons.h"
#include "gt_alignment_handling.h"
#include "gt_template_handling.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"

// Codes gt_status
#define GT_ISP_OK 1
#define GT_ISP_FAIL -1
#define GT_ISP_EOF 0

/*
 * Parsing error/state codes
 */
#define GT_ISP_PE_WRONG_FILE_FORMAT 10
#define GT_ISP_PE_PREMATURE_EOL 11
#define GT_ISP_PE_EXPECTED_NUMBER 12
#define GT_ISP_PE_BAD_CHARACTER 14
#define GT_ISP_PE_WRONG_READ_CONTENT 15
/* CIGAR */
#define GT_ISP_PE_CIGAR_PREMATURE_END 20
#define GT_ISP_PE_SAM_UNMAPPED_XA 21
/* PairedEnd Parsing */
#define GT_ISP_PE_WRONG_NUM_XA 30
#define GT_ISP_PE_WRONG_END_POS 31
#define GT_ISP_PE_UNSOLVED_PENDING_MAPS 32

/*
 * SAM file format constants
 */
#define GT_SAM_HEADER_BEGIN '@'
// SAM FLAGS
#define GT_SAM_FLAG_MULTIPLE_SEGMENTS 0x1
#define GT_SAM_FLAG_PROPERLY_ALIGNED 0x2
#define GT_SAM_FLAG_UNMAPPED 0x4
#define GT_SAM_FLAG_NEXT_UNMAPPED 0x8
#define GT_SAM_FLAG_REVERSE_COMPLEMENT 0x10
#define GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT 0x20
#define GT_SAM_FLAG_FIRST_SEGMENT 0x40
#define GT_SAM_FLAG_LAST_SEGMENT 0x80
#define GT_SAM_FLAG_SECONDARY_ALIGNMENT 0x100
#define GT_SAM_FLAG_NOT_PASSING_QC 0x200
#define GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE 0x400

/*
 * SAM File Format test
 */
GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors); // TODO: Install it in gt_input_file.h

/*
 * Parsing Building Blocks
 */
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_map_input);
GT_INLINE void gt_input_sam_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template(gt_buffered_input_file* const buffered_map_input,gt_template* const template);
GT_INLINE gt_status gt_input_sam_parser_get_alignment(gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment);

#endif /* GT_INPUT_SAM_PARSER_H_ */
