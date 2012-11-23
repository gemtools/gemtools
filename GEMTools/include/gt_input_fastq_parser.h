/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_fastq_parser.h
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */


#ifndef GT_INPUT_FASTQ_PARSER_H_
#define GT_INPUT_FASTQ_PARSER_H_

#include "gt_commons.h"
#include "gt_dna_read.h"
#include "gt_sequence_archive.h"

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
 * FASTQ File basics
 */
GT_INLINE bool gt_input_file_test_fastq(
    gt_input_file* const input_file,gt_fasta_file_format* const fasta_file_format,const bool show_errors);
GT_INLINE void gt_input_fastq_parser_prompt_error(
    gt_buffered_input_file* const buffered_fastq_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_fastq_parser_next_record(gt_buffered_input_file* const buffered_fastq_input);

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_fastq_parser_get_read(
    gt_buffered_input_file* const buffered_fastq_input,gt_dna_read* const dna_read);
GT_INLINE gt_status gt_input_fastq_parser_get_segmented_sequence(
    gt_buffered_input_file* const buffered_fastq_input,gt_segmented_sequence* const segmented_sequence);
GT_INLINE gt_status gt_input_fastq_parser_get_alignment(
    gt_buffered_input_file* const buffered_fastq_input,gt_alignment* const alignment);
GT_INLINE gt_status gt_input_fastq_parser_get_template(
    gt_buffered_input_file* const buffered_fastq_input,gt_template* const template);

GT_INLINE gt_status gt_input_fastq_parser_get_archive(
    gt_buffered_input_file* const buffered_fastq_input,gt_sequence_archive* const sequence_archive);

/*
 * FASTQ utils
 */
GT_INLINE uint64_t gt_fastq_tag_chomp_end_info(gt_string* const tag);


#endif /* GT_INPUT_FASTQ_PARSER_H_ */
