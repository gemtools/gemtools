/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_generic_parser.h
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic parser for {MAP,SAM,FASTQ}
 */


#ifndef GT_INPUT_GENERIC_PARSER_H_
#define GT_INPUT_GENERIC_PARSER_H_

#include "gt_commons.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"

#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"
#include "gt_input_map_parser.h"
#include "gt_input_sam_parser.h"

/*
 * Parsers Helpers
 */
#define GT_IGP_FAIL -1
#define GT_IGP_EOF 0
#define GT_IGP_OK 1

typedef struct {
  /* General attributes*/
  bool paired_read;
  /* SAM specific */
  gt_sam_parser_attr sam_parser_attr;
  /* MAP specific */
  uint64_t max_matches;
} gt_generic_parser_attr;


#define GENERIC_PARSER_ATTR_DEFAULT(parse_paired) { .sam_parser_attr.sam_soap_style=false, .max_matches=GT_ALL, .paired_read=parse_paired }

GT_INLINE gt_generic_parser_attr* gt_generic_parser_attr_new(bool const sam_soap_style, uint64_t const max_matches, bool const paired_read);

GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment,gt_generic_parser_attr* const attributes);
GT_INLINE gt_status gt_input_generic_parser_get_template(
    gt_buffered_input_file* const buffered_input,gt_template* const template,gt_generic_parser_attr* const attributes);


#endif /* GT_INPUT_GENERIC_PARSER_H_ */
