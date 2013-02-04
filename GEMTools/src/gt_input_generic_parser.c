/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_generic_parser.c
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic parser for {MAP,SAM,FASTQ}
 */

#include "gt_input_generic_parser.h"

GT_INLINE gt_generic_parser_attr* gt_generic_parser_attr_new(bool const sam_soap_style, uint64_t const max_matches, bool const paired_read){
  gt_generic_parser_attr* attr = malloc(sizeof(gt_generic_parser_attr));
  gt_cond_fatal_error(!attr,MEM_HANDLER);
  attr->max_matches = max_matches;
  attr->paired_read = paired_read;
  attr->sam_parser_attr = *(gt_sam_parser_attr_new(sam_soap_style));
  return attr;
}
/*
 * Accessors
 */
GT_INLINE void gt_input_generic_parser_attributes_reset_defaults(gt_generic_parser_attr* const attributes) {
  gt_input_map_parser_attributes_reset_defaults(&attributes->map_parser_attr);
  gt_input_sam_parser_attributes_reset_defaults(&attributes->sam_parser_attr);
}
GT_INLINE bool gt_input_generic_parser_attributes_is_paired(gt_generic_parser_attr* const attributes) {
  return attributes->map_parser_attr.read_paired;
}
GT_INLINE void gt_input_generic_parser_attributes_set_paired(gt_generic_parser_attr* const attributes,const bool is_paired) {
  gt_input_map_parser_attributes_set_paired(&attributes->map_parser_attr,is_paired);
}

/*
 * Parsers Helpers
 */
GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment,gt_generic_parser_attr* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      return gt_input_map_parser_get_alignment_g(buffered_input,alignment,&attributes->map_parser_attr);
      break;
    case SAM:
      return gt_input_sam_parser_get_alignment(buffered_input,alignment,&attributes->sam_parser_attr);
      break;
    case FASTA:
      return gt_input_fasta_parser_get_alignment(buffered_input,alignment);
      break;
    default:
      gt_fatal_error_msg("File type not supported");
      break;
  }
  return error_code;
}
GT_INLINE gt_status gt_input_generic_parser_get_template(
    gt_buffered_input_file* const buffered_input,gt_template* const template,gt_generic_parser_attr* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      return gt_input_map_parser_get_template_g(buffered_input,template,&attributes->map_parser_attr);
      break;
    case SAM:
      if (gt_input_generic_parser_attributes_is_paired(attributes)) {
        error_code = gt_input_sam_parser_get_template(buffered_input,template,&attributes->sam_parser_attr);
        gt_template_get_block_dyn(template,0);
        gt_template_get_block_dyn(template,1); // Make sure is a template
        return error_code;
      } else {
        return gt_input_sam_parser_get_alignment(
            buffered_input,gt_template_get_block_dyn(template,0),&attributes->sam_parser_attr);
      }
      break;
    case FASTA:
      return gt_input_fasta_parser_get_template(buffered_input,template,gt_input_generic_parser_attributes_is_paired(attributes));
      break;
    default:
      gt_fatal_error_msg("File type not supported");
      break;
  }
  return error_code;
}
