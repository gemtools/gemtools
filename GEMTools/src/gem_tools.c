/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

/*
 * Parsers Helpers
 */
GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment,gt_generic_parser_attr* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      return gt_input_map_parser_get_alignment_limited(buffered_input,alignment,attributes->max_matches);
      break;
    case SAM:
      return gt_input_sam_parser_get_alignment(buffered_input,alignment);
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
      if ((error_code=gt_input_map_parser_get_template_limited(buffered_input,template,attributes->max_matches))!=GT_IMP_OK) {
        return (error_code==GT_IMP_EOF) ? GT_IGP_EOF : GT_IGP_FAIL;
      }
      if (gt_template_get_num_blocks(template)==1 && attributes->paired_read) {
        if ((error_code=gt_input_map_parser_get_alignment_limited(
            buffered_input,gt_template_get_block_dyn(template,1),attributes->max_matches))!=GT_IMP_OK) {
          return GT_IGP_FAIL;
        }
      }
      // TODO: Tag check consistency
      break;
    case SAM:
      if (attributes->paired_read) {
        error_code = (!attributes->sam_soap_style) ?
            gt_input_sam_parser_get_template(buffered_input,template) :
            gt_input_sam_parser_get_soap_template(buffered_input,template);
        gt_template_get_block_dyn(template,0);
        gt_template_get_block_dyn(template,1); // Make sure is a template
        return error_code;
      } else {
        return gt_input_sam_parser_get_alignment(buffered_input,gt_template_get_block_dyn(template,0));
      }
      break;
    default:
      gt_fatal_error_msg("File type not supported");
      break;
  }
  return error_code;
}

