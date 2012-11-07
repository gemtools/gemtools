/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

/*
 * Parsers Helpers
 */
GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment) {
  gt_status error_code;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      if ((error_code=gt_input_map_parser_get_alignment(buffered_input,alignment))!=GT_IMP_OK) {
        return (error_code==GT_IMP_EOF) ? GT_IGP_EOF : GT_IGP_FAIL;
      }
      break;
    case SAM:
      // TODO
      break;
    default:
      return GT_IGP_FAIL;
      break;
  }
  return error_code;
}
GT_INLINE gt_status gt_input_generic_parser_get_template(
    gt_buffered_input_file* const buffered_input,gt_template* const template,const bool preserve_pairness) {
  gt_status error_code;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      if ((error_code=gt_input_map_parser_get_template(buffered_input,template))!=GT_IMP_OK) {
        return (error_code==GT_IMP_EOF) ? GT_IGP_EOF : GT_IGP_FAIL;
      }
      if (gt_template_get_num_blocks(template)==1 && preserve_pairness) {
        if ((error_code=gt_input_map_parser_get_alignment(
            buffered_input,gt_template_get_block_dyn(template,1)))!=GT_IMP_OK) {
          return GT_IGP_FAIL;
        }
      }
      break;
    case SAM:
      // TODO
      break;
    default:
      return GT_IGP_FAIL;
      break;
  }
  return error_code;
}

/*
 * Counters Helpers
 */
GT_INLINE uint64_t gt_calculate_num_maps(
    const uint64_t num_decoded_strata,const uint64_t num_decoded_matches,
    const uint64_t first_stratum_threshold) {
  // TODO
  return 0;
}
