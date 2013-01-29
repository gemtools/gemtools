/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

// OMP
#include "omp.h"

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

GT_INLINE gt_status gt_merge_map_read_template_sync(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file* const buffered_input_master,
    gt_buffered_input_file* const buffered_input_slave,const bool paired_end,
    gt_buffered_output_file* buffered_output,gt_template* const template_master,gt_template* const template_slave) {
  register gt_status error_code_master, error_code_slave;
  do {
    // Read Synch blocks
    error_code_master=gt_input_map_parser_synch_blocks_by_subset(
        input_mutex,buffered_input_master,buffered_input_slave,paired_end);
    if (error_code_master==GT_IMP_EOF) return GT_IMP_EOF;
    if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error synchronizing files");
    // Read master (always guaranteed)
    if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))==GT_IMP_FAIL) {
      gt_fatal_error_msg("Fatal error parsing file <<Master>>");
    }
    // Check slave
    if (gt_buffered_input_file_eob(buffered_input_slave)) { // Slave exhausted. Dump master & return EOF
      gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
      while ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))) {
        if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
      }
    } else {
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_slave,template_slave))==GT_IMP_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
      }
      // Synch loop
      while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
        // Print non correlative master's template
        gt_output_map_bofprint_gem_template(buffered_output,template_master,GT_ALL,true);
        // Fetch next master's template
        if (gt_buffered_input_file_eob(buffered_input_master)) gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
        if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))!=GT_IMP_OK) {
          gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        }
      }
      return GT_IMP_OK;
    }
  } while (true);
}

GT_INLINE void gt_merge_map_files(
    pthread_mutex_t* const input_mutex,gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,gt_output_file* const output_file) {
  // Reading+process
  gt_buffered_input_file* buffered_input_master = gt_buffered_input_file_new(input_map_master);
  gt_buffered_input_file* buffered_input_slave  = gt_buffered_input_file_new(input_map_slave);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_master,buffered_output);

  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  while (gt_merge_map_read_template_sync(input_mutex,buffered_input_master,buffered_input_slave,
      paired_end,buffered_output,template_1,template_2)) {
    // Merge maps
    register gt_template *ptemplate = gt_template_union_template_mmaps(template_1,template_2);
    // Print template
    gt_output_map_bofprint_gem_template(buffered_output,ptemplate,GT_ALL,true);
    // Delete template
    gt_template_delete(ptemplate);
  }

  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_master);
  gt_buffered_input_file_close(buffered_input_slave);
  gt_buffered_output_file_close(buffered_output);
}
