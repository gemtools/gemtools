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
GT_INLINE gt_status gt_merge_map_read_template_sync(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_map_parser_attr* const map_parser_attr,gt_template* const template_master,gt_template* const template_slave,
    gt_buffered_output_file* buffered_output) {
  register gt_status error_code_master, error_code_slave;
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  do {
    // Read Synch blocks
    error_code_master=gt_input_map_parser_synch_blocks_by_subset(
        input_mutex,map_parser_attr,buffered_input_master,buffered_input_slave);
    if (error_code_master==GT_IMP_EOF) return GT_IMP_EOF;
    if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error synchronizing files");
    // Read master (always guaranteed)
    if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))==GT_IMP_FAIL) {
      gt_fatal_error_msg("Fatal error parsing file <<Master>>");
    }
    // Check slave
    if (gt_buffered_input_file_eob(buffered_input_slave)) { // Slave exhausted. Dump master & return EOF
      gt_output_map_bofprint_gem_template(buffered_output,template_master,output_attributes);
      while ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))) {
        if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        gt_output_map_bofprint_gem_template(buffered_output,template_master,output_attributes);
      }
    } else {
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_slave,template_slave))==GT_IMP_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
      }
      // Synch loop
      while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
        // Print non correlative master's template
        gt_output_map_bofprint_gem_template(buffered_output,template_master,output_attributes);
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
    pthread_mutex_t* const input_mutex,
    gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,const bool files_contain_same_reads,gt_output_file* const output_file) {
  // Reading+process
  gt_buffered_input_file* buffered_input_master = gt_buffered_input_file_new(input_map_master);
  gt_buffered_input_file* buffered_input_slave  = gt_buffered_input_file_new(input_map_slave);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_master,buffered_output);

  gt_template *template_master = gt_template_new();
  gt_template *template_slave = gt_template_new();

  gt_map_parser_attr map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(paired_end);
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();

  if (files_contain_same_reads) {
    register gt_status error_code_master, error_code_slave;
    while (gt_input_map_parser_synch_blocks_va(input_mutex,&map_parser_attr,2,buffered_input_master,buffered_input_slave)) {
      // Read master (always guaranteed)
      if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master))==GT_IMP_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Master>>");
      }
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_slave,template_slave))==GT_IMP_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
      }
      if (error_code_master!=error_code_slave) {
        gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
      }
      // Check tags
      if (!gt_string_equals(template_master->tag,template_slave->tag)) {
        gt_fatal_error_msg("Files are not synchronized. Different TAGs found '"PRIgts"' '"PRIgts"' ",
            PRIgts_content(template_master->tag),PRIgts_content(template_slave->tag));
      }
      // Merge maps
      register gt_template *ptemplate = gt_template_union_template_mmaps(template_master,template_slave);
      // Print template
      gt_output_map_bofprint_gem_template(buffered_output,ptemplate,output_attributes);
      // Delete template
      gt_template_delete(ptemplate);
    }
  } else {
    while (gt_merge_map_read_template_sync(input_mutex,buffered_input_master,buffered_input_slave,&map_parser_attr,
        template_master,template_slave,buffered_output)) {
      // Merge maps
      register gt_template *ptemplate = gt_template_union_template_mmaps(template_master,template_slave);
      // Print template
      gt_output_map_bofprint_gem_template(buffered_output,ptemplate,output_attributes);
      // Delete template
      gt_template_delete(ptemplate);
    }
  }

  // Clean
  gt_template_delete(template_master);
  gt_template_delete(template_slave);
  gt_buffered_input_file_close(buffered_input_master);
  gt_buffered_input_file_close(buffered_input_slave);
  gt_buffered_output_file_close(buffered_output);
}
