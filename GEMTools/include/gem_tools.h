/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Essentials
#include "gt_essentials.h"
#include "gt_dna_string.h"

// Input handlers
#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
// Input parsers
#include "gt_input_parser.h"
#include "gt_input_map_parser.h"
#include "gt_input_sam_parser.h"
#include "gt_input_fasta_parser.h"
#include "gt_input_generic_parser.h"

// Output handlers
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
// Output printers (MAP,SAM,BAM,...)
#include "gt_output_fasta.h"
#include "gt_output_map.h"
#include "gt_output_sam.h"
#include "gt_output_generic_printer.h"

// GEM-Tools basic data structures: Template/Alignment/Maps/...
#include "gt_misms.h"
#include "gt_map.h"
#include "gt_dna_read.h"
#include "gt_attributes.h"
#include "gt_alignment.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"
#include "gt_template_utils.h"
#include "gt_counters_utils.h"
#include "gt_compact_dna_string.h"
#include "gt_sequence_archive.h"

// HighLevel Modules
#include "gt_stats.h"
#include "gt_gtf.h"

// GEM Idx Loader
#include "gt_gemIdx_loader.h"

// Merge functions (synch files)
#define gt_merge_synch_map_files(input_mutex,paired_end,output_file,input_map_master,input_map_slave) \
  gt_merge_synch_map_files_va(input_mutex,paired_end,output_file,input_map_master,1,input_map_slave)
GT_INLINE void gt_merge_synch_map_files_a(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file** const input_map_files,const uint64_t num_input_map_files);
GT_INLINE void gt_merge_synch_map_files_v(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,va_list v_args);
GT_INLINE void gt_merge_synch_map_files_va(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,...);

// Merge functions (unsynch files)
GT_INLINE void gt_merge_unsynch_map_files(
    pthread_mutex_t* const input_mutex,gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,gt_output_file* const output_file);

// General generic I/O loop (overkilling, but useful)
#define GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,paired_end,buffered_output,template) \
  /* Prepare IN/OUT buffers & printers */ \
  gt_status __error_code; \
  gt_buffered_input_file* __buffered_input = gt_buffered_input_file_new(input_file); \
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file); \
  gt_buffered_input_file_attach_buffered_output(__buffered_input,buffered_output); \
  /* Prepare Attributes for generic I/O */ \
  gt_generic_parser_attributes* __gparser_attr = gt_input_generic_parser_attributes_new(paired_end); \
  /* I/O Loop */ \
  gt_template* template = gt_template_new(); \
  while ((__error_code=gt_input_generic_parser_get_template(__buffered_input,template,__gparser_attr))) { \
    if (__error_code!=GT_IMP_OK) { \
      gt_error_msg("Fatal error parsing file '%s', line %"PRIu64"\n", \
          parameters.name_input_file,__buffered_input->current_line_num-1); \
      continue; \
    }
#define GT_END_READING_WRITING_LOOP(input_file,output_file,template) \
  } \
  /* Clean */ \
  gt_buffered_input_file_close(__buffered_input); \
  gt_buffered_output_file_close(buffered_output); \
  gt_input_generic_parser_attributes_delete(__gparser_attr); \
  gt_template_delete(template); \

#endif /* GEM_TOOLS_H_ */
