/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Common
#include "gt_commons.h"
#include "gt_string.h"
#include "gt_dna_string.h"

// Input handlers
#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
// Input parsers
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

// GEM-Tools basic data structures: Template/Alignment/Maps/...
#include "gt_misms.h"
#include "gt_map.h"
#include "gt_dna_read.h"
#include "gt_data_attributes.h"
#include "gt_alignment.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"
#include "gt_template_utils.h"
#include "gt_counters_utils.h"
#include "gt_compact_dna_string.h"
#include "gt_sequence_archive.h"

// HighLevel Modules
#include "gt_stats.h"

GT_INLINE void gt_merge_map_files(
    pthread_mutex_t* const input_mutex,gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,gt_output_file* const output_file);

#endif /* GEM_TOOLS_H_ */
