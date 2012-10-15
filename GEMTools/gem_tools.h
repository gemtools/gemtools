/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Common
#include "gt_commons.h"

// Input handlers
#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
// Input parsers
#include "gt_input_map_parser.h"
#include "gt_input_sam_parser.h"

// Output handlers
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
// Output printers (MAP,SAM,BAM,...)
#include "gt_output_map.h"

// GEM-Tools basic data structures: Template/Alignment/Maps/...
#include "gt_misms.h"
#include "gt_map.h"
#include "gt_alignment.h"
#include "gt_alignment_handling.h"
#include "gt_template.h"
#include "gt_template_handling.h"
// Iterators
#include "gt_iterators.h"

/*
 * Counters Helpers
 */
GT_INLINE uint64_t gt_calculate_num_maps(  // fixme
    const uint64_t num_decoded_strata,const uint64_t num_decoded_matches,
    const uint64_t first_stratum_threshold);


#endif /* GEM_TOOLS_H_ */
