/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_MAP_H_
#define GT_OUTPUT_MAP_H_

#include "gt_commons.h"
#include "gt_template.h"
#include "gt_iterators.h"
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"

/*
 * MAP building block printers
 */
// Generic Printer [G]
GT_INLINE gt_status gt_output_map_gprint_counters(gt_generic_printer* const gprinter,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_gprint_map(gt_generic_printer* const gprinter,gt_map* const map,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_template_maps(gt_generic_printer* const gprinter,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps(gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// Output Buffer [B]
GT_INLINE gt_status gt_output_map_bprint_counters(gt_output_buffer* const output_buffer,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_bprint_map(gt_output_buffer* const output_buffer,gt_map* const map,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_template_maps(gt_output_buffer* const output_buffer,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_alignment_maps(gt_output_buffer* const output_buffer,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// String [S]
GT_INLINE gt_status gt_output_map_sprint_counters(char **line_ptr,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_sprint_map(char **line_ptr,gt_map* const map,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_template_maps(char **line_ptr,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_alignment_maps(char **line_ptr,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// File [F]
GT_INLINE gt_status gt_output_map_fprint_counters(FILE* file,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_fprint_map(FILE* file,gt_map* const map,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_template_maps(FILE* file,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_alignment_maps(FILE* file,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);

/*
 * High-level MAP Printers
 */
// Generic Printer [G]
GT_INLINE gt_status gt_output_map_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// Output Buffer [B]
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer* const output_buffer,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer* const output_buffer,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// String [S]
GT_INLINE gt_status gt_output_map_sprint_template(char **line_ptr,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_alignment(char **line_ptr,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);
// File [F]
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);

#endif /* GT_OUTPUT_MAP_H_ */
