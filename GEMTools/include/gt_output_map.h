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
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"
#include "gt_input_map_parser.h"

/*
 * Error/state codes (Map Output Error)
 */
#define GT_MOE_INCONSISTENT_COUNTERS 10

/*
 * MAP building block printers
 */
// NOTE: Macro based definition for all printers
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_mismatch_string,gt_map* const map);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_counters,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map,gt_map* const map,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template_maps,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment_maps,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
// Sorted print by mismatch number
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template_maps_sorted,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment_maps_sorted,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);

/*
 * High-level MAP Printers
 */
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_template(gt_string* const string,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer* const output_buffer,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_ofprint_template(gt_output_file* const output_file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bofprint_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_alignment(gt_string* const string,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer* const output_buffer,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_ofprint_alignment(gt_output_file* const output_file,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bofprint_alignment(gt_buffered_output_file* const buffered_output_file,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);

/*
 * GEM printer
 */
// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_gem_template,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_gem_template(gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_fprint_gem_template(FILE* file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_sprint_gem_template(gt_string* const string,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bprint_gem_template(gt_output_buffer* const output_buffer,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_ofprint_gem_template(gt_output_file* const output_file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_bofprint_gem_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);

#endif /* GT_OUTPUT_MAP_H_ */
