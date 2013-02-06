/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
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
 * Output attributes
 */
typedef struct {
  bool print_scores; // Print alignment scores
  uint64_t max_printable_maps; // Maximum number of maps printed
  bool print_extra; // print extra information stored in attributes
  bool print_casava; // if available print casava ids, otherwise, appends /1 /2 for paird reads
} gt_output_map_attributes;

GT_INLINE gt_output_map_attributes* gt_output_map_attributes_new();
GT_INLINE void gt_output_map_attributes_delete(gt_output_map_attributes* attributes);
GT_INLINE void gt_output_map_attributes_reset_defaults(gt_output_map_attributes* const attributes);

GT_INLINE bool gt_output_map_attributes_is_print_scores(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_scores(gt_output_map_attributes* const attributes, bool print_scores);

GT_INLINE bool gt_output_map_attributes_is_print_extra(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_extra(gt_output_map_attributes* const attributes, bool print_extra);

GT_INLINE bool gt_output_map_attributes_is_print_casava(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_casava(gt_output_map_attributes* const attributes, bool print_casava);

GT_INLINE uint64_t gt_output_map_attributes_get_max_printable_maps(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_max_printable_maps(gt_output_map_attributes* const attributes, uint64_t max_printable_maps);


/*
 * TAG building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_tag,gt_string* const tag,gt_shash* const attributes,gt_output_map_attributes* const output_attributes);

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
GT_INLINE gt_status gt_output_map_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_sprint_template(gt_string* const string,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer* const output_buffer,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_ofprint_template(gt_output_file* const output_file,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bofprint_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_sprint_alignment(gt_string* const string,gt_alignment* const alignment,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer* const output_buffer,gt_alignment* const alignment,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_ofprint_alignment(gt_output_file* const output_file,gt_alignment* const alignment,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bofprint_alignment(gt_buffered_output_file* const buffered_output_file,gt_alignment* const alignment,gt_output_map_attributes* const attributes);

/*
 * GEM printer
 */
// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_gem_template,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_gprint_gem_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_fprint_gem_template(FILE* file,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_sprint_gem_template(gt_string* const string,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bprint_gem_template(gt_output_buffer* const output_buffer,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_ofprint_gem_template(gt_output_file* const output_file,gt_template* const template,gt_output_map_attributes* const attributes);
GT_INLINE gt_status gt_output_map_bofprint_gem_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,gt_output_map_attributes* const attributes);

#endif /* GT_OUTPUT_MAP_H_ */
