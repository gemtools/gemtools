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
 * MAP building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_mismatch_string,gt_map* const map);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_counters,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map,gt_map* const map,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template_maps,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment_maps,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);

/*
 * High-level MAP Printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template,gt_template* const template,const uint64_t num_maps,const bool print_scores);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores);

/*
 * GEM printer
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_gem_template,gt_template* const template,const uint64_t num_maps,const bool print_scores);

#endif /* GT_OUTPUT_MAP_H_ */
