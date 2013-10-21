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
// Input parsers/utils
#include "gt_input_parser.h"
#include "gt_input_map_parser.h"
#include "gt_input_map_utils.h"
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

// Utilities
#include "gt_json.h"

// GEM Idx Loader
#include "gt_gemIdx_loader.h"

/*
 * General generic I/O loop (overkilling, but useful)
 */
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
          input_file->file_name,__buffered_input->current_line_num-1); \
      continue; \
    }
#define GT_END_READING_WRITING_LOOP(input_file,output_file,template) \
  } \
  /* Clean */ \
  gt_buffered_input_file_close(__buffered_input); \
  gt_buffered_output_file_close(buffered_output); \
  gt_input_generic_parser_attributes_delete(__gparser_attr); \
  gt_template_delete(template)

/*
 * Options (Tools Menu)
 */
typedef enum { GT_OPT_NO_ARGUMENT=no_argument, GT_OPT_REQUIRED=required_argument, GT_OPT_OPTIONAL=optional_argument } gt_option_t;
typedef enum { GT_OPT_NONE, GT_OPT_INT, GT_OPT_FLOAT, GT_OPT_CHAR, GT_OPT_STRING, GT_OPT_BOOL } gt_option_argument_t;
typedef struct {
  int option_id;       // Integer ID or short character option
  char* long_option;   // Long string option
  gt_option_t option_type;            // Option type
  gt_option_argument_t argument_type; // Type of the argument
  uint64_t group_id;   // Label of the group it belongs to (zero if none)
  bool active;         // Enable/Disable option
  char* command_info;  // Extra command line syntax info
  char* description;   // Brief description
} gt_option;

extern gt_option gt_filter_options[];
extern char* gt_filter_groups[];

extern gt_option gt_stats_options[];
extern char* gt_stats_groups[];

extern gt_option gt_mapset_options[];
extern char* gt_mapset_groups[];

extern gt_option gt_scorereads_options[];
extern char* gt_scorereads_groups[];
extern gt_sam_attribute_option gt_scorereads_attribute_option_list[];

extern gt_option gt_map2sam_options[];
extern char* gt_map2sam_groups[];

extern gt_option gt_gtfcount_options[];
extern char* gt_gtfcount_groups[];

extern gt_option gt_region_options[];
extern char* gt_region_groups[];

GT_INLINE uint64_t gt_options_get_num_options(const gt_option* const options);
GT_INLINE struct option* gt_options_adaptor_getopt(const gt_option* const options);
GT_INLINE gt_string* gt_options_adaptor_getopt_short(const gt_option* const options);
GT_INLINE void gt_options_fprint_menu(
    FILE* const stream,const gt_option* const options,char* gt_filter_groups[],
    const bool print_description,const bool print_inactive);
GT_INLINE void gt_options_fprint_json_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive);

#endif /* GEM_TOOLS_H_ */
