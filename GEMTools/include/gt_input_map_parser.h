/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_MAP_PARSER_H_
#define GT_INPUT_MAP_PARSER_H_

#include "gt_commons.h"
#include "gt_dna_string.h"
#include "gt_dna_read.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"

#include "gt_buffered_output_file.h"

#include "gt_template.h"
#include "gt_template_utils.h"

// Codes gt_status
#define GT_IMP_OK   GT_STATUS_OK
#define GT_IMP_FAIL GT_STATUS_FAIL
#define GT_IMP_EOF  0

/*
 * Parsing error/state codes // TODO Move common errors to input_parser.h
 */
#define GT_IMP_PE_WRONG_FILE_FORMAT 10
#define GT_IMP_PE_NOT_IMPLEMENTED 11

#define GT_IMP_PE_PREMATURE_EOL 20
#define GT_IMP_PE_PENDING_BLOCKS 21
#define GT_IMP_PE_EOB 22
#define GT_IMP_PE_BAD_SEPARATOR 23
#define GT_IMP_PE_BAD_NUMBER_OF_BLOCKS 24
#define GT_IMP_PE_BAD_CHARACTER 25

#define GT_IMP_PE_READ_BAD_CHARACTER 30

#define GT_IMP_PE_QUAL_BAD_LENGTH 40
#define GT_IMP_PE_QUAL_BAD_CHARACTER 41
#define GT_IMP_PE_QUAL_BAD_SEPARATOR 42

#define GT_IMP_PE_COUNTERS_BAD_CHARACTER 50

#define GT_IMP_PE_MAP_PENDING_MAPS 60
#define GT_IMP_PE_MAP_ALREADY_PARSED 61
#define GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS 62
#define GT_IMP_PE_MAP_BAD_CHARACTER 63

#define GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER 65
#define GT_IMP_PE_SPLIT_MAP_BAD_NUM_ACCEPTORS 66
#define GT_IMP_PE_SPLIT_MAP_BAD_NUM_DONORS 67

#define GT_IMP_PE_MISMS_ALREADY_PARSED 70
#define GT_IMP_PE_MISMS_BAD_CHARACTER 71
#define GT_IMP_PE_MISMS_BAD_MISMS_POS 72

#define GT_IMP_PE_MAP_ATTR 100
#define GT_IMP_PE_MAP_GLOBAL_ATTR 101

/*
 * MAP file format constants
 */
#define GT_MAP_MCS '+'
#define GT_MAP_COUNTS_SEP ':'
#define GT_MAP_COUNTS_TIMES 'x'
#define GT_MAP_SEP ':'
#define GT_MAP_NONE '-'
#define GT_MAP_NEXT ','
#define GT_MAP_COUNTS_NOT_UNIQUE '!'
#define GT_MAP_SPLITMAP_OPEN_GEMv0 '['
#define GT_MAP_SPLITMAP_CLOSE_GEMv0 ']'
#define GT_MAP_SPLITMAP_DEF_GEMv0 '='
#define GT_MAP_SPLITMAP_SEP_GEMv0 '~'
#define GT_MAP_SPLITMAP_NEXT_GEMv0_0 '-'
#define GT_MAP_SPLITMAP_NEXT_GEMv0_1 ';'
#define GT_MAP_SCORE_GEMv0 '@'
#define GT_MAP_SCORE_SEP '/'

#define GT_MAP_SEP_S ":"
#define GT_MAP_NEXT_S ","
#define GT_MAP_MCS_S "+"
#define GT_MAP_COUNTS_TIMES_S "x"
#define GT_MAP_COUNTS_NOT_UNIQUE_S "!"
#define GT_MAP_TEMPLATE_SEP "::"
#define GT_MAP_TEMPLATE_SCORE ":::"

#define GT_MAP_STRAND_FORWARD_SYMBOL '+'
#define GT_MAP_STRAND_FORWARD_LETTER 'F'
#define GT_MAP_STRAND_REVERSE_SYMBOL '-'
#define GT_MAP_STRAND_REVERSE_LETTER 'R'

#define GT_MAP_INDEL_INSERTION '+'
#define GT_MAP_INDEL_DELETION '-'
#define GT_MAP_INDEL_SPLICE '*'

#define GT_MAP_SKIP_POSITIVE '+'
#define GT_MAP_SKIP_NEGATIVE '-'
#define GT_MAP_SKIP_SPLICE '*'

/*
 * Lazy parsing:
 *   PARSE_READ         Parses TAG,READ,QUALITY,COUNTERS
 *   PARSE_READ__MAPS   Parses TAG,READ,QUALITY,COUNTERS,MAPS
 *   PARSE_ALL          Parses TAG,READ,QUALITY,COUNTERS,MAPS,CIGAR
 */
typedef enum {PARSE_READ, PARSE_READ__MAPS, PARSE_ALL} gt_lazy_parse_mode;

/*
 * Parsing attributes
 */
typedef struct {
  bool force_read_paired; // Forces to read paired reads
  uint64_t max_parsed_maps; // Maximum number of maps to be parsed
  gt_string* src_text; // Src text line parsed
  bool skip_based_model; // Allows only mismatches & skips in the cigar string
  bool remove_duplicates; // TODO
  gt_lazy_parse_mode parse_mode;
} gt_map_parser_attr;
#define GT_MAP_PARSER_ATTR_DEFAULT(_force_read_paired) { \
  .force_read_paired=_force_read_paired,  \
  .max_parsed_maps=GT_ALL,  \
  .src_text=NULL, \
  .skip_based_model=false, \
  .remove_duplicates=false, \
  .parse_mode=PARSE_ALL \
}

GT_INLINE gt_map_parser_attr* gt_input_map_parser_attributes_new(const bool force_read_paired);
GT_INLINE void gt_input_map_parser_attributes_delete(gt_map_parser_attr* const attributes);

GT_INLINE void gt_input_map_parser_attributes_reset_defaults(gt_map_parser_attr* const attributes);
GT_INLINE bool gt_input_map_parser_attributes_is_paired(gt_map_parser_attr* const attributes);
GT_INLINE void gt_input_map_parser_attributes_set_paired(gt_map_parser_attr* const attributes,const bool force_read_paired);
GT_INLINE void gt_input_map_parser_attributes_set_max_parsed_maps(gt_map_parser_attr* const attributes,const uint64_t max_parsed_maps);
GT_INLINE void gt_input_map_parser_attributes_set_src_text(gt_map_parser_attr* const attributes,gt_string* const src_text);
GT_INLINE void gt_input_map_parser_attributes_set_skip_model(gt_map_parser_attr* const attributes,const bool skip_based_model);
GT_INLINE void gt_input_map_parser_attributes_set_duplicates_removal(gt_map_parser_attr* const attributes,const bool remove_duplicates);

/*
 * MAP File basics
 */
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors);
GT_INLINE void gt_input_map_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_map_parser_next_record(gt_buffered_input_file* const buffered_map_input);
GT_INLINE gt_status gt_input_map_parser_reload_buffer(
    gt_buffered_input_file* const buffered_map_input,const bool synchronized_map,const uint64_t num_lines);

/*
 * MAP string parsers
 */
GT_INLINE gt_status gt_input_map_parse_template(char* const string,gt_template* const template);
GT_INLINE gt_status gt_input_map_parse_alignment(char* const string,gt_alignment* const alignment);
GT_INLINE gt_status gt_input_map_parse_counters(char* const string,gt_vector* const counters,gt_shash* const attributes);
GT_INLINE gt_status gt_input_map_parse_map(char* const string,gt_map* const map);
GT_INLINE gt_status gt_input_map_parse_map_list(char* const string,gt_vector* const maps,const uint64_t num_maps);

/*
 * MAP High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_input_map_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template);
GT_INLINE gt_status gt_input_map_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment);

GT_INLINE gt_status gt_input_map_parser_get_template_g(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_map_parser_attr* const map_parser_attr);
GT_INLINE gt_status gt_input_map_parser_get_alignment_g(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_map_parser_attr* const map_parser_attr);

GT_INLINE gt_status gt_input_map_parser_synch_blocks(
    gt_buffered_input_file* const buffered_map_input1,gt_buffered_input_file* const buffered_map_input2,pthread_mutex_t* const input_mutex);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,gt_map_parser_attr* const map_parser_attr,
    const uint64_t num_map_inputs,gt_buffered_input_file* const buffered_map_input,...);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_map_input,
    const uint64_t total_files,gt_map_parser_attr* const map_parser_attr);

GT_INLINE gt_status gt_input_map_parser_synch_blocks_by_subset(
    pthread_mutex_t* const input_mutex,gt_map_parser_attr* const map_parser_attr,
    gt_buffered_input_file* const buffered_map_input_master,gt_buffered_input_file* const buffered_map_input_slave); // Used to merge files in parallel

#endif /* GT_INPUT_MAP_PARSER_H_ */
