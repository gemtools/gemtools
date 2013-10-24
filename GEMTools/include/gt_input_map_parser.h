/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_MAP_PARSER_H_
#define GT_INPUT_MAP_PARSER_H_

#include "gt_essentials.h"

#include "gt_dna_string.h"
#include "gt_dna_read.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"

#include "gt_buffered_output_file.h"

#include "gt_template.h"
#include "gt_template_utils.h"

/*
 * Codes gt_status
 */
#define GT_IMP_OK   GT_STATUS_OK
#define GT_IMP_FAIL GT_STATUS_FAIL
#define GT_IMP_EOF  0

/*
 * Parsing error/state codes
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
#define GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS 62
#define GT_IMP_PE_MAP_INCONSISTENT_BLOCKS 63
#define GT_IMP_PE_MAP_BAD_CHARACTER 64

#define GT_IMP_PE_MISMS_BAD_CHARACTER 70
#define GT_IMP_PE_MISMS_BAD_MISMS_POS 71

#define GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv0 100
#define GT_IMP_PE_MAP_ATTRIBUTE_SCORE_GEMv1 101
#define GT_IMP_PE_MMAP_ATTRIBUTE_SCORE 102

#define GT_IMP_NUM_LINES          GT_NUM_LINES_20K
#define GT_IMP_SUBSET_NUM_LINES   GT_IMP_NUM_LINES

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
#define GT_MAP_NONE_S "-"
#define GT_MAP_NEXT_S ","
#define GT_MAP_MCS_S "+"
#define GT_MAP_COUNTS_TIMES_S "x"
#define GT_MAP_COUNTS_NOT_UNIQUE_S "!"
#define GT_MAP_TEMPLATE_SEP "::"
#define GT_MAP_SCORE ":"
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

typedef enum {MISMATCH_STRING_GEMv0, MISMATCH_STRING_GEMv1} gt_misms_string_t; // Mismatch string format

/*
 * Map Parser Attributes
 */
typedef struct {
  /* PE/SE */
  bool force_read_paired; // Forces to read paired reads (2 lines of MAP format in case of unmapped)
  /* Map parsing */
  uint64_t max_parsed_maps; // Maximum number of maps to be parsed
  bool skip_based_model; // Allows only mismatches & skips in the cigar string
  bool remove_duplicates; // Instead of strictly parse the record, tries to merge duplicates (sort of cleanup in case of bugs ...)
  /* Auxiliary Buffers */
  gt_string* src_text; // Source text line parsed (parsing from file)
} gt_map_parser_attributes;
#define GT_MAP_PARSER_ATTR_DEFAULT(_force_read_paired) { \
  /* PE/SE */ \
  .force_read_paired=_force_read_paired,  \
  /* Map parsing */ \
  .max_parsed_maps=GT_ALL,  \
  .skip_based_model=false, \
  .remove_duplicates=false, \
  /* Auxiliary Buffers */ \
  .src_text=NULL, \
}
#define GT_MAP_PARSER_CHECK_ATTRIBUTES(attributes) \
  gt_map_parser_attributes __##attributes; \
  if (attributes==NULL) { /* Check null map_parser_attr */ \
    gt_input_map_parser_attributes_reset_defaults(&__##attributes); \
    attributes = &__##attributes; \
  }
#define GT_MAP_PARSER_MAP_BLOCKS_INITIAL_ELEMENTS 20

GT_INLINE gt_map_parser_attributes* gt_input_map_parser_attributes_new(const bool force_read_paired);
GT_INLINE void gt_input_map_parser_attributes_delete(gt_map_parser_attributes* const attributes);

GT_INLINE void gt_input_map_parser_attributes_reset_defaults(gt_map_parser_attributes* const attributes);
GT_INLINE bool gt_input_map_parser_attributes_is_paired(gt_map_parser_attributes* const attributes);
GT_INLINE void gt_input_map_parser_attributes_set_paired(gt_map_parser_attributes* const attributes,const bool force_read_paired);
GT_INLINE void gt_input_map_parser_attributes_set_max_parsed_maps(gt_map_parser_attributes* const attributes,const uint64_t max_parsed_maps);
GT_INLINE void gt_input_map_parser_attributes_set_src_text(gt_map_parser_attributes* const attributes,gt_string* const src_text);
GT_INLINE void gt_input_map_parser_attributes_set_skip_model(gt_map_parser_attributes* const attributes,const bool skip_based_model);
GT_INLINE void gt_input_map_parser_attributes_set_duplicates_removal(gt_map_parser_attributes* const attributes,const bool remove_duplicates);

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
GT_INLINE gt_status gt_input_map_parse_counters(const char* const string,gt_vector* const counters,gt_attributes* attributes);
GT_INLINE gt_status gt_input_map_parse_map(const char* const string,gt_map** const map,gt_map_parser_attributes* map_parser_attr);
GT_INLINE gt_status gt_input_map_parse_map_list(const char* const string,gt_vector* const maps,gt_map_parser_attributes* map_parser_attr);
GT_INLINE gt_status gt_input_map_parse_alignment(const char* const string,gt_alignment* const alignment);
GT_INLINE gt_status gt_input_map_parse_template(const char* const string,gt_template* const template);

/*
 * MAP High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_input_map_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_map_parser_attributes* map_parser_attr);
GT_INLINE gt_status gt_input_map_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_map_parser_attributes* map_parser_attr);

/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_map_parser_synch_blocks(
    gt_buffered_input_file* const buffered_input1,gt_buffered_input_file* const buffered_input2,pthread_mutex_t* const input_mutex);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    uint64_t num_inputs,gt_buffered_input_file* const buffered_input,va_list v_args);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,
    const uint64_t num_inputs,gt_map_parser_attributes* const map_parser_attr);
GT_INLINE gt_status gt_input_map_parser_synch_blocks_by_subset(
    pthread_mutex_t* const input_mutex,gt_map_parser_attributes* const map_parser_attr,
    gt_buffered_input_file* const buffered_map_input_master,gt_buffered_input_file* const buffered_map_input_slave); // Used to merge files in parallel
GT_INLINE gt_status gt_input_map_parser_synch_get_template(
		gt_buffered_input_file *buffered_input_file1,gt_buffered_input_file *buffered_input_file2,gt_template *template,pthread_mutex_t *mutex);

#endif /* GT_INPUT_MAP_PARSER_H_ */
