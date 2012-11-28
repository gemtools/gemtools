/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_fastq_parser.c
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_fastq_parser.h"

// Constants
#define GT_INPUT_SOAP_PARSER_NUM_LINES GT_NUM_LINES_10K
#define GT_INPUT_SOAP_PARSER_NUM_INITIAL_MAPS 5

/*
 * FASTQ File basics
 */
GT_INLINE bool gt_input_file_detect_fastq_format(char* const buffer,const uint64_t buffer_size,const bool show_errors) {
  // TODO
  return true;
}
GT_INLINE bool gt_input_file_test_fastq(
    gt_input_file* const input_file,gt_fasta_file_format* const fasta_file_format,const bool show_errors) {
  GT_INPUT_FILE_CHECK(input_file);
  //GT_INPUT_FILE_CHECK__FILL_BUFFER(input_file,NULL);
  if (gt_input_file_detect_fastq_format((char*)input_file->file_buffer,input_file->buffer_size,show_errors)) {
    // TODO
    return true;
  } else {
    return false;
  }
}
GT_INLINE gt_status gt_input_fastq_parser_check_fastq_file_format(gt_buffered_input_file* const buffered_map_input) {
//  register gt_input_file* const input_file = buffered_map_input->input_file;
//  if (gt_expect_false(input_file->file_format==UNKNOWN)) { // Unknown
//    gt_sam_headers sam_headers;
//    // Mutex format detection (because the first one must read the headers)
//    gt_input_file_lock(input_file);
//      register const bool is_sam_format =
//          gt_input_file_test_sam(input_file,&sam_headers,true);
//    gt_input_file_unlock(input_file);
//    if (!is_sam_format) return GT_ISP_PE_WRONG_FILE_FORMAT;
//    input_file->file_format = SAM;
//  } else if (gt_expect_false(input_file->file_format!=SAM)) {
//    return GT_ISP_PE_WRONG_FILE_FORMAT;
//  }
//  return 0;
  return 0;
}

/*
 * Parsing Building Blocks
 */
GT_INLINE void gt_input_fastq_parser_next_record(gt_buffered_input_file* const buffered_map_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  if (!gt_buffered_input_file_eob(buffered_map_input)) {
    GT_INPUT_FILE_SKIP_LINE(buffered_map_input);
  }
}
GT_INLINE void gt_input_fastq_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
//  // Display textual error msg
//  register const char* const file_name = (buffered_map_input != NULL) ?
//      buffered_map_input->input_file->file_name : "<<LazyParsing>>";
//  if ((buffered_map_input == NULL)) {
//    line_num = 0; column_pos = 0;
//  }
//  switch (error_code) {
//    case 0: /* No error */ break;
//    default:
//      gt_error(PARSE_SAM,buffered_map_input->input_file->file_name,line_num);
//      break;
//  }
}

GT_INLINE gt_status gt_ifp_read_tag(
    char** const text_line,char** tag,uint64_t* tag_length,uint64_t* end_position) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(tag_length);
  // Read tag
  *tag = *text_line;
  GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
  if (GT_IS_EOL(text_line)) return GT_IFP_PE_PREMATURE_EOL;
  register char* last_parsed_tag_char = *text_line-1;
  *tag_length = *text_line-*tag;
  if (**text_line==SPACE) {
    GT_READ_UNTIL(text_line,**text_line==TAB);
    if (GT_IS_EOL(text_line)) return GT_IFP_PE_PREMATURE_EOL;
  }
  GT_NEXT_CHAR(text_line);
  // Parse the end information {/1,/2,...}
  if (*tag_length>2 && *(last_parsed_tag_char-1)==SLASH) {
    *(last_parsed_tag_char-1)=EOS;
    if (*last_parsed_tag_char=='1') {
      *end_position = 0;
    } else if (*last_parsed_tag_char=='2') {
      *end_position = 1;
    }
    *tag_length -= 2;
  } else {
    *end_position = UINT64_MAX;
  }
  return 0;
}

#define GT_INPUT_SOAP_PARSER_CHECK_PREMATURE_EOL() \
  if (GT_IS_EOL(text_line)) return GT_IFP_PE_PREMATURE_EOL
#define GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT() \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL(); \
  GT_NEXT_CHAR(text_line)
#define GT_ISP_PARSE_SAM_ALG_SKIP_FIELD() \
  GT_READ_UNTIL(text_line,**text_line==TAB); \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT()
#define GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(number) \
  if (!gt_is_number(**text_line)) { gt_map_delete(map); return GT_ISP_PE_EXPECTED_NUMBER; } \
  GT_PARSE_NUMBER(text_line,number)

GT_INLINE gt_status gt_ifp_parse_fastq_record(
    char** const text_line,char* const sequence,char* const qualities) {
  // TODO
  return GT_IFP_FAIL;
}

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_fastq_parser_get_sequence(
    gt_buffered_input_file* const buffered_map_input,char* const sequence,char* const qualities) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
//  register const char* line_start = buffered_map_input->cursor;
//  register gt_status error_code;
//  // Check file format
//  register gt_input_file* input_file = buffered_map_input->input_file;
//  if (gt_input_sam_parser_check_sam_file_format(buffered_map_input)) {
//    gt_error(PARSE_SAM_BAD_FILE_FORMAT,input_file->file_name,1ul);
//    return GT_ISP_FAIL;
//  }
//  // Check the end_of_block. Reload buffer if needed
//  if (gt_buffered_input_file_eob(buffered_map_input)) {
//    register const uint64_t read_lines =
//        gt_buffered_input_file_get_block(buffered_map_input,GT_NUM_LINES_10K,true);
//    if (gt_expect_false(read_lines==0)) return GT_ISP_EOF;
//  }
//  // Allocate memory for the alignment
//  register const uint64_t line_num = buffered_map_input->current_line_num;
//  gt_alignment_clear(alignment);
//  alignment->alignment_id = line_num;
//  // Parse alignment
//  if ((error_code=gt_input_sam_parser_parse_alignment(buffered_map_input,alignment))) {
//    gt_input_sam_parser_prompt_error(buffered_map_input,line_num,
//        buffered_map_input->cursor-line_start,error_code);
//    gt_input_sam_parser_next_record(buffered_map_input);
//    return GT_ISP_FAIL;
//  }
//  gt_input_sam_parser_next_record(buffered_map_input);
  return GT_IFP_FAIL;
}
GT_INLINE gt_status gt_input_fastq_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment) {
  // TODO
  return GT_IFP_FAIL;
}
GT_INLINE gt_status gt_input_fastq_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template) {
  // TODO
  return GT_IFP_FAIL;
}

/*
 * FASTQ utils
 */
GT_INLINE uint64_t gt_fastq_tag_chomp_end_info(gt_string* const tag) {
  GT_STRING_CHECK(tag);
  // Parse the end information {/1,/2}
  register const uint64_t tag_length = gt_string_get_length(tag);
  if (tag_length>2 && *gt_string_char_at(tag,tag_length-2)==SLASH) {
    register const char tag_end = *gt_string_char_at(tag,tag_length-1);
    if (tag_end=='1') {
      gt_string_set_length(tag,tag_length-2);
      return 0;
    } else if (tag_end=='2') {
      gt_string_set_length(tag,tag_length-2);
      return 1;
    } else {
      return UINT64_MAX;
    }
  } else {
    return UINT64_MAX;
  }
}
