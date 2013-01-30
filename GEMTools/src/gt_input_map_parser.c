/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_map_parser.h"

#define GT_IMP_NUM_LINES          GT_NUM_LINES_20K
#define GT_IMP_SUBSET_NUM_LINES   GT_NUM_LINES_20K
#define GT_IMP_NUM_INIT_TAG_CHARS 200

/*
 * Useful macros
 */
// MAP/MAPQ related
#define gt_is_valid_template_separator(character) \
  (character==' ')
#define gt_is_valid_counter_separator(character) \
  (character=='+' || character==':' || character=='x')

/*
 * MAP File Format test
 */
GT_INLINE bool gt_input_map_parser_test_map(
    char* const file_name,const uint64_t line_num,char* const buffer,const uint64_t buffer_size,
    gt_map_file_format* const map_file_format,const bool show_errors) {
  // Count tabs
  register uint64_t buffer_pos=0, num_tabs=0;
  register int64_t begin_f2=-1, end_f2=-1;
  register int64_t begin_f3=-1, end_f3=-1;
  register int64_t begin_f4=-1, end_f4=-1;
  while (buffer_pos<buffer_size) {
    register const char c = buffer[buffer_pos];
    // Check TAB/EOL
    if (gt_expect_false(c==TAB)) {
      if (num_tabs==0) {begin_f2=buffer_pos+1;}
      if (num_tabs==1) {end_f2=buffer_pos; begin_f3=buffer_pos+1;}
      if (num_tabs==2) {end_f3=buffer_pos; begin_f4=buffer_pos+1;}
      if (num_tabs==3) {end_f4=buffer_pos;}
      ++num_tabs;
    } else if (gt_expect_false(c==EOL)) {
      break;
    }
    ++buffer_pos;
  }
  // Check MAP file
  //   MAP:   TAG\tREAD\tCOUNTERS\tMAPS
  //   MAPQ:  TAG\tREAD\tQUALS\tCOUNTERS\tMAPS
  // Required conditions:
  //   (1) 3|4 TABS
  //   (2) length(read)==length(quals)
  if (num_tabs!=3 && num_tabs!=4) {
    gt_cond_error(show_errors,PARSE_MAP_BAD_NUMBER_FIELDS,file_name,line_num,num_tabs);
    return false;
  } else if (gt_expect_false(num_tabs==4 && (end_f2-begin_f2)!=(end_f3-begin_f3))) {
    gt_cond_error(show_errors,PARSE_MAP_BAD_READ_QUAL_LENGTH,
        file_name,line_num,end_f2-begin_f2,end_f3-begin_f3);
    return false;
  }
  // Set MAP type {MAP,MAPQ}
  register const bool contains_qualities = (num_tabs==4);
  // Check counters format
  register const uint64_t begin_counters = (!contains_qualities) ? begin_f3 : begin_f4;
  register const uint64_t end_counters = (!contains_qualities) ? end_f3 : end_f4;
  register bool is_prev_number = false;
  for (buffer_pos=begin_counters;buffer_pos<end_counters;++buffer_pos) {
    register const char c = buffer[buffer_pos];
    if (gt_is_number(c)) {
      is_prev_number = true;
    } else {
      if (gt_expect_false(!is_prev_number || !gt_is_valid_counter_separator(c)) ) {
        gt_cond_error(show_errors,PARSE_MAP_COUNTERS,file_name,line_num,buffer_pos+1);
        return false;
      }
      is_prev_number = false;
    }
  }
  // Extract extra MAP format information
  register uint64_t num_blocks=1;
  register char block_separator = 0;
  // Detect read template blocks
  register const uint64_t f3_f4_diff = (begin_f4-begin_f3);
  for (buffer_pos=begin_f2;buffer_pos<end_f2;++buffer_pos) {
    register const char c = buffer[buffer_pos];
    if (!gt_is_dna(c)) { // Check template block separator
      if (gt_expect_false(!gt_is_valid_template_separator(c))) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not a valid separator");
        return false;
      } else if (gt_expect_false(block_separator!=0 && c!=block_separator)) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not consistent with previous separator");
        return false;
      } else if (gt_expect_false(contains_qualities && c!=buffer[buffer_pos+f3_f4_diff])) {
        gt_cond_error(show_errors,PARSE_MAP_BAD_TEMPLATE_SEP,
            file_name,line_num,buffer_pos+1,c,"Not synchronized with qualities separator");
        return false;
      }
      block_separator = c;
      ++num_blocks;
    }
  }
  // Set extra format attributes
  map_file_format->contains_qualities = contains_qualities;
  // TODO: IF qualities check (num_blocks(read==qualities)) ... and use DIFF_TEMPLATE_BLOCKS
  return true;
}
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors) {
  return gt_input_map_parser_test_map(
      input_file->file_name,input_file->processed_lines+1,(char*)input_file->file_buffer,
      input_file->buffer_size,map_file_format,show_errors);
}
GT_INLINE gt_status gt_input_map_parser_check_map_file_format(gt_buffered_input_file* const buffered_input_file) {
  register gt_input_file* const input_file = buffered_input_file->input_file;
  if (gt_expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    gt_map_file_format map_type;
    if (!gt_input_map_parser_test_map(
        input_file->file_name,buffered_input_file->current_line_num,
        gt_vector_get_mem(buffered_input_file->block_buffer,char),
        gt_vector_get_used(buffered_input_file->block_buffer),&map_type,true)) {
      return GT_IMP_PE_WRONG_FILE_FORMAT;
    }
    input_file->file_format = MAP;
    input_file->map_type = map_type;
  } else if (gt_expect_false(input_file->file_format!=MAP)) {
    return GT_IMP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * MAP File basics
 */
/* Error handler */
GT_INLINE void gt_input_map_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  register const char* const file_name = (buffered_map_input != NULL) ?
      buffered_map_input->input_file->file_name : "<<LazyParsing>>";
  if ((buffered_map_input == NULL)) {
    line_num = 0; column_pos = 0;
  }
  switch (error_code) {
    case 0: /* No error */ break;
    case GT_IMP_PE_WRONG_FILE_FORMAT: gt_error(PARSE_MAP_BAD_FILE_FORMAT,file_name,line_num); break;
    case GT_IMP_PE_NOT_IMPLEMENTED: gt_error(PARSE_MAP_NOT_IMPLEMENTED,file_name,line_num,column_pos); break;
    case GT_IMP_PE_PREMATURE_EOL: gt_error(PARSE_MAP_PREMATURE_EOL,file_name,line_num,column_pos); break;
    case GT_IMP_PE_BAD_NUMBER_OF_BLOCKS: /* (blocks(read)!=blocks(quals) || parse_alignment=>(num_blocks==1)) */
      gt_error(PARSE_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos);
      break;
    case GT_IMP_PE_BAD_CHARACTER: gt_error(PARSE_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_READ_BAD_CHARACTER: gt_error(PARSE_MAP_READ_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_SEPARATOR: gt_error(PARSE_MAP_QUAL_BAD_SEPARATOR,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_LENGTH: gt_error(PARSE_MAP_QUAL_BAD_LENGTH,file_name,line_num,column_pos); break;
    case GT_IMP_PE_QUAL_BAD_CHARACTER: gt_error(PARSE_MAP_QUAL_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_COUNTERS_BAD_CHARACTER: gt_error(PARSE_MAP_COUNTERS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MAP_ALREADY_PARSED: gt_error(PARSE_MAP_MAP_ALREADY_PARSED,file_name,line_num); break;
    case GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS: gt_error(PARSE_MAP_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MAP_BAD_CHARACTER: gt_error(PARSE_MAP_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER: gt_error(PARSE_MAP_SPLIT_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_SPLIT_MAP_BAD_NUM_ACCEPTORS: gt_error(PARSE_MAP_SPLIT_MAP_BAD_NUM_ACCEPTORS,file_name,line_num,column_pos); break;
    case GT_IMP_PE_SPLIT_MAP_BAD_NUM_DONORS: gt_error(PARSE_MAP_SPLIT_MAP_BAD_NUM_DONORS,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MISMS_ALREADY_PARSED: gt_error(PARSE_MAP_MISMS_ALREADY_PARSED,file_name,line_num); break;
    case GT_IMP_PE_MISMS_BAD_CHARACTER: gt_error(PARSE_MAP_MISMS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_IMP_PE_MISMS_BAD_MISMS_POS: gt_error(PARSE_MAP_MISMS_BAD_MISMS_POS,file_name,line_num,column_pos); break;
    default:
      gt_error(PARSE_MAP,file_name,line_num);
      break;
  }
}
/* MAP file. Skip record */
GT_INLINE void gt_input_map_parser_next_record(gt_buffered_input_file* const buffered_map_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  if (!gt_buffered_input_file_eob(buffered_map_input)) {
    GT_INPUT_FILE_SKIP_LINE(buffered_map_input);
  }
}
/* Read last record's tag (for block synchronization purposes ) */
GT_INLINE gt_status gt_input_map_parser_get_tag_last_read(gt_buffered_input_file* const buffered_map_input,gt_string* const last_tag) {
  register int64_t position = gt_vector_get_used(buffered_map_input->block_buffer);
  register char* text_line = gt_vector_get_last_elm(buffered_map_input->block_buffer,char);
  if (*text_line!=EOL) return GT_IMP_FAIL;
  // Skip EOL
  while (position>=0 && *text_line==EOL) {--text_line; --position;}
  // Skip Text and get previos EOL and TABS
  register uint64_t last_tab = 0;
  while (position>=0 && *text_line!=EOL) {
    if (*text_line==TAB) last_tab = position;
    --text_line; --position;
  }
  // Check end cases
  if (last_tab==0 || last_tab<=position) return GT_IMP_FAIL; // Sth is wrong
  // Set last record's tag
  if (position>=0) ++text_line;
  gt_string_set_nstring(last_tag,text_line,last_tab-position);
  return GT_IMP_OK;
}
/* */
GT_INLINE gt_status gt_input_map_parser_reload_buffer_matching_tag(
    gt_buffered_input_file* const buffered_map_input,gt_string* const reference_tag,const bool read_paired) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  // Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
  if (buffered_map_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_dump(buffered_map_input->buffered_output_file);
  }
  // Read new input block
  register gt_input_file* const input_file = buffered_map_input->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  gt_input_file_lock(input_file);
  if (input_file->eof) {
    gt_input_file_unlock(input_file);
    return GT_BMI_EOF;
  }
  buffered_map_input->block_id = gt_input_file_next_id(input_file) % UINT32_MAX;
  buffered_map_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clean(buffered_map_input->block_buffer); // Clear dst buffer
  // Read lines
  if (read_paired) gt_input_fastq_tag_chomp_end_info(reference_tag);
  register gt_string* const last_tag = gt_string_new(0);
  register uint64_t lines_read, total_lines_read = 0;
  uint64_t num_blocks = 0, num_tabs = 0;
  do {
    if ((lines_read=gt_input_file_next_record(input_file,
        buffered_map_input->block_buffer,last_tag,&num_blocks,&num_tabs))==0) break;
    if (read_paired) gt_input_fastq_tag_chomp_end_info(last_tag);
    ++total_lines_read;
  } while (!gt_string_equals(reference_tag,last_tag));
  if (read_paired && lines_read>0 && num_blocks%2!=0) { // Check paired read
    gt_input_file_next_line(input_file,buffered_map_input->block_buffer);
    ++total_lines_read;
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_map_input->block_buffer);
  if (total_lines_read > 0 && *gt_vector_get_last_elm(buffered_map_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_map_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=total_lines_read;
  buffered_map_input->lines_in_buffer = total_lines_read;
  gt_input_file_unlock(input_file);
  // Setup the block
  buffered_map_input->cursor = gt_vector_get_mem(buffered_map_input->block_buffer,char);
  // Assign block ID
  if (buffered_map_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_set_block_ids(
        buffered_map_input->buffered_output_file,buffered_map_input->block_id,0);
  }
  // Free
  gt_string_delete(last_tag);
  return GT_IMP_OK;
}
/* MAP file. Synchronized get block wrt to paired map records */
GT_INLINE gt_status gt_input_map_parser_get_block(
    gt_buffered_input_file* const buffered_map_input,const uint64_t num_records) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  register gt_input_file* const input_file = buffered_map_input->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  gt_input_file_lock(input_file);
  if (input_file->eof) {
    gt_input_file_unlock(input_file);
    return GT_BMI_EOF;
  }
  buffered_map_input->block_id = gt_input_file_next_id(input_file) % UINT32_MAX;
  buffered_map_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clean(buffered_map_input->block_buffer); // Clear dst buffer
  // Read lines
  uint64_t lines_read = 0, num_blocks = 0, num_tabs = 0;
  while ( (lines_read<num_records || num_blocks%2!=0) &&
      gt_input_file_next_record(input_file,buffered_map_input->block_buffer,NULL,&num_blocks,&num_tabs) ) ++lines_read;
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_map_input->block_buffer);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffered_map_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_map_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  buffered_map_input->lines_in_buffer = lines_read;
  gt_input_file_unlock(input_file);
  // Setup the block
  buffered_map_input->cursor = gt_vector_get_mem(buffered_map_input->block_buffer,char);
  return buffered_map_input->lines_in_buffer;
}
/* MAP file. Reload internal buffer */
GT_INLINE gt_status gt_input_map_parser_reload_buffer(
    gt_buffered_input_file* const buffered_map_input,const bool synchronized_map,const uint64_t num_lines) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  // Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
  if (buffered_map_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_dump(buffered_map_input->buffered_output_file);
  }
  // Read new input block
  register const uint64_t read_lines = (synchronized_map) ?
      gt_input_map_parser_get_block(buffered_map_input,num_lines):
      gt_buffered_input_file_get_block(buffered_map_input,num_lines);
  if (gt_expect_false(read_lines==0)) return GT_IMP_EOF;
  // Assign block ID
  if (buffered_map_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_set_block_ids(
        buffered_map_input->buffered_output_file,buffered_map_input->block_id,0);
  }
  return GT_IMP_OK;
}

/*
 * MAP format. Basic building block for parsing
 */
GT_INLINE gt_status gt_input_map_parse_tag(char** const text_line,gt_string* const tag) {
  // Delimit the tag
  register char* const tag_begin = *text_line;
  GT_READ_UNTIL(text_line,**text_line==TAB||**text_line==SPACE);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  register uint64_t tag_length = *text_line-tag_begin;
  // Copy string
  gt_string_set_nstring(tag,tag_begin,tag_length);
  // Place cursor at beginning of the next field
  if (gt_expect_false(**text_line==SPACE)) { // Drop everything beyond spaces
    GT_READ_UNTIL(text_line,**text_line==TAB);
    if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
    GT_NEXT_CHAR(text_line);
  } else {
    GT_NEXT_CHAR(text_line);
  }
  return 0;
}
GT_INLINE gt_status gt_input_map_parse_read_block(char** const text_line,gt_string* const read_block) {
  // Read READ_BLOCK
  register char* const read_block_begin = *text_line;
  while (gt_expect_true(**text_line!=TAB && !gt_is_valid_template_separator(**text_line) && !GT_IS_EOL(text_line))) {
    if (gt_expect_false(!gt_is_dna(**text_line))) return GT_IMP_PE_READ_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  // Copy string
  gt_string_set_nstring(read_block,read_block_begin,(*text_line-read_block_begin));
  // Place cursor at beginning of the next field
  register gt_status return_status;
  if (**text_line==TAB) {
    return_status = GT_IMP_PE_EOB;
    GT_NEXT_CHAR(text_line);
  } else if (gt_is_valid_template_separator(**text_line)) {
    return_status = GT_IMP_PE_PENDING_BLOCKS;
    GT_NEXT_CHAR(text_line);
  } else {
    return_status = GT_IMP_PE_READ_BAD_CHARACTER;
  }
  return return_status;
}
GT_INLINE gt_status gt_input_map_parse_qualities_block(char** const text_line,gt_string* const qualities_block) {
  // Read QUAL_BLOCK
  register char* const qualities_block_begin = *text_line;
  while (gt_expect_true(**text_line!=TAB && !gt_is_valid_template_separator(**text_line) && !GT_IS_EOL(text_line))) {
    if (gt_expect_false(!gt_is_valid_quality(**text_line))) return GT_IMP_PE_QUAL_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  // Copy string
  gt_string_set_nstring(qualities_block,qualities_block_begin,(*text_line-qualities_block_begin));
  // Place cursor at beginning of the next field
  register gt_status return_status;
  if (**text_line==TAB) {
    return_status = GT_IMP_PE_EOB;
    GT_NEXT_CHAR(text_line);
  } else if (gt_is_valid_template_separator(**text_line)) {
    return_status = GT_IMP_PE_PENDING_BLOCKS;
    GT_NEXT_CHAR(text_line);
  } else {
    return GT_IMP_PE_QUAL_BAD_SEPARATOR;
  }
  return return_status;
}
GT_INLINE gt_status gt_imp_counters(char** const text_line,gt_vector* const counters,gt_shash* const attributes) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_VECTOR_CHECK(counters);
  GT_HASH_CHECK(attributes);
  // Handle 'not-unique' flag
  if (**text_line==GT_MAP_COUNTS_NOT_UNIQUE) {
    bool not_unique = true;
    gt_attribute_set(attributes,GT_ATTR_NOT_UNIQUE,&not_unique,sizeof(bool));
    GT_NEXT_CHAR(text_line);
    if (gt_expect_false(**text_line!=TAB)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
    return 0;
  }
  // Parse counters
  register bool prev_char_was_sep = false, is_mcs_set = false;
  register uint64_t number;
  while (gt_expect_true(**text_line!=TAB && !GT_IS_EOL(text_line))) {
    if (gt_is_number(**text_line)) {
      GT_PARSE_NUMBER(text_line,number);
      gt_vector_insert(counters,number,uint64_t);
      prev_char_was_sep = false;
    } else if (**text_line==GT_MAP_COUNTS_TIMES) { // 0x10:1:1
      if (prev_char_was_sep) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      if (gt_expect_false(gt_vector_get_used(counters)==0)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      register uint64_t multiplier, i;
      // Parse multiplier
      GT_NEXT_CHAR(text_line);
      if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
      if (!gt_is_number(**text_line)) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,multiplier);
      number = *gt_vector_get_last_elm(counters,uint64_t);
      for (i=1;i<multiplier;++i) {
        gt_vector_insert(counters,number,uint64_t);
      }
    } else if (**text_line==GT_MAP_MCS) {
      if (prev_char_was_sep || is_mcs_set) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      is_mcs_set = true;
      gt_attribute_set(attributes,GT_ATTR_MAX_COMPLETE_STRATA,&(gt_vector_get_used(counters)),sizeof(uint64_t));
      GT_NEXT_CHAR(text_line);
      prev_char_was_sep = true;
    } else if (**text_line==GT_MAP_COUNTS_SEP) {
      if (prev_char_was_sep) return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
      if (*(*text_line+1)==GT_MAP_COUNTS_SEP) break;
      GT_NEXT_CHAR(text_line);
      prev_char_was_sep = true;
    } else {
      return GT_IMP_PE_COUNTERS_BAD_CHARACTER;
    }
  }
  // Set default MCS
  if (!is_mcs_set) {
    gt_attribute_set(attributes,GT_ATTR_MAX_COMPLETE_STRATA,&(gt_vector_get_used(counters)),sizeof(uint64_t));
  }
  // Parse attributes (if any)
  if (**text_line==GT_MAP_COUNTS_SEP && *(*text_line+1)==GT_MAP_COUNTS_SEP) { // 0:0::<Value>::<Value>
    return GT_IMP_PE_NOT_IMPLEMENTED; // Unlikely to be include into the grammar
  }
  return 0;
}
GT_INLINE gt_status gt_imp_parse_strand(char** const text_line,gt_strand* const strand) {
  switch ((**text_line)) {
    case GT_MAP_STRAND_FORWARD_SYMBOL:
    case GT_MAP_STRAND_FORWARD_LETTER:
      *strand=FORWARD;
    break;
    case GT_MAP_STRAND_REVERSE_SYMBOL:
    case GT_MAP_STRAND_REVERSE_LETTER:
      *strand=REVERSE;
    break;
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
      break;
  }
  return 0;
}
// OLD (v0): <+3>20A88C89C99
GT_INLINE gt_status gt_imp_parse_mismatch_string_v0(char** const text_line,gt_map* map) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  gt_map_clear_misms(map);
  // Parse Misms
  register uint64_t last_position = 0, last_cut_point = 0;  
  register const uint64_t global_length = gt_map_get_base_length(map);
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP &&
         !GT_IS_EOL(text_line) && (**text_line)!=GT_MAP_SCORE_GEMv0) {
    gt_misms misms;
    if (gt_is_dna((**text_line))) { // Mismatch
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      GT_NEXT_CHAR(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,misms.position);
      if (gt_expect_false(misms.position<=last_position)) {
        return GT_IMP_PE_MISMS_BAD_MISMS_POS;
      }
      --misms.position; // Zero based position
      last_position = misms.position;
      misms.position -= last_cut_point; // Split-offset correction
      // Add Mismatch
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='<') { // Indel
      register bool is_splice;
      GT_NEXT_CHAR(text_line);
      // Parse operation [+-*]
      switch ((**text_line)) {
        case GT_MAP_INDEL_INSERTION:
          misms.misms_type = INS;
          is_splice = false;
          break;
        case GT_MAP_INDEL_DELETION:
          misms.misms_type = DEL;
          is_splice = false;
          break;
        case GT_MAP_INDEL_SPLICE:
          is_splice = true;
          break;
        default:
          return GT_IMP_PE_MISMS_BAD_CHARACTER;
          break;
      }
      GT_NEXT_CHAR(text_line);
      // Parse size
      register uint64_t size, position;
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,size);
      // Parse Indel end ">"
      if (gt_expect_false((**text_line)!='>')) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_NEXT_CHAR(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,position);
      if (gt_expect_false(position<=last_position)) {
        return GT_IMP_PE_MISMS_BAD_MISMS_POS;
      }
      --position; // Zero based position
      last_position = position;
      // Add Indel
      if (gt_expect_true(!is_splice)) {
        misms.position = position-last_cut_point;
        misms.size = size;
        gt_map_add_misms(map,&misms);
      } else {
        // Close current map block
        gt_map_set_base_length(map,position-last_cut_point);
        last_cut_point = position;
        // Create a new map block
        gt_map* next_map = gt_map_new();
        gt_map_set_next_block(map,next_map,SPLICE);
        gt_map_set_seq_name(next_map,gt_string_get_string(map->seq_name),gt_string_get_length(map->seq_name));
        gt_map_set_position(next_map,gt_map_get_position(map)+gt_map_get_base_length(map)+size);
        gt_map_set_strand(next_map,gt_map_get_strand(map));
        gt_map_set_base_length(next_map,global_length-position);
        // Swap maps & Reset length,position
        map = next_map;
      }
    } else { // ?¿ Parsing error
      return GT_IMP_PE_MISMS_BAD_CHARACTER;
    }
  }
  return 0;
}
// NEW (v1): (5)43T46A9>24*  ||  33C9T30T24>1-(10)
GT_INLINE gt_status gt_imp_parse_mismatch_string_v1(char** const text_line,gt_map* map) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  gt_map_clear_misms(map);
  // Parse Misms
  register uint64_t position=0, length=0;
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP && !GT_IS_EOL(text_line)) {
    gt_misms misms;
    if (gt_is_number((**text_line))) { // Matching
      register uint64_t matching_characters;
      GT_PARSE_NUMBER(text_line,matching_characters);
      position+=matching_characters;
      length+=matching_characters;
    } else if (gt_is_dna((**text_line))) { // Mismatch
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      misms.position = position;
      ++position; ++length;
      GT_NEXT_CHAR(text_line);
      // Add Mismatch
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='(') { // Trim
      misms.misms_type = DEL;
      misms.position = position;
      GT_NEXT_CHAR(text_line);
      // Parse size
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,misms.size);
      position+=misms.size;
      // Parse Trim end ')'
      if (gt_expect_false((**text_line)!=')')) return GT_IMP_PE_MISMS_BAD_CHARACTER;
      GT_NEXT_CHAR(text_line);
      // Add Trim
      if (misms.size>0) gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='>') { // Indel/Skip
      GT_NEXT_CHAR(text_line);
      // Parse size
      register int64_t size;
      GT_PARSE_SIGNED_NUMBER_BLOCK(text_line,size) {
        if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_BAD_CHARACTER;
        GT_PARSE_NUMBER(text_line,size)
      } GT_PARSE_SIGNED_NUMBER_END_BLOCK(size);
      // Parse skip type
      if (size > 0 && ((**text_line)==GT_MAP_SKIP_POSITIVE || (**text_line)==GT_MAP_SKIP_NEGATIVE)) {  // INS/DEL
        misms.position = position;
        misms.size = size;
        if ((**text_line)==GT_MAP_SKIP_POSITIVE) {
          misms.misms_type = INS;
          length+=misms.size;
        } else {
          misms.misms_type = DEL;
          position+=misms.size;
        }
        GT_NEXT_CHAR(text_line);
        // Add Indel/Skip
        gt_map_add_misms(map,&misms);
      } else { // NSKIP/SPLICE
        register gt_junction_t junction;
        switch ((**text_line)) {
          case GT_MAP_SKIP_POSITIVE: junction=POSITIVE_SKIP; break;
          case GT_MAP_SKIP_NEGATIVE: junction=NEGATIVE_SKIP; break;
          case GT_MAP_SKIP_SPLICE: junction=SPLICE; break;
          default: return GT_IMP_PE_MISMS_BAD_CHARACTER; break;
        }
        GT_NEXT_CHAR(text_line);
        // Create a new map block
        gt_map* next_map = gt_map_new();
        gt_map_set_seq_name(next_map,gt_string_get_string(map->seq_name),gt_string_get_length(map->seq_name));
        gt_map_set_position(next_map,gt_map_get_position(map)+length+size);
        gt_map_set_strand(next_map,gt_map_get_strand(map));
        gt_map_set_base_length(next_map,gt_map_get_base_length(map)-position); // FIXED: position->length
        // Close current map block
        gt_map_set_base_length(map,position); // FIXED: position->length
        gt_map_set_next_block(map,next_map,junction);
        // Swap maps & Reset length,position
        map = next_map;
        position=0; length=0;
      }
    } else { // ?¿ Parsing error
     return GT_IMP_PE_MISMS_BAD_CHARACTER;
    }
  }
  return 0;
}
#define GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line) ((**text_line)==GT_MAP_SPLITMAP_NEXT_GEMv0_0 || (**text_line)==GT_MAP_SPLITMAP_NEXT_GEMv0_1)
GT_INLINE gt_status gt_imp_parse_split_map_v0(char** text_line,gt_vector* const split_maps,const uint64_t read_base_length) {
  GT_NULL_CHECK(text_line);
  GT_VECTOR_CHECK(split_maps);
  //
  // Split-map format (v0)
  //   [23]=chr6:R31322884~chr6:R31237276
  //   [70-71]=chr1:F188862944~chr19:F[53208292-53208293]
  //   [23-50]=chr1:F[188862944-188868041]~chr19:F53208292
  // But also we have ....
  //   [31;35]=chr16:R[2503415;2503411]~chr16:R2503271
  //   [30;34]=chr10:F74776624~chr10:F[74790025;74790029]
  // And also ...
  //   [25-26]=chr12:F95317898~chr12:F[95317924-95317925],[25-26]=chr11:R[100658012-100658011]~chr12:F[95317924-95317925]
  // Life is like a box of chocolates ...
  /*
   * Parse split-points
   */
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_OPEN_GEMv0)) return GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER;
  // Create and add the sm
  register gt_map* const donor_map = gt_map_new();
  gt_vector_insert(split_maps,donor_map,gt_map*);
  // Read split-points
  register uint64_t sm_position;
  register bool sm_elm_parsed = false;
  do {
    GT_NEXT_CHAR(text_line);
    if (!gt_is_number((**text_line))) return GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER;
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) {
      gt_map_set_base_length(donor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Read closing split-points and definition symbol
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) return GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_DEF_GEMv0)) return GT_IMP_PE_SPLIT_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  /*
   * Parse donor
   */
  // Read TAG
  register char* const donor_name = *text_line;
  GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  gt_map_set_seq_name(donor_map,donor_name,(*text_line-donor_name));
  GT_NEXT_CHAR(text_line);
  // Read Strand
  register gt_status error_code;
  if ((error_code=gt_imp_parse_strand(text_line,&(donor_map->strand)))) return error_code;
  GT_NEXT_CHAR(text_line);
  // Detect multiple donor position
  register bool multiple_sm;
  if ((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0) {
    GT_NEXT_CHAR(text_line); // [31;35]=chr16:R[2503415;2503411]~chr16:R2503271
    multiple_sm = true;
  } else {
    multiple_sm = false;
  }
  // Read Position
  sm_elm_parsed = false;
  do {
    if (sm_elm_parsed) GT_NEXT_CHAR(text_line);
    if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) { // Set donor
      gt_map_set_position(donor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Close multiple donor syntax
  if (multiple_sm) {
    if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  // Read separator (~)
  if (gt_expect_false(**text_line!=GT_MAP_SPLITMAP_SEP_GEMv0)) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  /*
   * Parse acceptor(s)
   */
  // Read acceptor's TAG
  register gt_map* const acceptor_map = gt_map_new();
  register char* const acceptor_name = *text_line;
  GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  gt_map_set_seq_name(acceptor_map,acceptor_name,(*text_line-acceptor_name));
  GT_NEXT_CHAR(text_line);
  // Read acceptor's Strand
  if ((error_code=gt_imp_parse_strand(text_line,&(acceptor_map->strand)))) return error_code;
  GT_NEXT_CHAR(text_line);
  // Link both donor & acceptor
  gt_map_set_base_length(acceptor_map,read_base_length-gt_map_get_base_length(donor_map));
  gt_map_set_next_block(donor_map,acceptor_map,SPLICE);
  // Detect multiple donor position
  if (gt_expect_false((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0)) { // [30;34]=chr10:F74776624~chr10:F[74790025;74790029]
    GT_NEXT_CHAR(text_line);
    multiple_sm = true;
  } else {
    multiple_sm = false;
  }
  // Parse acceptor position
  sm_elm_parsed = false;
  do {
    // Parse acceptor position
    if (sm_elm_parsed) GT_NEXT_CHAR(text_line);
    if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_PARSE_NUMBER(text_line,sm_position);
    if (!sm_elm_parsed) { // Store split-map's position
      gt_map_set_position(acceptor_map,sm_position);
      sm_elm_parsed=true;
    }
  } while (GT_IMP_PARSE_SPLITMAP_IS_SEP(text_line));
  // Close multiple acceptor syntax
  if (multiple_sm) {
    if (gt_expect_false((**text_line)!=GT_MAP_SPLITMAP_CLOSE_GEMv0)) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  // Return
  switch (**text_line) {
    case GT_MAP_NEXT: // ','
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_PENDING_MAPS;
    case EOL:
    case EOS:
      return GT_IMP_PE_EOB; // '\n'
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
}
GT_INLINE gt_status gt_imp_parse_map_score_v0(char** text_line,gt_map* const map) {
  // Parse score1
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_PARSE_NUMBER(text_line,map->score);
  // Skip score2 (no use)
  if (gt_expect_false(**text_line!=GT_MAP_SCORE_SEP)) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_READ_UNTIL(text_line,!gt_is_number((**text_line)));
  return 0;
}
GT_INLINE gt_status gt_imp_parse_map_score_v1(char** text_line,gt_map* const map) {
  // Parse score
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_PARSE_NUMBER(text_line,map->score);
  return 0;
}
GT_INLINE gt_status gt_imp_parse_map_attr_v0(char** text_line) {
  switch (**text_line) {
    case GT_MAP_NEXT: // ','
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_PENDING_MAPS;
    case GT_MAP_SCORE_GEMv0: // '@10/1'
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_ATTR;
    case EOL:
    case EOS:
      return GT_IMP_PE_EOB; // '\n'
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
}
GT_INLINE gt_status gt_imp_parse_map_attr_v1(char** text_line) {
  switch (**text_line) {
    case GT_MAP_NEXT: // ','
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_PENDING_MAPS;
    case EOL:
    case EOS:
      return GT_IMP_PE_EOB; // '\n'
    case GT_MAP_SEP:
      if ( (*(*text_line+1)) == GT_MAP_SEP) { // '::'
        if ( (*(*text_line+2)) == GT_MAP_SEP) { // ':::' (Attributes of the block group)
          (*text_line)+=3;
          return GT_IMP_PE_MAP_GLOBAL_ATTR;
        } else { // '::?'
          (*text_line)+=2;
          return GT_IMP_PE_PENDING_BLOCKS;
        }
      } else if (gt_is_number( (*(*text_line+1)) )) { // ':100'
        GT_NEXT_CHAR(text_line);
        return GT_IMP_PE_MAP_ATTR;
      } else { // ':?'
        return GT_IMP_PE_MAP_BAD_CHARACTER;
      }
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
}
GT_INLINE gt_status gt_imp_parse_map_global_attr(char** text_line) {
  switch (**text_line) {
    case GT_MAP_NEXT: // ','
      GT_NEXT_CHAR(text_line);
      return GT_IMP_PE_MAP_PENDING_MAPS;
    case EOL:
    case EOS:
      return GT_IMP_PE_EOB; // '\n'
    case GT_MAP_SEP:
      if ( (*(*text_line+1)) == GT_MAP_SEP &&
           (*(*text_line+2)) == GT_MAP_SEP) { // ':::'
          (*text_line)+=3;
          return GT_IMP_PE_MAP_GLOBAL_ATTR;
      } else { // ':?'
        return GT_IMP_PE_MAP_BAD_CHARACTER;
      }
    default:
      return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
}
GT_INLINE gt_status gt_imp_parse_map(char** text_line,gt_map* const map,const gt_lazy_parse_mode parse_mode) {
  /*
   * Parse MAP:
   *   OLD(v0)={chr7:F127708134G27T88}
   *   NEW(v2)={chr11:-:51590050:(5)43TTC5>5-46A9>24*}
   */
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  // Read TAG
  register char* const seq_name_start = *text_line;
  GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  gt_map_set_seq_name(map,seq_name_start,(*text_line-seq_name_start));
  if (gt_expect_false(**text_line!=GT_MAP_SEP)) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_NEXT_CHAR(text_line);
  // Read Strand
  register gt_status error_code;
  gt_strand map_strand;
  if ((error_code=gt_imp_parse_strand(text_line,&map_strand))) return error_code;
  gt_map_set_strand(map,map_strand);
  GT_NEXT_CHAR(text_line);
  // Determine format version
  register gt_misms_string_t misms_format;
  if ((**text_line)==GT_MAP_SEP) { // GEMv1
    misms_format = MISMATCH_STRING_GEMv1;
    GT_NEXT_CHAR(text_line);
  } else if (gt_is_number((**text_line))) { // GEMv0
    misms_format = MISMATCH_STRING_GEMv0;
  } else { // ?¿
    return GT_IMP_PE_MAP_BAD_CHARACTER;
  }
  // Position
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
  GT_PARSE_NUMBER(text_line,map->position);
  // Synch with mismatch string (GEMv1)
  if (misms_format==MISMATCH_STRING_GEMv1) {
    if (gt_expect_false((**text_line)!=GT_MAP_SEP)) return GT_IMP_PE_MAP_BAD_CHARACTER;
    GT_NEXT_CHAR(text_line);
  }
  // Parse Mismatch String
  if (parse_mode==PARSE_ALL) {
    gt_map_clear_misms_string(map);
    map->misms_txt_format = misms_format; // Just in case!
    if (misms_format==MISMATCH_STRING_GEMv1) {
      error_code=gt_imp_parse_mismatch_string_v1(text_line,map);
    } else {
      error_code=gt_imp_parse_mismatch_string_v0(text_line,map);
    }
    if (error_code) return error_code;
  } else {
    gt_map_set_misms_string(map,*text_line,misms_format);
    GT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP || (**text_line)==GT_MAP_NEXT || (**text_line)==GT_MAP_SCORE_GEMv0);
  }
  // Parse map's attributes (if any) and spot next item.
  register bool parsed_score = false;
  while (true) {
    error_code = (misms_format==MISMATCH_STRING_GEMv0) ?
        gt_imp_parse_map_attr_v0(text_line) : gt_imp_parse_map_attr_v1(text_line);
    switch (error_code) {
      case GT_IMP_PE_MAP_ATTR: /* Parse score */
        if (parsed_score) return GT_IMP_PE_MAP_BAD_CHARACTER;
        parsed_score = true;
        if (misms_format==MISMATCH_STRING_GEMv0) {
          if ((error_code=gt_imp_parse_map_score_v0(text_line,map))) return error_code;
        } else {
          if ((error_code=gt_imp_parse_map_score_v1(text_line,map))) return error_code;
        }
        break;
      default:
        return error_code;
        break;
    }
  }
}
#define GT_IMP_PARSE_MAP_ERROR(error_code) \
  (error_code!=GT_IMP_PE_PENDING_BLOCKS && \
   error_code!=GT_IMP_PE_MAP_PENDING_MAPS && \
   error_code!=GT_IMP_PE_EOB && \
   error_code!=GT_IMP_PE_MAP_GLOBAL_ATTR)
#define GT_IMP_PARSE_MAP_ATTR_ERROR(error_code) \
  (error_code!=GT_IMP_PE_MAP_PENDING_MAPS && \
   error_code!=GT_IMP_PE_EOB && \
   error_code!=GT_IMP_PE_MAP_GLOBAL_ATTR)

GT_INLINE gt_status gt_imp_parse_template_maps(
    char** text_line,gt_template* const template,
    const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  // Set as parsed (whatever the result is)
  template->maps_txt = NULL;
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_SKIP_LINE(text_line);
    return 0;
  }
  // Check max_num_maps to parse (stratum-wise)
  uint64_t max_num_maps = num_maps;
  if (num_maps<GT_ALL) {
    uint64_t strata;
    gt_counters_calculate_num_maps(gt_template_get_counters_vector(template),0,num_maps,&strata,&max_num_maps);
  }
  // Parse MAPS. Formats allowed:
  //   OLD (v0): chr7:F127708134G27T88::chr7:R127708509<+3>20A88C89C99
  //   NEW (v1): chr11:-:51590050:(5)43T46A9>24*::chr11:-:51579823:33C9T30T24>1-(10)
  register gt_status error_code = GT_IMP_PE_MAP_PENDING_MAPS;
  register const uint64_t num_blocks_template = gt_vector_get_used(template->blocks);
  register uint64_t num_maps_parsed = 0, num_blocks_parsed;
  register gt_vector* vector_maps = gt_vector_new(num_blocks_template,sizeof(gt_map*));
  gt_mmap_attributes mmap_attr;
  while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && num_maps_parsed<max_num_maps) {
    /*
     * Parse MAP
     */
    error_code = GT_IMP_PE_PENDING_BLOCKS;
    gt_vector_clear(vector_maps);
    num_blocks_parsed = 0;
    while (error_code==GT_IMP_PE_PENDING_BLOCKS) {
      register gt_alignment* const alignment = gt_template_get_block(template,num_blocks_parsed);
      // Allocate new map
      gt_map* map = gt_map_new();
      gt_map_set_base_length(map,gt_alignment_get_read_length(alignment)); // Set base length
      // Parse current MAP
      error_code = gt_imp_parse_map(text_line,map,parse_mode);
      if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
      // Add the map
      gt_alignment_add_map(alignment,map);
      gt_vector_insert(vector_maps,map,gt_map*);
      ++num_blocks_parsed;
    }
    /*
     * Parse attributes (if any) and calculate template attributes
     */
    gt_template_mmap_attr_new(&mmap_attr);
    mmap_attr.distance = gt_map_vector_get_distance(vector_maps);
    if (error_code==GT_IMP_PE_MAP_GLOBAL_ATTR) {
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_IMP_PE_MAP_BAD_CHARACTER;
      GT_PARSE_NUMBER(text_line,mmap_attr.score);
      error_code=gt_imp_parse_map_global_attr(text_line);
      if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
    } else {
      mmap_attr.score = GT_MAP_NO_SCORE;
    }
    /*
     * Check number of blocks parsed & Store MAP blocks parsed (only PE or more)
     */
    if (gt_expect_false(num_blocks_parsed>num_blocks_template)) {
      // Weird case of more blocks than blocks in the template (reorganize blocks)
      return GT_IMP_PE_NOT_IMPLEMENTED;
    }
    if (gt_expect_false(num_blocks_parsed<num_blocks_template || num_blocks_parsed==1)) {
      return GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS;
    }
    gt_template_add_mmap_gtvector(template,vector_maps,&mmap_attr);
    ++num_maps_parsed;
  }
  gt_vector_delete(vector_maps);
  return 0;
}
GT_INLINE gt_status gt_imp_parse_alignment_maps(
    char** text_line,gt_alignment* alignment,
    const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_ALIGNMENT_CHECK(alignment);
  // Set as parsed (whatever the result is)
  alignment->maps_txt = NULL;
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_SKIP_LINE(text_line);
    return 0;
  }
  // Check max_num_maps to parse (stratum-wise)
  uint64_t max_num_maps = num_maps;
  if (num_maps<GT_ALL) {
    uint64_t strata;
    gt_counters_calculate_num_maps(gt_alignment_get_counters_vector(alignment),0,num_maps,&strata,&max_num_maps);
  }
  // Parse MAPS. Formats allowed:
  //   OLD (v0): chr7:F127708134G27T88
  //   NEW (v1): chr11:-:51590050:(5)43T46A9>24*
  register const uint64_t alignment_base_length = gt_alignment_get_read_length(alignment);
  register uint64_t num_maps_parsed = 0;
  register gt_status error_code = GT_IMP_PE_MAP_PENDING_MAPS;
  register gt_vector* split_maps = NULL;
  while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && num_maps_parsed<max_num_maps) {
    // Parse current MAP
    if (gt_expect_false((**text_line)==GT_MAP_SPLITMAP_OPEN_GEMv0)) { // Old split-maps switch
      if (gt_expect_true(split_maps==NULL)) split_maps = gt_vector_new(2,sizeof(gt_map*));
      error_code = gt_imp_parse_split_map_v0(text_line,split_maps,alignment_base_length);
      if (error_code!=GT_IMP_PE_MAP_PENDING_MAPS && error_code!=GT_IMP_PE_EOB) {
        gt_vector_delete(split_maps);
        return error_code;
      }
      gt_alignment_add_map_gt_vector(alignment,split_maps);
      num_maps_parsed+=gt_vector_get_used(split_maps);
      gt_vector_clear(split_maps);
    } else {
      register gt_map* map = gt_map_new();
      // Set base length (needed to calculate the alignment's length in GEMv0)
      gt_map_set_base_length(map,alignment_base_length);
      error_code = gt_imp_parse_map(text_line,map,parse_mode);
      switch (error_code) {
        case GT_IMP_PE_PENDING_BLOCKS:
          if (gt_expect_false(split_maps!=NULL)) gt_vector_delete(split_maps);
          return GT_IMP_PE_MAP_BAD_NUMBER_OF_BLOCKS; /* Weird case of split blocks */
        case GT_IMP_PE_MAP_GLOBAL_ATTR:
          if (gt_expect_false(!gt_is_number((**text_line)))) {
            if (gt_expect_false(split_maps!=NULL)) gt_vector_delete(split_maps);
            return GT_IMP_PE_MAP_BAD_CHARACTER;
          }
          GT_PARSE_NUMBER(text_line,map->score);
          error_code=gt_imp_parse_map_attr_v1(text_line);
          break;
        case GT_IMP_PE_MAP_PENDING_MAPS:
        case GT_IMP_PE_EOB: break;
        default:
          if (gt_expect_false(split_maps!=NULL)) gt_vector_delete(split_maps);
          return error_code;
      }
      // Store map into alignment
      gt_alignment_add_map(alignment,map);
      ++num_maps_parsed;
    }
  }
  if (gt_expect_false(split_maps!=NULL)) gt_vector_delete(split_maps);
  return 0;
}
GT_INLINE gt_status gt_imp_map_blocks(char** const text_line,gt_map* const map) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_MAP_CHECK(map);
  // Clear map
  gt_map_clear(map);
  // Check null map
  if ((**text_line)==GT_MAP_NONE) return 0;
  // Read all map blocks
  register gt_status error_code = gt_imp_parse_map(text_line,map,PARSE_ALL);
  if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
  while (error_code==GT_IMP_PE_PENDING_BLOCKS) {
    register gt_map* next_map = gt_map_new();
    error_code = gt_imp_parse_map(text_line,next_map,PARSE_ALL);
    if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
    gt_map_append_block(map,next_map,JUNCTION_UNKNOWN);
  }
  // Skip attributes
  if (error_code==GT_IMP_PE_MAP_GLOBAL_ATTR) {
    while (**text_line!=GT_MAP_NEXT && !GT_IS_EOL(text_line)) {
      GT_NEXT_CHAR(text_line);
    }
    switch (**text_line) {
      case GT_MAP_NEXT: // ','
        GT_NEXT_CHAR(text_line);
        return GT_IMP_PE_MAP_PENDING_MAPS;
      case EOL:
      case EOS:
        return GT_IMP_PE_EOB; // '\n'
      default:
        return GT_IMP_PE_MAP_BAD_CHARACTER;
    }
  }
  return error_code;
}
/*
 * MAP string parsers
 */
GT_INLINE gt_status gt_imp_parse_template(
    char** const text_line,gt_template* const template,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line);
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  // TAG
  if ((error_code=gt_input_map_parse_tag(text_line,template->tag))) {
    return error_code;
  }
>>>>>>>>>>>>>>>>>>>> File 1
>>>>>>>>>>>>>>>>>>>> File 2

>>>>>>>>>>>>>>>>>>>> File 3

<<<<<<<<<<<<<<<<<<<<
  // READ
  register uint64_t num_blocks = 0;
  error_code=GT_IMP_PE_PENDING_BLOCKS;
  while (error_code==GT_IMP_PE_PENDING_BLOCKS) {
    gt_alignment* const alignment = gt_template_get_block_dyn(template,num_blocks);
    error_code=gt_input_map_parse_read_block(text_line,alignment->read);
    if (error_code!=GT_IMP_PE_PENDING_BLOCKS && error_code!=GT_IMP_PE_EOB) return error_code;
    ++num_blocks;
  }
>>>>>>>>>>>>>>>>>>>> File 1
>>>>>>>>>>>>>>>>>>>> File 2

>>>>>>>>>>>>>>>>>>>> File 3

<<<<<<<<<<<<<<<<<<<<
  // TAG Setup {Tag Splitting/Swapping}
  if (gt_expect_true(num_blocks>1)) {
    gt_template_deduce_alignments_tags(template);
  } else {
    gt_alignment* const alignment = gt_template_get_block(template,0);
    GT_SWAP(template->tag,alignment->tag);
    gt_template_deduce_template_tag(template,alignment);
  }
>>>>>>>>>>>>>>>>>>>> File 1
>>>>>>>>>>>>>>>>>>>> File 2

>>>>>>>>>>>>>>>>>>>> File 3

<<<<<<<<<<<<<<<<<<<<
  // QUALITIES
  if (has_map_quality) {
    register uint64_t i;
    error_code=GT_IMP_PE_PENDING_BLOCKS;
    for (i=0;i<num_blocks;++i) {
      if (error_code!=GT_IMP_PE_PENDING_BLOCKS) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
      gt_alignment* alignment = gt_template_get_block(template,i);
      error_code=gt_input_map_parse_qualities_block(text_line,alignment->qualities);
      if (gt_expect_false(gt_string_get_length(alignment->qualities)!=
                          gt_string_get_length(alignment->read))) return GT_IMP_PE_QUAL_BAD_LENGTH;
      if (error_code!=GT_IMP_PE_PENDING_BLOCKS && error_code!=GT_IMP_PE_EOB) return error_code;
    }
    if (error_code!=GT_IMP_PE_EOB) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
  }
>>>>>>>>>>>>>>>>>>>> File 1
>>>>>>>>>>>>>>>>>>>> File 2

>>>>>>>>>>>>>>>>>>>> File 3

<<<<<<<<<<<<<<<<<<<<
  // COUNTERS
  gt_vector_clear(gt_template_get_counters_vector(template));
  if (gt_expect_true(num_blocks>1)) {
    error_code=gt_imp_counters(text_line,gt_template_get_counters_vector(template),template->attributes);
  } else {
    register gt_alignment* const alignment = gt_template_get_block(template,0);
    error_code=gt_imp_counters(text_line,gt_alignment_get_counters_vector(alignment),alignment->attributes);
  }
  if (error_code) return error_code;
  if (gt_expect_false((**text_line)!=TAB)) return GT_IMP_PE_PREMATURE_EOL;
  GT_NEXT_CHAR(text_line);

  // MAPS
  if (parse_mode!=PARSE_READ) {
    template->maps_txt = NULL;
    if (gt_expect_true(num_blocks>1)) {
      error_code = gt_imp_parse_template_maps(text_line,template,parse_mode,num_maps);
    } else {
      error_code = gt_imp_parse_alignment_maps(text_line,gt_template_get_block(template,0),parse_mode,num_maps);
    }
  } else { // (lazy parsing)
    template->maps_txt = *text_line;
    error_code = 0;
  }
  return error_code;
}
GT_INLINE gt_status gt_imp_parse_alignment(
    char** const text_line,gt_alignment* alignment,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line);
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  // TAG
  if ((error_code=gt_input_map_parse_tag(text_line,alignment->tag))) return error_code;
  // READ
  error_code=gt_input_map_parse_read_block(text_line,alignment->read);
  if (gt_expect_false(error_code==GT_IMP_PE_PENDING_BLOCKS)) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
  if (gt_expect_false(error_code!=GT_IMP_PE_EOB)) return error_code;
  // QUALITIES
  if (has_map_quality) {
    error_code=gt_input_map_parse_qualities_block(text_line,alignment->qualities);
    if (gt_expect_false(gt_string_get_length(alignment->qualities)!=
                        gt_string_get_length(alignment->read))) return GT_IMP_PE_QUAL_BAD_LENGTH;
    if (gt_expect_false(error_code==GT_IMP_PE_PENDING_BLOCKS)) return GT_IMP_PE_BAD_NUMBER_OF_BLOCKS;
    if (gt_expect_false(error_code!=GT_IMP_PE_EOB)) return error_code;
  }
  // COUNTERS
  gt_vector_clear(gt_alignment_get_counters_vector(alignment));
  if ((error_code=gt_imp_counters(
      text_line,alignment->counters,alignment->attributes))) return error_code;
  if (GT_IS_EOL(text_line)) return GT_IMP_PE_PREMATURE_EOL;
  if (**text_line!=TAB) return GT_IMP_PE_BAD_SEPARATOR;
  GT_NEXT_CHAR(text_line);
  // MAPS
  if (parse_mode!=PARSE_READ) {
    alignment->maps_txt = NULL;
    error_code=gt_imp_parse_alignment_maps(text_line,alignment,parse_mode,num_maps);
  } else { // (lazy parsing)
    alignment->maps_txt = *text_line;
    error_code=0;
  }
  return error_code;
}
GT_INLINE gt_status gt_input_map_parse_template(char* const string,gt_template* const template) {
  GT_NULL_CHECK(string);
  GT_TEMPLATE_CHECK(template);
  char* _string = string; // Placeholder
  gt_template_clear(template,true); // Clear template
  // Count fields
  register uint64_t num_fields=0, i=0;
  while (gt_expect_true(_string[i]!=EOS)) {
    if (gt_expect_false(_string[i]==TAB)) ++num_fields;
    i++;
  }
  if (gt_expect_true(num_fields==4 || num_fields==3)) {
    return gt_imp_parse_template(&_string,template,num_fields==4,PARSE_ALL,UINT64_MAX);
  }
  return GT_IMP_PE_WRONG_FILE_FORMAT;
}
GT_INLINE gt_status gt_input_map_parse_alignment(char* const string,gt_alignment* const alignment) {
  GT_NULL_CHECK(string);
  GT_ALIGNMENT_CHECK(alignment);
  char* _string = string; // Placeholder
  gt_alignment_clear(alignment); // Clear alignment
  // Count fields
  register uint64_t num_fields=0, i=0;
  while (gt_expect_true(_string[i]!=EOS)) {
    if (gt_expect_false(_string[i]==TAB)) ++num_fields;
    i++;
  }
  if (gt_expect_true(num_fields==4 || num_fields==3)) {
    return gt_imp_parse_alignment(&_string,alignment,num_fields==4,PARSE_ALL,UINT64_MAX);
  }
  return GT_IMP_PE_WRONG_FILE_FORMAT;
}
GT_INLINE gt_status gt_input_map_parse_counters(char* const string,gt_vector* const counters,gt_shash* const attributes) {
  GT_NULL_CHECK(string);
  GT_VECTOR_CHECK(counters);
  GT_HASH_CHECK(attributes);
  char* _string = string; // Placeholder
  return gt_imp_counters(&_string,counters,attributes);
}
GT_INLINE gt_status gt_input_map_parse_map(char* const string,gt_map* const map) {
  GT_NULL_CHECK(string);
  GT_MAP_CHECK(map);
  char* _string = string; // Placeholder
  register gt_status error_code = gt_imp_map_blocks(&_string,map);
  if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
  return 0;
}
GT_INLINE gt_status gt_input_map_parse_map_list(char* const string,gt_vector* const maps,const uint64_t num_maps) {
  GT_NULL_CHECK(string);
  GT_VECTOR_CHECK(maps);
  char* _string = string; // Placeholder // GT_IMP_PE_MAP_PENDING_MAPS
  // Read maps
  register gt_status error_code;
  do {
    // Parse map list
    if (gt_expect_false((*_string)==GT_MAP_SPLITMAP_OPEN_GEMv0)) { // Old split-maps switch
      error_code = gt_imp_parse_split_map_v0(&_string,maps,UINT64_MAX);
      if (error_code!=GT_IMP_PE_MAP_PENDING_MAPS && error_code!=GT_IMP_PE_EOB) return error_code;
    } else {
      register gt_map* map = gt_map_new();
      error_code = gt_imp_map_blocks(&_string,map);
      if (GT_IMP_PARSE_MAP_ERROR(error_code)) return error_code;
      // Add map (NO check duplicates)
      gt_vector_insert(maps,map,gt_map*);
    }
  } while (error_code==GT_IMP_PE_MAP_PENDING_MAPS && gt_vector_get_used(maps)<num_maps);
  return 0;
}

/*
 * MAP Lazy Parsers
 */
GT_INLINE gt_status gt_input_map_parser_parse_template_maps(gt_template* template,uint64_t num_maps) {
  GT_TEMPLATE_CHECK(template);
  if (gt_expect_false(template->maps_txt==NULL)) return GT_IMP_PE_MAP_ALREADY_PARSED;
  register gt_status error_code = gt_imp_parse_template_maps(
      &template->maps_txt,template,PARSE_READ__MAPS,num_maps);
  template->maps_txt = NULL;
  return error_code;
}
GT_INLINE gt_status gt_input_map_parser_parse_alignment_maps(gt_alignment* alignment,uint64_t num_maps) {
  GT_ALIGNMENT_CHECK(alignment);
  if (gt_expect_false(alignment->maps_txt==NULL)) return GT_IMP_PE_MAP_ALREADY_PARSED;
  register gt_status error_code = gt_imp_parse_alignment_maps(
      &alignment->maps_txt,alignment,PARSE_READ__MAPS,num_maps);
  alignment->maps_txt = NULL;
  return error_code;
}
GT_INLINE gt_status gt_input_map_parse_template_mismatch_string(gt_template* template) {
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  gt_template_maps_iterator template_maps_iterator;
  gt_template_new_mmap_iterator(template,&template_maps_iterator);
  register const uint64_t num_blocks_template = gt_vector_get_used(template->blocks);
  gt_map** map_array;
  while (gt_template_next_mmap(&template_maps_iterator,&map_array,NULL)) {
    register uint64_t i;
    for (i=0;i<num_blocks_template;++i) {
      if (gt_expect_false(map_array[i]->misms_txt==NULL)) return GT_IMP_PE_MISMS_ALREADY_PARSED;
      if ((error_code = (gt_map_get_misms_string_format(map_array[i])==MISMATCH_STRING_GEMv1) ?
          gt_imp_parse_mismatch_string_v1(&(map_array[i]->misms_txt),map_array[i]):
          gt_imp_parse_mismatch_string_v0(&(map_array[i]->misms_txt),map_array[i]))) {
        return error_code;
      }
      gt_map_clear_misms_string(map_array[i]);
    }
  }
  free(map_array);
  return 0;
}
GT_INLINE gt_status gt_input_map_parse_alignment_mismatch_string(gt_alignment* alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  gt_map* map;
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    if (gt_expect_false(map->misms_txt==NULL)) return GT_IMP_PE_MISMS_ALREADY_PARSED;
    if ((error_code = (gt_map_get_misms_string_format(map)==MISMATCH_STRING_GEMv1) ?
        gt_imp_parse_mismatch_string_v1(&map->misms_txt,map):
        gt_imp_parse_mismatch_string_v0(&map->misms_txt,map))) {
      return error_code;
    }
    gt_map_clear_misms_string(map);
  }
  return 0;
}

/*
 * MAP High-level Parsers
 */
GT_INLINE gt_status gt_imp_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_string* const src_text,
    gt_template* const template,const gt_lazy_parse_mode parse_mode,const uint64_t num_mmaps) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  register gt_input_file* const input_file = buffered_map_input->input_file;
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_map_input)) {
    if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input,true,GT_IMP_NUM_LINES))!=GT_IMP_OK) return error_code;
  }
  // Check file format
  if (gt_input_map_parser_check_map_file_format(buffered_map_input)) {
    gt_error(PARSE_MAP_BAD_FILE_FORMAT,input_file->file_name,buffered_map_input->current_line_num);
    return GT_IMP_FAIL;
  }
  // Prepare the template
  register char* const line_start = buffered_map_input->cursor;
  register const uint64_t line_num = buffered_map_input->current_line_num;
  gt_template_clear(template,true);
  template->template_id = line_num;
  // Parse template
  if ((error_code=gt_imp_parse_template(&(buffered_map_input->cursor),
      template,input_file->map_type.contains_qualities,parse_mode,num_mmaps))) {
    gt_input_map_parser_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_input_map_parser_next_record(buffered_map_input);
    if (src_text) gt_string_set_nstring(src_text,line_start,buffered_map_input->cursor-line_start);
    return GT_IMP_FAIL;
  }
  // Store source record
  if (src_text) gt_string_set_nstring(src_text,line_start,buffered_map_input->cursor-line_start);
  // Next record
  gt_input_map_parser_next_record(buffered_map_input);
  return GT_IMP_OK;
}
GT_INLINE gt_status gt_imp_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_string* const src_text,
    gt_alignment* const alignment,const gt_lazy_parse_mode parse_mode,const uint64_t num_maps) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  register gt_input_file* const input_file = buffered_map_input->input_file;
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_map_input)) {
    if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input,false,GT_IMP_NUM_LINES))!=GT_IMP_OK) return error_code;
  }
  // Check file format
  if (gt_input_map_parser_check_map_file_format(buffered_map_input)) {
    gt_error(PARSE_MAP_BAD_FILE_FORMAT,input_file->file_name,buffered_map_input->current_line_num);
    return GT_IMP_FAIL;
  }
  // Allocate memory for the alignment
  register char* const line_start = buffered_map_input->cursor;
  register const uint64_t line_num = buffered_map_input->current_line_num;
  gt_alignment_clear(alignment);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_imp_parse_alignment(&(buffered_map_input->cursor),
      alignment,input_file->map_type.contains_qualities,parse_mode,num_maps))) {
    gt_input_map_parser_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_input_map_parser_next_record(buffered_map_input);
    if (src_text) gt_string_set_nstring(src_text,line_start,buffered_map_input->cursor-line_start);
    return GT_IMP_FAIL;
  }
  // Store source record
  if (src_text) gt_string_set_nstring(src_text,line_start,buffered_map_input->cursor-line_start);
  // Next record
  gt_input_map_parser_next_record(buffered_map_input);
  return GT_IMP_OK;
}
GT_INLINE gt_status gt_input_map_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  return gt_imp_get_template(buffered_map_input,NULL,template,PARSE_ALL,GT_ALL);
}
GT_INLINE gt_status gt_input_map_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  return gt_imp_get_alignment(buffered_map_input,NULL,alignment,PARSE_ALL,GT_ALL);
}

GT_INLINE gt_status gt_input_map_parser_get_template__src_text(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_string* const src_text) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  GT_STRING_CHECK(src_text);
  gt_string_cast_static(src_text);
  return gt_imp_get_template(buffered_map_input,src_text,template,PARSE_ALL,GT_ALL);
}
GT_INLINE gt_status gt_input_map_parser_get_alignment__src_text(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_string* const src_text) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  GT_STRING_CHECK(src_text);
  gt_string_cast_static(src_text);
  return gt_imp_get_alignment(buffered_map_input,src_text,alignment,PARSE_ALL,GT_ALL);
}

GT_INLINE gt_status gt_input_map_parser_get_template_limited(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,const uint64_t num_mmaps) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  return gt_imp_get_template(buffered_map_input,NULL,template,PARSE_ALL,num_mmaps);
}
GT_INLINE gt_status gt_input_map_parser_get_alignment_limited(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,const uint64_t num_maps) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  return gt_imp_get_alignment(buffered_map_input,NULL,alignment,PARSE_ALL,num_maps);
}

GT_INLINE gt_status gt_input_map_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,uint64_t num_map_inputs,
    gt_buffered_input_file* const buffered_map_input,va_list v_args) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed (synch)
  if (gt_buffered_input_file_eob(buffered_map_input)) {
    GT_BEGIN_MUTEX_SECTION(*input_mutex) {
      // Reload main 'buffered_map_input' file
      if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input,true))!=GT_IMP_OK) {
        GT_END_MUTEX_SECTION(*input_mutex);
        return error_code;
      }
      --num_map_inputs;
      // Reload the rest of the 'buffered_map_input' files
      while (num_map_inputs>0) {
        register gt_buffered_input_file* buffered_map_input_file = va_arg(v_args,gt_buffered_input_file*);
        if ((error_code=gt_input_map_parser_reload_buffer(buffered_map_input_file,true))!=GT_IMP_OK) {
          GT_END_MUTEX_SECTION(*input_mutex);
          return error_code;
        }
        --num_map_inputs;
      }
    } GT_END_MUTEX_SECTION(*input_mutex);
  }
  return GT_IMP_OK;
}

GT_INLINE gt_status gt_input_map_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,uint64_t num_map_inputs,
    gt_buffered_input_file* const buffered_map_input,...) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input);
  GT_ZERO_CHECK(num_map_inputs);
  va_list v_args;
  va_start(v_args,buffered_map_input);
  register gt_status error_code =
      gt_input_map_parser_synch_blocks_v(input_mutex,num_map_inputs,buffered_map_input,v_args);
  va_end(v_args);
  return error_code;

GT_INLINE gt_status gt_input_map_parser_synch_blocks_by_subset(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file* const buffered_map_input_master,
    gt_buffered_input_file* const buffered_map_input_slave,const bool read_paired) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input_master);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_map_input_slave);
  register gt_status error_code_master, error_code_slave;
  // Check the end_of_block. Reload buffer if needed (synch)
  register const bool eob_master = gt_buffered_input_file_eob(buffered_map_input_master);
  register const bool eob_slave = gt_buffered_input_file_eob(buffered_map_input_slave);
  if (!eob_master) return GT_IMP_OK;
  if (!eob_slave) return GT_IMP_FAIL;
  /*
   * Dump buffer if BOF it attached to Map-input, and get new out block (always FIRST)
   */
  if (buffered_map_input_master->buffered_output_file!=NULL) {
    gt_buffered_output_file_dump(buffered_map_input_master->buffered_output_file);
  }
  /*
   * Read synch blocks
   */
  GT_BEGIN_MUTEX_SECTION(*input_mutex) {
    /*
     * Read new input block for the slave
     */
    if ((error_code_slave=gt_input_map_parser_reload_buffer(buffered_map_input_slave,read_paired,GT_IMP_SUBSET_NUM_LINES+1))!=GT_IMP_OK) {
      // Read for the master and finish
      error_code_master=gt_input_map_parser_reload_buffer(buffered_map_input_master,read_paired,GT_IMP_SUBSET_NUM_LINES+1);
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_master;
    }
    /*
     * Read new input block for the master
     */
    // Get the last tag of the slave block
    register gt_string* const last_tag = gt_string_new(0);
    if ((error_code_slave=gt_input_map_parser_get_tag_last_read(buffered_map_input_slave,last_tag)!=GT_IMP_OK)) {
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_slave;
    }
    // Read new input block for master matching the last tag of the slave block
    if ((error_code_master=gt_input_map_parser_reload_buffer_matching_tag(
        buffered_map_input_master,last_tag,read_paired))!=GT_IMP_OK) {
      gt_string_delete(last_tag); // Free
      GT_END_MUTEX_SECTION(*input_mutex);
      return error_code_master;
    }
    gt_string_delete(last_tag); // Free
  } GT_END_MUTEX_SECTION(*input_mutex);
  return GT_IMP_OK;
}
