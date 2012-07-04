/*
 * PROJECT: GEM-Tools library
 * FILE: gt_buffered_map_input.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_buffered_map_input.h"

#define GT_BMI_BUFFER_SIZE GT_BUFFER_SIZE_4M
#define GT_BMI_NUM_LINES GT_NUM_LINES_5K

// Parsing error/state codes
#define GT_BMI_PE_WRONG_FILE_FORMAT 20
#define GT_BMI_PE_NOT_IMPLEMENTED 27

#define GT_BMI_PE_PREMATURE_EOL 21
#define GT_BMI_PE_PENDING_BLOCKS 22
#define GT_BMI_PE_EOB 23
#define GT_BMI_PE_BAD_SEPARATOR 24
#define GT_BMI_PE_BAD_NUMBER_OF_BLOCKS 25
#define GT_BMI_PE_BAD_CHARACTER 26

#define GT_BMI_PE_READ_BAD_CHARACTER 30

#define GT_BMI_PE_QUAL_BAD_PREMATURE_EOB 40
#define GT_BMI_PE_QUAL_BAD_CHARACTER 41
#define GT_BMI_PE_QUAL_BAD_SEPARATOR 42

#define GT_BMI_PE_COUNTERS_BAD_CHARACTER 50

#define GT_BMI_PE_MAP_PENDING_MAPS 61
#define GT_BMI_PE_MAP_ALREADY_PARSED 71
#define GT_BMI_PE_MAP_BAD_NUMBER_OF_BLOCKS 73
#define GT_BMI_PE_MAP_BAD_CHARACTER 74

#define GT_BMI_PE_MISMS_ALREADY_PARSED 72
#define GT_BMI_PE_MISMS_BAD_CHARACTER 75
#define GT_BMI_PE_MISMS_BAD_MISMS_POS 77

// Useful macros
#define GT_BMI_PARAMS_FILE__LINE(buffered_map_input) \
  buffered_map_input->input_file->file_name,buffered_map_input->current_line_num

/*
 * Buffered map file handlers
 */
gt_buffered_map_input* gt_buffered_map_input_new(gt_input_file* const input_file) {
  GT_NULL_CHECK(input_file);
  gt_buffered_map_input* buffered_map_file = malloc(sizeof(gt_buffered_map_input));
  gt_cond_fatal_error(!buffered_map_file,MEM_HANDLER);
  /* Input file */
  buffered_map_file->input_file = input_file;
  /* Block buffer and cursors */
  buffered_map_file->block_id = UINT64_MAX;
  buffered_map_file->block_buffer = gt_vector_new(GT_BMI_BUFFER_SIZE,sizeof(uint8_t));
  buffered_map_file->cursor = (char*) gt_vector_get_mem(buffered_map_file->block_buffer,uint8_t);
  buffered_map_file->current_line_num = UINT64_MAX;
  return buffered_map_file;
}
gt_status gt_buffered_map_input_close(gt_buffered_map_input* const buffered_map_input) {
  GT_BMI_CHECK(buffered_map_input);
  gt_vector_delete(buffered_map_input->block_buffer);
  return GT_BMI_OK;
}
GT_INLINE uint64_t gt_buffered_map_input_get_cursor_pos(gt_buffered_map_input* const buffered_map_input) {
  GT_BMI_CHECK(buffered_map_input);
  GT_NULL_CHECK(buffered_map_input->cursor);
  return buffered_map_input->cursor-gt_vector_get_mem(buffered_map_input->block_buffer,char);
}
GT_INLINE bool gt_buffered_map_input_eob(gt_buffered_map_input* const buffered_map_input) {
  GT_BMI_CHECK(buffered_map_input);
  return gt_buffered_map_input_get_cursor_pos(buffered_map_input) >= gt_vector_get_used(buffered_map_input->block_buffer);
}
GT_INLINE gt_status gt_buffered_map_input_get_block(gt_buffered_map_input* const buffered_map_input) {
  GT_BMI_CHECK(buffered_map_input);
  register gt_input_file* const input_file = buffered_map_input->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  gt_input_file_lock(input_file);
  if (input_file->eof) {
    gt_input_file_unlock(input_file);
    return GT_BMI_EOF;
  }
  buffered_map_input->block_id = gt_input_file_next_id(input_file);
  buffered_map_input->current_line_num = input_file->processed_lines+1;
  buffered_map_input->lines_in_buffer =
      gt_input_file_get_lines(input_file,buffered_map_input->block_buffer,GT_BMI_NUM_LINES);
  gt_input_file_unlock(input_file);
  // Setup the block
  buffered_map_input->cursor = gt_vector_get_mem(buffered_map_input->block_buffer,char);
  return buffered_map_input->lines_in_buffer;
}

/*
 * MAP File Format test
 */
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors) {
  return gt_buffered_map_input_test_map(
      input_file->file_name,input_file->processed_lines+1,(char*)input_file->file_buffer,
      input_file->buffer_size,map_file_format,show_errors);
}
GT_INLINE bool gt_buffered_map_input_test_map(
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
GT_INLINE gt_status gt_buffered_map_input_check_map_file_format(gt_buffered_map_input* const buffered_map_input) {
  register gt_input_file* const input_file = buffered_map_input->input_file;
  if (gt_expect_false(input_file->file_format==UNKNOWN)) { // Unknown
    gt_map_file_format map_type;
    if (!gt_buffered_map_input_test_map(
        input_file->file_name,buffered_map_input->current_line_num,
        gt_vector_get_mem(buffered_map_input->block_buffer,char),
        gt_vector_get_used(buffered_map_input->block_buffer),&map_type,true)) {
      return GT_BMI_PE_WRONG_FILE_FORMAT;
    }
    input_file->file_format = MAP;
  } else if (gt_expect_false(input_file->file_format!=MAP)) {
    return GT_BMI_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * Internal Building Blocks for parsing
 */
#define GT_BMI_TEXT_NEXT(text_line) ++(*text_line)
#define GT_BMI_TEXT_READ_UNTIL(text_line,test) \
  while (gt_expect_true(!(test) && (**text_line)!=EOL)) { \
    GT_BMI_TEXT_NEXT(text_line); \
  }
#define GT_BMI_TEXT_IS_EOL(text_line) gt_expect_false((**text_line)==EOL)
#define GT_BMI_TEXT_PARSE_NUMBER(text_line,number) \
  number = 0; \
  while (gt_expect_true(gt_is_number(**text_line))) { \
    number = (number*10) + gt_get_cipher(**text_line); \
    GT_BMI_TEXT_NEXT(text_line); \
  }
#define GT_BMI_TEXT_PARSE_SIGNED_NUMBER(text_line,number) { \
  register bool is_negative; \
  switch ((**text_line)) { \
    case PLUS: is_negative = false; break; \
    case MINUS: is_negative = true; break; \
    default: \
      is_negative = false; \
      if (gt_expect_false(!gt_is_number((**text_line)))) { \
        return GT_BMI_PE_BAD_CHARACTER; \
      } \
      break; \
  } \
  GT_BMI_TEXT_PARSE_NUMBER(text_line,number); \
  if (is_negative) number = -number; \
}
#define GT_BMI_TEXT_SET_EOS__NEXT(text_line) (**text_line) = EOS; ++(*text_line)
#define GT_BMI_TEXT_SKIP_LINE(text_line) \
  while (__builtin_expect((**text_line)!='\n',1)) { \
    ++(*text_line); \
  }

/*
 * TODO: Erase THIS!!!!
 */
#define GT_BMI_READ_UNTIL(buffered_map_input,test) /*TODO: Deprecated */ \
  while (gt_expect_true(!(test) && buffered_map_input->cursor[0]!=EOL)) { \
    ++buffered_map_input->cursor; \
  }
#define GT_BMI_IS_EOL(buffered_map_input) gt_expect_false(buffered_map_input->cursor[0]==EOL) /* TODO: Deprecated */
#define GT_BMI_PARSE_NUMBER(buffered_map_input,number) /* TODO: Deprecated */ \
  number = 0; \
  while (gt_expect_true(gt_is_number(buffered_map_input->cursor[0]))) { \
    number = (number*10) + gt_get_cipher(buffered_map_input->cursor[0]); \
    ++buffered_map_input->cursor; \
  }
#define GT_BMI_SET_EOS__NEXT(buffered_map_input) buffered_map_input->cursor[0] = EOS; ++buffered_map_input->cursor /* TODO: Deprecated */
#define GT_BMI_SKIP_LINE(buffered_map_input) { \
  while (buffered_map_input->cursor[0]!=EOL) { \
    ++buffered_map_input->cursor; \
  } \
  buffered_map_input->cursor[0]=EOS; \
  ++buffered_map_input->cursor; \
  if (gt_expect_false(!gt_buffered_map_input_eob(buffered_map_input) && buffered_map_input->cursor[0]==EOS)) { \
    ++buffered_map_input->cursor; \
  } \
  ++buffered_map_input->current_line_num; \
}
GT_INLINE void gt_buffered_map_input_prompt_error(
    gt_buffered_map_input* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  register const char* const file_name = (buffered_map_input != NULL) ?
      buffered_map_input->input_file->file_name : "<<LazyParsing>>";
  if ((buffered_map_input == NULL)) {
    line_num = 0; column_pos = 0;
  }
  switch (error_code) {
    case 0: /* No error */ break;
    case GT_BMI_PE_WRONG_FILE_FORMAT: gt_error(PARSE_MAP_BAD_FILE_FORMAT,file_name,line_num); break;
    case GT_BMI_PE_NOT_IMPLEMENTED: gt_error(PARSE_MAP_NOT_IMPLEMENTED,file_name,line_num,column_pos); break;
    case GT_BMI_PE_PREMATURE_EOL: gt_error(PARSE_MAP_PREMATURE_EOL,file_name,line_num,column_pos); break;
    case GT_BMI_PE_BAD_NUMBER_OF_BLOCKS: /* (blocks(read)!=blocks(quals) || parse_alignment=>(num_blocks==1)) */
      gt_error(PARSE_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos);
      break;
    case GT_BMI_PE_BAD_CHARACTER: gt_error(PARSE_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_READ_BAD_CHARACTER: gt_error(PARSE_MAP_READ_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_QUAL_BAD_SEPARATOR: gt_error(PARSE_MAP_QUAL_BAD_SEPARATOR,file_name,line_num,column_pos); break;
    case GT_BMI_PE_QUAL_BAD_PREMATURE_EOB: gt_error(PARSE_MAP_QUAL_BAD_PREMATURE_EOB,file_name,line_num,column_pos); break;
    case GT_BMI_PE_QUAL_BAD_CHARACTER: gt_error(PARSE_MAP_QUAL_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_COUNTERS_BAD_CHARACTER: gt_error(PARSE_MAP_COUNTERS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_MAP_ALREADY_PARSED: gt_error(PARSE_MAP_MAP_ALREADY_PARSED,file_name,line_num); break;
    case GT_BMI_PE_MAP_BAD_NUMBER_OF_BLOCKS: gt_error(PARSE_MAP_MAP_BAD_NUMBER_OF_BLOCKS,file_name,line_num,column_pos); break;
    case GT_BMI_PE_MAP_BAD_CHARACTER: gt_error(PARSE_MAP_MAP_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_MISMS_ALREADY_PARSED: gt_error(PARSE_MAP_MISMS_ALREADY_PARSED,file_name,line_num); break;
    case GT_BMI_PE_MISMS_BAD_CHARACTER: gt_error(PARSE_MAP_MISMS_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_BMI_PE_MISMS_BAD_MISMS_POS: gt_error(PARSE_MAP_MISMS_BAD_MISMS_POS,file_name,line_num,column_pos); break;
    default:
      gt_error(PARSE_MAP,buffered_map_input->input_file->file_name,line_num);
      break;
  }
}
GT_INLINE void gt_buffered_map_input_next_record(gt_buffered_map_input* const buffered_map_input) {
  GT_BMI_CHECK(buffered_map_input);
  if (!gt_buffered_map_input_eob(buffered_map_input)) {
    GT_BMI_SKIP_LINE(buffered_map_input);
  }
}
GT_INLINE gt_status gt_bmi_parse_tag(
    gt_buffered_map_input* const buffered_map_input,char** tag,uint64_t* tag_length) {
  *tag = buffered_map_input->cursor;
  GT_BMI_READ_UNTIL(buffered_map_input,buffered_map_input->cursor[0]==TAB);
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  *tag_length = buffered_map_input->cursor-*tag;
  GT_BMI_SET_EOS__NEXT(buffered_map_input);
  return 0;
}
GT_INLINE gt_status gt_bmi_parse_tag_block(char** tag,char** tag_block,const uint64_t num_expected_blocks) {
  /*
   * Tag Block Conventions::
   *  (1) line1: NAME => line1&2: NAME
   *      line2: NAME
   *  (2) line1: NAME<SEP>SOMETHING_ELSE_1 => line1&2: NAME<SEP>SOMETHING_ELSE_1|NAME<SEP>SOMETHING_ELSE_2
   *      line2: NAME<SEP>SOMETHING_ELSE_2
   *  (3) line1: NAME<SEP>1 => line1&2: NAME
   *      line2: NAME<SEP>2
   * Tag Decomposition::
   *  (<TAG_BLOCK><SEP>){num_blocks}
   *    -> tag.num_blocks == read.num_blocks
   *    -> <SEP> := ' ' | '/' | '#' | "|"  // TODO
   */

  // TODO
  return 0;
}
GT_INLINE gt_status gt_bmi_parse_read_block(
    gt_buffered_map_input* const buffered_map_input,char** read_block,uint64_t* read_block_length) {
  *read_block = buffered_map_input->cursor;
  // Read READ_BLOCK
  while (gt_expect_true(buffered_map_input->cursor[0]!=TAB &&
      !gt_is_valid_template_separator(buffered_map_input->cursor[0]) &&
      buffered_map_input->cursor[0]!=EOL)) {
    if (gt_expect_false(!gt_is_dna(buffered_map_input->cursor[0]))) return GT_BMI_PE_READ_BAD_CHARACTER;
    ++buffered_map_input->cursor;
  }
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  *read_block_length = buffered_map_input->cursor-*read_block;
  // Return READ_BLOCK
  register gt_status return_status;
  if (buffered_map_input->cursor[0]==TAB) {
    return_status = GT_BMI_PE_EOB;
  } else if (gt_is_valid_template_separator(buffered_map_input->cursor[0])) {
    return_status = GT_BMI_PE_PENDING_BLOCKS;
  } else {
    return_status = GT_BMI_PE_READ_BAD_CHARACTER;
  }
  GT_BMI_SET_EOS__NEXT(buffered_map_input);
  return return_status;
}
GT_INLINE gt_status gt_bmi_parse_qualities_block(
    gt_buffered_map_input* const buffered_map_input,
    uint64_t next_qualities_block_length,char** qualities_block) {
  // Read QUAL_BLOCK
  register uint64_t num_characters = 0;
  *qualities_block = buffered_map_input->cursor;
  while (gt_expect_true(num_characters<next_qualities_block_length &&
      buffered_map_input->cursor[0]!=TAB &&
      buffered_map_input->cursor[0]!=EOL)) {
    if (gt_expect_false(!gt_is_valid_quality(buffered_map_input->cursor[0]))) {
      return GT_BMI_PE_QUAL_BAD_CHARACTER;
    }
    ++buffered_map_input->cursor; ++num_characters;
  }
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  if (gt_expect_false(num_characters<next_qualities_block_length)) return GT_BMI_PE_QUAL_BAD_PREMATURE_EOB;
  // Return QUAL_BLOCK
  register gt_status return_status;
  if (gt_is_valid_template_separator(buffered_map_input->cursor[0])) {
    return_status = GT_BMI_PE_PENDING_BLOCKS;
  } else if (buffered_map_input->cursor[0]==TAB) {
    return_status = GT_BMI_PE_EOB;
  } else {
    return GT_BMI_PE_QUAL_BAD_SEPARATOR;
  }
  GT_BMI_SET_EOS__NEXT(buffered_map_input);
  return return_status;
}
GT_INLINE gt_status gt_bmi_parse_counters(
    gt_buffered_map_input* const buffered_map_input,gt_vector* const counters,uint64_t* const mcs) {
  register uint64_t number, strata=1;
  *mcs = 0;
  gt_vector_clean(counters);
  while (gt_expect_true(buffered_map_input->cursor[0]!=TAB && buffered_map_input->cursor[0]!=EOL)) {
    if (gt_is_number(buffered_map_input->cursor[0])) {
      GT_BMI_PARSE_NUMBER(buffered_map_input,number);
      gt_vector_insert(counters,number,uint64_t);
    } else if (buffered_map_input->cursor[0]==GT_MAP_COUNTS_TIMES) { // 0x10:1:1
      if (gt_expect_false(gt_vector_get_used(counters)==0)) return GT_BMI_PE_COUNTERS_BAD_CHARACTER;
      register uint64_t multiplier, i;
      // Parse multiplier
      ++buffered_map_input->cursor;
      if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
      if (!gt_is_number(buffered_map_input->cursor[0])) return GT_BMI_PE_COUNTERS_BAD_CHARACTER;
      GT_BMI_PARSE_NUMBER(buffered_map_input,multiplier);
      number = *gt_vector_get_last_elm(counters,uint64_t);
      for (i=0;i<multiplier;++i) {
        gt_vector_insert(counters,number,uint64_t);
      }
      strata+=(multiplier-1);
    } else if (buffered_map_input->cursor[0]==GT_MAP_MCS) {
      strata++;
      *mcs = strata;
      ++buffered_map_input->cursor;
    } else if (buffered_map_input->cursor[0]==GT_MAP_COUNTS_SEP) {
      strata++;
      ++buffered_map_input->cursor;
    } else {
      return GT_BMI_PE_COUNTERS_BAD_CHARACTER;
    }
    // TODO: Consider 0:0::<Value>::<Value>
  }
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  GT_BMI_SET_EOS__NEXT(buffered_map_input);
  return 0;
}
// OLD (v0): <+3>20A88C89C99
GT_INLINE gt_status gt_bmi_parse_mismatch_string_v0(char** text_line,gt_map* map) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  if (gt_expect_false(*text_line==NULL)) return GT_BMI_PE_MISMS_ALREADY_PARSED;
  gt_map_clear_misms(map);
  // Parse Misms
  register uint64_t last_position = 0, last_cut_point = 0;  
  register const uint64_t global_length = map->base_length;
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP && (**text_line)!=EOL) {      
    gt_misms misms;
    if (gt_is_dna((**text_line))) { // Mismatch
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      GT_BMI_TEXT_NEXT(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,misms.position);
      if (gt_expect_false(misms.position<=last_position)) {
        return GT_BMI_PE_MISMS_BAD_MISMS_POS;
      }
      --misms.position; // Zero based position
      last_position = misms.position;
      // Add Mismatch
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='<') { // Indel
      register bool is_splice;
      GT_BMI_TEXT_NEXT(text_line);
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
          return GT_BMI_PE_MISMS_BAD_CHARACTER;
          break;
      }
      GT_BMI_TEXT_NEXT(text_line);
      // Parse size
      register uint64_t size, position;
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,size);
      // Parse Indel end ">"
      if (gt_expect_false((**text_line)!='>')) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_NEXT(text_line);
      // Parse Position
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,position);
      if (gt_expect_false(position<=last_position)) {
        return GT_BMI_PE_MISMS_BAD_MISMS_POS;
      }
      --position; // Zero based position
      last_position = position;
      // Add Indel
      if (gt_expect_true(!is_splice)) {
        misms.position = position;
        misms.size = size;
        gt_map_add_misms(map,&misms);
      } else { // Create a new map block
        gt_map* next_map = gt_map_new();
        next_map->seq_name = map->seq_name;
        next_map->position = map->position+position+size;
        next_map->direction = map->direction;
        next_map->base_length = global_length-position;        
        // Close current map block
        map->base_length = position - last_cut_point;
        last_cut_point = position;
        gt_map_set_next_block(map,next_map,SPLICE);
        // Swap maps & Reset length,position
        map = next_map;
      }
    } else { // ?¿ Parsing error
      return GT_BMI_PE_MISMS_BAD_CHARACTER;
    }
  }
  return 0;
}
// NEW (v1): (5)43T46A9>24*  ||  33C9T30T24>1-(10)
GT_INLINE gt_status gt_bmi_parse_mismatch_string_v1(char** text_line,gt_map* map) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  if (gt_expect_false(*text_line==NULL)) return GT_BMI_PE_MISMS_ALREADY_PARSED;
  gt_map_clear_misms(map);
  // Parse Misms
  register uint64_t position=0, length=0;
  while ((**text_line)!=GT_MAP_NEXT && (**text_line)!=GT_MAP_SEP && (**text_line)!=EOL) {
    gt_misms misms;
    if (gt_is_number((**text_line))) { // Matching
      register uint64_t matching_characters;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,matching_characters);
      position+=matching_characters;
      length+=matching_characters;
    } else if (gt_is_dna((**text_line))) { // Mismatch
      misms.misms_type = MISMS;
      misms.base = (**text_line);
      misms.position = position;
      ++position; ++length;
      GT_BMI_TEXT_NEXT(text_line);
      // Add Mismatch
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='(') { // Trim // FIXME: Only at the ends
      misms.misms_type = DEL;
      misms.position = position;
      GT_BMI_TEXT_NEXT(text_line);
      // Parse size
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,misms.size);
      position+=misms.size;
      // Parse Trim end ')'
      if (gt_expect_false((**text_line)!=')')) return GT_BMI_PE_MISMS_BAD_CHARACTER;
      GT_BMI_TEXT_NEXT(text_line);
      // Add Trim
      gt_map_add_misms(map,&misms);
    } else if ((**text_line)=='>') { // Indel/Skip
      GT_BMI_TEXT_NEXT(text_line);
      // Parse size
      register int64_t size;
      GT_BMI_TEXT_PARSE_SIGNED_NUMBER(text_line,size);
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
        GT_BMI_TEXT_NEXT(text_line);
        // Add Indel/Skip
        gt_map_add_misms(map,&misms);
      } else { // NSKIP/SPLICE
        register gt_junction_t junction;
        switch ((**text_line)) {
          case GT_MAP_SKIP_POSITIVE: junction=POSITIVE_SKIP; break;
          case GT_MAP_SKIP_NEGATIVE: junction=NEGATIVE_SKIP; break;
          case GT_MAP_SKIP_SPLICE: junction=SPLICE; break;
          default: return GT_BMI_PE_MISMS_BAD_CHARACTER; break;
        }
        // Create a new map block
        gt_map* next_map = gt_map_new();
        next_map->seq_name = map->seq_name;
        next_map->position = map->position+length+size;
        next_map->direction = map->direction;
        next_map->base_length = map->base_length-length;
        // Close current map block
        map->base_length = length;
        gt_map_set_next_block(map,next_map,junction);
        // Swap maps & Reset length,position
        map = next_map;
        position=0; length=0;
      }
    } else {
     return GT_BMI_PE_MISMS_BAD_CHARACTER;
    }
  }
  return 0;
}
GT_INLINE gt_status gt_bmi_parse_map(char** text_line,gt_map* const map,const gt_lazy_parse_mode parse_mode) {
  GT_NULL_CHECK(text_line);
  GT_MAP_CHECK(map);
  if (gt_expect_false(*text_line==NULL)) return GT_BMI_PE_MAP_ALREADY_PARSED;
  // Read TAG
  map->seq_name = *text_line;
  GT_BMI_TEXT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP);
  if (GT_BMI_TEXT_IS_EOL(text_line)) return GT_BMI_PE_PREMATURE_EOL;
  GT_BMI_TEXT_SET_EOS__NEXT(text_line);
  // Read Strand
  switch ((**text_line)) {
    case GT_MAP_STRAND_FORWARD_SYMBOL:
    case GT_MAP_STRAND_FORWARD_LETTER:
      map->direction = FORWARD;
    break;
    case GT_MAP_STRAND_REVERSE_SYMBOL:
    case GT_MAP_STRAND_REVERSE_LETTER:
      map->direction = REVERSE;
    break;
    default:
      return GT_BMI_PE_MAP_BAD_CHARACTER;
      break;
  }
  GT_BMI_TEXT_NEXT(text_line);
  // Determine format version
  if ((**text_line)==GT_MAP_SEP) { // GEMv1
    map->map_misms_format = MISMATCH_STRING_GEMv1;
    GT_BMI_TEXT_NEXT(text_line);
  } else if (gt_is_number((**text_line))) { // GEMv0
    map->map_misms_format = MISMATCH_STRING_GEMv0;
  } else { // ?¿
    return GT_BMI_PE_MAP_BAD_CHARACTER;
  }
  // Position
  if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MAP_BAD_CHARACTER;
  GT_BMI_TEXT_PARSE_NUMBER(text_line,map->position);
  // Synch with mismatch string (GEMv1)
  if (map->map_misms_format==MISMATCH_STRING_GEMv1) {
    if (gt_expect_false((**text_line)!=GT_MAP_SEP)) return GT_BMI_PE_MAP_BAD_CHARACTER;
    GT_BMI_TEXT_NEXT(text_line);
  }
  // Parse Mismatch String
  if (parse_mode==PARSE_ALL) {
    register gt_status error_code;
    map->mismatches_txt = NULL;
    if (map->map_misms_format==MISMATCH_STRING_GEMv1) {
      error_code=gt_bmi_parse_mismatch_string_v1(text_line,map);
    } else {
      error_code=gt_bmi_parse_mismatch_string_v0(text_line,map);
    }
    if (error_code) return error_code;
  } else {
    map->mismatches_txt = *text_line;
    GT_BMI_TEXT_READ_UNTIL(text_line,(**text_line)==GT_MAP_SEP || (**text_line)==GT_MAP_NEXT);
  }
  // Parse Quality Score (if any)
  if ((**text_line)==GT_MAP_SEP && (*(*text_line+1))!=GT_MAP_SEP) { // ':'
    GT_BMI_TEXT_NEXT(text_line);
    if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MAP_BAD_CHARACTER;
    GT_BMI_TEXT_PARSE_NUMBER(text_line,map->score);
  }
  // Detect next character (MAP,BLOCK,EOL)
  if ((**text_line)==GT_MAP_NEXT) { // ','
    GT_BMI_TEXT_NEXT(text_line);
    return GT_BMI_PE_MAP_PENDING_MAPS;
  } else if ((**text_line)==GT_MAP_SEP) { // ':'
    if ((*(*text_line+1))==GT_MAP_SEP) { // '::'
      if ((*(*text_line+2))==GT_MAP_SEP) { // ':::' (Attributes of the block group)
        return GT_BMI_PE_EOB;
      } else { // '::?'
        (*text_line)+=2;
        return GT_BMI_PE_PENDING_BLOCKS;
      }
    } else { // ':?'
      return GT_BMI_PE_MAP_BAD_CHARACTER;
    }
  } else if ((**text_line)==EOL) { // '\n'
    return GT_BMI_PE_EOB;
  } else {
    return GT_BMI_PE_MAP_BAD_CHARACTER;
  }
}
#define GT_BMI_PARSE_MAP_ERROR(error_code) \
  (error_code!=GT_BMI_PE_PENDING_BLOCKS && \
   error_code!=GT_BMI_PE_MAP_PENDING_MAPS && \
   error_code!=GT_BMI_PE_EOB )
// Formats allowed:
//   OLD (v0): chr7:F127708134G27T88::chr7:R127708509<+3>20A88C89C99
//   NEW (v1): chr11:-:51590050:(5)43T46A9>24*::chr11:-:51579823:33C9T30T24>1-(10)
GT_INLINE gt_status gt_bmi_parse_template_maps(
    char** text_line,gt_template* const template,
    const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_TEMPLATE_CONSISTENCY_CHECK(template);
  // Set as parsed (whatever the result is)
  template->maps_txt = NULL;
  // Check null maps
  if ((**text_line)==GT_MAP_NONE) {
    GT_BMI_TEXT_SKIP_LINE(text_line);
    return 0;
  }
  // Parse MAPS
  register gt_status error_code = GT_BMI_PE_MAP_PENDING_MAPS;
  register const uint64_t num_blocks_template = gt_vector_get_used(template->blocks);
  register uint64_t num_maps_parsed = 0, num_blocks_parsed;
  register gt_vector* vector_maps = gt_vector_new(num_blocks_template,sizeof(gt_map*));
  gt_mmap_attributes mmap_attr;
  while (error_code==GT_BMI_PE_MAP_PENDING_MAPS && num_maps_parsed<num_maps) {
    // Parse MAP
    error_code = GT_BMI_PE_PENDING_BLOCKS;
    gt_vector_clean(vector_maps);
    num_blocks_parsed = 0;
    while (error_code==GT_BMI_PE_PENDING_BLOCKS) {
      // Allocate new map
      gt_map* map = gt_map_new();
      gt_vector_insert(vector_maps,map,gt_map*);
      // Set base length (needed to calculate the alignment's length in GEMv0)
      gt_map_set_base_length(map,gt_template_get_block(template,num_blocks_parsed)->read_length);
      // Parse current MAP
      error_code = gt_bmi_parse_map(text_line,map,parse_mode);
      if (GT_BMI_PARSE_MAP_ERROR(error_code)) return error_code;
      ++num_blocks_parsed;
    }
    // Check number of blocks parsed
    if (gt_expect_false(num_blocks_parsed<num_blocks_template)) return GT_BMI_PE_MAP_BAD_NUMBER_OF_BLOCKS;
    if (gt_expect_false(num_blocks_parsed>num_blocks_template)) { // Weird case of more blocks than blocks in the template (reorganize blocks)
      return GT_BMI_PE_NOT_IMPLEMENTED; // FIXME
    }
    // Add MAPs to corresponding alignments
    GT_VECTOR_ITERATE(vector_maps,map_ptr,map_pos,gt_map*) {
      register gt_alignment* const alignment = gt_template_get_block(template,map_pos);
      gt_alignment_insert_map(alignment,*map_ptr);
    }
    // Parse quality (if any) and calculate template attributes
    gt_template_clear_mmap_attributes(&mmap_attr);
    mmap_attr.distance = gt_map_vector_get_distance(vector_maps);
    if ((**text_line)==GT_MAP_SEP && (*(*text_line+1))==GT_MAP_SEP && (*(*text_line+2))==GT_MAP_SEP) { // ':::'
      (*text_line)+=3;
      if (gt_expect_false(!gt_is_number((**text_line)))) return GT_BMI_PE_MAP_BAD_CHARACTER;
      GT_BMI_TEXT_PARSE_NUMBER(text_line,mmap_attr.score);
    }
    // Store MAP blocks parsed
    gt_template_add_mmap_gtvector(template,vector_maps,&mmap_attr);
    ++num_maps_parsed;
  }
  gt_vector_delete(vector_maps);
  return 0;
}
// Formats allowed:
//   OLD (v0): chr7:F127708134G27T88
//   NEW (v1): chr11:-:51590050:(5)43T46A9>24*
GT_INLINE gt_status gt_bmi_parse_alignment_maps(
    char** text_line,gt_alignment* alignment,
    const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK((*text_line));
  GT_ALIGNMENT_CHECK(alignment);
  // Set as parsed (whatever the result is)
  alignment->maps_txt = NULL;
  // Check null maps
  if ((**text_line)!=GT_MAP_NONE) {
    GT_BMI_TEXT_SKIP_LINE(text_line);
    return 0;
  }
  // Parse MAPS
  register const uint64_t alignment_base_length = alignment->read_length;
  register uint64_t num_maps_parsed = 0;
  register gt_status error_code = GT_BMI_PE_MAP_PENDING_MAPS;
  while (error_code==GT_BMI_PE_MAP_PENDING_MAPS && num_maps_parsed<num_maps) {
    register gt_map* map = gt_map_new();
    // Set base length (needed to calculate the alignment's length in GEMv0)
    gt_map_set_base_length(map,alignment_base_length);
    // Parse current MAP
    error_code = gt_bmi_parse_map(text_line,map,parse_mode);
    if (error_code==GT_BMI_PE_PENDING_BLOCKS) return GT_BMI_PE_MAP_BAD_NUMBER_OF_BLOCKS; /* TODO: Weird case of split blocks (FIXME) */
    gt_alignment_add_map(alignment,map);
    ++num_maps_parsed;
  }
  return 0;
}

///*
// * TODO
// * Change all above to take char** so as to parse whatever (not only buffered files)
// */
// TODO:
//// Parse Template/Alignment
//GT_INLINE gt_status gt_buffered_map_input_parse_template(
//    gt_buffered_map_input* const buffered_map_input,gt_template* const template,
//    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps);
//GT_INLINE gt_status gt_buffered_map_input_parse_alignment(
//    gt_buffered_map_input* const buffered_map_input,gt_alignment* alignment,
//    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps);


/*
 * MAP/MAPQ/MMAP/MMAPQ Lazy Parsers
 *   Lazy parsing works in 3 steps
 *     (1) Parses TAG(s),READ(s),QUALITIES(s),COUNTERS -> gt_template/gt_alignment
 *     (2) Parses MAPS' SEQUENCE,STRAND,POSITION
 *     (3) Parses MAPS' CIGAR string
 */
GT_INLINE gt_status gt_buffered_map_input_parse_template(
    gt_buffered_map_input* const buffered_map_input,gt_template* const template,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_BMI_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  // TAG
  if ((error_code=gt_bmi_parse_tag(buffered_map_input,&template->tag,&template->tag_length))) {
    return error_code;
  }
  // READ
  register gt_input_file* const input_file = buffered_map_input->input_file;
  register uint64_t num_blocks = 0;
  error_code=GT_BMI_PE_PENDING_BLOCKS;
  while (error_code==GT_BMI_PE_PENDING_BLOCKS) {
    gt_alignment* const alignment = gt_template_dyn_get_block(template,num_blocks);
    error_code=gt_bmi_parse_read_block(buffered_map_input,&alignment->read,&alignment->read_length);
    if (error_code!=GT_BMI_PE_PENDING_BLOCKS && error_code!=GT_BMI_PE_EOB) return error_code;
    ++num_blocks;
  }

  // Tag Splitting (try to deduce alignments' tag out of the one template's tag) // TODo
  //gt_bmi_parse_tag_block(char** tag,char** tag_block,const uint64_t num_expected_blocks);

  // QUALITIES
  if (input_file->map_type.contains_qualities) {
    register uint64_t i;
    error_code=GT_BMI_PE_PENDING_BLOCKS;
    for (i=0;i<num_blocks;++i) {
      if (error_code!=GT_BMI_PE_PENDING_BLOCKS) return GT_BMI_PE_BAD_NUMBER_OF_BLOCKS;
      gt_alignment* alignment = gt_template_get_block(template,i);
      error_code=gt_bmi_parse_qualities_block(
          buffered_map_input,alignment->read_length,&alignment->qualities);
      if (error_code!=GT_BMI_PE_PENDING_BLOCKS && error_code!=GT_BMI_PE_EOB) return error_code;
    }
    if (error_code!=GT_BMI_PE_EOB) return GT_BMI_PE_BAD_NUMBER_OF_BLOCKS;
  }
  // COUNTERS
  if ((error_code=gt_bmi_parse_counters(buffered_map_input,
      template->counters,&template->max_complete_strata))) return error_code;
  // MAPS (lazy parsing)
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  if (parse_mode!=PARSE_READ) {
    template->maps_txt = NULL;
    error_code = gt_bmi_parse_template_maps(&(buffered_map_input->cursor),template,parse_mode,num_maps);
  } else {
    template->maps_txt = buffered_map_input->cursor;
    error_code = 0;
  }
  return error_code;
}
GT_INLINE gt_status gt_buffered_map_input_parse_alignment(
    gt_buffered_map_input* const buffered_map_input,gt_alignment* alignment,
    const bool has_map_quality,const gt_lazy_parse_mode parse_mode,uint64_t num_maps) {
  GT_BMI_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  // TAG
  if ((error_code=gt_bmi_parse_tag(buffered_map_input,&alignment->tag,&alignment->tag_length))) return error_code;
  // READ
  error_code=gt_bmi_parse_read_block(buffered_map_input,&alignment->read,&alignment->read_length);
  if (gt_expect_false(error_code==GT_BMI_PE_PENDING_BLOCKS)) return GT_BMI_PE_BAD_NUMBER_OF_BLOCKS;
  if (gt_expect_false(error_code!=GT_BMI_PE_EOB)) return error_code;
  // QUALITIES
  if (buffered_map_input->input_file->map_type.contains_qualities) {
    error_code=gt_bmi_parse_qualities_block(buffered_map_input,alignment->read_length,&alignment->qualities);
    if (gt_expect_false(error_code==GT_BMI_PE_PENDING_BLOCKS)) return GT_BMI_PE_BAD_NUMBER_OF_BLOCKS;
    if (gt_expect_false(error_code!=GT_BMI_PE_EOB)) return error_code;
  }
  // COUNTERS
  if ((error_code=gt_bmi_parse_counters(buffered_map_input,
      alignment->counters,&alignment->max_complete_strata))) return error_code;
  // MAPS (lazy parsing)
  if (GT_BMI_IS_EOL(buffered_map_input)) return GT_BMI_PE_PREMATURE_EOL;
  if (parse_mode!=PARSE_READ) {
    alignment->maps_txt = NULL;
    error_code=gt_bmi_parse_alignment_maps(&(buffered_map_input->cursor),alignment,parse_mode,num_maps);
  } else {
    alignment->maps_txt = buffered_map_input->cursor;
    error_code=0;
  }
  return error_code;
}
// Parse Maps
GT_INLINE gt_status gt_buffered_map_input_parse_template_maps(gt_template* template,uint64_t num_maps) {
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  error_code = gt_bmi_parse_template_maps(&template->maps_txt,template,PARSE_READ__MAPS,num_maps);
  template->maps_txt = NULL;
  return error_code;
}
GT_INLINE gt_status gt_buffered_map_input_parse_alignment_maps(gt_alignment* alignment,uint64_t num_maps) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  error_code = gt_bmi_parse_alignment_maps(&alignment->maps_txt,alignment,PARSE_READ__MAPS,num_maps);
  alignment->maps_txt = NULL;
  return error_code;
}
// Parse Mismatches
GT_INLINE gt_status gt_buffered_map_input_parse_template_mismatch_string(gt_template* template) {
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  gt_template_maps_iterator template_maps_iterator;
  gt_template_new_maps_iterator(template,&template_maps_iterator);
  register const uint64_t num_blocks_template = gt_vector_get_used(template->blocks);
  gt_map** map_array;
  while (gt_template_next_maps(&template_maps_iterator,&map_array)) {
    register uint64_t i;
    for (i=0;i<num_blocks_template;++i) {
      if ((error_code = (map_array[i]->map_misms_format==MISMATCH_STRING_GEMv1) ?
          gt_bmi_parse_mismatch_string_v1(&(map_array[i]->mismatches_txt),map_array[i]):
          gt_bmi_parse_mismatch_string_v0(&(map_array[i]->mismatches_txt),map_array[i]))) {
        return error_code;
      }
    }
  }
  free(map_array);
  return 0;
}
GT_INLINE gt_status gt_buffered_map_input_parse_alignment_mismatch_string(gt_alignment* alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  gt_map* map;
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    if ((error_code = (map->map_misms_format==MISMATCH_STRING_GEMv1) ?
        gt_bmi_parse_mismatch_string_v1(&map->mismatches_txt,map):
        gt_bmi_parse_mismatch_string_v0(&map->mismatches_txt,map))) {
      return error_code;
    }
  }
  return 0;
}
/*
 * MAP/MAPQ/MMAP/MMAPQ High-level Parsers
 *   - High-level parsing to extract one template/alignment from the buffered file (reads one line)
 *   - Syntax checking
 *   - Transparent buffer block reload
 *   - Template/Alignment transparent memory management
 */
GT_INLINE gt_status gt_bmi_get_template(
    gt_buffered_map_input* const buffered_map_input,gt_template* template,
    const gt_lazy_parse_mode parse_mode) {
  GT_BMI_CHECK(buffered_map_input);
  GT_TEMPLATE_CHECK(template);
  register const char* line_start = buffered_map_input->cursor;
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_map_input_eob(buffered_map_input)) {
    register const uint64_t read_lines = gt_buffered_map_input_get_block(buffered_map_input);
    if (gt_expect_false(read_lines==0)) return GT_BMI_EOF;
  }
  // Check file format
  register gt_input_file* input_file = buffered_map_input->input_file;
  if (gt_buffered_map_input_check_map_file_format(buffered_map_input)) {
    gt_error(PARSE_MAP_BAD_FILE_FORMAT,input_file->file_name,buffered_map_input->current_line_num);
    gt_buffered_map_input_next_record(buffered_map_input);
    return GT_BMI_FAIL;
  }
  // Prepare the template
  register const uint64_t line_num = buffered_map_input->current_line_num;
  gt_template_clear(template);
  template->template_id = line_num;
  // Parse template
  if ((error_code=gt_buffered_map_input_parse_template(buffered_map_input,
      template,input_file->map_type.contains_qualities,parse_mode,UINT64_MAX))) {
    gt_buffered_map_input_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_buffered_map_input_next_record(buffered_map_input);
    return GT_BMI_FAIL;
  }
  gt_buffered_map_input_next_record(buffered_map_input);
//  // Parse ALL template's maps
//  if (parse_mode==PARSE_READ) return GT_BMI_OK; // Lazy
//  if ((error_code=gt_buffered_map_input_parse_template_maps(template,UINT64_MAX))) {
//    gt_buffered_map_input_prompt_error(buffered_map_input,error_code);
//    return GT_BMI_FAIL;
//  }
//  // Parse ALL mismatch strings
//  if (parse_mode==PARSE_READ__MAPS) return GT_BMI_OK; // Lazy
//  if ((error_code=gt_buffered_map_input_parse_template_mismatch_string(template))) {
//    gt_buffered_map_input_prompt_error(buffered_map_input,error_code);
//    return GT_BMI_FAIL;
//  }
  return GT_BMI_OK;
}
GT_INLINE gt_status gt_bmi_get_alignment(
    gt_buffered_map_input* const buffered_map_input,gt_alignment* alignment,
    const gt_lazy_parse_mode parse_mode) {
  GT_BMI_CHECK(buffered_map_input);
  GT_ALIGNMENT_CHECK(alignment);
  register const char* line_start = buffered_map_input->cursor;
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_map_input_eob(buffered_map_input)) {
    register const uint64_t read_lines = gt_buffered_map_input_get_block(buffered_map_input);
    if (gt_expect_false(read_lines==0)) return GT_BMI_EOF;
  }
  // Check file format
  register gt_input_file* input_file = buffered_map_input->input_file;
  if (gt_buffered_map_input_check_map_file_format(buffered_map_input)) {
    gt_error(PARSE_MAP_BAD_FILE_FORMAT,input_file->file_name,buffered_map_input->current_line_num);
    gt_buffered_map_input_next_record(buffered_map_input);
    return GT_BMI_FAIL;
  }
  // Allocate memory for the alignment
  register const uint64_t line_num = buffered_map_input->current_line_num;
  gt_alignment_clear(alignment);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_buffered_map_input_parse_alignment(buffered_map_input,
      alignment,input_file->map_type.contains_qualities,parse_mode,UINT64_MAX))) {
    gt_buffered_map_input_prompt_error(buffered_map_input,line_num,
        buffered_map_input->cursor-line_start,error_code);
    gt_buffered_map_input_next_record(buffered_map_input);
    return GT_BMI_FAIL;
  }
  gt_buffered_map_input_next_record(buffered_map_input);
//  // Parse ALL alignment's maps
//  if (parse_mode==PARSE_READ) return GT_BMI_OK; // Lazy
//  if ((error_code=gt_buffered_map_input_parse_alignment_maps(alignment,UINT64_MAX))) {
//    gt_buffered_map_input_prompt_error(buffered_map_input,error_code);
//    return GT_BMI_FAIL;
//  }
//  // Parse ALL mismatch strings
//  if (parse_mode==PARSE_READ__MAPS) return GT_BMI_OK; // Lazy
//  if ((error_code=gt_buffered_map_input_parse_alignment_mismatch_string(alignment))) {
//    gt_buffered_map_input_prompt_error(buffered_map_input,error_code);
//    return GT_BMI_FAIL;
//  }
  return GT_BMI_OK;
}
GT_INLINE gt_status gt_buffered_map_input_get_template(
    gt_buffered_map_input* const buffered_map_input,gt_template* template) {
  return gt_bmi_get_template(buffered_map_input,template,PARSE_ALL);
}
GT_INLINE gt_status gt_buffered_map_input_get_alignment(
    gt_buffered_map_input* const buffered_map_input,gt_alignment* alignment) {
  return gt_bmi_get_alignment(buffered_map_input,alignment,PARSE_ALL);
}
