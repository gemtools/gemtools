/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.c
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_sam_parser.h"


// Constants
#define GT_ISP_NUM_LINES GT_NUM_LINES_10K
#define GT_ISP_NUM_INITIAL_MAPS 5

// Internal pair-pending
typedef struct {
  // Current map info
  char* map_seq_name;
  uint64_t map_position;
  uint64_t end_position; // 0/1
  // Next map info
  char* next_seq_name;
  uint64_t next_position;
  // Map location and span info
  uint64_t map_displacement; // In alignment's map vector
  uint64_t num_maps; // Maps in the vector coupled to the first one
} gt_sam_pending_end;

/*
 * SAM File Format test
 */
GT_INLINE bool gt_input_file_detect_sam_format(char* const buffer,const uint64_t buffer_size,const bool show_errors) {
  /*
   * After 2 header lines or 1 SAM line we are fine to say this is SAM-format
   * By default, like 16MB of buffered data should be fine to assess the SAM format
   *  (1) @SQ     SN:chr10        LN:135534747
   *      @SQ     SN:chr11        LN:135006516
   *  (2) 1/1     16  chr12  57338496  37  75M  *  0  0  TCTGGTT...TTTGN  _b....VNNQQB  XT:A:U  NM:i:1
   */
  // TODO
  return true;
}
GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors) {
  GT_INPUT_FILE_CHECK(input_file);
  if (gt_input_file_detect_sam_format((char*)input_file->file_buffer,input_file->buffer_size,show_errors)) {
    register uint64_t lines_read = 0;
    while (!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==GT_SAM_HEADER_BEGIN) {

      // TODO => sam_headers
      while (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)!=EOL)) { // FIXM: DOS_EOL
        GT_INPUT_FILE_NEXT_CHAR(input_file,NULL);
        GT_INPUT_FILE_CHECK__FILL_BUFFER(input_file,NULL);
      }
      if (!input_file->eof) {
        GT_INPUT_FILE_NEXT_CHAR(input_file,NULL); // Check DOS EOF
        if (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==DOS_EOL)) {
          GT_INPUT_FILE_NEXT_CHAR(input_file,NULL);
        }
      }
      ++lines_read;
      // TODO => sam_headers

    }
    // Skip header from input-file buffer
    input_file->buffer_begin = input_file->buffer_pos;
    input_file->processed_lines = lines_read;
    return true;
  } else {
    return false;
  }
}
GT_INLINE gt_status gt_input_sam_parser_check_sam_file_format(gt_buffered_input_file* const buffered_sam_input) {
  register gt_input_file* const input_file = buffered_sam_input->input_file;
  if (gt_expect_false(input_file->file_format==FILE_FORMAT_UNKNOWN)) { // Unknown
    gt_sam_headers sam_headers;
    // Mutex format detection (because the first one must read the headers)
    gt_input_file_lock(input_file);
      register const bool is_sam_format =
          gt_input_file_test_sam(input_file,&sam_headers,true);
    gt_input_file_unlock(input_file);
    if (!is_sam_format) return GT_ISP_PE_WRONG_FILE_FORMAT;
    input_file->file_format = SAM;
  } else if (gt_expect_false(input_file->file_format!=SAM)) {
    return GT_ISP_PE_WRONG_FILE_FORMAT;
  }
  return 0;
}

/*
 * SAM File basics
 */
/* Error handler */
GT_INLINE void gt_input_sam_parser_prompt_error(
    gt_buffered_input_file* const buffered_sam_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code) {
  // Display textual error msg
  register const char* const file_name = (buffered_sam_input != NULL) ?
      buffered_sam_input->input_file->file_name : "<<LazyParsing>>";
  if ((buffered_sam_input == NULL)) {
    line_num = 0; column_pos = 0;
  }
  switch (error_code) {
    case 0: /* No error */ break;
    case GT_ISP_PE_WRONG_FILE_FORMAT: gt_error(PARSE_SAM_BAD_FILE_FORMAT,file_name,line_num,column_pos); break;
    case GT_ISP_PE_PREMATURE_EOL: gt_error(PARSE_SAM_PREMATURE_EOL,file_name,line_num,column_pos); break;
    case GT_ISP_PE_EXPECTED_NUMBER: gt_error(PARSE_SAM_EXPECTED_NUMBER,file_name,line_num,column_pos); break;
    case GT_ISP_PE_BAD_CHARACTER: gt_error(PARSE_SAM_BAD_CHARACTER,file_name,line_num,column_pos); break;
    case GT_ISP_PE_WRONG_READ_CONTENT: gt_error(PARSE_SAM_WRONG_READ_CONTENT,file_name,line_num,column_pos); break;
    case GT_ISP_PE_CIGAR_PREMATURE_END: gt_error(PARSE_SAM_CIGAR_PREMATURE_END,file_name,line_num,column_pos); break;
    case GT_ISP_PE_SAM_UNMAPPED_XA: gt_error(PARSE_SAM_UNMAPPED_XA,file_name,line_num,column_pos); break;
    case GT_ISP_PE_WRONG_NUM_XA: gt_error(PARSE_SAM_WRONG_NUM_XA,file_name,line_num,column_pos); break;
    case GT_ISP_PE_UNSOLVED_PENDING_MAPS: gt_error(PARSE_SAM_UNSOLVED_PENDING_MAPS,file_name,line_num,column_pos); break;
    default:
      gt_error(PARSE_SAM,file_name,line_num,column_pos);
      break;
  }
}
/* MAP file. Skip record */
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_sam_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  if (!gt_buffered_input_file_eob(buffered_sam_input)) {
    GT_INPUT_FILE_SKIP_LINE(buffered_sam_input);
  }
}
/*
 * SAM file. Reload internal buffer
 */
/* SAM file. Synchronized get block wrt to sam records */
GT_INLINE gt_status gt_input_sam_parser_get_block(
    gt_buffered_input_file* const buffered_sam_input,const uint64_t num_records) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  register gt_input_file* const input_file = buffered_sam_input->input_file;
  // Read lines
  if (input_file->eof) return GT_BMI_EOF;
  gt_input_file_lock(input_file);
  if (input_file->eof) {
    gt_input_file_unlock(input_file);
    return GT_BMI_EOF;
  }
  buffered_sam_input->block_id = gt_input_file_next_id(input_file) % UINT32_MAX;
  buffered_sam_input->current_line_num = input_file->processed_lines+1;
  gt_vector_clean(buffered_sam_input->block_buffer); // Clear dst buffer
  // Read lines & synch SAM records
  uint64_t lines_read = 0;
  while (lines_read<num_records &&
      gt_input_file_next_sam_record(input_file,buffered_sam_input->block_buffer,NULL) ) ++lines_read;
  if (lines_read==num_records) { // !EOF, Synch wrt to tag content
    register gt_string* const reference_tag = gt_string_new(30);
    if (gt_input_file_next_sam_record(input_file,buffered_sam_input->block_buffer,reference_tag)) {
      gt_fastq_tag_chomp_end_info(reference_tag);
      while (gt_input_file_cmp_next_sam_record(input_file,reference_tag)) {
        if (!gt_input_file_next_sam_record(input_file,buffered_sam_input->block_buffer,NULL)) break;
        ++lines_read;
      }
    }
    gt_string_delete(reference_tag);
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffered_sam_input->block_buffer);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffered_sam_input->block_buffer,char) != EOL) {
    gt_vector_insert(buffered_sam_input->block_buffer,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  buffered_sam_input->lines_in_buffer = lines_read;
  gt_input_file_unlock(input_file);

  // Setup the block
  buffered_sam_input->cursor = gt_vector_get_mem(buffered_sam_input->block_buffer,char);
  return buffered_sam_input->lines_in_buffer;
}
/* SAM file. Reload internal buffer */
GT_INLINE gt_status gt_input_sam_parser_reload_buffer(gt_buffered_input_file* const buffered_sam_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  // Dump buffer if BOF it attached to SAM-input, and get new out block (always FIRST)
  if (buffered_sam_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_dump(buffered_sam_input->buffered_output_file);
  }
  // Read new input block
  register const uint64_t read_lines =
      gt_input_sam_parser_get_block(buffered_sam_input,GT_ISP_NUM_LINES);
  if (gt_expect_false(read_lines==0)) return GT_ISP_EOF;
  // Assign block ID
  if (buffered_sam_input->buffered_output_file!=NULL) {
    gt_buffered_output_file_set_block_ids(
        buffered_sam_input->buffered_output_file,buffered_sam_input->block_id,0);
  }
  return GT_ISP_OK;
}
/*
 * SAM format. Basic building block for parsing
 */
GT_INLINE char* gt_isp_read_tag(char** const text_line,gt_string* const tag) {
  // Save text_line state
  char* text_cp = *text_line;
  register char** const ptext_cp = &text_cp;
  // Read tag
  register char* const tag_begin = *ptext_cp;
  GT_READ_UNTIL(ptext_cp,**ptext_cp==TAB || **ptext_cp==SPACE);
  if (GT_IS_EOL(ptext_cp)) return NULL;
  // Set tag
  register uint64_t const tag_length = *ptext_cp-tag_begin;
  gt_string_set_nstring(tag,tag_begin,tag_length);
  // Read the rest till next field
  if (**ptext_cp==SPACE) {
    GT_READ_UNTIL(ptext_cp,**ptext_cp==TAB);
    if (GT_IS_EOL(ptext_cp)) return NULL;
  }
  GT_NEXT_CHAR(ptext_cp);
  return *ptext_cp;
}
GT_INLINE gt_status gt_isp_parse_sam_cigar(char** const text_line,gt_map* const map) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_MAP_CHECK(map);
  register uint64_t position = 0;
  if (**text_line==STAR) {
    GT_NEXT_CHAR(text_line);
    return 0;
  }
  // 5M1D95M3I40M
  while (**text_line!=TAB && **text_line!=EOL) {
    // Parse misms_op length
    register uint64_t length;
    if (!gt_is_number(**text_line)) return GT_ISP_PE_EXPECTED_NUMBER;
    GT_PARSE_NUMBER(text_line,length);
    // Parse misms_op
    if (gt_expect_false(**text_line==EOL || **text_line==TAB)) return GT_ISP_PE_CIGAR_PREMATURE_END;
    gt_misms misms;
    register const char cigar_op = **text_line;
    GT_NEXT_CHAR(text_line);
    switch (cigar_op) {
      case 'M':
      case '=':
      case 'X':
        position += length;
        break;
      case 'P': // Padding. Nothing specific implemented
      case 'S': // Soft clipping. Nothing specific implemented
      case 'H': // Hard clipping. Nothing specific implemented
      case 'I': // Insertion to the reference
        misms.misms_type = DEL;
        misms.position = position;
        misms.size = length;
        position += length;
        gt_map_add_misms(map,&misms);
        break;
      case 'D': // Deletion from the reference
        misms.misms_type = INS;
        misms.position = position;
        misms.size = length;
        gt_map_add_misms(map,&misms);
        break;
        break;
      case 'N':
        gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      default:
        return GT_ISP_PE_BAD_CHARACTER;
        break;
    }
  }
  gt_map_set_base_length(map,position);
  return 0;
}

#define GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL() \
  if (GT_IS_EOL(text_line)) { \
    gt_map_delete(map); return GT_ISP_PE_PREMATURE_EOL; \
  }
#define GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT() \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL(); \
  GT_NEXT_CHAR(text_line)
#define GT_ISP_PARSE_SAM_ALG_SKIP_FIELD() \
  GT_READ_UNTIL(text_line,**text_line==TAB); \
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT()
#define GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(number) \
  if (!gt_is_number(**text_line)) { gt_map_delete(map); return GT_ISP_PE_EXPECTED_NUMBER; } \
  GT_PARSE_NUMBER(text_line,number)
#define GT_ISP_IF_OPT_FIELD(text_line,char1,char2,type_char) { \
  if ((*text_line)[0]==char1 && (*text_line)[1]==char2 && (*text_line)[2]!=EOL && (*text_line)[3]==type_char) {
#define GT_ISP_END_OPT_FIELD \
    return 0; \
  }}


GT_INLINE gt_status gt_isp_parse_sam_optional_field(
    char** const text_line,gt_alignment* const alignment,
    gt_vector* const maps_vector,gt_sam_pending_end* const pending,
    const bool is_mapped) {
  /*
   * XA:Z:chr17,-34553512,125M,0;chr17,-34655077,125M,0;
   */
  GT_ISP_IF_OPT_FIELD(text_line,'X','A','Z') {
    if (!is_mapped) return GT_ISP_PE_SAM_UNMAPPED_XA;
    *text_line+=5;
    while (**text_line!=TAB && **text_line!=EOL) { // Read new attached maps
      gt_map* map = gt_map_new();
      gt_map_set_base_length(map,gt_alignment_get_read_length(alignment));
      // Sequence-name/Chromosome
      register char* const seq_name = *text_line;
      GT_READ_UNTIL(text_line,**text_line==COMA);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      gt_map_set_seq_name(map,seq_name,*text_line-seq_name);
      GT_NEXT_CHAR(text_line);
      // Position
      if (**text_line==MINUS) {
        gt_map_set_strand(map,REVERSE);
        GT_NEXT_CHAR(text_line);
      } else if (**text_line==PLUS) {
        gt_map_set_strand(map,FORWARD);
        GT_NEXT_CHAR(text_line);
      } else {
        gt_map_delete(map);
        return GT_ISP_PE_BAD_CHARACTER;
      }
      GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->position);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      GT_NEXT_CHAR(text_line);
      // CIGAR // TODO: Parse it !!
      GT_READ_UNTIL(text_line,**text_line==COMA);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      GT_NEXT_CHAR(text_line);
      // Edit distance
      GT_READ_UNTIL(text_line,**text_line==SEMICOLON);
      if (**text_line==SEMICOLON) GT_NEXT_CHAR(text_line);
      // Add it to the list
      gt_vector_insert(maps_vector,map,gt_map*);
      // Consider pending relation
      ++pending->num_maps;
    }
  } GT_ISP_END_OPT_FIELD;

  // TODO: MD field, and more ....

  // TODO: Store in attributes

  // Unknown field (skip it)
  GT_READ_UNTIL(text_line,**text_line==TAB);
  return 0;
}

// TODO: Increase the level of checking SAM consistency
GT_INLINE gt_status gt_isp_parse_sam_alignment(
    char** const text_line,gt_template* const _template,gt_alignment* const _alignment,
    uint64_t* const alignment_flag,gt_sam_pending_end* const pending,const bool override_pairing) {
  register gt_status error_code;
  register bool is_mapped = true, is_single_segment;
  register gt_map* map = gt_map_new();
  /*
   * Parse FLAG
   */
  if (!gt_is_number(**text_line)) {
    gt_map_delete(map);
    return GT_ISP_PE_EXPECTED_NUMBER;
  }
  GT_PARSE_NUMBER(text_line,*alignment_flag);
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  // Process flags
  is_mapped = !(*alignment_flag&GT_SAM_FLAG_UNMAPPED);
  is_single_segment = override_pairing || !(*alignment_flag&GT_SAM_FLAG_MULTIPLE_SEGMENTS);
  pending->end_position = (is_single_segment) ? 0 : ((*alignment_flag&GT_SAM_FLAG_FIRST_SEGMENT)?0:1);
  if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT)  {
    gt_map_set_strand(map,REVERSE);
  } else {
    gt_map_set_strand(map,FORWARD);
  }

  // Allocate template/alignment handlers
  register gt_alignment* alignment;
  if (_template) {
    alignment = gt_template_get_block_dyn(_template,0);
    if (pending->end_position==1) {
      alignment = gt_template_get_block_dyn(_template,1);
    }
  } else {
    GT_NULL_CHECK(_alignment);
    alignment = _alignment;
  }
  if (!gt_alignment_get_attr(alignment,GT_ATTR_SAM_FLAGS)) {
    gt_alignment_set_attr(alignment,GT_ATTR_SAM_FLAGS,alignment_flag,sizeof(uint64_t));
  }
  /*
   * Parse RNAME (Sequence-name/Chromosome)
   */
  if (gt_expect_false(**text_line==STAR)) {
    is_mapped=false; /* Unmapped */
    GT_NEXT_CHAR(text_line);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  } else {
    register char* const seq_name = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    gt_map_set_seq_name(map,seq_name,*text_line-seq_name);
    GT_NEXT_CHAR(text_line);
  }
  /*
   * Parse POS (Position 1-based)
   */
  GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->position);
  if (map->position==0) is_mapped=false; /* Unmapped */
  GT_NEXT_CHAR(text_line);
  /*
   * Parse MAPQ (Score)
   */
  GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->score);
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  /*
   * Parse CIGAR
   */
  if ((error_code=gt_isp_parse_sam_cigar(text_line,map))) {
    gt_map_delete(map); return error_code;
  }
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  /*
   * Parse RNEXT (Sequence-name of the next segment)
   */
  if (**text_line==STAR || is_single_segment || !is_mapped ||
      (*alignment_flag&GT_SAM_FLAG_NEXT_UNMAPPED)) {
    pending->next_seq_name = NULL;
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD(); // RNEXT
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD(); // PNEXT
  } else {
    // Parse RNEXT
    if (**text_line==EQUAL) {
      pending->next_seq_name = gt_map_get_seq_name(map);
      GT_NEXT_CHAR(text_line);
      if (**text_line!=TAB) {
        gt_map_delete(map); return GT_ISP_PE_BAD_CHARACTER;
      }
      GT_NEXT_CHAR(text_line);
    } else {
      pending->next_seq_name = *text_line;
      GT_READ_UNTIL(text_line,**text_line==TAB);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      **text_line=EOS; // Writing to buffer's info
      GT_NEXT_CHAR(text_line);
    }
    // Parse PNEXT (Position of the next segment)
    GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(pending->next_position);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
    if (pending->next_position==0) {
      pending->next_seq_name = NULL;
    } else {
      pending->num_maps = 1;
      pending->map_seq_name = gt_map_get_seq_name(map);
      pending->map_position = gt_map_get_position(map);
    }
  }
  /*
   * Parse TLEN (Template Length)
   */
  GT_ISP_PARSE_SAM_ALG_SKIP_FIELD();
  /*
   * Parse SEQ (READ)
   */
  if (gt_expect_false(**text_line==STAR)) {
    GT_NEXT_CHAR(text_line);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  } else {
    register char* const seq_read = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    register const uint64_t read_length = *text_line-seq_read;
    GT_NEXT_CHAR(text_line);
//    // FIXME: Check disabled due to wrong SAM outputs produced some mappers when trimming
//    if (!gt_string_is_null(alignment->read)) {
//      register gt_string* const check_read = gt_string_new(0); // Writing to buffer's info
//      gt_string_set_nstring(check_read,seq_read,read_length);
//      if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
//        gt_dna_string_reverse_complement(check_read);
//      }
//      register const bool equals = gt_string_equals(alignment->read,check_read);
//      gt_string_delete(check_read);
//
//      if (!equals) {
//        gt_map_delete(map);
//        return GT_ISP_PE_WRONG_READ_CONTENT;
//      }
//    }
    if (gt_string_is_null(alignment->read)) {
      gt_string_set_nstring(alignment->read,seq_read,read_length);
      if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
        gt_dna_string_reverse_complement(alignment->read);
      }
    }
    if (gt_map_get_base_length(map)==0) gt_map_set_base_length(map,gt_alignment_get_read_length(alignment));
  }
  /*
   * Parse QUAL (QUALITY STRING)
   */
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
  if (gt_expect_false(**text_line==STAR)) {
    GT_NEXT_CHAR(text_line);
  } else {
    register char* const seq_qual = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    register const uint64_t read_length = *text_line-seq_qual;
    if (gt_string_is_null(alignment->qualities)) {
      gt_string_set_nstring(alignment->qualities,seq_qual,read_length);
      if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
        gt_string_reverse(alignment->qualities);
      }
    }
    if (gt_map_get_base_length(map)==0) gt_map_set_base_length(map,gt_string_get_length(alignment->qualities));
  }
  if (!gt_string_is_null(alignment->read) && !gt_string_is_null(alignment->qualities)) {
    gt_fatal_check(gt_string_get_length(alignment->read)!=gt_string_get_length(alignment->qualities),ALIGN_READ_QUAL_LENGTH);
  }
  // Build a list of alignments
  register gt_vector *maps_vector = NULL;
  if (is_mapped) {
    maps_vector = gt_vector_new(10,sizeof(gt_map*));
    gt_vector_insert(maps_vector,map,gt_map*);
  } else {
    gt_map_delete(map);
  }
  /*
   * OPTIONAL FIELDS
   */
  while (**text_line==TAB) {
    GT_NEXT_CHAR(text_line);
    if ((error_code=gt_isp_parse_sam_optional_field(text_line,alignment,maps_vector,pending,is_mapped))) {
      if (maps_vector) {
        GT_VECTOR_ITERATE(maps_vector,map_elm,map_pos,gt_map*) gt_map_delete(*map_elm);
      }
      return error_code;
    }
  }
  // Add the main map
  if (**text_line==EOL && **text_line==EOS) return GT_ISP_PE_BAD_CHARACTER;
  if (is_mapped) {
    pending->map_displacement = gt_alignment_get_num_maps(alignment);
    gt_alignment_insert_map_gt_vector(alignment,maps_vector);
    gt_vector_delete(maps_vector);
  }
  return 0;
}

GT_INLINE bool gt_isp_fetch_next_line(
    gt_buffered_input_file* const buffered_sam_input,gt_string* const expected_tag,const bool chomp_tag) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_NULL_CHECK(expected_tag);
  // Check next record/line
  gt_input_sam_parser_next_record(buffered_sam_input);
  if (gt_buffered_input_file_eob(buffered_sam_input)) return false;
  // Fetch next tag
  register gt_string* const next_tag = gt_string_new(0);
  register char* const ptext_line = gt_isp_read_tag(&(buffered_sam_input->cursor),next_tag);
  gt_fastq_tag_chomp_end_info(next_tag);
  register const bool same_tag = ptext_line!=NULL && gt_string_equals(expected_tag,next_tag);
  gt_string_delete(next_tag);
  if (same_tag) {
    buffered_sam_input->cursor = ptext_line;
    return true;
  } else {
    return false;
  }
}

GT_INLINE void gt_isp_add_mmap(
    gt_template* const template,const uint64_t start_pos_end1,const uint64_t start_pos_end2,
    const uint64_t pending_maps_end1,const uint64_t pending_maps_end2) {
  register gt_map** mmap_end1 = gt_vector_get_elm(gt_template_get_block(template,0)->maps,start_pos_end1,gt_map*);
  register gt_map** mmap_end2 = gt_vector_get_elm(gt_template_get_block(template,1)->maps,start_pos_end2,gt_map*);
  register uint64_t i;
  gt_mmap_attributes attr;
  register const uint64_t pending_maps = GT_MAX(pending_maps_end1,pending_maps_end2);
  for (i=0;i<pending_maps;++i) {
    register gt_map* const map_end1 = (pending_maps_end1>1) ? mmap_end1[i] : mmap_end1[0];
    register gt_map* const map_end2 = (pending_maps_end2>1) ? mmap_end2[i] : mmap_end2[0];
    attr.distance = gt_map_get_global_distance(map_end1)+gt_map_get_global_distance(map_end2);
    gt_template_add_mmap_va(template,&attr,map_end1,map_end2);
  }
}

GT_INLINE bool gt_isp_check_pending_record__add_mmap(
    gt_template* const template,gt_sam_pending_end* const pending,
    const uint64_t end_position,char* const seq_name,const uint64_t position,
    const uint64_t map_displacement,const uint64_t num_maps) {
  if (pending->end_position!=end_position &&
      gt_streq(pending->next_seq_name,seq_name) &&
      pending->next_position==position) { // Found!
    // BWA_Compact. MMAPs paired against MMaps need to have the same cardinality (otherwise it's unpaired)
    if (pending->num_maps!=1 && num_maps!=1 && pending->num_maps!=num_maps) return false; //FIXME return true;
    // Insert mmap(s)
    if (pending->end_position==1) {
      gt_isp_add_mmap(template,map_displacement,pending->map_displacement,num_maps,pending->num_maps);
    } else {
      gt_isp_add_mmap(template,pending->map_displacement,map_displacement,pending->num_maps,num_maps);
    }
    return true;
  }
  return false;
}

GT_INLINE void gt_isp_solve_pending_maps(
    gt_vector* pending_v,gt_sam_pending_end* pending,gt_template* const template) {
  register bool found_match = false;
  // Look into pending records
  GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
    if ((found_match=gt_isp_check_pending_record__add_mmap(
            template,pending_elm,pending->end_position,pending->map_seq_name,
            pending->map_position,pending->map_displacement,pending->num_maps))) {
      pending_elm->next_seq_name = NULL; // Mark as solved
      break;
    }
  }
  // Look into alignment maps
  if (!found_match) {
    register const uint64_t map_end = (pending->end_position+1)%2;
    register uint64_t pos = 0;
    GT_ALIGNMENT_ITERATE(gt_template_get_block(template,map_end),map) {
      if ((found_match=gt_isp_check_pending_record__add_mmap(template,pending,map_end,
          gt_string_get_string(map->seq_name),map->position,pos,1))) {
        break;
      }
      ++pos;
    }
  }
  // Queue if not found
  if (!found_match) gt_vector_insert(pending_v,*pending,gt_sam_pending_end);
}

/* SAM general */
GT_INLINE gt_status gt_input_sam_parser_parse_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template) { // FIXME: On error read the rest of the tag SAM records
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  register char** text_line = &(buffered_sam_input->cursor);
  // Read initial TAG (QNAME := Query template)
  if (!((*text_line)=gt_isp_read_tag(text_line,template->tag))) return GT_ISP_PE_PREMATURE_EOL;
  // Read all maps related to this TAG
  register gt_status error_code;
  gt_vector* pending_v = gt_vector_new(GT_ISP_NUM_INITIAL_MAPS,sizeof(gt_sam_pending_end));
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending;
    uint64_t alignment_flag;
    if (gt_expect_false(error_code=gt_isp_parse_sam_alignment(
          text_line,template,NULL,&alignment_flag,&pending,false))) {
      gt_vector_delete(pending_v);
      return error_code;
    }
    // Solve pending ends
    if (pending.next_seq_name!=NULL) gt_isp_solve_pending_maps(pending_v,&pending,template);
  } while (gt_isp_fetch_next_line(buffered_sam_input,template->tag,true));
  // Check for unsolved pending maps
  error_code = 0;
  GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
    if (pending_elm->next_seq_name!=NULL) {
      error_code = GT_ISP_PE_UNSOLVED_PENDING_MAPS;
      break;
    }
  }
  gt_vector_delete(pending_v);
  return error_code;
}
/* SOAP2-SAM */
GT_INLINE gt_status gt_input_sam_parser_parse_soap_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  register char** text_line = &(buffered_sam_input->cursor);
  register gt_status error_code;
  // Read initial TAG (QNAME := Query template)
  if (!((*text_line)=gt_isp_read_tag(text_line,template->tag))) return GT_ISP_PE_PREMATURE_EOL; // TODO: Alignment tag builds
  gt_fastq_tag_chomp_end_info(template->tag);
  // Read all maps related to this TAG
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending;
    uint64_t alignment_flag;
    if (gt_expect_false(error_code=gt_isp_parse_sam_alignment(
          text_line,template,NULL,&alignment_flag,&pending,false))) {
      return error_code;
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,template->tag,true));
  // SOAP2 paired maps convention. Add maps
  register gt_alignment* const alignment_end0 = gt_template_get_block(template,0);
  register gt_alignment* const alignment_end1 = gt_template_get_block(template,1);
  if (gt_alignment_get_num_maps(alignment_end0) !=
      gt_alignment_get_num_maps(alignment_end1)) return GT_ISP_PE_UNSOLVED_PENDING_MAPS;
  register uint64_t pos_end_it=0;
  GT_ALIGNMENT_ITERATE(alignment_end0,map_end1) {
    gt_map* map_end2 = gt_alignment_get_map(alignment_end1,pos_end_it);
    gt_mmap_attributes attr;
    attr.distance = gt_map_get_global_distance(map_end1)+gt_map_get_global_distance(map_end2);
    gt_template_add_mmap_va(template,&attr,map_end1,map_end2);
    ++pos_end_it;
  }
  return 0;
}
/* SE-SAM */
GT_INLINE gt_status gt_input_sam_parser_parse_alignment(
    gt_buffered_input_file* const buffered_sam_input,gt_alignment* alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  register char** text_line = &(buffered_sam_input->cursor);
  // Read initial TAG (QNAME := Query template)
  if (!((*text_line)=gt_isp_read_tag(text_line,alignment->tag))) return GT_ISP_PE_PREMATURE_EOL;
  // Read all maps related to this TAG
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending;
    uint64_t alignment_flag;
    register gt_status error_code;
    if (gt_expect_false((error_code=gt_isp_parse_sam_alignment(
        text_line,NULL,alignment,&alignment_flag,&pending,true))!=0)) {
      return error_code;
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,alignment->tag,false));
  return 0;
}


/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template_(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template,const bool soap_sam) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    if ((error_code=gt_input_sam_parser_reload_buffer(buffered_sam_input))!=GT_ISP_OK) return error_code;
  }
  // Check file format
  register gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,input_file->file_name,buffered_sam_input->current_line_num,0ul);
    return GT_ISP_FAIL;
  }
  // Prepare the template
  register char* const line_start = buffered_sam_input->cursor;
  register const uint64_t line_num = buffered_sam_input->current_line_num;
  gt_template_clear(template,true);
  template->template_id = line_num;
  // Parse template
  if (gt_expect_false(soap_sam)) {
    error_code=gt_input_sam_parser_parse_soap_template(buffered_sam_input,template);
  } else {
    error_code=gt_input_sam_parser_parse_template(buffered_sam_input,template);
  }
  if (error_code) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  //gt_input_sam_parser_next_record(buffered_sam_input);
  return GT_ISP_OK;
}
GT_INLINE gt_status gt_input_sam_parser_get_template(gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  return gt_input_sam_parser_get_template_(buffered_sam_input,template,false);
}
GT_INLINE gt_status gt_input_sam_parser_get_soap_template(gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  return gt_input_sam_parser_get_template_(buffered_sam_input,template,true);
}
GT_INLINE gt_status gt_input_sam_parser_get_alignment(
    gt_buffered_input_file* const buffered_sam_input,gt_alignment* const alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    if ((error_code=gt_input_sam_parser_reload_buffer(buffered_sam_input))!=GT_ISP_OK) return error_code;
  }
  // Check file format
  register gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,input_file->file_name,buffered_sam_input->current_line_num,0ul);
    return GT_ISP_FAIL;
  }
  // Allocate memory for the alignment
  register char* const line_start = buffered_sam_input->cursor;
  register const uint64_t line_num = buffered_sam_input->current_line_num;
  gt_alignment_clear(alignment);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_input_sam_parser_parse_alignment(buffered_sam_input,alignment))) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  //gt_input_sam_parser_next_record(buffered_sam_input);
  return GT_ISP_OK;
}
