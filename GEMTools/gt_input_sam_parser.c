/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.c
 * DATE: 17/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_sam_parser.h"
#include "gt_input_parser.h"

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
  GT_INPUT_FILE_CHECK__FILL_BUFFER(input_file,NULL);
  if (gt_input_file_detect_sam_format((char*)input_file->file_buffer,input_file->buffer_size,show_errors)) {
    while (!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==GT_SAM_HEADER_BEGIN) {

      // TODO => sam_headers
      while (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)!=EOL)) {
        GT_INPUT_FILE_NEXT_CHAR(input_file,NULL);
        GT_INPUT_FILE_CHECK__FILL_BUFFER(input_file,NULL);
      }
      if (!input_file->eof) {
        GT_INPUT_FILE_NEXT_CHAR(input_file,NULL); // Check DOS EOF
        if (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==DOS_EOL)) {
          GT_INPUT_FILE_NEXT_CHAR(input_file,NULL);
        }
      }
      // TODO => sam_headers

    }
    // Skip header from input-file buffer
    input_file->buffer_begin = input_file->buffer_pos;
    return true;
  } else {
    return false;
  }
}
GT_INLINE gt_status gt_input_sam_parser_check_sam_file_format(gt_buffered_input_file* const buffered_sam_input) {
  register gt_input_file* const input_file = buffered_sam_input->input_file;
  if (gt_expect_false(input_file->file_format==UNKNOWN)) { // Unknown
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
 * Parsing Building Blocks
 */
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_sam_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  if (!gt_buffered_input_file_eob(buffered_sam_input)) {
    GT_INPUT_FILE_SKIP_LINE(buffered_sam_input);
  }
}
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
    default:
      gt_error(PARSE_SAM,buffered_sam_input->input_file->file_name,line_num);
      break;
  }
}
GT_INLINE gt_status gt_isp_read_tag(
    char** const text_line,char** tag,uint64_t* tag_length,uint64_t* end_position) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(tag_length);
  // Read tag
  *tag = *text_line;
  GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
  if (GT_IS_EOL(text_line)) return GT_ISP_PE_PREMATURE_EOL;
  register char* last_parsed_tag_char = *text_line-1;
  *tag_length = *text_line-*tag;
  if (**text_line==SPACE) {
    GT_READ_UNTIL(text_line,**text_line==TAB);
    if (GT_IS_EOL(text_line)) return GT_ISP_PE_PREMATURE_EOL;
  }
  GT_NEXT_CHAR(text_line);
  // Parse the end information {/1,/2,...}
  if (*tag_length>2 && *(last_parsed_tag_char-1)==SLASH) {
    // *(last_parsed_tag_char-1)=EOS;
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
GT_INLINE gt_status gt_isp_parse_sam_cigar(char** const text_line,gt_map* const map) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_MAP_CHECK(map);
  register uint64_t position = 0;
  register uint64_t map_length = 0;
  if (**text_line==STAR) GT_NEXT_CHAR(text_line);
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
      case 'I': // Insertion to the reference
      case 'P': // Padding. Nothing implemented
      case 'S': // Soft clipping. Nothing implemented
        misms.misms_type = INS;
        misms.position = position;
        misms.size = length;
        position += length;
        map_length -= length;
        gt_map_add_misms(map,&misms);
        break;
      case 'D':  // Deletion from the reference
        misms.misms_type = DEL;
        misms.position = position;
        misms.size = length;
        map_length += length;
        gt_map_add_misms(map,&misms);
        break;
      case 'H': // Hard clipping. Nothing implemented
        break;
      case 'N':
        gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      default:
        return GT_ISP_PE_BAD_CHARACTER;
        break;
    }
  }
  return 0;
}

#define GT_ISP_PARSE_SAM_EXTRA_ALG_CHECK_PREMATURE_EOL() \
  if (GT_IS_EOL(text_line)) { \
    gt_map_delete(map); \
    GT_VECTOR_ITERATE(maps_vector,map_it,map_pos,gt_map*) { \
      gt_map_delete(*map_it); \
    } \
    gt_vector_delete(maps_vector); \
    return GT_ISP_PE_PREMATURE_EOL; \
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
#define GT_ISP_IF_OPT_FIELD(text_line,char1,char2) { \
  if ((*text_line)[0]==char1 && (*text_line)[1]==char2) {
#define GT_ISP_END_OPT_FIELD \
    keep_parsing = (**text_line!=EOL); \
    continue; \
  }}
GT_INLINE gt_status gt_isp_parse_sam_alignment(
    char** const text_line,gt_alignment** const alignment_blocks,
    uint64_t* const alignment_flag,gt_sam_pending_end* const pending) {
  GT_NULL_CHECK(text_line); GT_NULL_CHECK(*text_line);
  GT_NULL_CHECK(alignment_flag);
  GT_NULL_CHECK(alignment_blocks);
  register gt_status error_code;
  register bool is_mapped = true, is_single_segment;
  register uint64_t read_length;
  register gt_map* map = gt_map_new();
  /*
   * Parse FLAG
   */
  if (!gt_is_number(**text_line)) return GT_ISP_PE_EXPECTED_NUMBER;
  GT_PARSE_NUMBER(text_line,*alignment_flag);
  GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  // Process flags
  is_mapped = !(*alignment_flag&GT_SAM_FLAG_UNMAPPED);
  is_single_segment = !(*alignment_flag&GT_SAM_FLAG_MULTIPLE_SEGMENTS);
  pending->end_position = (is_single_segment) ? 0 : ((*alignment_flag&GT_SAM_FLAG_FIRST_SEGMENT)?0:1);
  if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT)  {
    gt_map_set_direction(map,REVERSE);
  } else {
    gt_map_set_direction(map,FORWARD);
  }
  /*
   * Parse RNAME (Sequence-name/Chromosome)
   */
  if (gt_expect_false(**text_line==STAR)) {
    is_mapped=false; /* Unmapped */
    GT_NEXT_CHAR(text_line);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
  } else {
    gt_map_set_seq_name(map,*text_line);
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    **text_line=EOS; GT_NEXT_CHAR(text_line);
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
  if (**text_line==STAR) {
    pending->next_seq_name = NULL;
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD();
    GT_ISP_PARSE_SAM_ALG_SKIP_FIELD();
  } else {
    if (**text_line==EQUAL) {
      pending->next_seq_name = gt_map_get_seq_name(map);
      GT_NEXT_CHAR(text_line);
      if (**text_line!=TAB) return GT_ISP_PE_BAD_CHARACTER;
      GT_NEXT_CHAR(text_line);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    } else {
      pending->next_seq_name = *text_line;
      GT_READ_UNTIL(text_line,**text_line==TAB);
      GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
      **text_line=EOS; GT_NEXT_CHAR(text_line);
    }
    // Parse PNEXT (Position of the next segment)
    GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(pending->next_position);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL__NEXT();
    if (pending->next_position==0) {
      pending->next_seq_name = NULL;
    } else {
      pending->num_maps = 1;
      pending->map_seq_name = map->seq_name;
      pending->map_position = map->position;
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
    register char* seq_read = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    GT_ISP_PARSE_SAM_ALG_CHECK_PREMATURE_EOL();
    if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
      gt_reverse_complent(seq_read,*text_line-seq_read);
    }
    read_length = *text_line-seq_read;
    **text_line=EOS; GT_NEXT_CHAR(text_line);
    GT_NULL_CHECK(alignment_blocks[pending->end_position]);
    if (alignment_blocks[pending->end_position]->read!=NULL) {
      if (strcmp(alignment_blocks[pending->end_position]->read,seq_read)!=0) {
        gt_map_delete(map); return GT_ISP_PE_WRONG_READ_CONTENT;
      }
    } else {
      alignment_blocks[pending->end_position]->read_length = read_length;
      alignment_blocks[pending->end_position]->read = gt_string_cpy(seq_read,read_length);
      gt_map_set_base_length(map,read_length);
    }
  }
  /*
   * Parse QUAL (QUALITY STRING)
   */
  register bool keep_parsing;
  if (gt_expect_false(**text_line==STAR)) {
    GT_NEXT_CHAR(text_line);
    keep_parsing = (**text_line!=EOL);
  } else {
    register char* qual_read = *text_line;
    GT_READ_UNTIL(text_line,**text_line==TAB);
    if (*alignment_flag&GT_SAM_FLAG_REVERSE_COMPLEMENT) {
      gt_reverse(qual_read,*text_line-qual_read);
    }
    if (alignment_blocks[pending->end_position]->qualities==NULL) {
      alignment_blocks[pending->end_position]->qualities =
          gt_string_cpy(qual_read,*text_line-qual_read);
    }
    keep_parsing = (**text_line!=EOL);
    **text_line=EOS;
  }
  // Build a list of alignments
  register gt_map* const primary_map = map;
  register gt_vector *maps_vector;
  if (is_mapped) {
    maps_vector = gt_vector_new(10,sizeof(gt_map*));
    gt_vector_insert(maps_vector,map,gt_map*);
  }
  /*
   * OPTIONAL FIELDS
   */
  if (keep_parsing) {
    GT_NEXT_CHAR(text_line);
    while (keep_parsing) {
      // XA:Z:chr17,-34553512,125M,0;chr17,-34655077,125M,0;
      GT_ISP_IF_OPT_FIELD(text_line,'X','A') {
        if (!is_mapped) return GT_ISP_PE_SAM_UNMAPPED_XA;
        *text_line+=5;
        while (**text_line!=TAB && **text_line!=EOL) { // Read new attached maps
          map = gt_map_new();
          gt_map_set_base_length(map,read_length);
          // Sequence-name/Chromosome
          gt_map_set_seq_name(map,*text_line);
          GT_READ_UNTIL(text_line,**text_line==COMA);
          GT_ISP_PARSE_SAM_EXTRA_ALG_CHECK_PREMATURE_EOL();
          **text_line=EOS; GT_NEXT_CHAR(text_line);
          // Position
          if (**text_line==MINUS) {
            map->direction = REVERSE; GT_NEXT_CHAR(text_line);
          } else if (**text_line==PLUS) {
            map->direction = FORWARD; GT_NEXT_CHAR(text_line);
          } else {
            map->direction = FORWARD;
          }
          GT_ISP_PARSE_SAM_ALG_PARSE_NUMBER(map->position);
          // CIGAR
          GT_READ_UNTIL(text_line,**text_line==COMA);
          GT_ISP_PARSE_SAM_EXTRA_ALG_CHECK_PREMATURE_EOL();
          GT_NEXT_CHAR(text_line);
          // TODO Recalculate CIGAR (including mismatches)
          // Edit distance
          GT_READ_UNTIL(text_line,**text_line==COMA);
          GT_ISP_PARSE_SAM_EXTRA_ALG_CHECK_PREMATURE_EOL();
          GT_NEXT_CHAR(text_line);
          // Add it to the list
          gt_vector_insert(maps_vector,map,gt_map*);
          // Consider pending relation
          ++pending->num_maps;
        }
      } GT_ISP_END_OPT_FIELD;
      // TODO: MD field
      // Unknown field (skip it)
      GT_READ_UNTIL(text_line,**text_line==TAB);
      if ((keep_parsing=(**text_line==TAB))) GT_NEXT_CHAR(text_line);
    }
  }
  // Add the main map
  if (is_mapped) {
    pending->map_displacement = gt_alignment_get_num_maps(alignment_blocks[pending->end_position]);
    gt_alignment_insert_map_gt_vector(alignment_blocks[pending->end_position],maps_vector);
    gt_vector_delete(maps_vector);
  }
  return 0;
}

GT_INLINE bool gt_isp_fetch_next_line(
    gt_buffered_input_file* const buffered_sam_input,
    char* const expected_tag,const uint64_t expected_tag_length) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_NULL_CHECK(expected_tag);
  uint64_t end_position;
  // Check next record/line
  gt_input_sam_parser_next_record(buffered_sam_input);
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    if (gt_buffered_input_file_get_block(buffered_sam_input,GT_NUM_LINES_10K,true) == 0) {
      return false;
    }
  }
  // Fetch next tag
  register char* const initial_cursor = buffered_sam_input->cursor;
  char* next_tag;
  uint64_t next_tag_length;
  if (gt_isp_read_tag(&(buffered_sam_input->cursor),&next_tag,&next_tag_length,&end_position) ||
      next_tag_length!=expected_tag_length || (strncmp(next_tag,expected_tag,next_tag_length)!=0)) {
    buffered_sam_input->cursor = initial_cursor; // Rewind
    return false;
  }
  return true;
}

GT_INLINE gt_status gt_input_sam_parser_parse_template(
    gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  register char** text_line = &(buffered_sam_input->cursor);
  register gt_status error_code;
  // Read initial TAG (QNAME := Query template)
  uint64_t end_position;
  if ((error_code=gt_isp_read_tag(text_line,&template->tag,&template->tag_length,&end_position))) {
    return error_code;
  }
  template->tag = gt_string_cpy(template->tag,template->tag_length);
  // Read all maps related to this TAG
  gt_vector* pending_v = gt_vector_new(10,sizeof(gt_sam_pending_end));
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending;
    uint64_t alignment_flag;
    if (gt_expect_false((error_code=gt_isp_parse_sam_alignment(
        text_line,gt_vector_get_mem(template->blocks,gt_alignment*),&alignment_flag,&pending))!=0)) {
      gt_vector_delete(pending_v); return error_code;
    }
    // Keep track of pending ends
    if (pending.next_seq_name!=NULL) {
      register bool found_match = false;
      GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
        if (pending_elm->next_seq_name!=NULL &&
            gt_string_eq(pending_elm->next_seq_name,pending.next_seq_name) &&
            pending_elm->next_position==pending.map_position &&
            pending_elm->map_position==pending.next_position) { // Found!
          found_match = true;
          // Check full consistency
          if (pending_elm->num_maps != pending.num_maps) {
            gt_vector_delete(pending_v); return GT_ISP_PE_WRONG_NUM_XA;
          }
          if (pending_elm->end_position!=pending.end_position) {
            gt_vector_delete(pending_v); return GT_ISP_PE_WRONG_END_POS;
          }
          // Insert mmap(s)
          register const uint64_t start_pos_end1 = (pending.end_position==0)?pending.map_position:pending_elm->num_maps;
          register const uint64_t start_pos_end2 = (pending.end_position==1)?pending.map_position:pending_elm->num_maps;
          register gt_map** map_end1 = gt_vector_get_elm(gt_template_get_block(template,0)->maps,start_pos_end1,gt_map*);
          register gt_map** map_end2 = gt_vector_get_elm(gt_template_get_block(template,1)->maps,start_pos_end2,gt_map*);
          register uint64_t i;
          for (i=0;i<pending.num_maps;++i) {
            gt_mmap_attributes attr;
            attr.distance = gt_map_get_global_distance(map_end1[i])+gt_map_get_global_distance(map_end2[i]);
            gt_template_add_mmap_va(template,&attr,map_end1[i],map_end2[i]);
          }
          // Mark as solved
          pending_elm->next_seq_name = NULL;
        }
      }
      // Queue if not found
      if (!found_match) {

      }
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,template->tag,template->tag_length));
  // Check for unsolved pending maps
  GT_VECTOR_ITERATE(pending_v,pending_elm,pending_counter,gt_sam_pending_end) {
    if (pending_elm->next_seq_name!=NULL) {
      gt_vector_delete(pending_v); return GT_ISP_PE_UNSOLVED_PENDING_MAPS;
    }
  }
  gt_vector_delete(pending_v);
  return 0;
}

GT_INLINE gt_status gt_input_sam_parser_parse_alignment(
    gt_buffered_input_file* const buffered_sam_input,gt_alignment* alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  register char** text_line = &(buffered_sam_input->cursor);
  register gt_status error_code;
  uint64_t alignment_flag, end_position;
  // Read initial TAG (QNAME := Query template)
  if ((error_code=gt_isp_read_tag(text_line,&alignment->tag,&alignment->tag_length,&end_position))) {
    return error_code;
  }
  alignment->tag = gt_string_cpy(alignment->tag,alignment->tag_length);
  // Read all maps related to this TAG
  do {
    // Parse SAM Alignment
    gt_sam_pending_end pending;
    if (gt_expect_false((error_code=gt_isp_parse_sam_alignment(
        text_line,&alignment,&alignment_flag,&pending))!=0)) {
      return error_code;
    }
  } while (gt_isp_fetch_next_line(buffered_sam_input,alignment->tag,alignment->tag_length));
  return 0;
}


/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template(gt_buffered_input_file* const buffered_sam_input,gt_template* const template) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_TEMPLATE_CHECK(template);
  register const char* line_start = buffered_sam_input->cursor;
  register gt_status error_code;
  // Check file format
  register gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,input_file->file_name,1ul);
    return GT_ISP_FAIL;
  }
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    register const uint64_t read_lines =
        gt_buffered_input_file_get_block(buffered_sam_input,GT_NUM_LINES_10K,true);
    if (gt_expect_false(read_lines==0)) return GT_ISP_EOF;
  }
  // Prepare the template
  register const uint64_t line_num = buffered_sam_input->current_line_num;
  GT_ALIGNMENT_ITERATE(template,alignment) {
    gt_cfree(alignment->tag);
    gt_cfree(alignment->read);
    gt_cfree(alignment->qualities);
  }
  gt_template_clear(template);
  template->template_id = line_num;
  // Parse template
  if ((error_code=gt_input_sam_parser_parse_template(buffered_sam_input,template))) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  return GT_ISP_OK;
}
GT_INLINE gt_status gt_input_sam_parser_get_alignment(gt_buffered_input_file* const buffered_sam_input,gt_alignment* const alignment) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_ALIGNMENT_CHECK(alignment);
  register const char* line_start = buffered_sam_input->cursor;
  register gt_status error_code;
  // Check file format
  register gt_input_file* input_file = buffered_sam_input->input_file;
  if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
    gt_error(PARSE_SAM_BAD_FILE_FORMAT,input_file->file_name,1ul);
    return GT_ISP_FAIL;
  }
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_sam_input)) {
    register const uint64_t read_lines =
        gt_buffered_input_file_get_block(buffered_sam_input,GT_NUM_LINES_10K,true);
    if (gt_expect_false(read_lines==0)) return GT_ISP_EOF;
  }
  // Allocate memory for the alignment
  register const uint64_t line_num = buffered_sam_input->current_line_num;
  gt_alignment_clear(alignment); // FIXME: Persistent & Coherent model of Alg/Tem
  gt_cfree(alignment->tag);
  gt_cfree(alignment->read);
  gt_cfree(alignment->qualities);
  alignment->alignment_id = line_num;
  // Parse alignment
  if ((error_code=gt_input_sam_parser_parse_alignment(buffered_sam_input,alignment))) {
    gt_input_sam_parser_prompt_error(buffered_sam_input,line_num,
        buffered_sam_input->cursor-line_start,error_code);
    gt_input_sam_parser_next_record(buffered_sam_input);
    return GT_ISP_FAIL;
  }
  return GT_ISP_OK;
}
