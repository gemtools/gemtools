/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_tag_parser.c
 * DATE: 01/06/2012
 * AUTHOR(S): Thasso Griebel <thasso.griebel@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_parser.h"

// count number of ':' in field + 1
uint64_t count_casava_fields(char* const field, uint64_t length){
  register uint64_t i = 0;
  register uint64_t count = 1;
  for(;i<length;i++){
    if(field[i] == ':') count++;
  }
  return count;
}

GT_INLINE gt_status gt_input_parse_tag(char** const text_line,gt_string* const tag,gt_shash* const attributes){
  // Delimit the tag
  register char* const tag_begin = *text_line;
  // read until first SPACE or TAB
  GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
  // tag length
  register uint64_t tag_length = *text_line-tag_begin;
  gt_string_set_nstring(tag,tag_begin,tag_length);

  int64_t pair = GT_PAIR_SE;
  if (tag_length>2 && *gt_string_char_at(tag,tag_length-2)==SLASH) {
    register const char tag_end = *gt_string_char_at(tag,tag_length-1);
    if (tag_end=='1') {
      gt_string_set_length(tag,tag_length-2);
      pair = GT_PAIR_PE_1;
    } else if (tag_end=='2' || tag_end == '3') {
      gt_string_set_length(tag,tag_length-2);
      pair = GT_PAIR_PE_2;
    }
    gt_string_append_eos(tag);
  }

  if (gt_expect_true(**text_line==SPACE)) {
    GT_NEXT_CHAR(text_line);
    register char* attr_begin = *text_line;
    // additional attributes
    // check for casava 1.8 attributes
    GT_READ_UNTIL(text_line,**text_line==TAB || **text_line==SPACE);
    register uint64_t attr_length = *text_line-attr_begin;
    if (count_casava_fields(attr_begin,attr_length) == 4) {
      if (gt_expect_true(attr_begin[0]=='1')) {
        pair = GT_PAIR_PE_1;
      } else if (gt_expect_true(attr_begin[0]=='2' || attr_begin[0]=='3')) {
        pair = GT_PAIR_PE_2;
      }
      gt_string* casava_string = gt_string_new(attr_length);
      gt_string_set_nstring(casava_string,attr_begin,attr_length);
      // continue to next only if its not a tab
      if(!gt_expect_false((**text_line)==TAB) && !GT_IS_EOL(text_line)) {
        GT_NEXT_CHAR(text_line);
      }
      attr_begin = *text_line;
      gt_attribute_set_string(attributes,GT_TAG_CASAVA,casava_string);
    }
    GT_READ_UNTIL(text_line, **text_line==TAB);
    attr_length = *text_line-attr_begin;
    gt_string* extra_string = gt_string_new(attr_length);
    gt_string_set_nstring(extra_string,attr_begin,attr_length);

    // continue to next only if its not a tab
    if (!gt_expect_false((**text_line)==TAB) && ! GT_IS_EOL(text_line)) {
      GT_NEXT_CHAR(text_line);
    }
    attr_begin += attr_length;
    gt_attribute_set_string(attributes,GT_TAG_EXTRA,extra_string);
  }
  gt_attribute_set(attributes,GT_TAG_PAIR,&pair,int64_t);

  if (!GT_IS_EOL(text_line)) {
    GT_NEXT_CHAR(text_line);
  }
  return GT_STATUS_OK;
}

