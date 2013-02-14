/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_parser.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Thasso Griebel <thasso.griebel@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_PARSER_H_
#define GT_INPUT_PARSER_H_

#include "gt_commons.h"
#include "gt_data_attributes.h"

#define GT_TAG_PAIR "pair"
#define GT_TAG_CASAVA "casava"
#define GT_TAG_EXTRA "extra"

#define GT_PAIR_SE 0
#define GT_PAIR_BOTH 3
#define GT_PAIR_PE_1 1
#define GT_PAIR_PE_2 2

/*
 * Generic tag parser
 */
GT_INLINE gt_status gt_input_parse_tag(char** const text_line,gt_string* const tag,gt_shash* const attributes);

/*
 * Internal Building Blocks for parsing
 */
#define GT_NEXT_CHAR(text_line) ++(*text_line)
#define GT_IS_EOL(text_line) gt_expect_false((**text_line)==EOL || (**text_line)==EOS)
#define GT_READ_UNTIL(text_line,test) \
  while (gt_expect_true(!(test) && !GT_IS_EOL(text_line))) { \
    GT_NEXT_CHAR(text_line); \
  }
#define GT_PARSE_NUMBER(text_line,number) \
  number = 0; \
  while (gt_expect_true(gt_is_number(**text_line))) { \
    number = (number*10) + gt_get_cipher(**text_line); \
    GT_NEXT_CHAR(text_line); \
  }
#define GT_PARSE_SIGNED_NUMBER_BLOCK(text_line,number) { \
  register bool is_negative; \
  switch ((**text_line)) { \
    case PLUS: is_negative = false; break; \
    case MINUS: is_negative = true; break; \
    default: is_negative = false; break; \
  }
#define GT_PARSE_SIGNED_NUMBER_END_BLOCK(number) \
  if (is_negative) number = -number; \
}
#define GT_SKIP_LINE(text_line) \
  while (!GT_IS_EOL(text_line)) { \
    ++(*text_line); \
  }

#endif /* GT_INPUT_PARSER_H_ */
