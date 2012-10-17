/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_parser.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_PARSER_H_
#define GT_INPUT_PARSER_H_

#include "gt_commons.h"


/*
 * Internal Building Blocks for parsing
 */
#define GT_NEXT_CHAR(text_line) ++(*text_line)
#define GT_READ_UNTIL(text_line,test) \
  while (gt_expect_true(!(test) && (**text_line)!=EOL)) { \
    GT_NEXT_CHAR(text_line); \
  }
#define GT_IS_EOL(text_line) gt_expect_false((**text_line)==EOL)
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
  while (__builtin_expect((**text_line)!=EOL,1)) { \
    ++(*text_line); \
  }



#endif /* GT_INPUT_PARSER_H_ */
