/*
 * PROJECT: GEM-Tools library
 * FILE: gt_generic_printer.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_GENERIC_PRINTER_H_
#define GT_GENERIC_PRINTER_H_

#include "gt_commons.h"
#include "gt_output_buffer.h"

typedef enum { GT_FILE_PRINTER, GT_STRING_PRINTER, GT_BUFFER_PRINTER } gt_generic_printer_t;
typedef struct {
  gt_generic_printer_t printer_type;
  union {
    FILE* file;
    gt_string *string;
    gt_output_buffer* output_buffer;
  }; // Printer Object
} gt_generic_printer;

/*
 * Checkers
 */
#define GT_GENERIC_PRINTER_CHECK(generic_printer) \
  GT_NULL_CHECK(generic_printer); \
  switch ((generic_printer)->printer_type) { \
    case GT_FILE_PRINTER: GT_NULL_CHECK((generic_printer)->file); break; \
    case GT_STRING_PRINTER: GT_STRING_CHECK((generic_printer)->string); break; \
    case GT_BUFFER_PRINTER: GT_OUTPUT_BUFFER_CHECK((generic_printer)->output_buffer); break; \
    default: gt_fatal_check(true,SELECTION_NOT_VALID); break; \
  }


/*
 * Generic printer
 */
GT_INLINE void gt_generic_new_file_printer(gt_generic_printer* const generic_printer,FILE* const file);
GT_INLINE void gt_generic_new_string_printer(gt_generic_printer* const generic_printer,gt_string* const string);
GT_INLINE void gt_generic_new_buffer_printer(gt_generic_printer* const generic_printer,gt_output_buffer* const output_buffer);

GT_INLINE void gt_vgprintf(gt_generic_printer* const generic_printer,const char *template,va_list v_args);
GT_INLINE void gt_gprintf(gt_generic_printer* const generic_printer,const char *template,...);

#endif /* GT_OUTPUT_PRINTER_H_ */
