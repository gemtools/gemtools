/*
 * PROJECT: GEM-Tools library
 * FILE: gt_generic_printer.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_generic_printer.h"

#define GT_GEN_PRINTER_INITIAL_STRING_SIZE GT_BUFFER_SIZE_1K

/*
 * Generic printer
 */
GT_INLINE void gt_generic_new_file_printer(gt_generic_printer* const generic_printer,FILE* const file) {
  GT_NULL_CHECK(generic_printer);
  GT_NULL_CHECK(file);
  generic_printer->printer_type = GT_FILE_PRINTER;
  generic_printer->file = file;
}
GT_INLINE void gt_generic_new_string_printer(gt_generic_printer* const generic_printer,gt_string* const string) {
  GT_NULL_CHECK(generic_printer);
  GT_STRING_CHECK(string);
  generic_printer->printer_type = GT_STRING_PRINTER;
  generic_printer->string = string;
}
GT_INLINE void gt_generic_new_buffer_printer(gt_generic_printer* const generic_printer,gt_output_buffer* const output_buffer) {
  GT_NULL_CHECK(generic_printer);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  generic_printer->printer_type = GT_BUFFER_PRINTER;
  generic_printer->output_buffer = output_buffer;
}
GT_INLINE void gt_vgprintf(gt_generic_printer* const generic_printer,const char *template,va_list v_args) {
  GT_GENERIC_PRINTER_CHECK(generic_printer);
  GT_NULL_CHECK(template);
  switch (generic_printer->printer_type) {
    case GT_FILE_PRINTER:
      GT_NULL_CHECK(generic_printer->file);
      gt_cond_fatal_error(vfprintf(generic_printer->file,template,v_args)<0,FPRINTF);
      break;
    case GT_STRING_PRINTER:
      GT_STRING_CHECK(generic_printer->string);
      if (!gt_string_is_static(generic_printer->string)) { // Allocate memory
        register const uint64_t mem_required = gt_calculate_memory_required_v(template,v_args);
        gt_string_resize(generic_printer->string,mem_required);
      }
      register const int64_t mem_used = vsprintf(gt_string_get_string(generic_printer->string),template,v_args);
      gt_cond_fatal_error(mem_used<0,SPRINTF);
      gt_string_set_length(generic_printer->string,mem_used);
      break;
    case GT_BUFFER_PRINTER:
      GT_OUTPUT_BUFFER_CHECK(generic_printer->output_buffer);
      gt_cond_fatal_error(gt_vbprintf(&generic_printer->output_buffer,template,v_args)<0,BPRINTF);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      break;
  }
}
GT_INLINE void gt_gprintf(gt_generic_printer* const generic_printer,const char *template,...) {
  GT_GENERIC_PRINTER_CHECK(generic_printer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  gt_vgprintf(generic_printer,template,v_args);
}
