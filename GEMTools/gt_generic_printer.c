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
GT_INLINE void gt_generic_new_file_printer(gt_generic_printer* const generic_printer,FILE* file) {
  GT_NULL_CHECK(generic_printer); GT_NULL_CHECK(file);
  generic_printer->printer_type = GT_FILE_PRINTER;
  generic_printer->file = file;
}
GT_INLINE void gt_generic_new_string_printer(gt_generic_printer* const generic_printer,char** string) {
  GT_NULL_CHECK(generic_printer); GT_NULL_CHECK(string);
  generic_printer->printer_type = GT_STRING_PRINTER;
  if (*string==NULL) { // Dynamic Allocate memory for the string
    generic_printer->buffer = gt_vector_new(GT_GEN_PRINTER_INITIAL_STRING_SIZE,sizeof(char));
    generic_printer->src_string = string;
    *(generic_printer->src_string) = gt_vector_get_mem(generic_printer->buffer,char);
    generic_printer->is_dynamic = true;
  } else { // Use memory segment provided by the user (risky ... risky ...)
    generic_printer->buffer = malloc(sizeof(gt_vector));
    gt_cond_fatal_error(!generic_printer->buffer,MEM_HANDLER);
    generic_printer->buffer->memory=*string;
    generic_printer->buffer->used=0;
    generic_printer->is_dynamic = false;
  }
}
GT_INLINE void gt_generic_new_buffer_printer(gt_generic_printer* const generic_printer,gt_output_buffer* output_buffer) {
  GT_NULL_CHECK(generic_printer); GT_NULL_CHECK(output_buffer);
  generic_printer->printer_type = GT_BUFFER_PRINTER;
  generic_printer->output_buffer = output_buffer;
}
GT_INLINE void gt_generic_delete_printer(gt_generic_printer* const generic_printer) {
  GT_NULL_CHECK(generic_printer);
  switch (generic_printer->printer_type) {
    case GT_FILE_PRINTER:
      break;
    case GT_STRING_PRINTER:
      if (generic_printer->is_dynamic) {
        gt_vector_delete(generic_printer->buffer);
      } else {
        free(generic_printer->buffer);
      }
      break;
    case GT_BUFFER_PRINTER:
      break;
    default:
      gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      break;
  }
}
GT_INLINE void gt_vgprintf(gt_generic_printer* const generic_printer,const char *template,va_list v_args) {
  GT_NULL_CHECK(generic_printer); GT_NULL_CHECK(template);
  switch (generic_printer->printer_type) {
    case GT_FILE_PRINTER:
      GT_NULL_CHECK(generic_printer->file);
      gt_fatal_check(vfprintf(generic_printer->file,template,v_args)<0,FPRINTF);
      break;
    case GT_STRING_PRINTER:
      GT_NULL_CHECK(generic_printer->buffer);
      if (generic_printer->is_dynamic) { // Allocate memory
        register const uint64_t mem_required = gt_calculate_memory_required_v(template,v_args);
        gt_vector_reserve_additional(generic_printer->buffer,mem_required);
        *(generic_printer->src_string) = gt_vector_get_mem(generic_printer->buffer,char);
      }
      register const int64_t mem_used = vsprintf(
          gt_vector_get_free_elm(generic_printer->buffer,char),template,v_args);
      gt_fatal_check(mem_used<0,SPRINTF);
      gt_vector_add_used(generic_printer->buffer,mem_used);
      break;
    case GT_BUFFER_PRINTER:
      GT_NULL_CHECK(generic_printer->output_buffer);
      gt_fatal_check(gt_vbprintf(generic_printer->output_buffer,template,v_args)<0,BPRINTF);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      break;
  }
}
GT_INLINE void gt_gprintf(gt_generic_printer* const generic_printer,const char *template,...) {
  GT_NULL_CHECK(generic_printer); GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  gt_vgprintf(generic_printer,template,v_args);
}
