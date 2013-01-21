/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_error.h"

FILE* gt_error_stream=NULL;

#ifdef GT_DEBUG_INFO

#include <execinfo.h>
#include <stdlib.h>

#define GT_STACK_TRACE_SIZE 15

inline void gt_print_stack_trace() {
  void *stack[GT_STACK_TRACE_SIZE];
  size_t size = backtrace(stack,GT_STACK_TRACE_SIZE);
  fprintf(stderr,">>GT::StackTrace\n");
  backtrace_symbols_fd(stack,size,2); // FIXME
  fprintf(stderr,"<<\n");
}

#else

inline void gt_print_stack_trace() {}

#endif


