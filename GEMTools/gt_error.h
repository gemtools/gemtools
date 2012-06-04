/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_ERROR_H_
#define GT_ERROR_H_

#include <stdio.h>

// Base-name of the sources
#define GT_ERROR_BASENAME(S) \
  ({ register const char* const slash=strrchr((S),'/'); \
     slash ? slash + 1 : (S); })

// Gem-tools error/output streams
FILE* gt_error_stream=NULL;
// Getters/Setters
#define gt_error_get_stream() (gt_error_stream?gt_error_stream:stderr)
#define gt_error_set_stream(stream) gt_error_stream=(stream)

/*
 * Low-level error handling functions
 *
 * E.g. EXITING ERROR (FATAL)
 *  _gt_error_begin_block(gt_error_name,args...) {
 *    ..Code Block..
 *  } _gt_error_end_block(exit_error_code);
 *
 * E.g. NOT-EXITING
 *  _gt_error_begin_block(gt_error_name,args...) {
 *    ..Code Block..
 *  } _gt_error_end_block(exit_error_code);
 */
#define _gt_error_begin_block(gt_error_name,args...) \
  do { \
    register FILE* const error_stream=gt_error_get_stream(); \
    fprintf(error_stream, \
      "Fatal error (%s:%d,%s)\n "GT_ERROR_##gt_error_name"\n", \
       GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
    fflush(error_stream);
#define _gt_error_end_block() \
  } while (0)
#define _gt_error_end_block__exit(exit_error_code) \
    exit(exit_error_code); \
  } while (0)
/*
 * Error handlers
 */
#define _gt_fatal_error(error_handling_function,exit_error_code,gt_error_name,args...) \
  _gt_error_begin_block(gt_error_name,##args) { \
    error_handling_function; \
  } _gt_error_end_block__exit(exit_error_code)
#define _gt_error(error_handling_function,gt_error_name,args...) \
  _gt_error_begin_block(gt_error_name,##args) { \
    error_handling_function; \
  } _gt_error_end_block()
/*
 * Succinct error handlers
 */
#define gt_fatal_error(gt_error_name,args...) \
  _gt_error_begin_block(gt_error_name,##args) \
  _gt_error_end_block__exit(1)
#define gt_error(gt_error_name,args...) \
  _gt_error_begin_block(gt_error_name,##args) \
  _gt_error_end_block()
/*
 * Exception handlers (conditional error handlers)
 */
#define gt_cond_fatal_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error(gt_error_name,##args); \
    } \
  } while (0)
#define gt_cond_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_error(gt_error_name,##args); \
    } \
  } while (0)
/*
 * ERROR CODES/MSG
 *   #define GT_ERROR_<CODE> "<MSG>"
 */

// Memory errors
#define GT_ERROR_MEM_HANDLER "Could not allocate handler"
#define GT_ERROR_MEM_ALLOC "Could not allocate memory"

// System errors
#define GT_ERROR_SYS_MMAP "Could not map file '%s' to memory"
#define GT_ERROR_SYS_UNMAP "Could not unmap memory"
#define GT_ERROR_SYS_THREAD "Could not create thread"
#define GT_ERROR_SYS_MUTEX "Mutex call error"

// File errors
#define GT_ERROR_FILE_STAT "Could not stat file '%s'"
#define GT_ERROR_FILE_OPEN "Could not open file '%s'"
#define GT_ERROR_FILE_READ "Could not read from file '%s'"
#define GT_ERROR_FILE_WRITE "Could not write to file '%s'"
#define GT_ERROR_FILE_CLOSE "Could not close file '%s'"
#define GT_ERROR_FILE_FORMAT "Could not determine file format"

#endif /* GT_ERROR_H_ */
