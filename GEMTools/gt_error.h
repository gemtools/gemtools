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
extern FILE* gt_error_stream;
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
 * Robust checkers
 */
#ifndef GT_NO_CONSISTENCY_CHECKS
  #define gt_fatal_check(condition,gt_error_name,args...) gt_cond_fatal_error(condition,gt_error_name,##args)
  #define gt_check(condition,gt_error_name,args...) gt_cond_error(condition,gt_error_name,##args)
#else
  #define gt_fatal_check(condition,gt_error_name,args...)
  #define gt_check(condition,gt_error_name,args...)
#endif

/*
 * ERROR CODES/MSG
 *   #define GT_ERROR_<CODE> "<MSG>"
 */

// Memory errors
#define GT_ERROR_MEM_HANDLER "Could not allocate handler"
#define GT_ERROR_MEM_ALLOC "Could not allocate memory"
#define GT_ERROR_NULL_HANDLER "Null handler or fields not properly allocated"
#define GT_ERROR_NULL_HANDLER_INFO "Null handler %s "
#define GT_ERROR_NOT_ZERO "Value Zero. Variable %s must be non-zero"

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

// Template/Alignment/Map/Misms errors
#define GT_ERROR_POSITION_OUT_OF_RANGE "Requested position out of range"
#define GT_ERROR_POSITION_OUT_OF_RANGE_INFO "Requested position (%"PRIu64") out of range [%"PRId64",%"PRId64"]"
#define GT_ERROR_MISMS_TYPE "Misms incorrect type"
#define GT_ERROR_MISMS_SPLICE_POS "Splicing distance must be positive (non-zero)"
#define GT_ERROR_COUNTERS_POS_STRATUM "Stratum must be strictly positive (stratum>0)"
#define GT_ERROR_MAP_MISMS_NOT_PARSED "Map's mismatches not parsed yet"
#define GT_ERROR_ALIGN_READ_QUAL_LENGTH "Read and quality length differs"
#define GT_ERROR_ALIGN_MAPS_NOT_PARSED "Alignment's maps not parsed yet"
#define GT_ERROR_TEMPLATE_MAPS_NOT_PARSED "Template's maps not parsed yet"
#define GT_ERROR_TEMPLATE_ZERO_BLOCKS "Zero alignment blocks (num_blocks_template>0)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION "Template inconsistency. Incorrect number of matches' elements (check num_blocks_template)"
#define GT_ERROR_TEMPLATE_ADD_BAD_NUM_BLOCKS "Trying to add wrong number of blocks to the template"
#define GT_ERROR_PALIGN_BAD_NUM_BLOCKS "Invalid Paired-alignment. Wrong number of alignment blocks (%"PRIu64")"

// Parsing MAP File format errors
#define GT_ERROR_PARSE_MAP "Parsing MAP error(%s:%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_FILE_FORMAT "Parsing MAP error(%s:%"PRIu64"). Not a MAP file"
#define GT_ERROR_PARSE_MAP_BAD_NUMBER_FIELDS "Parsing MAP error(%s:%"PRIu64"). Wrong number of TAB separated fields (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_READ_QUAL_LENGTH "Parsing MAP error(%s:%"PRIu64"). Mismatching Read length (%"PRIu64") and Quality length (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_COUNTERS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Error parsing counters"
#define GT_ERROR_PARSE_MAP_BAD_TEMPLATE_SEP "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Read character '%c' not valid (%s)"
#define GT_ERROR_PARSE_MAP_MISMS_TEMPLATE_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Different number of template blocks {read(%"PRIu64"),qualities(%"PRIu64")}"
#define GT_ERROR_PARSE_MAP_NOT_AN_ALIGNMENT "Parsing MAP error(%s:%"PRIu64"). File doesn't contains simple alignments (use template)"

#define GT_ERROR_PARSE_MAP_MAP_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Maps already parsed or null lazy-parsing handler"
#define GT_ERROR_PARSE_MAP_MISMS_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Mismatch string already parsed or null lazy-parsing handler"


/*
 * General purpose checkers
 */
#define GT_NULL_CHECK(object) gt_fatal_check(object==NULL,NULL_HANDLER_INFO,((char*)GT_QUOTE(object)))
#define GT_ZERO_CHECK(object) gt_fatal_check(object==0,NOT_ZERO,((char*)GT_QUOTE(object)))


#endif /* GT_ERROR_H_ */
