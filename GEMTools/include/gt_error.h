/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_ERROR_H_
#define GT_ERROR_H_

#include <stdio.h>

inline void gt_print_stack_trace();

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
    gt_print_stack_trace(); \
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
 * Error msg
 */
#define gt_fatal_error_msg(gt_error_msg,args...) \
  do { \
  register FILE* const error_stream=gt_error_get_stream(); \
  fprintf(error_stream, \
    "Fatal error (%s:%d,%s)\n "gt_error_msg"\n", \
     GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
  fflush(error_stream); \
  gt_print_stack_trace(); \
  exit(1); \
  } while (0)
#define gt_error_msg(gt_error_msg,args...) \
  do { \
  register FILE* const error_stream=gt_error_get_stream(); \
  fprintf(error_stream, \
    "Fatal error (%s:%d,%s)\n "gt_error_msg"\n", \
     GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
  fflush(error_stream); \
  } while (0)
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
  #define gt_check_block(condition) if (condition)
#else
  #define gt_fatal_check(condition,gt_error_name,args...)
  #define gt_check(condition,gt_error_name,args...)
  #define gt_check_block(condition) if (false)
#endif

/*
 * ERROR CODES/MSG
 *   #define GT_ERROR_<CODE> "<MSG>"
 */

// Library/Program errors
#define GT_ERROR_NOT_ZERO "Value Zero. Variable %s must be non-zero"
#define GT_ERROR_POSITION_OUT_OF_RANGE "Requested position out of range"
#define GT_ERROR_POSITION_OUT_OF_RANGE_INFO "Requested position (%"PRIu64") out of range [%"PRIu64",%"PRId64"]"
#define GT_ERROR_SELECTION_NOT_IMPLEMENTED "Library error. Selection not implemented or corrupted value"
#define GT_ERROR_SELECTION_NOT_VALID "Library error. Selection not valid"
#define GT_ERROR_ALG_INCONSISNTENCY "Library error. Algorithmic inconsistency, check your program"

// Memory errors
#define GT_ERROR_MEM_HANDLER "Could not allocate handler"
#define GT_ERROR_MEM_ALLOC "Could not allocate memory"
#define GT_ERROR_MEM_ALLOC_INFO "Could not allocate memory (%"PRIu64" requested)"
#define GT_ERROR_MEM_REALLOC "Could not reallocate memory"
#define GT_ERROR_NULL_HANDLER "Null handler or fields not properly allocated"
#define GT_ERROR_NULL_HANDLER_INFO "Null handler %s "

// System errors
#define GT_ERROR_SYS_MMAP "Could not map file '%s' to memory"
#define GT_ERROR_SYS_UNMAP "Could not unmap memory"
#define GT_ERROR_SYS_THREAD "Could not create thread"
#define GT_ERROR_SYS_MUTEX "Mutex call error"
#define GT_ERROR_SYS_MUTEX_INIT "Mutex initialization error"
#define GT_ERROR_SYS_MUTEX_DESTROY "Mutex destroy call error"
#define GT_ERROR_SYS_COND_VAR "Conditional variable call error"
#define GT_ERROR_SYS_COND_VAR_INIT "Conditional variable initialization error"
#define GT_ERROR_SYS_COND_VAR_DESTROY "Conditional variable destroy call error"

// String errors
#define GT_ERROR_STRING_STATIC "Could not perform operation on static string"

// File errors
#define GT_ERROR_FILE_STAT "Could not stat file '%s'"
#define GT_ERROR_FILE_OPEN "Could not open file '%s'"
#define GT_ERROR_FILE_READ "Could not read from file '%s'"
#define GT_ERROR_FILE_WRITE "Could not write to file '%s'"
#define GT_ERROR_FILE_CLOSE "Could not close file '%s'"
#define GT_ERROR_FILE_FORMAT "Could not determine file format"
#define GT_ERROR_FILE_GZIP_OPEN "Could not open GZIPPED file '%s'"
#define GT_ERROR_FILE_BZIP_OPEN "Could not open BZIPPED file '%s'"

// Output errors
#define GT_ERROR_FPRINTF "Printing output. 'fprintf' call failed"
#define GT_ERROR_SPRINTF "Printing output. 'sprintf' call failed"
#define GT_ERROR_BPRINTF "Printing output. Buffer print formated 'gt_bprintf' call failed"
#define GT_ERROR_OFPRINTF "Printing output. Output File print formated 'gt_ofprintf' call failed"
#define GT_ERROR_BOFPRINTF "Printing output. Buffered Output file print formated 'gt_bofprintf' call failed"

#define GT_ERROR_PRINT_FORMAT "Incorrect print format. Expected format character"

// Template/Alignment/Map/Misms errors
#define GT_ERROR_MISMS_TYPE "Misms incorrect type"
#define GT_ERROR_MISMS_SPLICE_POS "Splicing distance must be positive (non-zero)"
#define GT_ERROR_COUNTERS_POS_STRATUM "Stratum must be strictly positive (stratum>0)"
#define GT_ERROR_MAP_MISMS_NOT_PARSED "Map's mismatches not parsed yet"
#define GT_ERROR_MAP_NEG_LENGTH "Negative Map total length"
#define GT_ERROR_ALIGNMENT_READ_QUAL_LENGTH "Read and quality length differs"
#define GT_ERROR_ALIGNMENT_MAPS_NOT_PARSED "Alignment's maps not parsed yet"
#define GT_ERROR_ALIGNMENT_INCONSISTENT_COUNTERS "Alignment inconsistency. Maps inconsistent with counters values"
#define GT_ERROR_TEMPLATE_MAPS_NOT_PARSED "Template's maps not parsed yet"
#define GT_ERROR_TEMPLATE_ZERO_BLOCKS "Zero alignment blocks (num_blocks_template>0)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT "Template inconsistency. Multimaps' members must be contained by single alignments"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ATTRB_RELATION "Template inconsistency. Incorrect number of mmaps and mmaps' attributes"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION "Template inconsistency. Incorrect number of matches' elements (check num_blocks_template)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_BLOCKS "Template inconsistency. Number of blocks must be the same across templates"
#define GT_ERROR_TEMPLATE_INCONSISTENT_COUNTERS "Template inconsistency. MMaps inconsistent with counters values"
#define GT_ERROR_TEMPLATE_ADD_BAD_NUM_BLOCKS "Trying to add wrong number of blocks to the template"
#define GT_ERROR_PALIGN_BAD_NUM_BLOCKS "Invalid Paired-alignment. Wrong number of alignment blocks (%"PRIu64")"

// Sequence Archive/Segmented Sequence errors
#define GT_ERROR_SEGMENTED_SEQ_IDX_OUT_OF_RANGE "Error accessing segmented sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
#define GT_ERROR_CDNA_IT_OUT_OF_RANGE "Error seeking sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
#define GT_ERROR_SEQ_ARCHIVE_NOT_FOUND "Sequence '%s' not found in reference archive"
#define GT_ERROR_SEQ_ARCHIVE_POS_OUT_OF_RANGE "Requested position '%"PRIu64"' out of sequence boundaries"
#define GT_ERROR_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE "Requested sequence string [%"PRIu64",%"PRIu64") out of sequence boundaries"

/*
 * Parsing FASTQ File format errors
 */
// IFP (Input FASTA Parser). General
#define GT_ERROR_PARSE_FASTA "Parsing FASTA/FASTQ error(%s:%"PRIu64")"

/*
 * Parsing MAP File format errors
 */
// IMP (Input MAP Parser). General
#define GT_ERROR_PARSE_MAP "Parsing MAP error(%s:%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_FILE_FORMAT "Parsing MAP error(%s:%"PRIu64"). Not a MAP file"
#define GT_ERROR_PARSE_MAP_BAD_NUMBER_FIELDS "Parsing MAP error(%s:%"PRIu64"). Wrong number of TAB separated fields (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_BAD_READ_QUAL_LENGTH "Parsing MAP error(%s:%"PRIu64"). Mismatching Read length (%"PRIu64") and Quality length (%"PRIu64")"
#define GT_ERROR_PARSE_MAP_COUNTERS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Error parsing counters"
#define GT_ERROR_PARSE_MAP_BAD_TEMPLATE_SEP "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Read character '%c' not valid (%s)"
#define GT_ERROR_PARSE_MAP_DIFF_TEMPLATE_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Different number of template blocks {read(%"PRIu64"),qualities(%"PRIu64")}"
#define GT_ERROR_PARSE_MAP_NOT_AN_ALIGNMENT "Parsing MAP error(%s:%"PRIu64"). File doesn't contains simple alignments (use template)"
#define GT_ERROR_PARSE_MAP_MAP_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Maps already parsed or null lazy-parsing handler"
#define GT_ERROR_PARSE_MAP_MISMS_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Mismatch string already parsed or null lazy-parsing handler"
#define GT_ERROR_PARSE_MAP_NOT_IMPLEMENTED "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Feature not implemented yet (sorry)"
#define GT_ERROR_PARSE_MAP_PREMATURE_EOL "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
#define GT_ERROR_PARSE_MAP_BAD_NUMBER_OF_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Wrong number of blocks"
#define GT_ERROR_PARSE_MAP_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Bad character found"
// IMP (Input MAP Parser). Parsing Read Errors
#define GT_ERROR_PARSE_MAP_READ_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing read, bad character found"
// IMP (Input MAP Parser). Parsing Qualities Errors
#define GT_ERROR_PARSE_MAP_QUAL_BAD_SEPARATOR "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, bad block-separator found"
#define GT_ERROR_PARSE_MAP_QUAL_BAD_LENGTH "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong length (w.r.t. read length)"
#define GT_ERROR_PARSE_MAP_QUAL_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong quality value (bad character)"
// IMP (Input MAP Parser). Parsing Counters Errors
#define GT_ERROR_PARSE_MAP_COUNTERS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing counters, bad character found"
// IMP (Input MAP Parser). Parsing Maps Errors
#define GT_ERROR_PARSE_MAP_MAP_BAD_NUMBER_OF_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, wrong number of blocks"
#define GT_ERROR_PARSE_MAP_MAP_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, bad character found"
#define GT_ERROR_PARSE_MAP_SPLIT_MAP_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad character found"
#define GT_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_ACCEPTORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of acceptors"
#define GT_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_DONORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of donors"
// IMP (Input MAP Parser). Parsing Mismatch String Errors
#define GT_ERROR_PARSE_MAP_MISMS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, bad character found"
#define GT_ERROR_PARSE_MAP_MISMS_BAD_MISMS_POS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, unsorted mismatches"

/*
 * Parsing SAM File format errors
 */
// ISP (Input SAM Parser). General
#define GT_ERROR_PARSE_SAM "Parsing SAM error(%s:%"PRIu64":%"PRIu64")"
#define GT_ERROR_PARSE_SAM_BAD_FILE_FORMAT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Not a SAM file"
#define GT_ERROR_PARSE_SAM_BAD_CHARACTER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Bad character found"
#define GT_ERROR_PARSE_SAM_UNMAPPED_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Unmapped read contains XA field (inconsistency)"
#define GT_ERROR_PARSE_SAM_PREMATURE_EOL "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
#define GT_ERROR_PARSE_SAM_EXPECTED_NUMBER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Expected number"
#define GT_ERROR_PARSE_SAM_WRONG_READ_CONTENT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Read in template doesn't match previously parse reads with same tag"
#define GT_ERROR_PARSE_SAM_CIGAR_PREMATURE_END "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature end of CIGAR string"
#define GT_ERROR_PARSE_SAM_WRONG_NUM_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Wrong number of eXtra mAps (as to pair them)"
#define GT_ERROR_PARSE_SAM_UNSOLVED_PENDING_MAPS "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Failed to pair maps"

// Output File
#define GT_ERROR_OUTPUT_FILE_INCONSISTENCY "Output file state inconsistent"
#define GT_ERROR_OUTPUT_FILE_FAIL_WRITE "Output file. Error writing to to file"
#define GT_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump"

/*
 * General purpose checkers
 */
#define GT_NULL_CHECK(object) gt_fatal_check(object==NULL,NULL_HANDLER_INFO,((char*)GT_QUOTE(object)))
#define GT_ZERO_CHECK(object) gt_fatal_check(object==0,NOT_ZERO,((char*)GT_QUOTE(object)))


#endif /* GT_ERROR_H_ */
