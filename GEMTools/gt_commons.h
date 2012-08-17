/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_COMMONS_H_
#define GT_COMMONS_H_

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include <string.h>
#include <math.h>
#include <stdarg.h>

#include <ctype.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>

#include <errno.h>
#include <err.h>
#include <assert.h>

#include <pthread.h>

// Internally to Gem-tools error codes are returned as gt_status
typedef int32_t gt_status;

// Special characters
#define EOS '\0'
#define EOL '\n'
#define TAB '\t'
#define DOS_EOL '\r'
#define PLUS '+'
#define MINUS '-'
#define FORMAT '%'
#define SPACE ' '
#define SLASH '/'
#define STAR '*'
#define EQUAL '='
#define COMA ','
#define SEMICOLON ';'

// Buffer sizes
#define GT_BUFFER_SIZE_1K   ((1<<10)-64)
#define GT_BUFFER_SIZE_2K   ((1<<11)-64)
#define GT_BUFFER_SIZE_4K   ((1<<12)-64)
#define GT_BUFFER_SIZE_8K   ((1<<13)-64)
#define GT_BUFFER_SIZE_16K  ((1<<14)-64)
#define GT_BUFFER_SIZE_32K  ((1<<15)-64)
#define GT_BUFFER_SIZE_64K  ((1<<16)-64)
#define GT_BUFFER_SIZE_128K ((1<<17)-64)
#define GT_BUFFER_SIZE_256K ((1<<18)-64)
#define GT_BUFFER_SIZE_512K ((1<<19)-64)
#define GT_BUFFER_SIZE_1M  ((1<<20)-64)
#define GT_BUFFER_SIZE_2M  ((1<<21)-64)
#define GT_BUFFER_SIZE_4M  ((1<<22)-64)
#define GT_BUFFER_SIZE_8M  ((1<<23)-64)
#define GT_BUFFER_SIZE_16M ((1<<24)-64)
#define GT_BUFFER_SIZE_32M ((1<<25)-64)
#define GT_BUFFER_SIZE_64M ((1<<26)-64)
#define GT_BUFFER_SIZE_128M ((1<<27)-64)
#define GT_BUFFER_SIZE_256M ((1<<28)-64)
#define GT_BUFFER_SIZE_512M ((1<<29)-64)

// Number of lines
#define GT_NUM_LINES_1K      (1000)
#define GT_NUM_LINES_2K      (2000)
#define GT_NUM_LINES_5K      (5000)
#define GT_NUM_LINES_10K    (10000)
#define GT_NUM_LINES_20K    (20000)
#define GT_NUM_LINES_50K    (50000)
#define GT_NUM_LINES_100K  (100000)
#define GT_NUM_LINES_200K  (200000)
#define GT_NUM_LINES_500K  (500000)
#define GT_NUM_LINES_1M   (1000000)
#define GT_NUM_LINES_2M   (2000000)
#define GT_NUM_LINES_5M   (5000000)
#define GT_NUM_LINES_10M (10000000)
#define GT_NUM_LINES_20M (20000000)
#define GT_NUM_LINES_50M (50000000)

// Conditional expect
#define gt_expect_true(condition) __builtin_expect(condition,1)
#define gt_expect_false(condition) __builtin_expect(condition,0)

// GemTools Inline
#define GT_INLINE inline

// Macro Stringify
#define GT_QUOTE(value) #value

/*
 * Is functions
 */
extern bool gt_dna[256];
extern char gt_complement_table[256];
#define gt_is_dna(character) gt_dna[(int)character]
#define gt_get_complement(character) (gt_complement_table[(int)character])
#define gt_is_valid_quality(character) (33 <= (character) && (character) <= 127)
#define gt_is_number(character) ('0' <= (character) && (character) <= '9')
#define gt_get_cipher(character) ((character) - '0')
// MAP/MAPQ related
#define gt_is_valid_template_separator(character) \
  (character==' ' || character=='|' || character=='#')
#define gt_is_valid_counter_separator(character) \
  (character=='+' || character==':' || character=='x')

/*
 * Helper functions (OPERATIVE)
 */
#define gt_cfree(handler) if (handler!=NULL) free(handler);
#define gt_string_eq(a,b) (a!=NULL && b!=NULL && strcmp(a,b)==0)
GT_INLINE void gt_reverse_complent(char* const sequence,const uint64_t length);
GT_INLINE void gt_cpy_reverse_complent(char* const sequence_dst,char* const sequence_ori,const uint64_t length);
GT_INLINE void gt_reverse(char* const sequence,const uint64_t length);
GT_INLINE void gt_cpy_reverse(char* const sequence_dst,char* const sequence_ori,const uint64_t length);
GT_INLINE char* gt_string_cpy(char* sequence,const uint64_t length);
#define GT_MIN(a,b) ((a)<=(b)?(a):(b))
#define GT_MAX(a,b) ((a)>=(b)?(a):(b))
#define GT_ABS(a) ((a)>=0?(a):-(a))

/*
 * Error value return wrapper
 */
#define GT_DELEGATE_ERROR(funtion_call,error_code) if ((error_code=funtion_call)) { return error_code; }

/*
 * Mutex Helpers
 */
#define GT_BEGIN_MUTEX_SECTION(mutex) \
  gt_cond_fatal_error(pthread_mutex_lock(&(mutex)),SYS_MUTEX);
#define GT_END_MUTEX_SECTION(mutex) \
  gt_cond_fatal_error(pthread_mutex_unlock(&(mutex)),SYS_MUTEX);

/*
 * GEM-Tools basic includes
 */
#include "gt_error.h"
#include "gt_vector.h"
#include "gt_hash.h"

/*
 * Common constants
 */
#define GT_STREAM_FILE_NAME "<<STREAM>>"
#define GT_ALL UINT64_MAX

#endif /* GT_COMMONS_H_ */
