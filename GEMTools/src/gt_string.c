/*
 * PROJECT: GEM-Tools library
 * FILE: gt_string.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple string implementation.
 *   Static stings gt_string_new(0), which share memory across instances (stores mem ptr)
 *   Dynamic strings gt_string_new(n>0), which handle their own memory and hold copy of the string
 */

#include "gt_string.h"

#define GT_STRING_STATIC 0
#define GT_STRING_DEFAULT_BUFFER_SIZE 200

/*
 * Constructor & Accessors
 */
GT_INLINE gt_string* gt_string_new(const uint64_t initial_buffer_size) {
  gt_string* string = malloc(sizeof(gt_string));
  gt_cond_fatal_error(!string,MEM_HANDLER);
  // Initialize string
  if (gt_expect_true(initial_buffer_size>0)) {
    string->buffer = malloc(initial_buffer_size);
    gt_cond_fatal_error(!string->buffer,MEM_ALLOC_INFO,initial_buffer_size);
    string->buffer[0] = EOS;
  } else {
    string->buffer = NULL;
  }
  string->allocated = initial_buffer_size;
  string->length = 0;
  return string;
}
GT_INLINE void gt_string_resize(gt_string* const string,const uint64_t new_buffer_size) {
  GT_STRING_CHECK_BUFFER(string);
  if (string->allocated > 0 && string->allocated < new_buffer_size) {
    string->buffer = realloc(string->buffer,new_buffer_size);
    gt_cond_fatal_error(!string->buffer,MEM_REALLOC);
    string->allocated = new_buffer_size;
  }
}
GT_INLINE void gt_string_clear(gt_string* const string) {
  GT_STRING_CHECK(string);
  if (string->allocated) string->buffer[0] = EOS;
  else string->buffer = NULL;
  string->length = 0;
}
GT_INLINE void gt_string_delete(gt_string* const string) {
  GT_STRING_CHECK(string);
  if (string->allocated) free(string->buffer);
  free(string);
}

GT_INLINE bool gt_string_is_static(gt_string* const string) {
  GT_STRING_CHECK(string);
  return string->allocated==0;
}
GT_INLINE void gt_string_cast_static(gt_string* const string) {
  GT_STRING_CHECK(string);
  if (string->allocated > 0) {
    free(string->buffer);
    string->allocated = 0;
  }
  string->buffer = NULL;
  string->length = 0;
}
GT_INLINE void gt_string_cast_dynamic(gt_string* const string,const uint64_t initial_buffer_size) {
  GT_STRING_CHECK(string);
  if (gt_expect_false(initial_buffer_size==0)) {
    gt_string_cast_static(string);
  } else {
    if (string->buffer!=NULL) {
      string->buffer = gt_strndup(string->buffer,string->length);
      string->allocated = string->length+1;
    } else {
      string->buffer = malloc(initial_buffer_size);
      gt_cond_fatal_error(!string->buffer,MEM_ALLOC);
      string->buffer[0] = EOS;
    }
  }
}

GT_INLINE void gt_string_set_string(gt_string* const string,char* const string_src) {
  GT_NULL_CHECK(string_src);
  register const uint64_t length = strlen(string_src);
  gt_string_set_nstring(string,string_src,length);
}
GT_INLINE void gt_string_set_nstring(gt_string* const string,char* const string_src,const uint64_t length) {
  GT_STRING_CHECK(string);
  GT_NULL_CHECK(string_src);
  if (gt_expect_true(string->allocated>0)) {
    gt_string_resize(string,length+1);
    gt_strncpy(string->buffer,string_src,length);
  } else {
    string->buffer = string_src;
  }
  string->length = length;
}

GT_INLINE char* gt_string_get_string(gt_string* const string) {
  GT_STRING_CHECK(string);
  return string->buffer;
}
GT_INLINE uint64_t gt_string_get_length(gt_string* const string) {
  GT_STRING_CHECK(string);
  return string->length;
}
GT_INLINE void gt_string_set_length(gt_string* const string,const uint64_t length) {
  GT_STRING_CHECK(string);
  string->length = length;
}

GT_INLINE char* gt_string_char_at(gt_string* const string,const uint64_t pos) {
  GT_STRING_CHECK_BUFFER(string);
  gt_fatal_check(pos>string->length,POSITION_OUT_OF_RANGE_INFO,pos,(uint64_t)0,string->length);
  return string->buffer+pos;
}

GT_INLINE void gt_string_append_char(gt_string* const string_dst,char const character) {
  GT_STRING_CHECK_NO_STATIC(string_dst);
  gt_string_resize(string_dst,string_dst->length+1);
  string_dst->buffer[string_dst->length] = character; // NOTE: No EOS appended
  ++string_dst->length;
}
GT_INLINE void gt_string_append_eos(gt_string* const string_dst) {
  GT_STRING_CHECK_NO_STATIC(string_dst);
  gt_string_resize(string_dst,string_dst->length+1);
  string_dst->buffer[string_dst->length] = EOS;
}

GT_INLINE void gt_string_append_string(gt_string* const string_dst,char* const string_src,const uint64_t length) {
  GT_STRING_CHECK_NO_STATIC(string_dst);
  GT_NULL_CHECK(string_src);
  register const uint64_t final_length = string_dst->length+length;
  gt_string_resize(string_dst,final_length);
  gt_strncpy(string_dst->buffer+string_dst->length,string_src,length);
  string_dst->length = final_length;
}
GT_INLINE void gt_string_append_gt_string(gt_string* const string_dst,gt_string* const string_src) {
  GT_STRING_CHECK_NO_STATIC(string_dst);
  GT_STRING_CHECK(string_src);
  register const uint64_t final_length = string_dst->length+string_src->length;
  gt_string_resize(string_dst,final_length);
  gt_strncpy(string_dst->buffer+string_dst->length,string_src->buffer,string_src->length);
  string_dst->length = final_length;
}

/*
 * Cmp functions
 */
GT_INLINE bool gt_string_is_null(gt_string* const string) {
  GT_STRING_CHECK(string);
  return (gt_expect_true(string->allocated>0)) ? string->buffer[0]==EOS : string->buffer==NULL;
}
GT_INLINE int64_t gt_string_cmp(gt_string* const string_a,gt_string* const string_b) {
  GT_STRING_CHECK(string_a);
  GT_STRING_CHECK(string_b);
  register char* const buffer_a = string_a->buffer;
  register char* const buffer_b = string_b->buffer;
  register const uint64_t min_length = GT_MIN(string_a->length,string_b->length);
  register uint64_t i;
  register int8_t diff;
  for (i=0;i<min_length;++i) {
    if ((diff=(buffer_a[i]-buffer_b[i]))) {
      return (diff>0)?(i+1):(-i-1);
    }
  }
  if (string_a->length==string_b->length) {
    return 0;
  } else if (string_a->length<string_b->length) {
    return min_length+1;
  } else {
    return -min_length-1;
  }
}
GT_INLINE int64_t gt_string_ncmp(gt_string* const string_a,gt_string* const string_b,const uint64_t length) {
  GT_STRING_CHECK(string_a);
  GT_STRING_CHECK(string_b);
  register char* const buffer_a = string_a->buffer;
  register char* const buffer_b = string_b->buffer;
  register const uint64_t min_length = GT_MIN(GT_MIN(string_a->length,string_b->length),length);
  register uint64_t i;
  register int8_t diff;
  for (i=0;i<min_length;++i) {
    if ((diff=(buffer_a[i]-buffer_b[i]))) {
      return (diff>0)?(i+1):(-i-1);
    }
  }
  return 0;
}
GT_INLINE bool gt_string_equals(gt_string* const string_a,gt_string* const string_b) {
  GT_STRING_CHECK(string_a);
  GT_STRING_CHECK(string_b);
  return gt_string_cmp(string_a,string_b)==0;
}
GT_INLINE bool gt_string_nequals(gt_string* const string_a,gt_string* const string_b,const uint64_t length) {
  GT_STRING_CHECK(string_a);
  GT_STRING_CHECK(string_b);
  return gt_string_ncmp(string_a,string_b,length)==0;
}

/*
 * Handlers
 */
GT_INLINE gt_string* gt_string_dup(gt_string* const sequence) {
  GT_STRING_CHECK(sequence);
  register gt_string* sequence_cpy = gt_string_new(sequence->length+1);
  gt_strncpy(sequence_cpy->buffer,sequence->buffer,sequence->length);
  sequence_cpy->length = sequence->length;
  return sequence_cpy;
}
GT_INLINE void gt_string_copy(gt_string* const sequence_dst,gt_string* const sequence_src) {
  GT_STRING_CHECK_NO_STATIC(sequence_dst);
  GT_STRING_CHECK(sequence_src);
  gt_string_resize(sequence_dst,sequence_src->length+1);
  gt_strncpy(sequence_dst->buffer,sequence_src->buffer,sequence_src->length);
  sequence_dst->length = sequence_src->length;
}
GT_INLINE void gt_string_reverse_copy(gt_string* const sequence_dst,gt_string* const sequence_src) {
  GT_STRING_CHECK_NO_STATIC(sequence_dst);
  GT_STRING_CHECK(sequence_src);
  register const uint64_t string_length = sequence_src->length;
  gt_string_resize(sequence_dst,string_length+1);
  register char* const buffer_src = sequence_src->buffer;
  register char* const buffer_dst = sequence_dst->buffer;
  register uint64_t i;
  for (i=0;i<string_length;++i) {
    buffer_dst[i] = buffer_src[string_length-1-i];
  }
  buffer_dst[string_length] = EOS;
  sequence_dst->length = string_length;
}
GT_INLINE void gt_string_reverse(gt_string* const sequence) {
  GT_STRING_CHECK(sequence);
  register const uint64_t string_length = sequence->length;
  register const uint64_t middle = string_length/2;
  register char* const buffer = sequence->buffer;
  register uint64_t i;
  for (i=0;i<middle;++i) {
    register const char aux = buffer[i];
    buffer[i] = buffer[string_length-i-1];
    buffer[string_length-i-1] = aux;
  }
}

/*
 * String Printers
 */
GT_INLINE gt_status gt_vsprintf(gt_string* const sequence,const char *template,va_list v_args) {
  GT_STRING_CHECK(sequence);
  register gt_status chars_printed;
  if (!gt_string_is_static(sequence)) { // Allocate memory
    register const uint64_t mem_required = gt_calculate_memory_required_v(template,v_args);
    gt_string_resize(sequence,mem_required+1);
  }
  chars_printed=vsprintf(gt_string_get_string(sequence),template,v_args);
  gt_string_set_length(sequence,chars_printed);
  gt_string_get_string(sequence)[chars_printed] = EOS;
  return chars_printed;
}
GT_INLINE gt_status gt_sprintf(gt_string* const sequence,const char *template,...) {
  GT_STRING_CHECK(sequence);
  va_list v_args;
  va_start(v_args,template);
  register const gt_status chars_printed = gt_vsprintf(sequence,template,v_args);
  va_end(v_args);
  return chars_printed;
}
GT_INLINE gt_status gt_vsprintf_append(gt_string* const sequence,const char *template,va_list v_args) {
  GT_STRING_CHECK(sequence);
  register gt_status chars_printed = gt_string_get_length(sequence);
  if (!gt_string_is_static(sequence)) { // Allocate memory
    register const uint64_t mem_required = gt_calculate_memory_required_v(template,v_args);
    gt_string_resize(sequence,mem_required+chars_printed+1);
  }
  chars_printed+=vsprintf(gt_string_get_string(sequence)+chars_printed,template,v_args);
  gt_string_set_length(sequence,chars_printed);
  gt_string_get_string(sequence)[chars_printed] = EOS;
  return chars_printed;
}
GT_INLINE gt_status gt_sprintf_append(gt_string* const sequence,const char *template,...) {
  GT_STRING_CHECK(sequence);
  va_list v_args;
  va_start(v_args,template);
  register const gt_status chars_printed = gt_vsprintf_append(sequence,template,v_args);
  va_end(v_args);
  return chars_printed;
}
