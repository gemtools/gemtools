/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * DESCRIPTION: Base module containing general purpose functions
 */

#include "gt_commons.h"

/*
 * String/Buffer functions
 */
GT_INLINE void gt_strncpy(char* const buffer_dst,char* const buffer_src,const uint64_t length) {
  GT_NULL_CHECK(buffer_dst); GT_NULL_CHECK(buffer_src);
  memcpy(buffer_dst,buffer_src,length);
  buffer_dst[length] = EOS;
}
GT_INLINE char* gt_strndup(char* const buffer,const uint64_t length) {
  GT_NULL_CHECK(buffer);
  register char* const buffer_cpy = malloc(length+1);
  gt_cond_fatal_error(!buffer_cpy,MEM_ALLOC);
  gt_strncpy(buffer_cpy,buffer,length);
  return buffer_cpy;
}
GT_INLINE int gt_strcmp(char* const buffer_a,char* const buffer_b) {
  GT_NULL_CHECK(buffer_a); GT_NULL_CHECK(buffer_b);
  return strcmp(buffer_a,buffer_b);
}
GT_INLINE bool gt_streq(char* const buffer_a,char* const buffer_b) { // TODO: Add safe mem str_cmp
  GT_NULL_CHECK(buffer_a); GT_NULL_CHECK(buffer_b);
  return strcmp(buffer_a,buffer_b)==0;
}
GT_INLINE int gt_strncmp(char* const buffer_a,char* const buffer_b,const uint64_t length) {
  GT_NULL_CHECK(buffer_a); GT_NULL_CHECK(buffer_b);
  return strncmp(buffer_a,buffer_b,length);
}
GT_INLINE bool gt_strneq(char* const buffer_a,char* const buffer_b,const uint64_t length) {
  GT_NULL_CHECK(buffer_a); GT_NULL_CHECK(buffer_b);
  return strncmp(buffer_a,buffer_b,length)==0;
}

/*
 * Memory usage helper functions
 */
GT_INLINE uint64_t gt_calculate_memory_required_v(const char *template,va_list v_args) {
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Copy to avoid spoiling v_args
  va_list v_args_cpy;
  va_copy(v_args_cpy,v_args);
  // Calculate memory required to print the template{v_args}
  register uint64_t mem_required = 0, precision=0;
  register const char* centinel;
  for (centinel=template;*centinel!=EOS;++centinel,precision=0) {
    if (*centinel==FORMAT) {
      ++centinel;
      // Read modifiers
      while (gt_is_number(*centinel)) ++centinel;
      if (*centinel==DOT){
        ++centinel;
        if (*centinel==STAR) {
          ++centinel;
          precision = va_arg(v_args_cpy,int);
        } else {
          while (gt_is_number(*centinel)) ++centinel;
        }
      }
      gt_check(centinel==EOS,PRINT_FORMAT);
      // Check format
      switch (*centinel) {
        case 's': { // String requires fetching the argument length // FIXME: %.*s
          register char* const string = va_arg(v_args_cpy,char*);
          mem_required += (precision>0) ? precision : strlen(string);
          break;
        }
        default:
          // As for the rest, we estimate the memory usage
          // Also we assume an upper bound over the possible formats (int, chars, floats, ...)
          va_arg(v_args_cpy,int);
          mem_required+=20;
          break;
      }
    } else {
      ++mem_required;
    }
  }
  return mem_required;
}
GT_INLINE uint64_t gt_calculate_memory_required_va(const char *template,...) {
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  return gt_calculate_memory_required_v(template,v_args);
}
