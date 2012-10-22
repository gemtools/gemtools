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
GT_INLINE bool gt_streq(char* const buffer_a,char* const buffer_b) {
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
