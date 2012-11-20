/*
 * PROJECT: GEM-Tools library
 * FILE: gt_dna_read.c
 * DATE: 20/08/2012
 * DESCRIPTION: // TODO
 */

#include "gt_dna_read.h"

#define GT_DNA_READ_TAG_INITIAL_LENGTH 40
#define GT_DNA_READ_INITIAL_LENGTH 100

/*
 * Constructor
 */
GT_INLINE gt_dna_read* gt_dna_read_new(void) {
  gt_dna_read* read = malloc(sizeof(gt_dna_read));
  gt_cond_fatal_error(!read,MEM_HANDLER);
  read->tag = gt_string_new(GT_DNA_READ_TAG_INITIAL_LENGTH);
  read->read = gt_string_new(GT_DNA_READ_INITIAL_LENGTH);
  read->qualities = gt_string_new(GT_DNA_READ_INITIAL_LENGTH);
  return read;
}
GT_INLINE void gt_dna_read_clear(gt_dna_read* read) {
  GT_DNA_READ_CHECK(read);
  gt_string_clear(read->tag);
  gt_string_clear(read->read);
  gt_string_clear(read->read);
}
GT_INLINE void gt_dna_read_delete(gt_dna_read* read) {
  GT_DNA_READ_CHECK(read);
  gt_string_delete(read->tag);
  gt_string_delete(read->read);
  gt_string_delete(read->read);
}

/*
 * Accessors
 */
GT_INLINE void gt_dna_read_set_ntag(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->tag,text,length);
}
GT_INLINE void gt_dna_read_set_tag(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->tag,text);
}
GT_INLINE char* gt_dna_read_get_tag(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->tag);
}

GT_INLINE void gt_dna_read_set_nread(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->read,text,length);
}
GT_INLINE void gt_dna_read_set_read(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->read,text);
}
GT_INLINE char* gt_dna_read_get_read(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->read);
}

GT_INLINE void gt_dna_read_set_nqualities(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->qualities,text,length);
}
GT_INLINE void gt_dna_read_set_qualities(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->qualities,text);
}
GT_INLINE char* gt_dna_read_get_qualities(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->qualities);
}

/*
 * Handlers
 */
GT_INLINE gt_qualities_offset_t gt_dna_read_get_qualities_offset(gt_dna_read* const read) {
  register uint8_t min = UINT8_MAX;
  GT_STRING_ITERATE(read->qualities,string,pos) {
    if (string[pos]<min) min = string[pos];
  }
  if (min >= 64) return GT_QUALS_OFFSET_64;
  if (min >= 33) return GT_QUALS_OFFSET_33;
  return GT_QUALS_OFFSET_INVALID;
}
