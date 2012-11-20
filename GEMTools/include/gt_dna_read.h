/*
 * PROJECT: GEM-Tools library
 * FILE: gt_dna_read.h
 * DATE: 20/08/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_DNA_READ_H_
#define GT_DNA_READ_H_

#include "gt_commons.h"

typedef struct {
  gt_string* tag;
  gt_string* read;
  gt_string* qualities;
} gt_dna_read;

typedef enum { GT_QUALS_OFFSET_33, GT_QUALS_OFFSET_64, GT_QUALS_OFFSET_INVALID } gt_qualities_offset_t;

/*
 * Checkers
 */
#define GT_DNA_READ_CHECK(read) \
  GT_NULL_CHECK(read); \
  GT_STRING_CHECK(read->tag); \
  GT_STRING_CHECK(read->read); \
  GT_STRING_CHECK(read->qualities)

/*
 * Constructor
 */
GT_INLINE gt_dna_read* gt_dna_read_new(void);
GT_INLINE void gt_dna_read_clear(gt_dna_read* read);
GT_INLINE void gt_dna_read_delete(gt_dna_read* read);

/*
 * Accessors
 */
GT_INLINE void gt_dna_read_set_ntag(gt_dna_read* read,char* const text,const uint64_t length);
GT_INLINE void gt_dna_read_set_tag(gt_dna_read* read,char* const text);
GT_INLINE char* gt_dna_read_get_tag(gt_dna_read* const read);

GT_INLINE void gt_dna_read_set_nread(gt_dna_read* read,char* const text,const uint64_t length);
GT_INLINE void gt_dna_read_set_read(gt_dna_read* read,char* const text);
GT_INLINE char* gt_dna_read_get_read(gt_dna_read* const read);

GT_INLINE void gt_dna_read_set_nqualities(gt_dna_read* read,char* const text,const uint64_t length);
GT_INLINE void gt_dna_read_set_qualities(gt_dna_read* read,char* const text);
GT_INLINE char* gt_dna_read_get_qualities(gt_dna_read* const read);

/*
 * Handlers
 */
GT_INLINE gt_qualities_offset_t gt_dna_read_get_qualities_offset(gt_dna_read* const read);

#endif /* GT_DNA_READ_H_ */
