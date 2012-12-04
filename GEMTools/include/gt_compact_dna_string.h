/*
 * PROJECT: GEM-Tools library
 * FILE: gt_compact_dna_string.h
 * DATE: 20/08/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_COMPACT_DNA_STRING_H_
#define GT_COMPACT_DNA_STRING_H_

#include "gt_commons.h"

typedef struct {
  uint64_t* bitmaps;
  uint64_t allocated;
  uint64_t length;
} gt_compact_dna_string;

typedef struct {
  gt_compact_dna_string* cdna_string;
  uint64_t current_pos;
  uint64_t* current_bitmap;
  uint64_t bp_0;
  uint64_t bp_1;
  uint64_t bp_2;
} gt_compact_dna_sequence_iterator;

/*
 * Checkers
 */
#define GT_COMPACT_DNA_STRING_CHECK(cdna_string) gt_fatal_check(cdna_string==NULL||cdna_string->bitmaps==NULL,NULL_HANDLER)
#define GT_COMPACT_DNA_STRING_ITERATOR_CHECK(cdna_string_iterator) \
  GT_NULL_CHECK(cdna_string_iterator) \
  GT_COMPACT_DNA_STRING_CHECK(cdna_string_iterator->cdna_string) \
  GT_NULL_CHECK(cdna_string_iterator->current_bitmap)

/*
 * Constructor
 */
GT_INLINE gt_compact_dna_string* gt_cdna_string_new(const uint64_t initial_buffer_size);
GT_INLINE void gt_cdna_string_resize(gt_compact_dna_string* const cdna_string,const uint64_t new_buffer_size);
GT_INLINE void gt_cdna_string_clear(gt_compact_dna_string* const cdna_string);
GT_INLINE void gt_cdna_string_delete(gt_compact_dna_string* const cdna_string);

/*
 * Handlers
 */
GT_INLINE char gt_cdna_string_get_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos);
GT_INLINE void gt_cdna_string_set_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos,const char character);

GT_INLINE uint64_t gt_cdna_string_get_length(gt_compact_dna_string* const cdna_string);

GT_INLINE void gt_cdna_string_append_string(gt_compact_dna_string* const cdna_string,char* const string,const uint64_t length);

/*
 * Compact DNA String Sequence Iterator
 */
GT_INLINE void gt_cdna_string_new_iterator(
    gt_compact_dna_string* const cdna_string,const uint64_t pos,gt_strand const strand,
    gt_compact_dna_sequence_iterator* const cdna_string_iterator);
GT_INLINE void gt_cdna_string_iterator_seek(gt_compact_dna_sequence_iterator* const cdna_string_iterator,const uint64_t pos);
GT_INLINE bool gt_cdna_string_iterator_eos(gt_compact_dna_sequence_iterator* const cdna_string_iterator);
GT_INLINE char gt_cdna_string_iterator_next(gt_compact_dna_sequence_iterator* const cdna_string_iterator);
GT_INLINE char gt_cdna_string_iterator_previous(gt_compact_dna_sequence_iterator* const cdna_string_iterator);

#endif /* GT_COMPACT_DNA_STRING_H_ */
