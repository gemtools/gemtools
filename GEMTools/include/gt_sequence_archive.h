/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.h
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Data structures needed to store a dictionary of DNA-sequences(chromosomes,contigs,etc) indexed by tag
 *   Internal sequence representation is based on a memory-segmented DNA-string
 */

#ifndef GT_SEQUENCE_ARCHIVE_H_
#define GT_SEQUENCE_ARCHIVE_H_

#include "gt_commons.h"
#include "gt_map.h"
#include "gt_compact_dna_string.h"

typedef struct {
  gt_string* seq_name;
  gt_vector* blocks; /* (gt_compact_dna_string*) */
  uint64_t sequence_total_length;
} gt_segmented_sequence;

typedef struct {
  gt_segmented_sequence* sequence;
  gt_string_traversal direction;
  uint64_t global_pos;
  int64_t local_pos;
  gt_compact_dna_sequence_iterator cdna_string_iterator;
} gt_segmented_sequence_iterator;

typedef struct {
  gt_shash* sequences; /* (gt_segmented_sequence*) */
} gt_sequence_archive;

typedef struct {
  gt_sequence_archive* sequence_archive;
  gt_shash_element *shash_it;
} gt_sequence_archive_iterator;

/*
 * Checkers
 */
#define GT_SEQUENCE_ARCHIVE_CHECK(seq_archive) \
  GT_NULL_CHECK(seq_archive); \
  GT_NULL_CHECK(seq_archive->sequences)
#define GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator) \
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive_iterator->sequence_archive)

#define GT_SEGMENTED_SEQ_CHECK(segmented_sequence) \
  GT_NULL_CHECK(segmented_sequence); \
  GT_NULL_CHECK(segmented_sequence->blocks); \
  GT_STRING_CHECK(segmented_sequence->seq_name)
#define GT_SEGMENTED_SEQ_POSITION_CHECK(segmented_sequence,position) \
  gt_fatal_check(position>=segmented_sequence->sequence_total_length, \
      SEGMENTED_SEQ_IDX_OUT_OF_RANGE,position,segmented_sequence->sequence_total_length);
#define GT_SEGMENTED_SEQ_ITERATOR_CHECK(segmented_sequence_iterator) \
  GT_SEGMENTED_SEQ_CHECK(segmented_sequence_iterator->sequence)

/*
 * SequenceARCHIVE Constructor
 */
GT_INLINE gt_sequence_archive* gt_sequence_archive_new(void);
GT_INLINE void gt_sequence_archive_clear(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive);
/*
 * SequenceARCHIVE handler
 */
GT_INLINE void gt_sequence_archive_add_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence);
GT_INLINE void gt_sequence_archive_remove_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_get_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);

GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int (*gt_cmp_string)(char*,char*));
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive);
/*
 * SequenceARCHIVE Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_previous(gt_sequence_archive_iterator* const seq_archive_iterator);

/*
 * SegmentedSEQ Constructor
 */
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(void);
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence);
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence);
/*
 * SegmentedSEQ Sequence handler
 */
GT_INLINE void gt_segmented_sequence_set_name(gt_segmented_sequence* const sequence,char* const seq_name,const uint64_t seq_name_length);
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence);

GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t position);
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t position,const char character);
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,char* const string,const uint64_t length);

GT_INLINE void gt_segmented_sequence_get_sequence(
    gt_segmented_sequence* const sequence,const uint64_t position,const uint64_t length,gt_string* const string);
/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t position,gt_string_traversal const direction,
    gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE void gt_segmented_sequence_iterator_seek(gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t position,gt_string_traversal const direction);
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator);

#endif /* GT_SEQUENCE_ARCHIVE_H_ */
