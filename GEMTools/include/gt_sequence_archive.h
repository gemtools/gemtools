/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.h
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
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
  gt_strand strand;
  uint64_t global_pos;
  uint64_t local_pos;
  gt_compact_dna_string* dna_string;
} gt_segmented_sequence_iterator;

typedef struct {
  gt_shash* sequences; /* (gt_segmented_sequence*) */
} gt_sequence_archive;

typedef struct {
  gt_sequence_archive* sequence_archive;
  // ... TODO
} gt_sequence_archive_iterator;

/*
 * Checkers
 */
#define GT_SEQUENCE_ARCHIVE_CHECK(seq_archive) \
  GT_NULL_CHECK(seq_archive); \
  GT_NULL_CHECK(seq_archive->sequences)
#define GT_SEGMENTED_SEQ_CHECK(segmented_sequence) \
  GT_NULL_CHECK(segmented_sequence); \
  GT_NULL_CHECK(segmented_sequence->blocks); \
  GT_STRING_CHECK(segmented_sequence->seq_name)

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
GT_INLINE void gt_sequence_archive_get_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);

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
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(char* const seq_name,const uint64_t seq_name_length);
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence);
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence);
/*
 * SegmentedSEQ Sequence handler
 */
GT_INLINE void gt_segmented_sequence_set_name(char* const seq_name,const uint64_t seq_name_length);
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence);

GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t pos);
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t pos,const char character);
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,char* const string,const uint64_t length);
/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t pos,gt_strand const strand,
    gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE void gt_segmented_sequence_iterator_seek(gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t pos);
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_previous(gt_segmented_sequence_iterator* const sequence_iterator);


#endif /* GT_SEQUENCE_ARCHIVE_H_ */
