/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.h
 * DATE: 3/09/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_SEQUENCE_ARCHIVE_H_
#define GT_SEQUENCE_ARCHIVE_H_

#include "gt_commons.h"

#define GT_SEQ_ARCHIVE_BLOCK_SIZE GT_BUFFER_SIZE_2M

#define gt_string_impl gt_dna_string

typedef struct {
  gt_vector* blocks;
  uint64_t sequence_total_length;
} gt_segmented_sequence;

typedef struct {
  gt_segmented_sequence* sequence;
  gt_strand strand;
  uint64_t global_pos;
  uint64_t local_pos;
  gt_string_impl* dna_string;
} gt_segmented_sequence_iterator;

typedef struct {
  gt_shash* sequences;
} gt_sequence_archive;

typedef struct {
  gt_sequence_archive* sequence_archive;
  // ... TODO
} gt_sequence_archive_iterator;

/*
 * Checkers
 */
#define GT_SEQ_ARCHIVE_CHECK(reference) gt_fatal_check(reference==NULL||reference->dna_sequences==NULL,NULL_HANDLER)
#define GT_SEG_DNA_SEQ_CHECK(dna_sequence) \
  gt_fatal_check(dna_sequence==NULL,NULL_HANDLER); \
  GT_VECTOR_CHECK(dna_sequence->blocks)

/*
 * Constructor
 */
GT_INLINE gt_sequence_archive* gt_sequence_archive_new(void);
GT_INLINE void gt_sequence_archive_clear(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive);
/*
 * Multi-Sequence handler
 */
GT_INLINE void gt_sequence_archive_add_sequence(gt_sequence_archive* const seq_archive,char* const seq_id,gt_segmented_sequence* const sequence);
GT_INLINE void gt_sequence_archive_remove_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);
GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int64_t (*gt_cmp_string)(char*,char*));
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive);
/*
 * Constructor
 */
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(void);
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence);
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence);
/*
 * Sequence handler
 */
GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t pos);
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t pos,const char character);
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,char* const string,const uint64_t length);
/*
 * Segmented Sequence Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t pos,gt_strand const strand,
    gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE void gt_segmented_sequence_iterator_seek(gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t pos);
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator);
GT_INLINE char gt_segmented_sequence_iterator_previous(gt_segmented_sequence_iterator* const sequence_iterator);
/*
 * Sequence Archive Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_previous(gt_sequence_archive_iterator* const seq_archive_iterator);


#endif /* GT_SEQUENCE_ARCHIVE_H_ */
