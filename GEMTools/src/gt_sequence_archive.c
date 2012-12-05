/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.c
 * DATE: 3/09/2012
 * DESCRIPTION: // TODO
 */

#include "gt_sequence_archive.h"

#define GT_SEQ_ARCHIVE_NUM_BLOCKS 1000
#define GT_SEQ_ARCHIVE_BLOCK_SIZE GT_BUFFER_SIZE_2M

/*
 * SequenceARCHIVE Constructor
 */
GT_INLINE gt_sequence_archive* gt_sequence_archive_new(void) {
  gt_sequence_archive* seq_archive = malloc(sizeof(gt_sequence_archive));
  gt_cond_fatal_error(!seq_archive,MEM_HANDLER);
  seq_archive->sequences = gt_shash_new();
  return seq_archive;
}
GT_INLINE void gt_sequence_archive_clear(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  gt_shash_clean(seq_archive);
}
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  gt_shash_delete(seq_archive,true,true);
}

/*
 * SequenceARCHIVE handler
 */
GT_INLINE void gt_sequence_archive_add_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_shash_insert(seq_archive->sequences,gt_string_get_string(sequence->seq_name),sequence,gt_segmented_sequence*);
}
GT_INLINE void gt_sequence_archive_remove_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  gt_shash_remove_element(seq_archive->sequences,seq_id);
}
GT_INLINE void gt_sequence_archive_get_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  return gt_shash_get(seq_archive->sequences,seq_id,gt_segmented_sequence*);
}


int gt_sequence_archive_lexicographical_sort_fx(char *a,char *b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return gt_strcmp(a,b);
}
int gt_sequence_archive_karyotypic_sort_fx(char *a,char *b) {
//
}
GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int64_t (*gt_cmp_string)(char*,char*)) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);

}
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences,gt_sequence_archive_lexicographical_sort_fx);
}
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);

}

/*
 * SequenceARCHIVE Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator) {
  // TODO
}
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator) {
  // TODO
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator) {
  // TODO
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_previous(gt_sequence_archive_iterator* const seq_archive_iterator) {
  // TODO
}



/*
 * SegmentedSEQ Constructor
 */
GT_INLINE gt_segmented_sequence* gt_segmented_sequence_new(void) {
  gt_segmented_sequence* sequence = malloc(sizeof(gt_segmented_sequence));
  gt_cond_fatal_error(!sequence,MEM_HANDLER);
  sequence->blocks = gt_vector_new(GT_SEQ_ARCHIVE_NUM_BLOCKS,sizeof(gt_compact_dna_string*));
  sequence->sequence_total_length = 0;
  sequence->seq_name = gt_string_new(10);
  return sequence;
}
GT_INLINE void gt_segmented_sequence_clear(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_VECTOR_ITERATE(sequence->blocks,block,block_num,uint64_t*) {
    free(*block);
  }
  gt_vector_clean(sequence->blocks);
  gt_string_clear(sequence->seq_name);
}
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_VECTOR_ITERATE(sequence->blocks,block,block_num,uint64_t*) {
    free(*block);
  }
  gt_vector_delete(sequence->blocks);
  gt_string_delete(sequence->seq_name);
}
/*
 * SegmentedSEQ handler
 */
GT_INLINE void gt_segmented_sequence_set_name(char* const seq_name,const uint64_t seq_length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_string_set_nstring(sequence->seq_name,seq_name,seq_length);
}
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  return gt_string_get_string(sequence->seq_name);
}

GT_INLINE gt_compact_dna_string* gt_segmented_sequence_get_block(gt_segmented_sequence* const sequence,const uint64_t pos) {
  register const uint64_t num_block = pos/GT_SEQ_ARCHIVE_BLOCK_SIZE;
  register const uint64_t blocks_used = gt_vector_get_used(sequence->blocks);
  // Allocate new blocks (if needed)
  if (num_block>=blocks_used) {
    register uint64_t i;
    for (i=blocks_used;i<=num_block;++i) {
      gt_vector_insert(sequence->blocks,gt_cdna_string_new(GT_SEQ_ARCHIVE_BLOCK_SIZE),gt_compact_dna_string*);
    }
  }
  return *gt_vector_get_elm(sequence->blocks,num_block,gt_compact_dna_string*);
}
GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t pos) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_fatal_check(pos>=sequence->sequence_total_length,SEGMENTED_SEQ_IDX_OUT_OF_RANGE,pos,sequence->sequence_total_length);
  register const uint64_t pos_in_block = pos%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  return gt_cdna_string_get_char_at(gt_segmented_sequence_get_block(sequence,pos),pos_in_block);
}
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t pos,const char character) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_SEGMENTED_SEQ_CHECK(sequence);
  register const uint64_t pos_in_block = pos%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  // Adjust sequence total length
  if (pos>=sequence->sequence_total_length) sequence->sequence_total_length = pos+1;
  // Set character in compact dna string
  gt_cdna_string_set_char_at(gt_segmented_sequence_get_block(sequence,pos),pos_in_block,character);
}
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,char* const string,const uint64_t length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  register uint64_t chars_written = 0;
  while (chars_written < length) {
    register const uint64_t chunk_size = ((length-chars_written)<GT_SEQ_ARCHIVE_BLOCK_SIZE) ?
        length-chars_written : GT_SEQ_ARCHIVE_BLOCK_SIZE;
    gt_cdna_string_append_string(gt_segmented_sequence_get_block(sequence,pos),string+chars_written,chunk_size);
    chars_written += chunk_size;
  }
}

/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t pos,gt_strand const strand,
    gt_segmented_sequence_iterator* const sequence_iterator) {
  // TODO
}
GT_INLINE void gt_segmented_sequence_iterator_seek(gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t pos) {
  // TODO
}
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator) {
  // TODO
}
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator) {
  // TODO
}
GT_INLINE char gt_segmented_sequence_iterator_previous(gt_segmented_sequence_iterator* const sequence_iterator) {
  // TODO
}
