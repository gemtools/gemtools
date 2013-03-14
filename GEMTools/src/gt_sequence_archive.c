/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.c
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Data structures needed to store a dictionary of DNA-sequences(chromosomes,contigs,etc) indexed by tag
 *   Internal sequence representation is based on a memory-segmented DNA-string
 */

#include "gt_sequence_archive.h"

#define GT_SEQ_ARCHIVE_NUM_BLOCKS 15000
#define GT_SEQ_ARCHIVE_BLOCK_SIZE GT_BUFFER_SIZE_256K

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
  GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->sequences,sequence,gt_segmented_sequence) {
    gt_segmented_sequence_delete(sequence);
  } GT_SHASH_END_ITERATE;
  gt_shash_clear(seq_archive->sequences,false);
}
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->sequences,sequence,gt_segmented_sequence) {
    gt_segmented_sequence_delete(sequence);
  } GT_SHASH_END_ITERATE;
  gt_shash_delete(seq_archive->sequences,false);
  free(seq_archive);
}

/*
 * SequenceARCHIVE handler
 */
GT_INLINE void gt_sequence_archive_add_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_shash_insert(seq_archive->sequences,gt_string_get_string(sequence->seq_name),sequence,gt_segmented_sequence);
}
GT_INLINE void gt_sequence_archive_remove_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  // TODO: Retrieve seq and deallocate by hand
  gt_shash_remove(seq_archive->sequences,seq_id,false);
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_get_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  return gt_shash_get(seq_archive->sequences,seq_id,gt_segmented_sequence);
}

/*
 * SequenceARCHIVE High-level Retriever
 */
GT_INLINE gt_status gt_sequence_archive_get_sequence_string(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,gt_string* const string) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_NULL_CHECK(seq_id);
  GT_ZERO_CHECK(length);
  GT_STRING_CHECK_NO_STATIC(string);
  register gt_status error_code;
  // Retrieve the sequence
  register gt_segmented_sequence* seg_seq = gt_sequence_archive_get_sequence(seq_archive,seq_id);
  if (seg_seq==NULL) return GT_SEQ_ARCHIVE_NOT_FOUND;
  // Get the actual chunk
  if ((error_code=gt_segmented_sequence_get_sequence(seg_seq,position,length,string))) {
    return error_code;
  }
  // RC (if needed)
  if (strand==REVERSE) {
    gt_dna_string_reverse_complement(string);
  }
  return 0;
}
GT_INLINE gt_status gt_sequence_archive_retrieve_sequence_chunk(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,const uint64_t extra_length,gt_string* const string) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_ZERO_CHECK(position);
  GT_NULL_CHECK(seq_id);
  GT_STRING_CHECK_NO_STATIC(string);
  // Retrieve the sequence
  register gt_segmented_sequence* seg_seq = gt_sequence_archive_get_sequence(seq_archive,seq_id);
  if (seg_seq==NULL) {
    gt_error(SEQ_ARCHIVE_NOT_FOUND,seq_id);
    return GT_SEQ_ARCHIVE_NOT_FOUND;
  }
  // Calculate sequence's boundaries
  register const uint64_t sequence_total_length = seg_seq->sequence_total_length;
  register uint64_t init_position, total_length;
  if (position-1 >= seg_seq->sequence_total_length) { // Check position
    gt_error(SEQ_ARCHIVE_POS_OUT_OF_RANGE,position-1);
    return GT_SEQ_ARCHIVE_POS_OUT_OF_RANGE;
  }
  // Adjust init_position,total_length wrt strand
  init_position = position-1;
  if (strand==REVERSE) {
    init_position = (extra_length>init_position) ? 0 : init_position-extra_length;
  }
  total_length = length+extra_length;
  if (total_length >= sequence_total_length) total_length = seg_seq->sequence_total_length-1;
  // Get the actual chunk
  //  register gt_status error_code; /* Error checking & reporting version */
  //  if ((error_code=gt_segmented_sequence_get_sequence(seg_seq,init_position,total_length,string))) {
  //    gt_error(SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE,init_position,init_position+total_length,seq_id);
  //    return GT_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE;
  //  }
  gt_segmented_sequence_get_sequence(seg_seq,init_position,total_length,string);
  // RC (if needed)
  if (strand==REVERSE) {
    gt_dna_string_reverse_complement(string);
  }
  return 0;
}


/*
 * SequenceARCHIVE sorting functions
 */
int gt_sequence_archive_lexicographical_sort_fx(char *a,char *b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return gt_strcmp(a,b);
}
int gt_sequence_archive_karyotypic_sort_fx(char *a,char *b) {
  /*
   * Karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs
   */
  register const uint64_t alen = strlen(a);
  register const uint64_t alast = alen>1 ? alen-1 : 0;
  register const uint64_t blen = strlen(b);
  register const uint64_t blast = blen>1 ? blen-1 : 0;
  register const int str_cmp_ab = gt_strcmp(a,b);
  if (str_cmp_ab==0) return 0;
  // Chromosome M
  if (strncmp(a,"chrM",4)==0 || a[alast]=='M' || a[alast]=='m') return INT16_MAX;
  if (strncmp(b,"chrM",4)==0 || b[blast]=='M' || b[blast]=='m') return -INT16_MAX;
  // Chromosome Y
  if (strncmp(a,"chrY",4)==0 || a[alast]=='Y' || a[alast]=='y') return INT16_MAX-1;
  if (strncmp(b,"chrY",4)==0 || b[blast]=='Y' || b[blast]=='y') return -(INT16_MAX-1);
  // Chromosome X
  if (strncmp(a,"chrX",4)==0 || a[alast]=='X' || a[alast]=='x') return INT16_MAX-2;
  if (strncmp(b,"chrX",4)==0 || b[blast]=='X' || b[blast]=='x') return -(INT16_MAX-2);
  // Other Chromosome
  return str_cmp_ab;
}

#define gt_cmp_string_wrapper(arg1,arg2) gt_cmp_string((char*)arg1,(char*)arg2)
#define gt_sequence_archive_lexicographical_sort_fx_wrapper(arg1,arg2) gt_sequence_archive_lexicographical_sort_fx((char*)arg1,(char*)arg2)
#define gt_sequence_archive_karyotypic_sort_fx_wrapper(arg1,arg2) gt_sequence_archive_karyotypic_sort_fx((char*)arg1,(char*)arg2)
GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int (*gt_cmp_string)(char*,char*)) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_cmp_string_wrapper);
}
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_sequence_archive_lexicographical_sort_fx_wrapper);
}
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_sequence_archive_karyotypic_sort_fx_wrapper);
}

/*
 * SequenceARCHIVE Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_NULL_CHECK(seq_archive_iterator);
  seq_archive_iterator->sequence_archive = seq_archive;
  seq_archive_iterator->shash_it = seq_archive->sequences->shash_head;
}
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator);
  return seq_archive_iterator->shash_it==NULL;
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator);
  if (seq_archive_iterator->shash_it) {
    register gt_segmented_sequence* elm =  seq_archive_iterator->shash_it->element;
    seq_archive_iterator->shash_it = seq_archive_iterator->shash_it->hh.next;
    return elm;
  } else {
    return NULL;
  }
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
    if (*block) free(*block);
  }
  gt_vector_clear(sequence->blocks);
  gt_string_clear(sequence->seq_name);
}
GT_INLINE void gt_segmented_sequence_delete(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_VECTOR_ITERATE(sequence->blocks,block,block_num,gt_compact_dna_string*) {
    if (*block) gt_cdna_string_delete(*block);
  }
  gt_vector_delete(sequence->blocks);
  gt_string_delete(sequence->seq_name);
  free(sequence);
}
/*
 * SegmentedSEQ handler
 */
GT_INLINE void gt_segmented_sequence_set_name(gt_segmented_sequence* const sequence,char* const seq_name,const uint64_t seq_name_length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_string_set_nstring(sequence->seq_name,seq_name,seq_name_length);
}
GT_INLINE char* gt_segmented_sequence_get_name(gt_segmented_sequence* const sequence) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  return gt_string_get_string(sequence->seq_name);
}

GT_INLINE gt_compact_dna_string* gt_segmented_sequence_get_block(gt_segmented_sequence* const sequence,const uint64_t position) {
  register const uint64_t num_block = position/GT_SEQ_ARCHIVE_BLOCK_SIZE;
  register const uint64_t blocks_used = gt_vector_get_used(sequence->blocks);
  // Allocate new blocks (if needed)
  if (num_block>=blocks_used) {
    register uint64_t i;
    for (i=blocks_used;i<num_block;++i) {
      gt_vector_insert(sequence->blocks,NULL,gt_compact_dna_string*);
    }
    register gt_compact_dna_string* const block = gt_cdna_string_new(GT_SEQ_ARCHIVE_BLOCK_SIZE);
    gt_vector_insert(sequence->blocks,block,gt_compact_dna_string*);
    return block;
  } else {
    register gt_compact_dna_string* block = *gt_vector_get_elm(sequence->blocks,num_block,gt_compact_dna_string*);
    if (!block) {
      block = gt_cdna_string_new(GT_SEQ_ARCHIVE_BLOCK_SIZE);
      gt_vector_set_elm(sequence->blocks,num_block,gt_compact_dna_string*,block);
    }
    return block;
  }
}
GT_INLINE char gt_segmented_sequence_get_char_at(gt_segmented_sequence* const sequence,const uint64_t position) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence,position);
  register const uint64_t pos_in_block = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  return gt_cdna_string_get_char_at(gt_segmented_sequence_get_block(sequence,position),pos_in_block);
}
GT_INLINE void gt_segmented_sequence_set_char_at(gt_segmented_sequence* const sequence,const uint64_t position,const char character) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  register const uint64_t pos_in_block = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  // Adjust sequence total length
  if (position>=sequence->sequence_total_length) sequence->sequence_total_length = position+1;
  // Set character in compact dna string
  gt_cdna_string_set_char_at(gt_segmented_sequence_get_block(sequence,position),pos_in_block,character);
}
GT_INLINE void gt_segmented_sequence_append_string(gt_segmented_sequence* const sequence,char* const string,const uint64_t length) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  register uint64_t block_free_space = GT_SEQ_ARCHIVE_BLOCK_SIZE-(sequence->sequence_total_length%GT_SEQ_ARCHIVE_BLOCK_SIZE);
  register uint64_t current_length = sequence->sequence_total_length;
  register uint64_t chars_written = 0;
  while (chars_written < length) {
    register const uint64_t chunk_size = ((length-chars_written)<block_free_space) ?
        length-chars_written : block_free_space;
    gt_cdna_string_append_string(
        gt_segmented_sequence_get_block(sequence,current_length),string+chars_written,chunk_size);
    chars_written += chunk_size;
    current_length += chunk_size;
    block_free_space = GT_SEQ_ARCHIVE_BLOCK_SIZE;
  }
  sequence->sequence_total_length = current_length;
}

GT_INLINE gt_status gt_segmented_sequence_get_sequence(
    gt_segmented_sequence* const sequence,const uint64_t position,const uint64_t length,gt_string* const string) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence,position);
  GT_STRING_CHECK(string);
  GT_ZERO_CHECK(length);
  // Clear string
  gt_string_clear(string);
  // Check position
  if (gt_expect_false(position >= sequence->sequence_total_length)) return GT_SEQ_ARCHIVE_POS_OUT_OF_RANGE;
  // Retrieve String
  register uint64_t i=0;
  gt_segmented_sequence_iterator sequence_iterator;
  gt_segmented_sequence_new_iterator(sequence,position,GT_ST_FORWARD,&sequence_iterator);
  while (i<length && !gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
    gt_string_append_char(string,gt_segmented_sequence_iterator_next(&sequence_iterator));
    ++i;
  }
  gt_string_append_eos(string);
  return (i==length) ? GT_SEQ_ARCHIVE_OK : GT_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE;
}

/*
 * SegmentedSEQ Iterator
 */
GT_INLINE void gt_segmented_sequence_new_iterator(
    gt_segmented_sequence* const sequence,const uint64_t position,gt_string_traversal const direction,
    gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_CHECK(sequence);
  GT_NULL_CHECK(sequence_iterator);
  // Set iterator
  sequence_iterator->sequence = sequence;
  if (gt_expect_true(position<sequence->sequence_total_length)) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,position,direction);
  } else {
    sequence_iterator->global_pos = position; // EOS
  }
}
GT_INLINE void gt_segmented_sequence_iterator_seek(
    gt_segmented_sequence_iterator* const sequence_iterator,const uint64_t position,gt_string_traversal const direction) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  GT_SEGMENTED_SEQ_POSITION_CHECK(sequence_iterator->sequence,position);
  // Set iterator direction
  sequence_iterator->direction = direction;
  // Set sequence location fields
  sequence_iterator->global_pos = position;
  sequence_iterator->local_pos = position%GT_SEQ_ARCHIVE_BLOCK_SIZE;
  // Init sequence locator
  register gt_compact_dna_string* const cdna_string = gt_segmented_sequence_get_block(sequence_iterator->sequence,position);
  gt_cdna_string_new_iterator(cdna_string,sequence_iterator->local_pos,direction,&sequence_iterator->cdna_string_iterator);
}
GT_INLINE bool gt_segmented_sequence_iterator_eos(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  return (sequence_iterator->global_pos>=sequence_iterator->sequence->sequence_total_length);
}
GT_INLINE char gt_segmented_sequence_iterator_following(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  // Seek to proper position (load block if necessary)
  if (sequence_iterator->local_pos==GT_SEQ_ARCHIVE_BLOCK_SIZE) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,sequence_iterator->global_pos,sequence_iterator->direction);
    gt_check(sequence_iterator->global_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE != 0,ALG_INCONSISNTENCY);
  }
  // Update position
  ++sequence_iterator->global_pos;
  ++sequence_iterator->local_pos;
  // Return next!
  return gt_cdna_string_iterator_next(&sequence_iterator->cdna_string_iterator);
}
GT_INLINE char gt_segmented_sequence_iterator_previous(gt_segmented_sequence_iterator* const sequence_iterator) {
  GT_SEGMENTED_SEQ_ITERATOR_CHECK(sequence_iterator);
  // Seek to proper position (load block if necessary)
  if (sequence_iterator->local_pos==-1) {
    gt_segmented_sequence_iterator_seek(sequence_iterator,sequence_iterator->global_pos,sequence_iterator->direction);
    gt_check(sequence_iterator->local_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE ==
        sequence_iterator->global_pos%GT_SEQ_ARCHIVE_BLOCK_SIZE,ALG_INCONSISNTENCY);
  }
  // Update position
  --sequence_iterator->global_pos;
  --sequence_iterator->local_pos;
  // Return next!
  return gt_cdna_string_iterator_next(&sequence_iterator->cdna_string_iterator);
}
GT_INLINE char gt_segmented_sequence_iterator_next(gt_segmented_sequence_iterator* const sequence_iterator) {
  return (sequence_iterator->direction==GT_ST_FORWARD) ?
    gt_segmented_sequence_iterator_following(sequence_iterator) :
    gt_segmented_sequence_iterator_previous(sequence_iterator);
}
