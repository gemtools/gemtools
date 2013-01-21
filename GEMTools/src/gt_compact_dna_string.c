/*
 * PROJECT: GEM-Tools library
 * FILE: ggt_compact_dna_string.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_compact_dna_string.h"

/*
 * Blocks encoding dimensions
 */
#define GT_CDNA_BLOCK_CHARS 64
#define GT_CDNA_BLOCK_BITMAP_SIZE 8
#define GT_CDNA_BLOCK_BITMAPS 3
#define GT_CDNA_BLOCK_SIZE (GT_CDNA_BLOCK_BITMAP_SIZE*GT_CDNA_BLOCK_BITMAPS)

/*
 * CDNA Bitmaps Encoding
 */
#define GT_CDNA_ENC_CHAR_A 0
#define GT_CDNA_ENC_CHAR_C 1
#define GT_CDNA_ENC_CHAR_G 2
#define GT_CDNA_ENC_CHAR_T 3
#define GT_CDNA_ENC_CHAR_N 4

#define GT_CDNA_ZERO_MASK 0xFFFFFFFFFFFFFFFEull
#define GT_CDNA_ONE_MASK  0x0000000000000001ull

#define GT_CDNA_EXTRACT_MASK GT_CDNA_ONE_MASK

#define GT_CDNA_ENCODED_CHAR_BM0 ((uint8_t)1)
#define GT_CDNA_ENCODED_CHAR_BM1 ((uint8_t)2)
#define GT_CDNA_ENCODED_CHAR_BM2 ((uint8_t)4)

const uint8_t gt_cdna_encoded_char_bm[3] = {
  GT_CDNA_ENCODED_CHAR_BM0, GT_CDNA_ENCODED_CHAR_BM1, GT_CDNA_ENCODED_CHAR_BM2
};

/*
 * CDNA translation tables
 */
const char gt_cdna_decode[8] = {
  GT_DNA_CHAR_A, GT_DNA_CHAR_C, GT_DNA_CHAR_G, GT_DNA_CHAR_T,
  GT_DNA_CHAR_N, GT_DNA_CHAR_N, GT_DNA_CHAR_N
};
const uint8_t gt_cdna_encode[256] = {
    [0 ... 255] = GT_CDNA_ENC_CHAR_N,
    ['A'] = GT_CDNA_ENC_CHAR_A,['C'] = GT_CDNA_ENC_CHAR_C,['G'] = GT_CDNA_ENC_CHAR_G,['T'] = GT_CDNA_ENC_CHAR_T,
    ['a'] = GT_CDNA_ENC_CHAR_A,['c'] = GT_CDNA_ENC_CHAR_C,['g'] = GT_CDNA_ENC_CHAR_G,['t'] = GT_CDNA_ENC_CHAR_T,
};

/*
 * Block locator functions
 */
#define GT_CDNA_GET_NUM_BLOCKS(num_chars) ((num_chars+(GT_CDNA_BLOCK_CHARS-1))/GT_CDNA_BLOCK_CHARS)
#define GT_CDNA_GET_NUM_CHARS(num_blocks) (GT_CDNA_BLOCK_CHARS*num_blocks)
#define GT_CDNA_GET_BLOCKS_MEM(num_blocks) (GT_CDNA_BLOCK_SIZE*num_blocks)

/*
 * CDNA Internal Bitmap Handlers
 */

#define GT_CDNA_GET_BLOCK_POS(global_pos,block_num,block_pos) \
  block_num = global_pos/GT_CDNA_BLOCK_CHARS; \
  block_pos = global_pos % GT_CDNA_BLOCK_CHARS

#define GT_CDNA_INIT_BLOCK(block_mem) \
  ((uint64_t*)block_mem)[0]=UINT64_ZEROS; \
  ((uint64_t*)block_mem)[1]=UINT64_ZEROS; \
  ((uint64_t*)block_mem)[2]=UINT64_ONES

#define GT_CDNA_GET_BLOCKS(bitmaps_mem,block_num,bm_0,bm_1,bm_2) { \
  register uint64_t* const block_mem = ((uint64_t*)bitmaps_mem)+((block_num)*GT_CDNA_BLOCK_BITMAPS); \
  bm_0 = *block_mem; \
  bm_1 = *block_mem+1; \
  bm_2 = *block_mem+2; \
}

#define GT_CDNA_SHIFT_CHARS(block_pos,bm_0,bm_1,bm_2) \
  bm_0 >>= block_pos; bm_1 >>= block_pos; bm_2 >>= block_pos

#define GT_CDNA_EXTRACT_CHAR(bm_0,bm_1,bm_2) \
  (((bm_0&GT_CDNA_EXTRACT_MASK))    |        \
   ((bm_1&GT_CDNA_EXTRACT_MASK)<<1) |        \
   ((bm_2&GT_CDNA_EXTRACT_MASK)<<2))

#define GT_CDNA_PROYECT_CHAR(block_mem,block_pos,enc_char,bm_pos) \
  if ((enc_char) & gt_cdna_encoded_char_bm[bm_pos]) { \
    register const uint64_t bm_mask = GT_CDNA_ONE_MASK<<(block_pos); \
    block_mem[bm_pos] |= bm_mask; \
  } else { \
    register const uint64_t bm_mask = GT_CDNA_ZERO_MASK<<(block_pos); \
    block_mem[bm_pos] &= bm_mask; \
  }

#define GT_CDNA_SET_CHAR(block_mem,block_pos,enc_char)  \
  GT_CDNA_PROYECT_CHAR(block_mem,block_pos,enc_char,0); \
  GT_CDNA_PROYECT_CHAR(block_mem,block_pos,enc_char,1); \
  GT_CDNA_PROYECT_CHAR(block_mem,block_pos,enc_char,2)

/*
 * Constructor
 */
GT_INLINE gt_compact_dna_string* gt_cdna_string_new(const uint64_t initial_chars) {
  gt_compact_dna_string* cdna_string = malloc(sizeof(gt_compact_dna_string));
  gt_cond_fatal_error(!cdna_string,MEM_HANDLER);
  register const uint64_t initial_blocks = GT_CDNA_GET_NUM_BLOCKS(initial_chars);
  cdna_string->bitmaps = malloc(GT_CDNA_GET_BLOCKS_MEM(initial_blocks));
  gt_cond_fatal_error(!cdna_string->bitmaps,MEM_ALLOC);
  cdna_string->allocated = GT_CDNA_GET_NUM_CHARS(initial_blocks);
  cdna_string->length = 0;
  GT_CDNA_INIT_BLOCK(cdna_string->bitmaps); // Init 0-block // FIXME: Block not zero
  return cdna_string;
}
GT_INLINE void gt_cdna_string_resize(gt_compact_dna_string* const cdna_string,const uint64_t num_chars) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  if (num_chars > cdna_string->allocated) {
    register const uint64_t num_blocks = GT_CDNA_GET_NUM_BLOCKS(num_chars);
    realloc(cdna_string->bitmaps,GT_CDNA_GET_BLOCKS_MEM(num_blocks));
    gt_cond_fatal_error(!cdna_string->bitmaps,MEM_REALLOC);
    cdna_string->allocated = GT_CDNA_GET_NUM_CHARS(num_blocks);
  }
}
GT_INLINE void gt_cdna_string_clear(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  cdna_string->length = 0;
  GT_CDNA_INIT_BLOCK(cdna_string->bitmaps); // Init 0-block // FIXME: Block not zero
}
GT_INLINE void gt_cdna_string_delete(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  free(cdna_string->bitmaps);
  free(cdna_string);
}

/*
 * Handlers
 */
GT_INLINE char gt_cdna_string_get_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  if (pos >= cdna_string->length) return GT_DNA_CHAR_N;
  register uint64_t block_num, block_pos;
  register uint64_t bm_0, bm_1, bm_2;
  GT_CDNA_GET_BLOCK_POS(pos,block_num,block_pos);
  GT_CDNA_GET_BLOCKS(cdna_string->bitmaps,block_num,bm_0,bm_1,bm_2);
  GT_CDNA_SHIFT_CHARS(block_pos,bm_0,bm_1,bm_2);
  return cdna_decode[GT_CDNA_EXTRACT_CHAR(bm_0,bm_1,bm_2)];
}
GT_INLINE void gt_cdna_allocate__init_blocks(gt_compact_dna_string* const cdna_string,const uint64_t pos) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  if (pos >= cdna_string->length) {
    // Check allocated blocks
    gt_cdna_string_resize(cdna_string,pos+1);
    // Initialize remaining new blocks
    register const uint64_t next_block_num = (gt_expect_true(cdna_string->length>0) ?
        (cdna_string->length-1)/GT_CDNA_BLOCK_CHARS : 0) + 1;
    register const uint64_t top_block_num = pos/GT_CDNA_BLOCK_CHARS;
    if (next_block_num <= top_block_num) {
      register uint64_t i;
      register uint64_t* block_mem = cdna_string->bitmaps+(next_block_num*GT_CDNA_BLOCK_BITMAPS);
      for (i=next_block_num; i<=top_block_num; ++i) {
        GT_CDNA_INIT_BLOCK(block_mem);
        block_mem+=GT_CDNA_BLOCK_BITMAPS;
      }
    }
  }
}
GT_INLINE void gt_cdna_string_set_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos,const char character) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  // Check allocated bitmaps
  gt_cdna_allocate__init_blocks(cdna_string,pos);
  // Encode char
  register uint64_t block_num, block_pos;
  GT_CDNA_GET_BLOCK_POS(pos,block_num,block_pos);
  register uint64_t* const block_mem = (cdna_string->bitmaps)+(block_num*GT_CDNA_BLOCK_BITMAPS);
  register const uint8_t enc_char = gt_cdna_encode[character];
  GT_CDNA_SET_CHAR(block_mem,block_pos,enc_char);
}
GT_INLINE uint64_t gt_cdna_string_get_length(gt_compact_dna_string* const cdna_string) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  return cdna_string->length;
}

GT_INLINE void gt_cdna_string_append_string(gt_compact_dna_string* const cdna_string,char* const string,const uint64_t length) {
  GT_COMPACT_DNA_STRING_CHECK(cdna_string);
  // Check allocated bitmaps
  gt_cdna_allocate__init_blocks(cdna_string,pos);
  // TODO
}

/*
 * Compact DNA String Sequence Iterator
 */
GT_INLINE void gt_cdna_new_iterator(
    gt_compact_dna_string* const cdna_string,const uint64_t pos,gt_strand const strand,
    gt_compact_dna_sequence_iterator* const cdna_sequence_iterator) {
  // TODO
}
GT_INLINE void gt_cdna_iterator_seek(gt_compact_dna_sequence_iterator* const cdna_sequence_iterator,const uint64_t pos) {
  // TODO
}
GT_INLINE bool gt_cdna_iterator_eos(gt_compact_dna_sequence_iterator* const cdna_sequence_iterator) {
  // TODO
}
GT_INLINE char gt_cdna_iterator_next(gt_compact_dna_sequence_iterator* const cdna_sequence_iterator) {
  // TODO
}
GT_INLINE char gt_cdna_iterator_previous(gt_compact_dna_sequence_iterator* const cdna_sequence_iterator) {
  // TODO
}
