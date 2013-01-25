/*
 * PROJECT: GEM-Tools library
 * FILE: gt_dna_string.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_dna_string.h"

const bool gt_dna[256] =
{
    [0 ... 255] = false,
    ['A'] = true,['C'] = true,['G'] = true,['T'] = true,
    /*['a'] = true,['c'] = true,['g'] = true,['t'] = true, *//* TODO: a,c,g,t,n compatibility via gt_dna_string */
    ['N'] = true
};
const char gt_dna_normalized[256] =
{
    [0 ... 255] = 'N',
    ['A'] = 'A',['C'] = 'C',['G'] = 'G',['T'] = 'T',
    ['a'] = 'A',['c'] = 'C',['g'] = 'G',['t'] = 'T',
};
const char gt_complement_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C',  ['T'] = 'A', ['N'] = 'N'
};
const bool gt_iupac_code[256] =
{
  [0 ... 255] = false,
  /* Upper case */
  ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true,
  ['N'] = true, ['R'] = true, ['N'] = true, ['D'] = true,
  ['E'] = true, ['H'] = true, ['I'] = true, ['L'] = true,
  ['K'] = true, ['M'] = true, ['F'] = true, ['P'] = true,
  ['S'] = true, ['W'] = true, ['Y'] = true, ['V'] = true, ['Q'] = true,
  ['B'] = true, ['Z'] = true, ['X'] = true, ['U'] = true, ['R'] = true,
  /* Lower case*/
  ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true,
  ['n'] = true, ['r'] = true, ['n'] = true, ['d'] = true,
  ['e'] = true, ['h'] = true, ['i'] = true, ['l'] = true,
  ['k'] = true, ['m'] = true, ['f'] = true, ['p'] = true,
  ['s'] = true, ['w'] = true, ['y'] = true, ['v'] = true, ['q'] = true,
  ['b'] = true, ['z'] = true, ['x'] = true, ['u'] = true, ['r'] = true
};


/*
 * DNA String handler
 */
GT_INLINE bool gt_dna_string_is_dna_string(gt_dna_string* const dna_string) {
  GT_STRING_CHECK(dna_string);
  register const uint64_t length = dna_string->length;
  register const char* const buffer = dna_string->buffer;
  register uint64_t i;
  for (i=0;i<length;++i) {
    if (gt_expect_false(!gt_is_dna(buffer[i]))) return false;
  }
  return true;
}

GT_INLINE char gt_dna_string_get_char_at(gt_dna_string* const dna_string,const uint64_t pos) {
  // TODO
  return 'a';
}
GT_INLINE void gt_dna_string_set_char_at(gt_dna_string* const dna_string,const uint64_t pos,const char character) {
  // TODO
}

GT_INLINE void gt_dna_string_reverse_complement(gt_dna_string* const dna_string) {
  GT_STRING_CHECK(dna_string);
  register const uint64_t length = dna_string->length;
  register const uint64_t middle = length/2;
  register char* const buffer = dna_string->buffer;
  register uint64_t i;
  for (i=0;i<middle;++i) {
    register const char aux = buffer[i];
    buffer[i] = gt_get_complement(buffer[length-i-1]);
    buffer[length-i-1] = gt_get_complement(aux);
  }
  if (length%2==1) {
    buffer[middle] = gt_get_complement(buffer[middle]);
  }
}
GT_INLINE void gt_dna_string_reverse_complement_copy(gt_dna_string* const dna_string_dst,gt_dna_string* const dna_string_src) {
  GT_STRING_CHECK(dna_string_dst);
  GT_STRING_CHECK(dna_string_src);
  register const uint64_t length = dna_string_src->length;
  gt_string_resize(dna_string_dst,length+1);
  register char* const buffer_src = dna_string_src->buffer;
  register char* const buffer_dst = dna_string_dst->buffer;
  register uint64_t i;
  for (i=0;i<length;++i) {
    buffer_dst[i] = gt_get_complement(buffer_src[length-1-i]);
  }
  buffer_dst[length] = EOS;
  dna_string_dst->length = length;
}

/*
 * DNA String Iterator
 */
GT_INLINE void gt_dna_string_new_iterator(
    gt_dna_string* const dna_string,const uint64_t pos,gt_string_traversal const direction,
    gt_dna_string_iterator* const dna_string_iterator) {
  // TODO
}
GT_INLINE void gt_dna_string_iterator_seek(gt_dna_string_iterator* const dna_string_iterator,const uint64_t pos) {
  // TODO
}
GT_INLINE bool gt_dna_string_iterator_eos(gt_dna_string_iterator* const dna_string_iterator) {
  // TODO
  return true;
}
GT_INLINE char gt_dna_string_iterator_next(gt_dna_string_iterator* const dna_string_iterator) {
  // TODO
  return 'a';
}
GT_INLINE char gt_dna_string_iterator_previous(gt_dna_string_iterator* const dna_string_iterator) {
  // TODO
  return 'a';
}
