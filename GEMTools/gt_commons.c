/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_commons.h"

bool gt_dna[256] =
{
    [0 ... 255] = false,
    ['A'] = true,['C'] = true,['G'] = true,['T'] = true,
    ['a'] = true,['c'] = true,['g'] = true,['t'] = true,
    ['N'] = true
};
char gt_complement_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C',  ['T'] = 'A', ['N'] = 'N'
};

GT_INLINE void gt_reverse_complent(char* const sequence,const uint64_t length) {
  GT_NULL_CHECK(sequence);
  register const uint64_t middle = length/2;
  register uint64_t i;
  for (i=0;i<middle;++i) {
    register const char aux = sequence[i];
    sequence[i] = gt_get_complement(sequence[length-i-1]);
    sequence[length-i-1] = gt_get_complement(aux);
  }
  if (length%2) {
    sequence[middle] = gt_get_complement(sequence[middle]);
  }
}
GT_INLINE void gt_cpy_reverse_complent(char* const sequence_dst,char* const sequence_ori,const uint64_t length) {
  GT_NULL_CHECK(sequence_dst); GT_NULL_CHECK(sequence_ori);
  register uint64_t i;
  for (i=0;i<length;++i) {
    sequence_dst[i] = gt_get_complement(sequence_ori[i]);
  }
}
GT_INLINE void gt_reverse(char* const sequence,const uint64_t length) {
  GT_NULL_CHECK(sequence);
  register const uint64_t middle = length/2;
  register uint64_t i;
  for (i=0;i<middle;++i) {
    register const char aux = sequence[i];
    sequence[i] = sequence[length-i-1];
    sequence[length-i-1] = aux;
  }
  if (length%2) {
    sequence[middle] = sequence[middle];
  }
}
GT_INLINE void gt_cpy_reverse(char* const sequence_dst,char* const sequence_ori,const uint64_t length) {
  GT_NULL_CHECK(sequence_dst); GT_NULL_CHECK(sequence_ori);
  register uint64_t i;
  for (i=0;i<length;++i) {
    sequence_dst[i] = sequence_ori[i];
  }
}

GT_INLINE char* gt_string_cpy(char* sequence,const uint64_t length) {
  GT_NULL_CHECK(sequence);
  register char* sequence_cpy = malloc(length+1);
  register uint64_t i;
  for (i=0;i<length;++i) {
    sequence_cpy[i] = sequence[i];
  }
  sequence_cpy[length] =EOS;
  return sequence_cpy;
}
