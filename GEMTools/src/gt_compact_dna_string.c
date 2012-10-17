/*
 * PROJECT: GEM-Tools library
 * FILE: ggt_compact_dna_string.c
 * DATE: 20/08/2012
 * DESCRIPTION: // TODO
 */

#include "gt_compact_dna_string.h"

/*
 * Constructor
 */
GT_INLINE gt_compact_dna_string* gt_cdna_string_new(const uint64_t initial_buffer_size) {
  // TODO
}
GT_INLINE void gt_cdna_string_resize(gt_compact_dna_string* const cdna_string) {
  // TODO
}
GT_INLINE void gt_cdna_string_clear(gt_compact_dna_string* const cdna_string) {
  // TODO
}
GT_INLINE void gt_cdna_string_delete(gt_compact_dna_string* const cdna_string) {
  // TODO
}

/*
 * Handlers
 */
GT_INLINE char gt_cdna_string_get_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos) {
  // TODO
}
GT_INLINE void gt_cdna_string_set_char_at(gt_compact_dna_string* const cdna_string,const uint64_t pos,const char character) {
  // TODO
}

/*
 * Compact DNA String Sequence Iterator
 */
GT_INLINE void gt_cdna_new_iterator(
    gt_compact_dna_string* const cdna_string,const uint64_t pos,,gt_strand const strand,
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
