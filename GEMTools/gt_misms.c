/*
 * PROJECT: GEM-Tools library
 * FILE: gt_misms.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_misms.h"

/*
 * Constructors
 */
GT_INLINE void gt_misms_set_mismatch(gt_misms* const misms,const uint64_t position,const char base) {
  misms->misms_type = MISMS;
  misms->position = position;
  misms->base = base;
}
GT_INLINE void gt_misms_set_indel(gt_misms* const misms,const uint64_t position,const int64_t size) {
  misms->misms_type = (size>0) ? INS : DEL;
  misms->position = position;
  misms->size = (size>0) ? size : -size;
}

/*
 * Accessors
 */
GT_INLINE gt_misms_t gt_misms_get_type(gt_misms* const misms) {
  return misms->misms_type;
}
GT_INLINE uint64_t gt_misms_get_position(gt_misms* const misms) {
  return misms->position;
}
// Mismatches
GT_INLINE char gt_misms_get_base(gt_misms* const misms) {
  gt_check(misms->misms_type!=MISMS,MISMS_TYPE);
  return misms->base;
}
// Insertion/Deletion
GT_INLINE uint64_t gt_misms_get_size(gt_misms* const misms) {
  gt_check(misms->misms_type==MISMS,MISMS_TYPE);
  return misms->size;
}
