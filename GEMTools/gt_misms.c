/*
 * PROJECT: GEM-Tools library
 * FILE: gt_misms.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_misms.h"

GT_INLINE void gt_misms_set_mismatch(gt_misms misms,const uint64_t position,const char base) {
  misms.misms_type = MISMS;
  misms.position = position;
  misms.base = base;
}
GT_INLINE void gt_misms_set_indel(gt_misms misms,const uint64_t position,const int64_t size) {
  misms.misms_type = (size>0) ? INS : DEL;
  misms.position = position;
  misms.size = size;
}
GT_INLINE gt_misms_t gt_misms_get_type(const gt_misms misms) {
  return misms.misms_type;
}
GT_INLINE uint64_t gt_misms_get_position(const gt_misms misms) {
  return misms.position;
}
GT_INLINE uint64_t gt_misms_get_size(const gt_misms misms) {
  gt_check(misms.misms_type==MISMS,MISMS_TYPE);
  return misms.size;
}
GT_INLINE char gt_misms_get_base(const gt_misms misms) {
  gt_check(misms.misms_type!=MISMS,MISMS_TYPE);
  return misms.base;
}
