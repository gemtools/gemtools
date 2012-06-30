/*
 * PROJECT: GEM-Tools library
 * FILE: gt_misms.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_MISMS_H_
#define GT_MISMS_H_

#include "gt_commons.h"

/*
 * Type of Mismatches
 *  - MISMS: Substitution -> (position,base)
 *  - INS/DEL: Insertion/Deletion -> (position,size)
 */
typedef enum { MISMS, INS, DEL } gt_misms_t;
typedef struct {
  gt_misms_t misms_type;
  uint64_t position;
  union {
    int64_t size;
    char base;
  };
} gt_misms;

/*
 * Constructors
 */
GT_INLINE void gt_misms_set_mismatch(gt_misms* const misms,const uint64_t position,const char base);
GT_INLINE void gt_misms_set_indel(gt_misms* const misms,const uint64_t position,const int64_t size);

/*
 * Accessors
 */
GT_INLINE gt_misms_t gt_misms_get_type(gt_misms* const misms);
GT_INLINE uint64_t gt_misms_get_position(gt_misms* const misms);
// Mismatches
GT_INLINE char gt_misms_get_base(gt_misms* const misms);
// Insertion/Deletion
GT_INLINE uint64_t gt_misms_get_size(gt_misms* const misms);

#endif /* GT_MISMS_H_ */
