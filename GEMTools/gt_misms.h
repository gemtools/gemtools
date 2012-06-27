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
 *  - SPLICE: Splicing event (following the splice site consensus) -> (position,size)
 *  - SKIP: More general approach to describe INS/DEL/SPLICE/...
 *      It's currently supported but its use is not recommended.
 *      If you want to model INS/DEL don't use SKIPS
 *      MAP Parsing modules (gt_buffered_map_input) try to use INS/DEL instead of SKIP
 */
typedef enum { MISMS, INS, DEL, SPLICE, SKIP } gt_misms_t;
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
GT_INLINE void gt_misms_set_mismatch(gt_misms misms,const uint64_t position,const char base);
GT_INLINE void gt_misms_set_indel(gt_misms misms,const uint64_t position,const int64_t size);
GT_INLINE void gt_misms_set_splice(gt_misms misms,const uint64_t position,const int64_t size);

/*
 * Accessors
 */
GT_INLINE gt_misms_t gt_misms_get_type(const gt_misms misms);
GT_INLINE uint64_t gt_misms_get_position(const gt_misms misms);
// Mismatches
GT_INLINE char gt_misms_get_base(const gt_misms misms);
// Insertion/Deletion
GT_INLINE uint64_t gt_misms_get_size(const gt_misms misms);

#endif /* GT_MISMS_H_ */
