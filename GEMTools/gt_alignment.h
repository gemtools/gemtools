/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_ALIGNMENT_H_
#define GT_ALIGNMENT_H_

typedef struct {
  char* tag;
  uint64_t tag_length;
  char* read;
  uint64_t read_length;
  char* qualities;
  gt_vector* counters;
  uint64_t max_complete_strata;
  gt_vector* maps; /* (gt_map) */
} gt_alignment;

#endif /* GT_ALIGNMENT_H_ */
