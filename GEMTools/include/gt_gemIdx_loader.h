/*
 * PROJECT: GEM-Tools library
 * FILE: gt_gemIdx_loader.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_GEM_INDEX_LOADER_H_
#define GT_GEM_INDEX_LOADER_H_

#include "gt_commons.h"
#include "gt_compact_dna_string.h"
#include "gt_sequence_archive.h"

/*
 * Auxiliary Data Structures (as to read the idx)
 */
typedef struct {
  uint64_t bot;
  uint64_t top;
  int64_t sequence_offset;
  int64_t tag_offset;
} gem_loc_t;

/*
 * Setup
 */
GT_INLINE void gt_gemIdx_load_archive(
    char* const index_file_name,gt_sequence_archive* const sequence_archive,const bool load_sequences);

/*
 * Retrieve sequences from GEMindex
 */
GT_INLINE int64_t gt_sequence_archive_get_bed_sequence_string(
  gt_sequence_archive* const sequence_archive,char* const seq_id,
  const uint64_t position,const uint64_t length,gt_string* const string);

#endif /* GT_GEM_INDEX_LOADER_H_ */
