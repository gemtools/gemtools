/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_score.h
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_SCORE_H_
#define GT_MAP_SCORE_H_

#include "gt_essentials.h"

typedef enum { GT_XT_UNKNOWN, GT_XT_UNIQUE, GT_XT_REPEAT, GT_XT_UNMAPPED, GT_XT_MATE_SW } gt_xt_value;

#include "gt_alignment_utils.h"
#include "gt_template_utils.h"
#include "gt_dna_read.h"

/*
 * Parameters required for the scoring
 */
typedef struct {
	int64_t minimum_insert;
	int64_t maximum_insert;
	uint8_t *insert_phred; // Size of array: maximum_insert-minimum_insert+1
	gt_qualities_offset_t quality_format;
  uint32_t mapping_cutoff;
  uint32_t max_strata_searched;
  uint64_t max_pair_maps; // Max pair maps to store;
  int indel_penalty;
  int split_penalty;
} gt_map_score_attributes;

#define GT_MAP_SCORE_MAX_QUALITY 50
#define GT_MAP_SCORE_MISSING_QUALITY 30 // Value to use in calculations if qualities not available
#define GT_MAP_SCORE_MAX_GT_SCORE 0xffff
#define GT_MAP_SCORE_MAX_PAIR_MAPS 128
#define GT_MAP_SCORE_ATTRIBUTES_CHECK(ms_attr) GT_NULL_CHECK(ms_attr);

/*
 * Map Score Accessors
 */
GT_INLINE uint64_t gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score);
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map);
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score);

uint64_t gt_map_calculate_gt_score(gt_alignment *al, gt_map *map, gt_map_score_attributes *ms_attr);
void gt_map_pair_template(gt_template *template,gt_map_score_attributes *ms_attr);
void gt_map_calculate_template_mapq_score(gt_template *template,gt_map_score_attributes *ms_attr);
void gt_map_calculate_alignment_mapq_score(gt_alignment *alignment,gt_map_score_attributes *ms_attr);


/*
 * TODO
 * fx(MAP) Scoring functions
 * fx(ALIGNMENT,MAP) Scoring functions
 * fx(TEMPLATE) Scoring functions
 */

#endif /* GT_MAP_SCORE_H_ */
