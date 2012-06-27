/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_H_
#define GT_MAP_H_

#include "gt_commons.h"
#include "gt_misms.h"

#define MISMATCH_STRING_UNKNOWN UINT64_MAX-3
#define MISMATCH_STRING_GEMv0 UINT64_MAX-2
#define MISMATCH_STRING_GEMv1 UINT64_MAX-1 // OLD(v0)={chr7:F127708134G27T88} NEW(v2)={chr11:-:51590050:(5)43T46A9>24*}
typedef enum { FORWARD, REVERSE } gt_strand;

#define map_mismatch_string_format distance  /* Overload field */
typedef struct {
  /* Sequence-name(Chromosome), position and strand */
  char *seq_name;
  uint64_t position;
  uint64_t length;
  gt_strand direction;
  /* Metrics */
  uint64_t distance;
  float score;
  /* Mismatches/Indels/... */
  gt_vector *mismatches; /* (misms_t) */
  char* mismatches_txt;
  /* Multimap (splice-map across chromosomes, etc) */
  struct gt_map* next_map;
} gt_map;

// Checkers
#define GT_MAP_CHECK(map) gt_fatal_check(map==NULL||map->mismatches==NULL,NULL_HANDLER)
#define GT_MAP_EDITABLE_CHECK(map) \
  GT_MAP_CHECK(map); \
  gt_fatal_check(map->mismatches_txt!=NULL,MAP_MISMS_NOT_PARSED)

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new();
GT_INLINE void gt_map_clear(gt_map* const map);
GT_INLINE void gt_map_delete(gt_map* map);

/*
 * Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map);
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* seq_name);

GT_INLINE uint64_t gt_map_get_position(gt_map* const map);
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position);

GT_INLINE gt_strand gt_map_get_direction(gt_map* const map);
GT_INLINE void gt_map_set_direction(gt_map* const map,const gt_strand strand);

GT_INLINE uint64_t gt_map_get_distance(gt_map* const map);
GT_INLINE void gt_map_set_distance(gt_map* const map,const uint64_t distance);

GT_INLINE float gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const float score);

/*
 * Nested Maps Handlers
 */
GT_INLINE gt_map* gt_map_get_next_map_block(gt_map* const map);
GT_INLINE void gt_map_set_next_map_block(gt_map* const map,gt_map* const next_map);

GT_INLINE uint64_t gt_map_get_map_block_length(gt_map* const map);
GT_INLINE void gt_map_set_map_block_length(gt_map* const map,const uint64_t length);
GT_INLINE uint64_t gt_map_get_length(gt_map* const map);

/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms);
GT_INLINE void gt_map_clear_misms(gt_map* const map);
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,uint64_t offset);
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map);

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map);
// Macro generic iterator
//  GT_ALIGNMENT_MISMS_ITERATOR(map,misms_it,misms_pos) {
//    ..code..
//  }
#define GT_ALIGNMENT_MISMS_ITERATOR(map,misms_it,misms_pos) \
  GT_VECTOR_ITERATE(map->mismatches,misms_it,misms_pos,gt_misms)


#endif /* GT_MAP_H_ */
