/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */


#include "gt_map.h"

#define GT_MAP_NUM_INITIAL_MISMS 4

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new() {
  gt_map* map = malloc(sizeof(gt_map));
  gt_cond_fatal_error(!map,MEM_HANDLER);
  map->seq_name = NULL;
  map->position = 0;
  map->length = MISMATCH_STRING_UNKNOWN;
  map->distance = UINT64_MAX;
  map->score = 0.0;
  map->mismatches = gt_vector_new(GT_MAP_NUM_INITIAL_MISMS,sizeof(gt_misms));
  map->mismatches_txt = NULL;
  map->next_map = NULL;
  return map;
}
GT_INLINE void gt_map_clear(gt_map* const map) {
  GT_MAP_CHECK(map);
  map->seq_name = NULL;
  map->position = 0;
  map->length = MISMATCH_STRING_UNKNOWN;
  map->distance = UINT64_MAX;
  map->score = 0.0;
  gt_map_clear_misms(map);
  map->mismatches_txt = NULL;
  map->next_map = NULL;
}
GT_INLINE void gt_map_delete(gt_map* map) {
  GT_MAP_CHECK(map);
  gt_vector_delete(map->mismatches);
  free(map);
}

/*
 * Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->seq_name;
}
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* seq_name) {
  GT_MAP_CHECK(map);
  gt_check(seq_name==NULL,NULL_HANDLER);
  map->seq_name = seq_name;
}
GT_INLINE uint64_t gt_map_get_position(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->position;
}
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position) {
  GT_MAP_CHECK(map);
  map->position = position;
}
GT_INLINE gt_strand gt_map_get_direction(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->direction;
}
GT_INLINE void gt_map_set_direction(gt_map* const map,const gt_strand strand) {
  GT_MAP_CHECK(map);
  map->direction = strand;
}
/*
 * Map.distance and Map.score are hierarchical attributes.
 * Each map block contains the distance/score for all maps
 * from that point and on (including itself)
 *   The leading map block contains the global distance/score
 */
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->distance;
}
GT_INLINE void gt_map_set_distance(gt_map* const map,const uint64_t distance) {
  GT_MAP_CHECK(map);
  map->distance = distance;
}
GT_INLINE float gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const float score) {
  GT_MAP_CHECK(map);
  map->score = score;
}
GT_INLINE gt_map* gt_map_get_next_map_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (gt_map*) map->next_map;
}
GT_INLINE void gt_map_set_next_map_block(gt_map* const map,gt_map* next_map) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  map->next_map = (struct gt_map*) next_map;
}
/*
 * Map.length is a global attribute sum of all Map.blocks.length
 */
GT_INLINE uint64_t gt_map_get_map_block_length(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  return map->length;
}
GT_INLINE void gt_map_set_map_block_length(gt_map* const map,const uint64_t length) {
  GT_MAP_EDITABLE_CHECK(map);
  map->length = length;
}
GT_INLINE uint64_t gt_map_get_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  // Single Block
  if (gt_expect_true(map->next_map==NULL)) return gt_map_get_map_block_length(map);
  // Multiple Block
  register uint64_t total_length = 0;
  register gt_map* aux_map = map;
  do {
    total_length += gt_map_get_map_block_length(aux_map);
    aux_map = (gt_map*) aux_map->next_map;
  } while (aux_map!=NULL);
  return total_length;
}

/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms) {
  GT_MAP_EDITABLE_CHECK(map);
  gt_vector_insert(map->mismatches,misms,gt_misms*);
}
GT_INLINE void gt_map_clear_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_vector_clean(map->mismatches);
}
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,uint64_t offset) {
  GT_MAP_EDITABLE_CHECK(map);
  return gt_vector_get_elm(map->mismatches,offset,gt_misms);
}
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map) {
  GT_MAP_CHECK(map);
  gt_map* map_cpy = gt_map_new();
  map_cpy->seq_name = map->seq_name;
  map_cpy->position = map->position;
  map_cpy->length = map->length;
  map_cpy->direction = map->direction;
  map_cpy->distance = map->distance;
  map_cpy->score = map->score;
  map_cpy->mismatches = gt_vector_new(gt_vector_get_used(map->mismatches),sizeof(gt_misms));
  gt_vector_copy(map_cpy->mismatches,map->mismatches);
  map_cpy->mismatches_txt = map->mismatches_txt;
  map_cpy->next_map = (map->next_map==NULL) ? NULL : (struct gt_map*) gt_map_copy((gt_map*)map->next_map);
  return map_cpy;
}
