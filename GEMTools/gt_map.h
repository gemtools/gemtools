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

// Map format
#define map_misms_format score  /* Overload field */
#define MISMATCH_STRING_UNKNOWN UINT64_MAX-3
#define MISMATCH_STRING_GEMv0 UINT64_MAX-2
#define MISMATCH_STRING_GEMv1 UINT64_MAX-1 // OLD(v0)={chr7:F127708134G27T88} NEW(v2)={chr11:-:51590050:(5)43T46A9>24*}

// Orientation (strand)
typedef enum { FORWARD, REVERSE } gt_strand;
// Types of junctions between map blocks
typedef enum { NO_JUNCTION, SPLICE, POSITIVE_SKIP, NEGATIVE_SKIP, INSERT } gt_junction_t;
// Forward declaration of gt_map
typedef struct _gt_map gt_map;
// Junction
typedef struct {
  gt_map* map;
  gt_junction_t junction;
} gt_map_junction;
/*
 * Map type
 *   NOTE: Use gt_map (not _gt_map)
 */
struct _gt_map {
  /* Sequence-name(Chromosome/Contig/...), position and strand */
  char* seq_name;
  uint64_t position;
  uint64_t base_length; // Length not including indels
  gt_strand direction;
  /* Metrics */
  uint64_t score;
  /* Mismatches/Indels/... */
  gt_vector* mismatches; /* (misms_t) */
  char* mismatches_txt;
  /* Multiple Block Map (splice-map, local alignments, ...) */
  gt_map_junction* next_block;
};

// Iterators
typedef struct {
  gt_map* map;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_misms* next_misms;
} gt_map_mism_iterator;
typedef struct {
  gt_map* map;
  gt_map* next_map;
} gt_map_block_iterator;

// Checkers
#define GT_MAP_CHECK(map) gt_fatal_check(map==NULL||map->mismatches==NULL,NULL_HANDLER)
#define GT_MAP_EDITABLE_CHECK(map) \
  GT_MAP_CHECK(map); \
  gt_fatal_check(map->mismatches_txt!=NULL,MAP_MISMS_NOT_PARSED)
#define GT_MAP_NEXT_BLOCK_CHECK(map) \
  GT_NULL_CHECK(map->next_block); \
  GT_MAP_CHECK(map->next_block->map)

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new(void);
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
// Single block attributes (NOTE: To retrieve length/distance/score over all block use global methods)
GT_INLINE uint64_t gt_map_get_length(gt_map* const map); // Indel taken into account to calculate the length
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map); // Length of the base read (no indels)
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length);

GT_INLINE uint64_t gt_map_get_distance(gt_map* const map);

GT_INLINE uint64_t gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score);

/*
 * Multiple Block Maps Handlers
 */
GT_INLINE bool gt_map_has_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map);
GT_INLINE gt_junction_t gt_map_get_next_block_junction(gt_map* const map);
GT_INLINE int64_t gt_map_get_next_block_distance(gt_map* const map);
GT_INLINE void gt_map_set_next_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction);
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map);
GT_INLINE void gt_map_append_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction);
/*
 * Map blocks iterator
 *  GT_BEGIN_MAP_BLOCKS_ITERATOR(initial_map,map_it) {
 *    ..code(map_it)..
 *  } GT_END_MAP_BLOCKS_ITERATOR;
 */
#define GT_BEGIN_MAP_BLOCKS_ITERATOR(initial_map,map_it) { \
  register gt_map* map_it = initial_map; \
  while (map_it!=NULL) {
#define GT_END_MAP_BLOCKS_ITERATOR \
    map_it = gt_map_get_next_block(map_it); \
  } \
}

/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms);
GT_INLINE void gt_map_clear_misms(gt_map* const map);
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,uint64_t offset);
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map);

/*
 * High-level Procedures
 */
// Global metrics (over all blocks)
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_score(gt_map* const map);
// Vector based ( Metrics out of a set of maps )
GT_INLINE uint64_t gt_map_vector_get_length(gt_vector* const maps);
GT_INLINE uint64_t gt_map_vector_get_distance(gt_vector* const maps);
GT_INLINE uint64_t gt_map_vector_get_score(gt_vector* const maps);
// Distance procedures
GT_INLINE uint64_t gt_map_get_levenshtein_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_levenshtein_distance(gt_map* const map);

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map);
// Macro generic iterator
//  GT_ALIGNMENT_MISMS_ITERATOR(map,misms_it,misms_pos) {
//    ..code..
//  }
#define GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) \
  GT_VECTOR_ITERATE(map->mismatches,misms_it,misms_pos,gt_misms)

// Map's Blocks iterator
GT_INLINE void gt_map_new_block_iterator(gt_map* const map,gt_map_block_iterator* const map_block_iterator);
GT_INLINE gt_map* gt_map_next_block(gt_map_block_iterator* const map_block_iterator);

// Map's Mismatches iterator
GT_INLINE void gt_map_new_misms_iterator(gt_map* const map,gt_map_mism_iterator* const map_mism_iterator);
GT_INLINE gt_misms* gt_map_next_misms(gt_map_mism_iterator* const map_mism_iterator);
GT_INLINE gt_misms* gt_map_dinamic_next_misms(gt_map_mism_iterator* const map_mism_iterator);
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator);


#endif /* GT_MAP_H_ */
