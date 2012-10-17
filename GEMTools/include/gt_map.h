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

#define GT_MAP_NO_SCORE (-1)

// Orientation (strand)
typedef enum { FORWARD, REVERSE } gt_strand;
// Types of junctions between map blocks
typedef enum { NO_JUNCTION, SPLICE, POSITIVE_SKIP, NEGATIVE_SKIP, INSERT } gt_junction_t;
// Mismatch string format (lazy parsing)
typedef enum { MISMATCH_STRING_GEMv0, MISMATCH_STRING_GEMv1 } gt_misms_string_t;

typedef struct _gt_map gt_map; // Forward declaration of gt_map
/*
 * Junction (gt_map_junction)
 */
typedef struct {
  gt_map* map;
  gt_junction_t junction;
} gt_map_junction;
/*
 * Map (gt_map)
 */
struct _gt_map {
  /* Sequence-name(Chromosome/Contig/...), position and strand */
  gt_string* seq_name;
  uint64_t position;
  uint64_t base_length; // Length not including indels
  gt_strand strand;
  /* Metrics */
  int64_t score;
  /* Mismatches/Indels */
  gt_vector* mismatches; /* (misms_t) */
  char* misms_txt;
  gt_misms_string_t misms_txt_format;
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
#define GT_MAP_CHECK(map) gt_fatal_check((map)==NULL||(map)->mismatches==NULL,NULL_HANDLER)
// TODO: Scheduled for v2.0 (all lazy parsing)
//#define GT_MAP_EDITABLE_CHECK(map) \
//  GT_MAP_CHECK(map); \
//  gt_fatal_check(map->mismatches_txt!=NULL,MAP_MISMS_NOT_PARSED)
#define GT_MAP_NEXT_BLOCK_CHECK(map) \
  GT_NULL_CHECK(map->next_block); \
  GT_MAP_CHECK(map->next_block->map)

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new(void);
GT_INLINE gt_map* gt_map_new_(const bool static_seq_name);
GT_INLINE void gt_map_clear(gt_map* const map);
GT_INLINE void gt_map_delete(gt_map* const map);

/*
 * Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map);
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* const seq_name,const uint64_t length);

GT_INLINE uint64_t gt_map_get_position(gt_map* const map);
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position);

GT_INLINE gt_strand gt_map_get_strand(gt_map* const map);
GT_INLINE void gt_map_set_strand(gt_map* const map,const gt_strand strand);
// Single block attributes (NOTE: To retrieve length/distance/score over all block use global methods)
GT_INLINE uint64_t gt_map_get_length(gt_map* const map); // Indel taken into account to calculate the length
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map); // Length of the base read (no indels)
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length);

GT_INLINE uint64_t gt_map_get_distance(gt_map* const map);

GT_INLINE int64_t gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const int64_t score);

/*
 * Multiple Block Maps Handlers
 */
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map);
GT_INLINE bool gt_map_has_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map);
GT_INLINE void gt_map_set_next_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction);
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map);
GT_INLINE void gt_map_append_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction);

GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map);
GT_INLINE uint64_t gt_map_get_junction_distance(gt_map* const map);
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

GT_INLINE void gt_map_clear_misms_string(gt_map* const map);
GT_INLINE char* gt_map_get_misms_string(gt_map* const map);
GT_INLINE gt_misms_string_t gt_map_get_misms_string_format(gt_map* const map);
GT_INLINE void gt_map_set_misms_string(gt_map* const map,char* misms_string,const gt_misms_string_t format);

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
// Map compare
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2);
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated);
GT_INLINE int64_t gt_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps);
GT_INLINE int64_t gt_mmap_range_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps,const uint64_t range_tolerated);

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map);
// Macro generic iterator
//  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
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
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator);


#endif /* GT_MAP_H_ */
