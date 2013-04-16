/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_H_
#define GT_MAP_H_

#include "gt_essentials.h"
#include "gt_data_attributes.h"

#include "gt_misms.h"
#include "gt_dna_string.h"

/*
 * Constants
 */
#define GT_MAP_NO_GT_SCORE UINT64_MAX
#define GT_MAP_NO_PHRED_SCORE 255

/*
 * Junction (gt_map_junction)
 */
// Types of junctions between map blocks
typedef enum { NO_JUNCTION, SPLICE, POSITIVE_SKIP, NEGATIVE_SKIP, INSERT, JUNCTION_UNKNOWN } gt_junction_t;
typedef struct _gt_map gt_map; // Forward declaration of gt_map
typedef struct {
  gt_map* map;
  gt_junction_t junction;
  int64_t junction_size;
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
  int64_t gt_score;
  uint8_t phred_score;
  /* Mismatches/Indels */
  gt_vector* mismatches; /* (misms_t) */
  /* Multiple Block Map (splice-map, local alignments, ...) */
  gt_map_junction* next_block;
  /* Attributes */
  gt_shash* attributes;
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

/*
 * Checkers
 */
#define GT_MAP_CHECK(map) \
  GT_NULL_CHECK(map); \
  GT_NULL_CHECK((map)->mismatches)
#define GT_MAP_NEXT_BLOCK_CHECK(map) \
  GT_NULL_CHECK(map->next_block); \
  GT_MAP_CHECK(map->next_block->map)
/* TODO: Scheduled for v2.0 (all lazy parsing)
 * #define GT_MAP_EDITABLE_CHECK(map) \
 *   GT_MAP_CHECK(map); \
 *  gt_fatal_check(map->mismatches_txt!=NULL,MAP_MISMS_NOT_PARSED)
 */

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new(void);
GT_INLINE void gt_map_clear(gt_map* const map);
GT_INLINE void gt_map_delete(gt_map* const map);

/*
 * Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map);
GT_INLINE uint64_t gt_map_get_seq_name_length(gt_map* const map);
GT_INLINE gt_string* gt_map_get_string_seq_name(gt_map* const map);
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* const seq_name,const uint64_t length);
GT_INLINE void gt_map_set_string_seq_name(gt_map* const map,gt_string* const seq_name);

GT_INLINE uint64_t gt_map_get_global_position(gt_map* const map);
GT_INLINE uint64_t gt_map_get_position(gt_map* const map);
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position);
GT_INLINE gt_strand gt_map_get_strand(gt_map* const map);
GT_INLINE void gt_map_set_strand(gt_map* const map,const gt_strand strand);
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map); // Length of the base read (no indels)
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length);

// Attributes (lazy allocation of the attributes field)
GT_INLINE void* gt_map_attribute_get(gt_map* const map,char* const attribute_id);
GT_INLINE void gt_map_attribute_set_(gt_map* const map,char* const attribute_id,void* const attribute,const size_t element_size);
#define gt_map_attribute_set(map,attribute_id,attribute,element_type) \
    gt_map_attribute_set_(map,attribute_id,(void*)attribute,sizeof(element_type))

// Local metrics (over first blocks)
GT_INLINE uint64_t gt_map_get_length(gt_map* const map); // Indel taken into account to calculate the length
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map);
GT_INLINE int64_t gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const int64_t score);
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map);
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score);

/*
 * Multiple Block Maps Handlers
 */
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map);
GT_INLINE bool gt_map_has_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map);
GT_INLINE void gt_map_set_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);
GT_INLINE void gt_map_insert_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);
GT_INLINE void gt_map_append_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);

GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map);
GT_INLINE int64_t gt_map_get_junction_size(gt_map* const map);

GT_INLINE void gt_map_reverse_blocks_positions(gt_map* const head_map,const uint64_t start_position);
GT_INLINE void gt_map_reverse_misms(gt_map* const map);

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
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,const uint64_t offset);
GT_INLINE void gt_map_set_misms(gt_map* const map,gt_misms* misms,const uint64_t offset);
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map);
GT_INLINE void gt_map_set_num_misms(gt_map* const map,const uint64_t num_misms);

// Map check/recover operators
#define GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms_ptr,num_misms) { \
  if (misms_offset<num_misms) { \
    misms_ptr=gt_map_get_misms(map,misms_offset); \
  } else { \
    misms_ptr=NULL; \
  } \
}
#define GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_misms) { \
  if (misms_offset<num_misms) { \
    misms=*gt_map_get_misms(map,misms_offset); \
  } else { \
    misms.position=UINT64_MAX; \
  } \
}

// Counting
GT_INLINE uint64_t gt_map_get_num_mismatch_bases(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_indels(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_insertions(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_deletions(gt_map* const map);

/*
 * High-level Procedures
 */
// Global metrics (over all blocks)
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map);
// Begin/End Position
GT_INLINE uint64_t gt_map_get_begin_position(gt_map* const map);
GT_INLINE uint64_t gt_map_get_end_position(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_end_position(gt_map* const map);
// Trim helpers
GT_INLINE uint64_t gt_map_get_left_trim_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_right_trim_length(gt_map* const map);
// Bases aligned
GT_INLINE uint64_t gt_map_get_bases_aligned(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_bases_aligned(gt_map* const map);
// Distance
GT_INLINE uint64_t gt_map_get_levenshtein_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_levenshtein_distance(gt_map* const map);
// Vector based ( Metrics out of a set of maps )
GT_INLINE uint64_t gt_map_vector_get_length(gt_vector* const maps);
GT_INLINE uint64_t gt_map_vector_get_distance(gt_vector* const maps);
// Map compare functions
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2);
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated);
GT_INLINE int64_t gt_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps);
GT_INLINE int64_t gt_mmap_range_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps,const uint64_t range_tolerated);
GT_INLINE bool gt_map_less_than(gt_map* const map_1,gt_map* const map_2); // As to resolve ties

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map);
GT_INLINE gt_map** gt_mmap_array_copy(gt_map** mmap,const uint64_t num_blocks);

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

/*
 * Iterate over the map blocks of a map
 *   Map = (map.block1,map.block2)
 *   GT_MAP_BLOCKS_ITERATE := {map.block1,map.block2}
 */
#define GT_MAP_ITERATE(map,map_block) \
  /* Map. Iterate over map blocks */ \
  gt_map_block_iterator __##map_block##_iterator; \
  register gt_map* map_block; \
  gt_map_new_block_iterator(map,&(__##map_block##_iterator)); \
  while ((map_block=gt_map_next_block(&(__##map_block##_iterator))))

/*
 * Iterate over the mismatches(M/I/D) of a map
 */
#define GT_MISMS_ITERATE(map_block,misms) \
  /* Map Block. Iterate over all mismatches */ \
  gt_map_mism_iterator __##misms##_iterator; \
  register gt_misms* misms; \
  gt_map_new_misms_iterator(map_block,&(__##misms##_iterator)); \
  while ((misms=gt_map_next_misms(&(__##misms##_iterator))))

#endif /* GT_MAP_H_ */
