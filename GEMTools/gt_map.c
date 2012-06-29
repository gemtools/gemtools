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
  map->base_length = UINT32_MAX;
  map->map_misms_format = MISMATCH_STRING_UNKNOWN;
  map->mismatches = gt_vector_new(GT_MAP_NUM_INITIAL_MISMS,sizeof(gt_misms));
  map->mismatches_txt = NULL;
  map->next_map = NULL;
  return map;
}
GT_INLINE void gt_map_clear(gt_map* const map) {
  GT_MAP_CHECK(map);
  map->seq_name = NULL;
  map->position = 0;
  map->base_length = UINT32_MAX;
  map->map_misms_format = MISMATCH_STRING_UNKNOWN;
  gt_map_clear_misms(map);
  map->mismatches_txt = NULL;
  map->next_map = NULL;
}
GT_INLINE void gt_map_block_delete(gt_map* map) {
  GT_MAP_CHECK(map);
  gt_vector_delete(map->mismatches);
  free(map);
}
GT_INLINE void gt_map_delete(gt_map* map) {
  GT_MAP_CHECK(map);
  if (map->next_map != NULL) {
    gt_map_delete(map->next_map->map);
    free(map->next_map);
  }
  gt_map_block_delete(map);
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
GT_INLINE uint64_t gt_map_get_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t length = map->base_length; // FIXME: What if neg length
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
        length += gt_misms_get_size(misms_it);
      case DEL:
        length -= gt_misms_get_size(misms_it);
        break;
    }
  }
  return length;
}
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->base_length;
}
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  map->base_length = length;
}
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}
GT_INLINE uint64_t gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score) {
  GT_MAP_CHECK(map);
  map->score = score;
}

/*
 * Multiple Block Maps Handlers
 */
GT_INLINE bool gt_map_has_next_block(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  return (map->next_map!=NULL);
}
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  return (gt_expect_false(map->next_map==NULL)) ? NULL : map->next_map->map;
}
GT_INLINE gt_junction_t gt_map_get_next_block_junction(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  GT_MAP_NEXT_BLOCK_CHECK(map);
  return map->next_map->junction;
}
GT_INLINE int64_t gt_map_get_next_block_distance(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  GT_MAP_NEXT_BLOCK_CHECK(map);
  return (map->next_map->map->position-(map->position+gt_map_get_length(map)));
}
GT_INLINE void gt_map_set_next_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  register gt_next_map *aux_next_map = map->next_map;
  map->next_map = malloc(sizeof(gt_next_map));
  map->next_map->junction = junction;
  map->next_map->map = next_map;
  next_map->next_map = aux_next_map;
}
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  register gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
  }
  return aux_map;
}
GT_INLINE void gt_map_append_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction) {
  GT_MAP_EDITABLE_CHECK(map);
  GT_MAP_EDITABLE_CHECK(next_map);
  register gt_map* aux_map = gt_map_get_last_block(map);
  gt_map_set_next_block(aux_map,next_map,junction);
}

/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms) {
  GT_MAP_CHECK(map);
  gt_vector_insert(map->mismatches,*misms,gt_misms);
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
 * High-level Procedures
 */

// Global metrics
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  register uint64_t length = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    length += gt_map_get_length(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return length;
}
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  register uint64_t distance = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    distance += gt_map_get_distance(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return distance;
}
GT_INLINE uint64_t gt_map_get_global_score(gt_map* const map) {
  GT_MAP_EDITABLE_CHECK(map);
  register uint64_t score = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    score += gt_map_get_score(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return score;
}
// Vector based ( Metrics out of a set of maps )
GT_INLINE uint64_t gt_map_vector_get_length(gt_vector* const maps) {
  GT_NULL_CHECK(maps);
  register uint64_t length = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    length += gt_map_get_global_length(*map);
  }
  return length;
}
GT_INLINE uint64_t gt_map_vector_get_distance(gt_vector* const maps) {
  GT_NULL_CHECK(maps);
  register uint64_t distance = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    distance += gt_map_get_global_distance(*map);
  }
  return distance;
}
GT_INLINE uint64_t gt_map_vector_get_score(gt_vector* const maps) {
  GT_NULL_CHECK(maps);
  register uint64_t score = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    score += gt_map_get_global_score(*map);
  }
  return score;
}
// Distance procedures
GT_INLINE uint64_t gt_map_get_levenshtein_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t lev_distance = 0;
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS:
        ++lev_distance;
        break;
      case INS:
      case DEL:
        lev_distance += gt_misms_get_size(misms_it);
        break;
    }
  }
  return lev_distance;
}
GT_INLINE uint64_t gt_map_get_global_levenshtein_distance(gt_map* const map) {
  register uint64_t lev_distance = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    lev_distance += gt_map_get_levenshtein_distance(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return lev_distance;
}

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map) {
  GT_MAP_CHECK(map);
  gt_map* map_cpy = gt_map_new();
  map_cpy->seq_name = map->seq_name;
  map_cpy->position = map->position;
  map_cpy->base_length = map->base_length;
  map_cpy->direction = map->direction;
  map_cpy->score = map->score;
  gt_vector_copy(map_cpy->mismatches,map->mismatches);
  if (map->mismatches_txt!=NULL) {
    map_cpy->mismatches_txt = map->mismatches_txt;
    map_cpy->next_map = NULL;
  } else {
    map_cpy->mismatches_txt = NULL;
    if (map->next_map==NULL) {
      map_cpy->next_map = NULL;
    } else {
      map_cpy->next_map = malloc(sizeof(gt_next_map));
      map_cpy->next_map->junction = map->next_map->junction;
      map_cpy->next_map->map = gt_map_copy(map->next_map->map);
    }
  }
  return map_cpy;
}

/*
 * Map's Blocks iterator
 */
GT_INLINE void gt_map_new_block_iterator(gt_map* const map,gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  GT_MAP_EDITABLE_CHECK(map);
  map_block_iterator->map = map;
  map_block_iterator->next_map = map;
}
GT_INLINE gt_map* gt_map_next_block(gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  register gt_map* returned_map = map_block_iterator->next_map;
  map_block_iterator->next_map = (returned_map!=NULL && returned_map->next_map!=NULL) ?
      returned_map->next_map->map : NULL;
  return returned_map;
}

/*
 * Map's Mismatches iterator
 */
GT_INLINE void gt_map_new_misms_iterator(gt_map* const map,gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_EDITABLE_CHECK(map);
  map_mism_iterator->map = map;
  map_mism_iterator->next_pos = 0;
  map_mism_iterator->total_pos = gt_vector_get_used(map->mismatches);
  map_mism_iterator->next_misms = gt_vector_get_mem(map->mismatches,gt_misms*);
}
GT_INLINE gt_misms* gt_map_next_misms(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_EDITABLE_CHECK(map_mism_iterator->map);
  if (gt_expect_true(map_mism_iterator->next_pos<map_mism_iterator->total_pos)) {
    register gt_misms* const misms = *(map_mism_iterator->next_misms);
    ++map_mism_iterator->next_misms;
    ++map_mism_iterator->next_pos;
    return misms;
  } else {
    return NULL;
  }
}
GT_INLINE gt_misms* gt_map_dinamic_next_misms(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_EDITABLE_CHECK(map_mism_iterator->map);
  register gt_map* const map = map_mism_iterator->map;
  if (gt_expect_true(map_mism_iterator->next_pos<gt_vector_get_used(map->mismatches))) {
    register gt_misms* const misms = *gt_vector_get_elm(map->mismatches,map_mism_iterator->next_pos,gt_misms*);
    ++map_mism_iterator->next_pos;
    return misms;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_EDITABLE_CHECK(map_mism_iterator->map);
  return map_mism_iterator->next_pos;
}
