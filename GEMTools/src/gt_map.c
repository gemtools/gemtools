/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#include "gt_map.h"

#define GT_MAP_NUM_INITIAL_MISMS 4
#define GT_MAP_INITIAL_SEQ_NAME_SIZE 10

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new() {
  gt_map* map = malloc(sizeof(gt_map));
  gt_cond_fatal_error(!map,MEM_HANDLER);
  map->seq_name = gt_string_new(GT_MAP_INITIAL_SEQ_NAME_SIZE);
  map->position = 0;
  map->base_length = 0;
  map->score = GT_MAP_NO_SCORE;
  map->mismatches = gt_vector_new(GT_MAP_NUM_INITIAL_MISMS,sizeof(gt_misms));
  map->misms_txt = NULL;
  map->next_block = NULL;
  return map;
}
GT_INLINE void gt_map_clear(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_clear(map->seq_name);
  map->position = 0;
  map->base_length = 0;
  map->score = GT_MAP_NO_SCORE;
  gt_map_clear_misms(map);
  map->misms_txt = NULL;
  map->next_block = NULL;
}
GT_INLINE void gt_map_block_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_delete(map->seq_name);
  gt_vector_delete(map->mismatches);
  free(map);
}
GT_INLINE void gt_map_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  if (map->next_block != NULL) {
    gt_map_delete(map->next_block->map);
    free(map->next_block);
  }
  gt_map_block_delete(map);
}

/*
 * Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_string_get_string(map->seq_name);
}
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* const seq_name,const uint64_t length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(seq_name);
  gt_string_set_nstring(map->seq_name,seq_name,length);
}
GT_INLINE uint64_t gt_map_get_position(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->position;
}
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position) {
  GT_MAP_CHECK(map);
  map->position = position;
}
GT_INLINE gt_strand gt_map_get_strand(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->strand;
}
GT_INLINE void gt_map_set_strand(gt_map* const map,const gt_strand strand) {
  GT_MAP_CHECK(map);
  map->strand = strand;
}
GT_INLINE uint64_t gt_map_get_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  register int64_t length = map->base_length; // FIXME: What if neg length
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
        length += gt_misms_get_size(misms_it);
        break;
      case DEL:
        length -= gt_misms_get_size(misms_it);
        break;
    }
  }
  gt_cond_fatal_error(length<0,MAP_NEG_LENGTH);
  return (uint64_t)length;
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
GT_INLINE int64_t gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const int64_t score) {
  GT_MAP_CHECK(map);
  map->score = score;
}

/*
 * Multiple Block Maps Handlers
 */
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t num_blocks = 1;
  register gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
    ++num_blocks;
  }
  return num_blocks;
}
GT_INLINE bool gt_map_has_next_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->next_block!=NULL);
}
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (gt_expect_false(map->next_block==NULL)) ? NULL : map->next_block->map;
}
GT_INLINE void gt_map_set_next_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  register gt_map_junction *aux_next_block = map->next_block;
  map->next_block = malloc(sizeof(gt_map_junction));
  gt_cond_fatal_error(!map->next_block,MEM_HANDLER);
  map->next_block->junction = junction;
  map->next_block->map = next_map;
  next_map->next_block = aux_next_block;
}
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  register gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
  }
  return aux_map;
}
GT_INLINE void gt_map_append_block(gt_map* const map,gt_map* const next_map,gt_junction_t junction) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  register gt_map* aux_map = gt_map_get_last_block(map);
  gt_map_set_next_block(aux_map,next_map,junction);
}
GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map) {
  GT_MAP_CHECK(map);
  GT_MAP_NEXT_BLOCK_CHECK(map);
  return (map->next_block == NULL) ? NO_JUNCTION : map->next_block->junction;
}
GT_INLINE uint64_t gt_map_get_junction_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->next_block == NULL) ? 0 : (map->next_block->map->position-(map->position+gt_map_get_length(map)));
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
  gt_vector_clear(map->mismatches);
}
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,uint64_t offset) {
  GT_MAP_CHECK(map);
  return gt_vector_get_elm(map->mismatches,offset,gt_misms);
}
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}

GT_INLINE void gt_map_clear_misms_string(gt_map* const map) {
  GT_MAP_CHECK(map);
  map->misms_txt = NULL;
}
GT_INLINE char* gt_map_get_misms_string(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->misms_txt;
}
GT_INLINE gt_misms_string_t gt_map_get_misms_string_format(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->misms_txt_format;
}
GT_INLINE void gt_map_set_misms_string(gt_map* const map,char* misms_string,const gt_misms_string_t format) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(misms_string);
  map->misms_txt = misms_string;
  map->misms_txt_format = format;
}

/*
 * High-level Procedures
 */
// Trim helpers
GT_INLINE uint64_t gt_map_get_left_trim_length(gt_map* const map) {
  if (gt_map_get_num_misms(map)>0) {
    register gt_misms* const first_misms = gt_map_get_misms(map,0);
    if (first_misms->position==0 && first_misms->misms_type==DEL) return first_misms->size;
  }
  return 0;
}
GT_INLINE uint64_t gt_map_get_right_trim_length(gt_map* const map) {
  register const uint64_t num_misms = gt_map_get_num_misms(map);
  if (num_misms>0) {
    register gt_misms* const last_misms = gt_map_get_misms(map,num_misms-1);
    if (last_misms->position+last_misms->size==map->base_length && last_misms->misms_type==DEL) return last_misms->size;
  }
  return 0;
}
// Bases aligned
GT_INLINE uint64_t gt_map_get_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  register int64_t bases_aligned = map->base_length; // FIXME: What if neg length
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case INS:
        break;
      case MISMS:
        --bases_aligned;
        break;
      case DEL:
        bases_aligned -= gt_misms_get_size(misms_it);
        break;
    }
  }
  gt_cond_fatal_error(bases_aligned<0,MAP_NEG_LENGTH); // FIXME: Error msg
  return (uint64_t)bases_aligned;
}
GT_INLINE uint64_t gt_map_get_global_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t bases_aligned = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    bases_aligned += gt_map_get_bases_aligned(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return bases_aligned;
}
// Global metrics
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t length = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    length += gt_map_get_length(map_it) + gt_map_get_junction_distance(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return length;
}
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t distance = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    distance += gt_map_get_distance(map_it) + (gt_map_has_next_block(map_it)?1:0);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return distance;
}
GT_INLINE uint64_t gt_map_get_global_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  // FIXME: Don't know what to do with this
  return map->score;
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
  GT_MAP_CHECK(map);
  register uint64_t lev_distance = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    lev_distance += gt_map_get_levenshtein_distance(map_it);
  } GT_END_MAP_BLOCKS_ITERATOR;
  return lev_distance;
}
// Map compare
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  return gt_map_range_cmp(map_1,map_2,0);
}
GT_INLINE int64_t gt_map_cmp_strict(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  if (!(map_1->position==map_2->position && gt_map_get_distance(map_1)==gt_map_get_distance(map_2))) {
    return 1;
  } else {
    if (map_1->next_block==NULL || map_2->next_block==NULL) return 1;
    return gt_map_cmp_strict(map_1->next_block->map,map_2->next_block->map);
  }
}
GT_INLINE int64_t gt_map_cmp_false(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  return 1;
}
GT_INLINE int64_t gt_map_cmp_true(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  return 0;
}
#define GT_MAP_RANGE_CMP_NEXT_MAPS(map_1,map_2,range_tolerated) \
  if (map_1->next_block==NULL && map_2->next_block!=NULL) return INT64_MAX; \
  if (map_1->next_block!=NULL && map_2->next_block==NULL) return INT64_MIN; \
  if (map_1->next_block==NULL && map_2->next_block==NULL) return 0; \
  return gt_map_range_cmp(map_1->next_block->map,map_2->next_block->map,range_tolerated)
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2); // TODO: Should be corrected new/old format trim-based position
  if (gt_string_equals(map_1->seq_name,map_2->seq_name) && map_1->strand==map_2->strand) {
    // Cmp BEGIN position
    register const int64_t begin_distance = ((int64_t)map_1->position-(int64_t)gt_map_get_left_trim_length(map_1)) -
        ((int64_t)map_2->position-(int64_t)gt_map_get_left_trim_length(map_2));
    if (GT_ABS(begin_distance)<=range_tolerated) {
      GT_MAP_RANGE_CMP_NEXT_MAPS(map_1,map_2,range_tolerated);
    }
    // Cmp END position // TODO: Consider when the read had been trimmed then there is no base_length reliable
    register const int64_t end_distance = (int64_t)(map_1->position+gt_map_get_length(map_1)) -
        (int64_t)(map_2->position+gt_map_get_length(map_2));
    if (GT_ABS(end_distance)<=range_tolerated) {
      GT_MAP_RANGE_CMP_NEXT_MAPS(map_1,map_2,range_tolerated);
    }
    // Maps are different, return the one that begins first
    return begin_distance;
  } else {
    return INT64_MAX; // TODO: GT_string_cmp();
  }
}
GT_INLINE int64_t gt_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  GT_NULL_CHECK(map_1); GT_NULL_CHECK(map_2);
  register uint64_t i;
  for (i=0;i<num_maps;++i) {
    register const int64_t map_cmp = gt_map_cmp(map_1[i],map_2[i]);
    if (map_cmp!=0) return map_cmp;
  }
  return 0;
}
GT_INLINE int64_t gt_mmap_range_cmp(
    gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps,const uint64_t range_tolerated) {
  GT_NULL_CHECK(map_1); GT_NULL_CHECK(map_2);
  register uint64_t i;
  for (i=0;i<num_maps;++i) {
    register const int64_t map_cmp = gt_map_range_cmp(map_1[i],map_2[i],range_tolerated);
    if (map_cmp!=0) return map_cmp;
  }
  return 0;
}

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* map_cpy = gt_map_new();
  gt_string_copy(map_cpy->seq_name,map->seq_name);
  map_cpy->position = map->position;
  map_cpy->base_length = map->base_length;
  map_cpy->strand = map->strand;
  map_cpy->score = map->score;
  gt_vector_copy(map_cpy->mismatches,map->mismatches);
  map_cpy->misms_txt = map->misms_txt;
  map_cpy->misms_txt_format = map->misms_txt_format;
  if (map->next_block==NULL || map->misms_txt!=NULL) {
    map_cpy->next_block = NULL;
  } else {
    map_cpy->next_block = malloc(sizeof(gt_map_junction));
    map_cpy->next_block->junction = map->next_block->junction;
    map_cpy->next_block->map = gt_map_copy(map->next_block->map);
  }
  return map_cpy;
}
GT_INLINE gt_map** gt_mmap_array_copy(gt_map** mmap,const uint64_t num_blocks) {
  GT_ZERO_CHECK(num_blocks);
  gt_map** mmap_copy = malloc(num_blocks*sizeof(gt_map*));
  gt_cond_fatal_error(!mmap_copy,MEM_HANDLER);
  register uint64_t i;
  for (i=0;i<num_blocks;++i) {
    GT_MAP_CHECK(mmap[i]);
    mmap_copy[i] = gt_map_copy(mmap[i]);
  }
  return mmap_copy;
}

/*
 * Map's Blocks iterator
 */
GT_INLINE void gt_map_new_block_iterator(gt_map* const map,gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  GT_MAP_CHECK(map);
  map_block_iterator->map = map;
  map_block_iterator->next_map = map;
}
GT_INLINE gt_map* gt_map_next_block(gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  register gt_map* returned_map = map_block_iterator->next_map;
  map_block_iterator->next_map = (returned_map!=NULL && returned_map->next_block!=NULL) ?
      returned_map->next_block->map : NULL;
  return returned_map;
}

/*
 * Map's Mismatches iterator
 */
GT_INLINE void gt_map_new_misms_iterator(gt_map* const map,gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map);
  map_mism_iterator->map = map;
  map_mism_iterator->next_pos = 0;
  map_mism_iterator->total_pos = gt_vector_get_used(map->mismatches);
  map_mism_iterator->next_misms = gt_vector_get_mem(map->mismatches,gt_misms);
}
GT_INLINE gt_misms* gt_map_next_misms(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map_mism_iterator->map);
  register gt_map* const map = map_mism_iterator->map;
  if (gt_expect_true(map_mism_iterator->next_pos<gt_vector_get_used(map->mismatches))) {
    register gt_misms* const misms = gt_vector_get_elm(map->mismatches,map_mism_iterator->next_pos,gt_misms);
    ++map_mism_iterator->next_pos;
    return misms;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map_mism_iterator->map);
  return map_mism_iterator->next_pos;
}
