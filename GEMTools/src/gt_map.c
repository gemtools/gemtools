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
  gt_map* map = gt_alloc(gt_map);
  map->seq_name = gt_string_new(GT_MAP_INITIAL_SEQ_NAME_SIZE);
  map->position = 0;
  map->base_length = 0;
  map->gt_score = GT_MAP_NO_GT_SCORE;
  map->phred_score = GT_MAP_NO_PHRED_SCORE;
  map->mismatches = gt_vector_new(GT_MAP_NUM_INITIAL_MISMS,sizeof(gt_misms));
  map->next_block = NULL;
  map->attributes = NULL;
  return map;
}
GT_INLINE void gt_map_clear(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_clear(map->seq_name);
  map->position = 0;
  map->base_length = 0;
  map->gt_score = GT_MAP_NO_GT_SCORE;
  map->phred_score = GT_MAP_NO_PHRED_SCORE;
  gt_map_clear_misms(map);
  map->next_block = NULL;
  if (map->attributes!=NULL) gt_attribute_clear(map->attributes);
}
GT_INLINE void gt_map_block_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_delete(map->seq_name);
  gt_vector_delete(map->mismatches);
  if (map->attributes!=NULL) gt_attribute_delete(map->attributes);
  gt_free(map);
}
GT_INLINE void gt_map_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  if (map->next_block != NULL) {
    gt_map_delete(map->next_block->map);
    gt_free(map->next_block);
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
GT_INLINE uint64_t gt_map_get_seq_name_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_string_get_length(map->seq_name);
}
GT_INLINE gt_string* gt_map_get_string_seq_name(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->seq_name;
}
GT_INLINE void gt_map_set_seq_name(gt_map* const map,char* const seq_name,const uint64_t length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(seq_name);
  gt_string_set_nstring(map->seq_name,seq_name,length);
}
GT_INLINE void gt_map_set_string_seq_name(gt_map* const map,gt_string* const seq_name) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(seq_name);
  gt_string_set_nstring(map->seq_name,gt_string_get_string(seq_name),gt_string_get_length(seq_name));
}
GT_INLINE uint64_t gt_map_get_global_position(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->strand==FORWARD) ? gt_map_get_position(map): gt_map_get_position(gt_map_get_last_block(map));
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
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->base_length;
}
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  map->base_length = length;
}

/*
 * Attributes (lazy allocation of the attributes field)
 */
GT_INLINE void* gt_map_attribute_get(gt_map* const map,char* const attribute_id) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(attribute_id);
  return (map->attributes==NULL) ? NULL : gt_attribute_get(map->attributes,attribute_id);
}
GT_INLINE void gt_map_attribute_set_(gt_map* const map,char* const attribute_id,void* const attribute,const size_t element_size) {
  GT_MAP_CHECK(map);
  if (map->attributes==NULL) {
    map->attributes = gt_attribute_new();
  }
  gt_attribute_add_primitive(map->attributes,attribute_id,attribute,element_size);
}

/*
 * Local metrics (over first blocks)
 */
GT_INLINE uint64_t gt_map_get_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  register int64_t length = map->base_length;
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
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}

GT_INLINE int64_t gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->gt_score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const int64_t score) {
  GT_MAP_CHECK(map);
  map->gt_score = score;
}
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->phred_score;
}
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score) {
  GT_MAP_CHECK(map);
  map->phred_score = phred_score;
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
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  register gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
  }
  return aux_map;
}
GT_INLINE void gt_map_set_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  if (gt_expect_true(next_map!=NULL)) {
    GT_MAP_CHECK(next_map);
    if (map->next_block==NULL) map->next_block = gt_alloc(gt_map_junction);
    map->next_block->junction = junction;
    map->next_block->junction_size = junction_size;
    map->next_block->map = next_map;
  } else {
    if (map->next_block!=NULL) {
      gt_free(map->next_block);
      map->next_block = NULL;
    }
  }
}
GT_INLINE void gt_map_insert_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  register gt_map_junction *aux_next_block = map->next_block;
  if (map->next_block==NULL) map->next_block = gt_alloc(gt_map_junction);
  map->next_block->junction = junction;
  map->next_block->map = next_map;
  next_map->next_block = aux_next_block;
}
GT_INLINE void gt_map_append_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  register gt_map* aux_map = gt_map_get_last_block(map);
  gt_map_set_next_block(aux_map,next_map,junction,junction_size);
}
GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map) {
  GT_MAP_CHECK(map);
  GT_MAP_NEXT_BLOCK_CHECK(map);
  return (map->next_block == NULL) ? NO_JUNCTION : map->next_block->junction;
}
GT_INLINE int64_t gt_map_get_junction_size(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->next_block == NULL) ? 0 : (map->next_block->junction_size);
}

GT_INLINE uint64_t gt_map_reverse_blocks_positions_(gt_map* const map,const uint64_t start_position) {
  GT_MAP_CHECK(map);
  if (gt_map_has_next_block(map)) {
    register const uint64_t position =
        gt_map_reverse_blocks_positions_(gt_map_get_next_block(map),start_position) + gt_map_get_junction_size(map);
    gt_map_set_position(map,position);
    return position+gt_map_get_length(map);
  } else {
    gt_map_set_position(map,start_position);
    return start_position+gt_map_get_length(map);
  }
}
GT_INLINE void gt_map_reverse_blocks_positions(gt_map* const head_map,const uint64_t start_position) {
  GT_MAP_CHECK(head_map);
  gt_map_reverse_blocks_positions_(head_map,start_position);
}
#define GT_MAP_REVERSE_MISMS_ADJUST_POS(misms,base_length) \
  misms->position = base_length - misms->position; \
  switch (misms->misms_type) { \
    case DEL: misms->position-=misms->size; break; \
    case INS: misms->position++; break; \
    default: break; \
  }
GT_INLINE void gt_map_reverse_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  // Flip all mismatches
  register uint64_t z;
  register const uint64_t base_length = gt_map_get_base_length(map);
  register const uint64_t num_misms = gt_map_get_num_misms(map);
  register const uint64_t mid_point = num_misms/2;
  for (z=0;z<mid_point;++z) {
    gt_misms* misms_a = gt_map_get_misms(map,z);
    gt_misms* misms_b = gt_map_get_misms(map,num_misms-1-z);
    // Correct position
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms_a,base_length);
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms_b,base_length);
    // Flip mismatches
    gt_misms misms = *misms_a;
    gt_map_set_misms(map,misms_b,z);
    gt_map_set_misms(map,&misms,num_misms-1-z);
  }
  if (num_misms%2) {
    gt_misms* misms = gt_map_get_misms(map,mid_point);
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms,base_length);
  }
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
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,const uint64_t offset) {
  GT_MAP_CHECK(map);
  return gt_vector_get_elm(map->mismatches,offset,gt_misms);
}
GT_INLINE void gt_map_set_misms(gt_map* const map,gt_misms* misms,const uint64_t offset) {
  GT_MAP_CHECK(map);
  gt_vector_set_elm(map->mismatches,offset,gt_misms,*misms);
}
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}
GT_INLINE void gt_map_set_num_misms(gt_map* const map,const uint64_t num_misms) {
  GT_MAP_CHECK(map);
  gt_vector_set_used(map->mismatches,num_misms);
}

// Counting
GT_INLINE uint64_t gt_map_get_num_mismatch_bases(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t count = 0;
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS:
        ++count;
        break;
      case INS:
      case DEL:
        break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_indels(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t count = 0;
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
      case DEL:
        ++count;
        break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_insertions(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t count = 0;
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
        ++count;
        break;
      case DEL:   break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_deletions(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t count = 0;
  GT_MAP_MISMS_ITERATOR(map,misms_it,misms_pos) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:   break;
      case DEL:
        ++count;
        break;
    }
  }
  return count;
}

/*
 *
 * High-level Procedures
 *
 */
/*
 * Global metrics
 */
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  register uint64_t length = 0;
  GT_BEGIN_MAP_BLOCKS_ITERATOR(map,map_it) {
    length += gt_map_get_length(map_it) + gt_map_get_junction_size(map_it);
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
/*
 * Begin/End Position
 */
GT_INLINE uint64_t gt_map_get_begin_position(gt_map* const map) {
  return map->position-gt_map_get_left_trim_length(map);
}
GT_INLINE uint64_t gt_map_get_end_position(gt_map* const map) {
  return map->position+gt_map_get_length(map);
}
GT_INLINE uint64_t gt_map_get_global_end_position(gt_map* const map) {
  return gt_map_get_end_position(gt_map_get_last_block(map));
}
/*
 * Trim helpers
 */
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
/*
 * Bases aligned
 */
GT_INLINE uint64_t gt_map_get_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  register int64_t bases_aligned = map->base_length;
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
  gt_cond_fatal_error(bases_aligned<0,MAP_NEG_MAPPED_BASES);
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
/*
 * Distance procedures
 */
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
/*
 * Vector based ( Metrics out of a set of maps )
 */
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
    score += gt_map_get_score(*map);
  }
  return score;
}


/*
 * Strict Map compare functions
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- for mapBlock in map {
 *         4.1 Same begin position
 *         4.2 Same end position
 *       }
 */
GT_INLINE int64_t gt_map_sm_cmp(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  register const int64_t begin_distance = ((int64_t)gt_map_get_begin_position(map_1)) - ((int64_t)gt_map_get_begin_position(map_2));
  if (begin_distance != 0) return 1;
  register const int64_t end_distance = ((int64_t)gt_map_get_end_position(map_1)) - (int64_t)(gt_map_get_end_position(map_2));
  if (end_distance != 0) return 1;
  if (map_1->next_block==NULL && map_2->next_block==NULL) return 0;
  if (map_1->next_block!=NULL && map_2->next_block!=NULL) {
    return gt_map_sm_cmp(map_1->next_block->map,map_2->next_block->map);
  }
  return 1;
}
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  if (gt_string_cmp(map_1->seq_name,map_2->seq_name)!=0) {
    return 1;
  } else {
    if (map_1->strand==map_2->strand) {
      return gt_map_sm_cmp(map_1,map_2);
    } else {
      return 1;
    }
  }
}
/*
 * General Map compare function (differs from the standard one at @gt_map_cmp)
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- If (numBlocks == 1) {
 *           4.1 |begin_a-begin_b|<=range_tolerated || |end_a-end_b|<=range_tolerated
 *       } else {
 *         for mapBlock in map {
 *           5.1 |begin_a-begin_b|<=range_tolerated
 *           5.2 |end_a-end_b|<=range_tolerated
 *         }
 *       }
 */
GT_INLINE int64_t gt_map_range_sm_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated,const uint64_t num_maps_left) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  register const int64_t begin_distance = ((int64_t)gt_map_get_begin_position(map_1)) - ((int64_t)gt_map_get_begin_position(map_2));
  register const int64_t end_distance = ((int64_t)gt_map_get_end_position(map_1)) - (int64_t)(gt_map_get_end_position(map_2));
  if (GT_ABS(begin_distance)<=range_tolerated && GT_ABS(end_distance)<=range_tolerated) {
    return (num_maps_left==1) ? 0 : gt_map_range_sm_cmp(map_1->next_block->map,map_2->next_block->map,range_tolerated,num_maps_left-1);
  } else {
    return GT_ABS(begin_distance)+GT_ABS(end_distance);
  }
}
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  register int64_t cmp_tags = gt_string_cmp(map_1->seq_name,map_2->seq_name);
  if (cmp_tags!=0) {
    return cmp_tags;
  } else {
    if (map_1->strand==map_2->strand) {
      register const uint64_t num_blocks_map_1 = gt_map_get_num_blocks(map_1);
      register const uint64_t num_blocks_map_2 = gt_map_get_num_blocks(map_2);
      if (num_blocks_map_1==num_blocks_map_2) {
        if (num_blocks_map_1==1) { // Standard Mapping
          register const int64_t begin_distance = ((int64_t)gt_map_get_begin_position(map_1)) - ((int64_t)gt_map_get_begin_position(map_2));
          register const int64_t end_distance = ((int64_t)gt_map_get_end_position(map_1)) - (int64_t)(gt_map_get_end_position(map_2));
          if (GT_ABS(begin_distance)<=range_tolerated || GT_ABS(end_distance)<=range_tolerated) return 0;
          return GT_ABS(begin_distance)+GT_ABS(end_distance);
        } else { // Split Maps Involved
          return gt_map_range_sm_cmp(map_1,map_2,range_tolerated,num_blocks_map_1);
        }
      } else { // Different splits
        return (num_blocks_map_1-num_blocks_map_2);
      }
    } else {
      return map_1->strand==FORWARD ? 1 : -1;
    }
  }
}
/*
 * MMap compare functions (based on Map compare functions)
 */
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
GT_INLINE bool gt_map_less_than(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  if (gt_map_get_global_bases_aligned(map_1) < gt_map_get_global_bases_aligned(map_2)) return true;
  if (gt_map_get_num_indels(map_1) < gt_map_get_num_indels(map_2)) return true;
  if (gt_map_get_global_levenshtein_distance(map_1) < gt_map_get_global_levenshtein_distance(map_2)) return true;
  return false;
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
  map_cpy->phred_score = map->phred_score;
  map_cpy->gt_score = map->gt_score;
  gt_vector_copy(map_cpy->mismatches,map->mismatches);
  if (map->next_block==NULL) {
    map_cpy->next_block = NULL;
  } else {
    map_cpy->next_block = gt_alloc(gt_map_junction);
    map_cpy->next_block->junction = map->next_block->junction;
    map_cpy->next_block->junction_size = map->next_block->junction_size;
    map_cpy->next_block->map = gt_map_copy(map->next_block->map);
  }
  return map_cpy;
}
GT_INLINE gt_map** gt_mmap_array_copy(gt_map** mmap,const uint64_t num_blocks) {
  GT_ZERO_CHECK(num_blocks);
  gt_map** mmap_copy = gt_calloc(num_blocks,gt_map*,false);
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
