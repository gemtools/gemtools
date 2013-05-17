/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_align.h"
#include "gt_output_map.h"

// Factor to multiply the read length as to allow expansion at realignment
#define GT_MAP_REALIGN_EXPANSION_FACTOR (0.20)

// Types of misms (as to backtrack the mismatches/indels)
#define GT_MAP_ALG_MISMS_NONE 0
#define GT_MAP_ALG_MISMS_MISMS 1
#define GT_MAP_ALG_MISMS_INS 2
#define GT_MAP_ALG_MISMS_DEL 3

/*
 * Map check/recover operators
 */
GT_INLINE gt_status gt_map_block_check_alignment(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  // Check Alignment
  uint64_t pattern_centinel=0, sequence_centinel=0;
  // Init misms
  const uint64_t num_misms = gt_map_get_num_misms(map);
  uint64_t misms_offset=0;
  gt_misms* misms;
  GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
  // Traverse the sequence
  while (pattern_centinel<pattern_length || sequence_centinel<sequence_length) {
    if (misms!=NULL && misms->position==pattern_centinel) { // Misms
      switch (misms->misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) return GT_MAP_CHECK_ALG_MISMS_OUT_OF_SEQ;
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_NO_MISMS;
          if (misms->base!=sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_BAD_MISMS;
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms->size>sequence_length) return GT_MAP_CHECK_ALG_INS_OUT_OF_SEQ;
          sequence_centinel+=misms->size;
          break;
        case DEL:
          if (pattern_centinel+misms->size>pattern_length) return GT_MAP_CHECK_ALG_DEL_OUT_OF_SEQ;
          pattern_centinel+=misms->size;
          break;
      }
      ++misms_offset;
      GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
    } else { // Match
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_MISMATCH;
      ++pattern_centinel;
      ++sequence_centinel;
    }
  }
  // Ok, no complains
  return 0;
}
GT_INLINE gt_status gt_map_block_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  const uint64_t sequence_length = gt_map_get_length(map);
  gt_string* const sequence = gt_string_new(sequence_length+1);
  gt_status error_code;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Check Alignment
  error_code = gt_map_block_check_alignment(map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),sequence_length);
  gt_string_delete(sequence);
  return error_code;
}
GT_INLINE gt_status gt_map_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    return gt_map_block_check_alignment_sa(map,pattern,sequence_archive);
  } else { // Realigning Slit-Maps (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_check_alignment_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}

GT_INLINE gt_status gt_map_block_recover_mismatches(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  GT_ZERO_CHECK(pattern_length); GT_ZERO_CHECK(sequence_length);
  // Set misms pattern
  gt_misms new_misms;
  new_misms.misms_type = MISMS;
  const uint64_t num_base_misms = gt_map_get_num_misms(map);
  // Traverse pattern & annotate mismatches
  uint64_t sequence_centinel=0, pattern_centinel=0;
  uint64_t misms_offset=0;
  gt_misms misms = {.position = UINT64_MAX};
  GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_base_misms);
  while (sequence_centinel<sequence_length || pattern_centinel<pattern_length) {
    if (misms.position==pattern_centinel) { // Misms annotated
      switch (misms.misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
            return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
          }
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_NO_MISMS;
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms.size>sequence_length) return GT_MAP_CHECK_ALG_INS_OUT_OF_SEQ;
          sequence_centinel+=misms.size;
          break;
        case DEL:
          if (pattern_centinel+misms.size>pattern_length) return GT_MAP_CHECK_ALG_DEL_OUT_OF_SEQ;
          pattern_centinel+=misms.size;
          break;
      }
      // Add the misms + reload old misms
      gt_map_add_misms(map,&misms);
      ++misms_offset;
      GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_base_misms);
    } else { // Nothing annotated
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
        return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
      }
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) {
        new_misms.position = pattern_centinel;
        new_misms.base = sequence[sequence_centinel];
        gt_map_add_misms(map,&new_misms);
      }
      ++pattern_centinel; ++sequence_centinel;
    }
  }
  // Set new misms vector
  const uint64_t total_misms = gt_map_get_num_misms(map);
  uint64_t i, j;
  for (i=0,j=num_base_misms;j<total_misms;++j,++i) {
    gt_map_set_misms(map,gt_map_get_misms(map,j),i);
  }
  gt_map_set_num_misms(map,total_misms-num_base_misms);
  return 0;
}
GT_INLINE gt_status gt_map_block_recover_mismatches_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  gt_status error_code;
  const uint64_t pattern_length = gt_string_get_length(pattern);
  const uint64_t sequence_length = gt_map_get_length(map);
  gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Recover mismatches
  const uint64_t num_misms = gt_map_get_num_misms(map);
  if ((error_code=gt_map_block_recover_mismatches(map,gt_string_get_string(pattern),
      pattern_length,gt_string_get_string(sequence),sequence_length))) {
    gt_map_set_num_misms(map,num_misms); // Restore state
    gt_error(MAP_RECOVER_MISMS_WRONG_BASE_ALG);
    gt_string_delete(sequence);
    return error_code;
  }
  gt_string_delete(sequence);
  return 0;
}
GT_INLINE gt_status gt_map_recover_mismatches_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    if ((error_code=gt_map_block_recover_mismatches_sa(map,pattern,sequence_archive))) {
      return error_code;
    }
    return 0;
  } else { // Split-Map (let's try not to spoil the splice-site consensus)
    uint64_t offset = 0;
    gt_string* read_chunk = gt_string_new(0);
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_recover_mismatches_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}

/*
 * Map (Re)alignment operators
 */
GT_INLINE gt_status gt_map_block_realign_hamming(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_NULL_CHECK(sequence);
  GT_ZERO_CHECK(length);
  // Set misms pattern & clear map mismatches
  gt_misms misms;
  misms.misms_type = MISMS;
  gt_map_clear_misms(map);
  // Traverse pattern & annotate mismatches
  uint64_t i;
  for (i=0;i<length;++i) {
    if (pattern[i]!=sequence[i]) {
      misms.position = i;
      misms.base = sequence[i];
      gt_map_add_misms(map,&misms);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  gt_status error_code;
  const uint64_t pattern_length = gt_string_get_length(pattern);
  gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,0,sequence))) {
    gt_string_delete(sequence);
    return error_code;
  }
  // Realign Hamming
  error_code=gt_map_block_realign_hamming(map,gt_string_get_string(pattern),gt_string_get_string(sequence),pattern_length);
  gt_string_delete(sequence);
  return error_code;
}
GT_INLINE gt_status gt_map_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    return gt_map_block_realign_hamming_sa(map,pattern,sequence_archive);
  } else { // Realigning Slit-Maps (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_realign_hamming_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}

#define GT_DP(i,j) dp_array[(i)*pattern_len+(j)]
#define GT_DP_SET_MISMS(misms,position_pattern,position_sequence,prev_misms,num_misms) { \
  misms.misms_type = MISMS; \
  misms.position = position_pattern; \
  misms.base = sequence[position_sequence]; \
  gt_map_add_misms(map,&misms); ++num_misms; \
  prev_misms = GT_MAP_ALG_MISMS_MISMS; \
}
#define GT_DP_SET_INS(map,misms,pos,length,prev_misms,num_misms) { \
  if (prev_misms != GT_MAP_ALG_MISMS_INS) { \
    misms.position = pos-length+2; \
    misms.misms_type = INS; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
    _misms->size+=length; \
  } \
  prev_misms = GT_MAP_ALG_MISMS_INS; \
}
#define GT_DP_SET_DEL(map,misms,pos,length,prev_misms,num_misms) { \
    if (prev_misms != GT_MAP_ALG_MISMS_DEL) { \
    misms.position = pos-length+1; \
    misms.misms_type = DEL; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
    _misms->position-=length; \
    _misms->size+=length; \
  } \
  prev_misms = GT_MAP_ALG_MISMS_DEL; \
}
GT_INLINE void gt_map_realign_dp_matrix_print(
    uint64_t* const dp_array,const uint64_t pattern_len,const uint64_t sequence_len,
    const uint64_t pattern_limit,const uint64_t sequence_limit) {
  uint64_t i, j;
  for (j=0;j<pattern_limit;++j) {
    for (i=0;i<sequence_limit;++i) {
      fprintf(stderr,"%02"PRIu64" ",GT_DP(i,j));
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}
GT_INLINE gt_status gt_map_block_realign_levenshtein(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,const bool ends_free) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_ZERO_CHECK(pattern_length);
  GT_NULL_CHECK(sequence); GT_ZERO_CHECK(sequence_length);
  // Clear map misms
  gt_map_clear_misms(map);
  // Allocate DP matrix
  const uint64_t pattern_len = pattern_length+1;
  const uint64_t sequence_len = sequence_length+1;
  uint64_t* dp_array = gt_calloc(pattern_len*sequence_len,uint64_t,false);
  uint64_t min_val = UINT32_MAX;
  uint64_t i, j, i_pos = sequence_len-1;
  // Calculate DP-Matrix
  for (i=0;i<sequence_len;++i) GT_DP(i,0)=(ends_free)?0:i;
  for (j=0;j<pattern_len;++j) GT_DP(0,j)=j;
  for (i=1;i<sequence_len;++i) {
    for (j=1;j<pattern_len;++j) {
      const uint64_t ins = GT_DP(i-1,j) + 1;
      const uint64_t del = GT_DP(i,j-1) + 1;
      const uint64_t sub = GT_DP(i-1,j-1) + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
      GT_DP(i,j) = GT_MIN(sub,GT_MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && GT_DP(i,pattern_length) < min_val) {
      min_val = GT_DP(i,pattern_length);
      i_pos = i;
    }
  }
  // DEBUG gt_map_realign_dp_matrix_print(dp_array,pattern_len,sequence_len,30,30);
  // Backtrack all edit operations
  uint64_t num_misms = 0, prev_misms = GT_MAP_ALG_MISMS_NONE;
  gt_misms misms;
  for (i=i_pos,j=pattern_len-1;i>0 && j>0;) {
    const uint32_t current_cell = GT_DP(i,j);
    if (sequence[i-1]==pattern[j-1]) { // Match
      prev_misms = GT_MAP_ALG_MISMS_NONE;
      --i; --j;
    } else {
      if (GT_DP(i-1,j)+1 == current_cell) { // Ins
        GT_DP_SET_INS(map,misms,j-1,1,prev_misms,num_misms);
        --i;
      } else if (GT_DP(i,j-1)+1 == current_cell) { // Del
        GT_DP_SET_DEL(map,misms,j-1,1,prev_misms,num_misms);
        --j;
      } else if (GT_DP(i-1,j-1)+1 == current_cell) { // Misms
        GT_DP_SET_MISMS(misms,j-1,i-1,prev_misms,num_misms);
        --i; --j;
      }
    }
  }
  if (i>0) {
    if (!ends_free) { // Insert the rest of the pattern
      GT_DP_SET_INS(map,misms,i-2,i,prev_misms,num_misms);
    } else {
      map->position+=i;
    }
  }
  if (j>0) { // Delete the rest of the sequence
    GT_DP_SET_DEL(map,misms,j-1,j,prev_misms,num_misms);
  }
  // Flip all mismatches
  uint64_t z;
  const uint64_t mid_point = num_misms/2;
  for (z=0;z<mid_point;++z) {
    misms = *gt_map_get_misms(map,z);
    gt_map_set_misms(map,gt_map_get_misms(map,num_misms-1-z),z);
    gt_map_set_misms(map,&misms,num_misms-1-z);
  }
  // Set map base length
  gt_map_set_base_length(map,pattern_length);
  // Safe check
  gt_debug_block(gt_map_block_check_alignment(map,pattern,pattern_length,
      sequence+((ends_free)?i:0),gt_map_get_length(map))!=0) {
    gt_output_map_fprint_map_block_pretty(stderr,map,pattern,pattern_length,
       sequence+((ends_free)?i:0),gt_map_get_length(map));
//    gt_cond_fatal_error(gt_map_block_check_alignment(map,pattern,pattern_length,
//      sequence+((ends_free)?i:0),gt_map_get_length(map))!=0,MAP_ALG_WRONG_ALG);
  }
  // Free
  gt_free(dp_array);
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive,const uint64_t extra_length,const bool ends_free) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  gt_status error_code;
  // Retrieve the sequence
  const uint64_t decode_length = (ends_free) ? gt_string_get_length(pattern) : gt_map_get_length(map);
  const uint64_t extra_decode_length = (ends_free) ? extra_length : 0;
  gt_string* const sequence = gt_string_new(decode_length+extra_decode_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      decode_length,extra_decode_length,sequence))) {
    gt_string_delete(sequence); // Free
    return error_code;
  }
  // Realign Levenshtein
  error_code = gt_map_block_realign_levenshtein(map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),gt_string_get_length(sequence),ends_free);
  gt_string_delete(sequence); // Free
  return error_code;
}
GT_INLINE gt_status gt_map_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) {
    return gt_map_block_realign_levenshtein_sa(map,pattern,sequence_archive,GT_MAP_REALIGN_EXPANSION_FACTOR,true);
  } else { // Realigning SM (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_realign_levenshtein_sa(map_block,read_chunk,sequence_archive,0,false))) {
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}

GT_INLINE gt_status gt_map_block_realign_weighted(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,int32_t (*gt_weigh_fx)(char*,char*)) {
  // TODO
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive,const uint64_t extra_length,
    int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  gt_status error_code;
  // Retrieve the sequence
  const uint64_t pattern_length = gt_string_get_length(pattern);
  gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,extra_length,sequence))) return error_code;
  // Realign Weighted
  return gt_map_block_realign_weighted(map,
      gt_string_get_string(pattern),pattern_length,
      gt_string_get_string(sequence),gt_string_get_length(sequence),gt_weigh_fx);
}
GT_INLINE gt_status gt_map_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) {
    return gt_map_block_realign_weighted_sa(map,pattern,sequence_archive,GT_MAP_REALIGN_EXPANSION_FACTOR,gt_weigh_fx);
  } else { // Realigning SM (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_realign_weighted_sa(
          map_block,read_chunk,sequence_archive,GT_MAP_REALIGN_EXPANSION_FACTOR,gt_weigh_fx))) {
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}

