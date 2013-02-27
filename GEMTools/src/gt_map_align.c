/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_align.h"

// Types of misms (as to backtrack the mismatches/indels)
#define GT_MAP_ALG_MISMS_NONE 0
#define GT_MAP_ALG_MISMS_MISMS 1
#define GT_MAP_ALG_MISMS_INS 2
#define GT_MAP_ALG_MISMS_DEL 3

#define GT_MAP_CHECK_RELOAD_MISMS(map,misms_offset,misms_ptr) { \
  if (misms_offset<gt_map_get_num_misms(map)) { \
    misms=gt_map_get_misms(map,misms_offset); \
  } else { \
    misms=NULL; \
  } \
}
GT_INLINE gt_status gt_map_check_alignment(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  // Check Alignment
  register uint64_t pattern_centinel=0, sequence_centinel=0;
  // Init misms
  register uint64_t misms_offset=0;
  register gt_misms* misms;
  GT_MAP_CHECK_RELOAD_MISMS(map,misms_offset,misms);
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
      GT_MAP_CHECK_RELOAD_MISMS(map,misms_offset,misms);
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
GT_INLINE gt_status gt_map_check_alignment_sa(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  register const uint64_t sequence_length = gt_map_get_length(map);
  register gt_string* const sequence = gt_string_new(sequence_length+1);
  register gt_status error_code;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Check Alignment
  return gt_map_check_alignment(map,pattern,pattern_length,gt_string_get_string(sequence),sequence_length);
}
GT_INLINE gt_status gt_map_realign_hamming(
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
  register uint64_t i;
  for (i=0;i<length;++i) {
    if (pattern[i]!=sequence[i]) {
      misms.position = i;
      misms.base = sequence[i];
      gt_map_add_misms(map,&misms);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_map_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  register gt_status error_code;
  register const uint64_t pattern_length = gt_string_get_length(pattern);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,0,sequence))) return error_code;
  // Realign Hamming
  return gt_map_realign_hamming(map,gt_string_get_string(pattern),gt_string_get_string(sequence),pattern_length);
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
    register gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
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
    register gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
    _misms->position-=length; \
    _misms->size+=length; \
  } \
  prev_misms = GT_MAP_ALG_MISMS_DEL; \
}

GT_INLINE void gt_map_realign_dp_matrix_print(
    uint64_t* const dp_array,const uint64_t pattern_len,const uint64_t sequence_len,
    const uint64_t pattern_limit,const uint64_t sequence_limit) {
  register uint64_t i, j;
  for (j=0;j<pattern_limit;++j) {
    for (i=0;i<sequence_limit;++i) {
      fprintf(stderr,"%02lu ",GT_DP(i,j));
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}
GT_INLINE gt_status gt_map_realign_levenshtein(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,const bool ends_free) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_ZERO_CHECK(pattern_length);
  GT_NULL_CHECK(sequence); GT_ZERO_CHECK(sequence_length);
  // Clear map misms
  gt_map_clear_misms(map);
  // Allocate DP matrix
  register const uint64_t pattern_len = pattern_length+1;
  register const uint64_t sequence_len = sequence_length+1;
  register uint64_t* dp_array = gt_calloc(pattern_len*sequence_len,uint64_t);
  register uint64_t min_val = UINT32_MAX;
  register uint64_t i, j, i_pos = sequence_len-1;
  // Calculate DP-Matrix
  for (i=0;i<sequence_len;++i) GT_DP(i,0)=(ends_free)?0:i;
  for (j=0;j<pattern_len;++j) GT_DP(0,j)=j;
  for (i=1;i<sequence_len;++i) {
    for (j=1;j<pattern_len;++j) {
      register const uint64_t ins = GT_DP(i-1,j) + 1;
      register const uint64_t del = GT_DP(i,j-1) + 1;
      register const uint64_t sub = GT_DP(i-1,j-1) + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
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
  register uint64_t num_misms = 0, prev_misms = GT_MAP_ALG_MISMS_NONE;
  gt_misms misms;
  for (i=i_pos,j=pattern_len-1;i>0 && j>0;) {
    register const uint32_t current_cell = GT_DP(i,j);
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
  register uint64_t z;
  register const uint64_t mid_point = num_misms/2;
  for (z=0;z<mid_point;++z) {
    misms = *gt_map_get_misms(map,z);
    gt_map_set_misms(map,gt_map_get_misms(map,num_misms-1-z),z);
    gt_map_set_misms(map,&misms,num_misms-1-z);
  }
  // Set map base length
  gt_map_set_base_length(map,pattern_length);
  // Safe check
  gt_cond_fatal_error(gt_map_check_alignment(map,pattern,pattern_length,
    sequence+((ends_free)?i:0),gt_map_get_length(map))!=0,MAP_ALG_WRONG_ALG);
  // Free
  free(dp_array);
  return 0;
}
GT_INLINE gt_status gt_map_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive,const uint64_t extra_length,const bool ends_free) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  register gt_status error_code;
  // Retrieve the sequence
  register const uint64_t decode_length = (ends_free) ? gt_string_get_length(pattern) : gt_map_get_length(map);
  register const uint64_t extra_decode_length = (ends_free) ? extra_length : 0;
  register gt_string* const sequence = gt_string_new(decode_length+extra_decode_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      decode_length,extra_decode_length,sequence))) {
    gt_string_delete(sequence); // Free
    return error_code;
  }
  // Realign Levenshtein
  error_code = gt_map_realign_levenshtein(map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),gt_string_get_length(sequence),ends_free);
  gt_string_delete(sequence); // Free
  return error_code;
}
GT_INLINE gt_status gt_map_realign_weighted(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,int32_t (*gt_weigh_fx)(char*,char*)) {
  // TODO
  return 0;
}
GT_INLINE gt_status gt_map_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive,const uint64_t extra_length,
    int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  register gt_status error_code;
  // Retrieve the sequence
  register const uint64_t pattern_length = gt_string_get_length(pattern);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,extra_length,sequence))) return error_code;
  // Realign Weighted
  return gt_map_realign_weighted(map,
      gt_string_get_string(pattern),pattern_length,
      gt_string_get_string(sequence),gt_string_get_length(sequence),gt_weigh_fx);
}

