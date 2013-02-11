/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_align.h"

#define GT_MAP_REALIGN_EXPANSION_FACTOR (0.20)

GT_INLINE gt_status gt_map_check_alignment(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  // TODO
  return 0;
}
GT_INLINE gt_status gt_map_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  register const uint64_t pattern_length = gt_map_get_length(map);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  register gt_status error_code;
  if ((error_code=gt_sequence_archive_get_sequence_string(
      sequence_archive,gt_map_get_seq_name(map),gt_map_get_strand(map),
      gt_map_get_position(map),pattern_length,sequence))) return error_code;
  // Realign Hamming
  return gt_map_realign_hamming(map,
      gt_string_get_string(pattern),gt_string_get_string(sequence),pattern_length);
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
  register const uint64_t pattern_length = gt_string_get_length(pattern);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  register gt_status error_code;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,0,sequence))) return error_code;
  // Realign Hamming
  return gt_map_realign_hamming(map,gt_string_get_string(pattern),gt_string_get_string(sequence),pattern_length);
}

#define GT_DP(i,j) dp_array[(i)*pattern_len+(j)]
#define GT_DP_SET_MISMS(misms,position_pattern,position_sequence,num_misms) { \
  misms.misms_type = MISMS; \
  misms.position = position_pattern; \
  misms.base = sequence[position_sequence]; \
  gt_map_add_misms(map,&misms); ++num_misms; \
}
#define GT_DP_SET_INS(map,misms,pos,length,num_misms) { \
  if (num_misms==0 || misms.misms_type != INS) { \
    misms.misms_type = INS; \
    misms.position = pos; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_map_get_misms(map,num_misms-1)->size+=length; \
  } \
}
#define GT_DP_SET_DEL(map,misms,pos,length,num_misms) { \
  if (num_misms==0 || misms.misms_type != DEL) { \
    misms.misms_type = DEL; \
    misms.position = pos; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_map_get_misms(map,num_misms-1)->size+=length; \
  } \
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
    char* const sequence,const uint64_t sequence_length) {
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
  register uint64_t i, j, i_pos = UINT64_MAX;
  // Calculate DP-Matrix
  for (i=0;i<sequence_len;++i) GT_DP(i,0)=i;
  for (j=0;j<pattern_len;++j) GT_DP(0,j)=j;
  for (i=1;i<sequence_len;++i) {
    for (j=1;j<pattern_len;++j) {
      register const uint64_t ins = GT_DP(i-1,j) + 1;
      register const uint64_t del = GT_DP(i,j-1) + 1;
      register const uint64_t sub = GT_DP(i-1,j-1) + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
      GT_DP(i,j) = GT_MIN(sub,GT_MIN(ins,del));
    }
    // Check last cell value
    if (GT_DP(i,pattern_length) < min_val) {
      min_val = GT_DP(i,pattern_length);
      i_pos = i;
    }
  }
  // DEBUG gt_dp_matrix_display(dp_array,pattern_len,sequence_len,30,30);
  // Backtrack all edit operations
  register uint64_t num_misms = 0;
  gt_misms misms;
  for (i=i_pos,j=pattern_len-1;i>0 && j>0;) {
    register const uint32_t current_cell = GT_DP(i,j);
    if (sequence[i-1]==pattern[j-1]) { // Match
      --i; --j;
    } else {
      if (GT_DP(i-1,j)+1 == current_cell) { // Ins
        GT_DP_SET_INS(map,misms,j-1,1,num_misms);
        --i;
      } else if (GT_DP(i,j-1)+1 == current_cell) { // Del
        GT_DP_SET_DEL(map,misms,j-1,1,num_misms);
        --j;
      } else if (GT_DP(i-1,j-1)+1 == current_cell) { // Misms
        GT_DP_SET_MISMS(misms,j-1,i-1,num_misms);
        --i; --j;
      }
    }
  }
  if (i>0) { // Insert the rest of the sequence
    GT_DP_SET_INS(map,misms,0,i,num_misms);
  } else if (j>0) { // Delete the rest of the pattern
    GT_DP_SET_DEL(map,misms,j-1,j,num_misms);
  }
  // Flip all mismatches
  register const uint64_t mid_point = num_misms/2;
  for (i=0;i<mid_point;++i) {
    misms = *gt_map_get_misms(map,i);
    gt_map_set_misms(map,gt_map_get_misms(map,num_misms-1-i),i);
    gt_map_set_misms(map,&misms,num_misms-1-i);
  }
  // Free
  free(dp_array);
  return 0;
}
GT_INLINE gt_status gt_map_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  register gt_status error_code;
  // Retrieve the sequence
  register const uint64_t pattern_length = gt_string_get_length(pattern);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  register const uint64_t extra_length = GT_MAP_REALIGN_EXPANSION_FACTOR*pattern_length;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,extra_length,sequence))) return error_code;
  // Realign Levenshtein
  return gt_map_realign_levenshtein(map,
      gt_string_get_string(pattern),pattern_length,
      gt_string_get_string(sequence),gt_string_get_length(sequence));
}
GT_INLINE gt_status gt_map_realign_weighted(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,int32_t (*gt_weigh_fx)(char*,char*)) {
  // TODO
  return 0;
}
GT_INLINE gt_status gt_map_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  register gt_status error_code;
  // Retrieve the sequence
  register const uint64_t pattern_length = gt_string_get_length(pattern);
  register gt_string* const sequence = gt_string_new(pattern_length+1);
  register const uint64_t extra_length = GT_MAP_REALIGN_EXPANSION_FACTOR*pattern_length;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,extra_length,sequence))) return error_code;
  // Realign Weighted
  return gt_map_realign_weighted(map,
      gt_string_get_string(pattern),pattern_length,
      gt_string_get_string(sequence),gt_string_get_length(sequence),gt_weigh_fx);
}

