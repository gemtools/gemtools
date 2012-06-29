/*
 * PROJECT: GEM-Tools library
 * FILE: gt_paired_alignment.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

//#include "gt_commons.h"
//#include "gt_paired_alignment.h"
//
///*
// * Paired Alignment
// */
//GT_INLINE gt_paired_alignment* gt_paired_alignment_new() {
//  gt_template* template = gt_template_new();
////  gt_template_set_num_blocks_template(template,2);
//  return (gt_paired_alignment*) template;
//}
//GT_INLINE void gt_paired_alignment_delete(gt_paired_alignment* const paired_alignment) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_delete(paired_alignment);
//}
//GT_INLINE void gt_paired_alignment_clear(gt_paired_alignment* const paired_alignment) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_clear(paired_alignment);
//}
//
///*
// * Accessors
// */
//GT_INLINE void gt_paired_alignment_set_ends(
//    gt_paired_alignment* const paired_alignment,
//    gt_alignment* const alignment_end1,gt_alignment* const alignment_end2) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_add_block(paired_alignment,alignment_end1);
//  gt_template_add_block(paired_alignment,alignment_end2);
//}
//GT_INLINE void gt_paired_alignment_clear_ends(gt_paired_alignment* const paired_alignment) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_clear_blocks(paired_alignment);
//}
//
///*
// * Matches handlers
// */
//GT_INLINE void gt_paired_alignment_add_match(
//    gt_paired_alignment* const paired_alignment,const uint64_t total_distance,
//    gt_map* map_end1,gt_map* map_end2) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
////  gt_template_add_match(paired_alignment,total_distance,map_end1,map_end2);
//}
//GT_INLINE void gt_paired_alignment_clear_matches(gt_paired_alignment* const paired_alignment) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_clear_matches(paired_alignment);
//}
//GT_INLINE uint64_t gt_paired_alignment_get_num_matches(gt_paired_alignment* const paired_alignment) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_template_get_num_matches(paired_alignment);
//}
//GT_INLINE void gt_paired_alignment_get_match(
//    gt_paired_alignment* const paired_alignment,const uint64_t position,
//    gt_map** map_end1,gt_map** map_end2) {
//  GT_PAIRED_ALIGNMENT_CHECK(paired_alignment);
//  gt_map* maps_array[2];
//  gt_template_get_match(paired_alignment,position,maps_array);
//  *map_end1 = maps_array[0];
//  *map_end2 = maps_array[1];
//}
//
///*
// * Template added functions
// */
//GT_INLINE bool gt_template_is_paired(gt_template* const template) {
//  gt_fatal_check(template==NULL,NULL_HANDLER);
//  return gt_template_get_num_blocks(template)==2;
//}

