/*
 * PROJECT: GEM-Tools library
 * FILE: gt_iterators.h
 * DATE: 02/07/2012
 * DESCRIPTION: // TODO: Explode this
 */

#ifndef GT_ITERATORS_H_
#define GT_ITERATORS_H_

/*
 * Iterate over the alignment of a template (individual blocks)
 *   Template = Alignment_end1 + Alignment_end2 + {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
 *   GT_ALIGNMENT_ITERATE(template) := {Alignment_end1,Alignment_end2}
 */
#define GT_ALIGNMENT_ITERATE(template,alignment) \
  /* Template. Iterate over all alignments */ \
  gt_template_alignment_iterator alignment##_iterator; \
  register gt_alignment* alignment; \
  gt_template_new_alignment_iterator(template,&(alignment##_iterator)); \
  while ((alignment=gt_template_next_alignment(&(alignment##_iterator))))
/*
 * Iterate over the map(s) of the template
 *   (Eg. Single End => maps)
 *   (Eg. Paired Alignment => pairs of maps (map_end1,map_end2) )
 *   Template = Alignment_end1 + Alignment_end2 + {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
 *   GT_TEMPLATE_ITERATE(template) := {(end1.map1,end2.map1),(end1.map2,end2.map2)}
 */
#define GT_TEMPLATE_ITERATE_(template,map_array) \
  gt_map** map_array; \
  gt_template_maps_iterator __##template##_maps_iterator; \
  gt_template_new_mmap_iterator(template,&(__##template##_maps_iterator)); \
  while (gt_template_next_mmap(&(__##template##_maps_iterator),&map_array,NULL))
#define GT_TEMPLATE_ITERATE(template,map_array) \
  register const uint64_t __##map_array##_num_blocks = gt_template_get_num_blocks(template); \
  GT_TEMPLATE_ITERATE_(template,map_array)
// Querying also attributes {distance, score, ...}
#define GT_TEMPLATE__ATTR_ITERATE(template,map_array,map_array_attr) \
  register const uint64_t __map_array##_num_blocks = gt_template_get_num_blocks(template); \
  gt_map** map_array; \
  gt_mmap_attributes *map_array_attr; \
  gt_template_maps_iterator __##template##_maps_iterator; \
  gt_template_new_mmap_iterator(template,&(__##template##_maps_iterator)); \
  while (gt_template_next_mmap(&(__##template##_maps_iterator),&map_array,&map_array_attr))

/*
 * Iterate over array of maps provided by GT_TEMPLATE_ITERATE
 *   map_array = (end1.map1,end2.map1)
 *   GT_MAP_ARRAY_ITERATE(map_array) := {end1.map1,end2.map1}
 */
#define GT_MAP_ARRAY_ITERATE(map_array,map,end_position) \
  register uint64_t end_position; \
  register gt_map* map; \
  for (end_position=0,map=*map_array; \
       end_position<(__map_array##_num_blocks); \
       ++end_position,map=*(map_array+end_position))
/*
 * Iterate over the map of an alignment
 *   Alignment = {(map1),(map2.block1,map2.block2),(map3),(map4)}
 *   GT_MAPS_ITERATE := {(map1),(map2.block1,map2.block2),(map3),(map4)}
 */
#define GT_MAPS_ITERATE(alignment,map) \
  /* Alignment. Iterate over all maps */ \
  gt_alignment_map_iterator __##map##_iterator; \
  register gt_map* map; \
  gt_alignment_new_map_iterator(alignment,&(__##map##_iterator)); \
  while ((map=gt_alignment_next_map(&(__##map##_iterator))))
/*
 * Iterate over the map blocks of a map
 *   Map = (map.block1,map.block2)
 *   GT_MAP_BLOCKS_ITERATE := {map.block1,map.block2}
 */
#define GT_MAP_BLOCKS_ITERATE(map,map_block) \
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

#endif /* GT_ITERATORS_H_ */
