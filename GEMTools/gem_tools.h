/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Common
#include "gt_commons.h"

// Input handlers
#include "gt_input_file.h"
#include "gt_buffered_map_input.h"

// Output handlers
#include "gt_buffered_output_file.h"
//#include "gt_output_map.h"

// GEM-Tools basic data structures: Template/Alignment/Maps/...
#include "gt_misms.h"
#include "gt_map.h"
#include "gt_alignment.h"
#include "gt_template.h"

#define GT_ALIGNMENT_ITERATE(template,alignment) \
  /* Template. Iterate over all alignments */ \
  gt_template_alignment_iterator alignment##_iterator; \
  register gt_alignment* alignment; \
  gt_template_new_alignment_iterator(template,&(alignment##_iterator)); \
  while ((alignment=gt_template_next_alignment(&(alignment##_iterator))))

#define GT_MAPS_ITERATE(alignment,map) \
  /* Alignment. Iterate over all maps */ \
  gt_alignment_map_iterator map##_iterator; \
  register gt_map* map; \
  gt_alignment_new_map_iterator(alignment,&(map##_iterator)); \
  while ((map=gt_alignment_next_map(&(map##_iterator))))

#define GT_MAP_BLOCKS_ITERATE(map,map_block) \
  /* Map. Iterate over map blocks */ \
  gt_map_block_iterator map_block##_iterator; \
  register gt_map* map_block; \
  gt_map_new_block_iterator(map,&(map_block##_iterator)); \
  while ((map_block=gt_map_next_block(&(map_block##_iterator))))

#define GT_MISMS_ITERATE(map_block,misms) \
  /* Map Block. Iterate over all mismatches */ \
  gt_map_mism_iterator misms##_iterator; \
  register gt_misms* misms; \
  gt_map_new_misms_iterator(map_block,&(misms##_iterator)); \
  while ((misms=gt_map_next_misms(&(misms##_iterator))))


#endif /* GEM_TOOLS_H_ */
