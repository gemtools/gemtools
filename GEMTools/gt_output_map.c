/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_commons.h"
#include "gt_output_map.h"

#define GT_BMO_COMPACT_COUNTERS_ZEROS_TH 5

/*
 * Setup
 */
gt_buffered_map_output* gt_buffered_map_output_new(gt_buffered_output_file* const buffered_output_file) {
  //TODO
}
gt_status gt_buffered_map_output_close(gt_buffered_map_output* const buffered_map_output) {
  //TODO
}

/*
 * MAP building block printers
 */
GT_INLINE gt_status gt_output_map_fprint_counters(
    FILE* file,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact) {
  register const uint64_t num_counters = gt_vector_get_used(counters);
  register uint64_t i;
  for (i=0;i<num_counters;) {
    if (i>0) fprintf(file,"%c",gt_expect_false(i==max_complete_strata)?GT_MAP_MCS:GT_MAP_COUNTS_SEP);
    register const uint64_t counter = *gt_vector_get_elm(counters,i,uint64_t);
    if (gt_expect_false(compact && counter==0)) {
      register uint64_t j=i+1;
      while (j<num_counters && *gt_vector_get_elm(counters,j,uint64_t)==0) ++j;
      if (gt_expect_false((j-i)>=GT_BMO_COMPACT_COUNTERS_ZEROS_TH)) {
        fprintf(file,"0" GT_MAP_COUNTS_STIMES "%"PRIu64,(j-i)); i=j;
      } else {
        fprintf(file,"0"); ++i;
      }
    } else {
      fprintf(file,"%"PRIu64,counter); ++i;
    }
  }
  return 0;
}
GT_INLINE gt_status gt_output_map_fprint_map(FILE* file,gt_map* const map,const bool print_scores) {
  // chr11:-:51590050:(5)43T46A9>24*

}
GT_INLINE gt_status gt_output_map_fprint_template_maps(
    FILE* file,gt_template* const template,const uint64_t num_maps,const bool print_scores) {
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sort beforehand
  register uint64_t i = 0;
  GT_TEMPLATE_ITERATE(template,map_array) {
    if (i>=num_maps) break;
    if ((i++)>0) fprintf(file,GT_MAP_SNEXT);
    GT_MAP_ARRAY_ITERATE(map_array,map,end_position) {
      if (end_position>0) fprintf(file,GT_MAP_TEMPLATE_SEP);
      gt_output_map_fprint_map(file,map,print_scores);
    }
  }
  fprintf(file,"\n");
  return 0;
}
GT_INLINE gt_status gt_output_map_fprint_alignment_maps(
    FILE* file,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores) {
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sort beforehand
  register uint64_t i = 0;
  GT_MAPS_ITERATE(alignment,map) {
    if (i>=num_maps) break;
    if ((i++)>0) fprintf(file,GT_MAP_SNEXT);
    gt_output_map_fprint_map(file,map,print_scores);
  }
  return 0;
}

/*
 * High-level MAP Printers
 */
GT_INLINE gt_status gt_buffered_map_output_print_template(gt_buffered_map_output* const buffered_map_output,gt_template* const template) {
  //TODO
}
GT_INLINE gt_status gt_buffered_map_output_print_alignment(gt_buffered_map_output* const buffered_map_output,gt_alignment* const alignment) {
  //TODO
}
/* */
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer *output_buffer,gt_template* const template,const uint64_t num_maps) {
  //TODO
}
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer *output_buffer,gt_alignment* const alignment,const uint64_t num_maps) {
  //TODO
}
GT_INLINE gt_status gt_output_map_sprint_template(char **line_ptr,gt_template* const template,const uint64_t num_maps) {
  //TODO
}
GT_INLINE gt_status gt_output_map_sprint_alignment(char **line_ptr,gt_alignment* const alignment,const uint64_t num_maps) {
  //TODO
}

GT_INLINE gt_status gt_output_map_fprint_template(
    FILE* file,gt_template* const template,const uint64_t num_maps,const bool print_scores) {
  // Print TAG
  fprintf(file,"%s",gt_template_get_tag(template));
  // Print READ(s)
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i = 0;
  fprintf(file,"\t%s",gt_alignment_get_read(gt_template_get_block(template,i)));
  while (++i<num_blocks) {
    fprintf(file," %s",gt_alignment_get_read(gt_template_get_block(template,i)));
  }
  // Print QUALITY
  i = 0;
  if (gt_alignment_get_qualities(gt_template_get_block(template,i))!=NULL) {
    fprintf(file,"\t%s",gt_alignment_get_qualities(gt_template_get_block(template,i)));
    while (++i<num_blocks) {
      fprintf(file," %s",gt_alignment_get_qualities(gt_template_get_block(template,i)));
    }
  }
  // Print COUNTERS
  fprintf(file,"\t");
  gt_output_map_fprint_counters(file,gt_template_get_counters_vector(template),gt_template_get_mcs(template),false);
  // Print MAPS
  fprintf(file,"\t");
  gt_output_map_fprint_template_maps(file,template,num_maps,print_scores);
  return 0;
}
GT_INLINE gt_status gt_output_map_fprint_alignment(
    FILE* file,gt_alignment* const alignment,const uint64_t num_maps,const bool print_scores) {
  // Print TAG
  fprintf(file,"%s",gt_alignment_get_tag(alignment));
  // Print READ(s)
  fprintf(file,"\t%s",gt_alignment_get_read(alignment));
  // Print QUALITY
  if (gt_alignment_get_qualities(alignment)!=NULL) {
    fprintf(file,"\t%s",gt_alignment_get_qualities(alignment));
  }
  // Print COUNTERS
  fprintf(file,"\t");
  gt_output_map_fprint_counters(file,gt_alignment_get_counters_vector(alignment),gt_alignment_get_mcs(alignment),false);
  // Print MAPS
  fprintf(file,"\t");
  gt_output_map_fprint_alignment_maps(file,alignment,num_maps,print_scores);
  return 0;
}








