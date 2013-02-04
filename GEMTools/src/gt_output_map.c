/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_map.h"

#define GT_OUTPUT_MAP_COMPACT_COUNTERS_ZEROS_TH 5

/*
 * MAP building block printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS counters,max_complete_strata,compact
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_counters,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact);
GT_INLINE gt_status gt_output_map_gprint_counters(
    gt_generic_printer* const gprinter,gt_vector* const counters,const uint64_t max_complete_strata,const bool compact) {
  GT_NULL_CHECK(gprinter); GT_NULL_CHECK(counters);
  register const uint64_t num_counters = gt_vector_get_used(counters);
  register uint64_t i;
  if (num_counters==0) {
    gt_gprintf(gprinter,"0");
    return 0;
  }
  for (i=0;i<num_counters;) {
    if (i>0) gt_gprintf(gprinter,"%c",gt_expect_false(i==max_complete_strata)?GT_MAP_MCS:GT_MAP_COUNTS_SEP);
    register const uint64_t counter = *gt_vector_get_elm(counters,i,uint64_t);
    if (gt_expect_false(compact && counter==0)) {
      register uint64_t j=i+1;
      while (j<num_counters && *gt_vector_get_elm(counters,j,uint64_t)==0) ++j;
      if (gt_expect_false((j-i)>=GT_OUTPUT_MAP_COMPACT_COUNTERS_ZEROS_TH)) {
        gt_gprintf(gprinter,"0" GT_MAP_COUNTS_TIMES_S "%"PRIu64,(j-i)); i=j;
      } else {
        gt_gprintf(gprinter,"0"); ++i;
      }
    } else {
      gt_gprintf(gprinter,"%"PRIu64,counter); ++i;
    }
  }
  return 0;
}

GT_INLINE gt_status gt_output_map_gprint_mismatch_string_(
    gt_generic_printer* const gprinter,gt_map* const map,const bool begin_trim,const bool end_trim) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  register const uint64_t map_length = gt_map_get_base_length(map);
  register uint64_t centinel = 0;
  GT_MISMS_ITERATE(map,misms) {
    register const uint64_t misms_pos = gt_misms_get_position(misms);
    if (misms_pos!=centinel) {
      gt_gprintf(gprinter,"%"PRIu64,misms_pos-centinel);
      centinel = misms_pos;
    }
    switch (gt_misms_get_type(misms)) {
      case MISMS:
        gt_gprintf(gprinter,"%c",gt_misms_get_base(misms));
        centinel=misms_pos+1;
        break;
      case INS:
        gt_gprintf(gprinter,">%"PRIu64"+",gt_misms_get_size(misms));
        break;
      case DEL: {
        register const uint64_t init_centinel = centinel;
        centinel+=gt_misms_get_size(misms);
        if (gt_expect_false((init_centinel==0 && begin_trim) || (centinel==map_length && end_trim))) { // Trim
          gt_gprintf(gprinter,"(%"PRIu64")",gt_misms_get_size(misms));
        } else {
          gt_gprintf(gprinter,">%"PRIu64"-",gt_misms_get_size(misms));
        }
        break;
      }
      default:
        gt_error(SELECTION_NOT_VALID);
        break;
    }
  }
  if (centinel < map_length) {
    gt_gprintf(gprinter,"%"PRIu64,map_length-centinel);
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_string,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_mismatch_string(gt_generic_printer* const gprinter,gt_map* const map) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  return gt_output_map_gprint_mismatch_string_(gprinter,map,true,true);
}

GT_INLINE gt_status gt_output_map_gprint_map_(gt_generic_printer* const gprinter,gt_map* const map,const bool print_scores,const bool print_trims) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  // FORMAT => chr11:-:51590050:(5)43T46A9>24*
  // Print sequence name
  gt_gprintf(gprinter,PRIgts,PRIgts_content(map->seq_name));
  // Print strand
  gt_gprintf(gprinter,GT_MAP_SEP_S"%c",gt_map_get_strand(map)==FORWARD?GT_MAP_STRAND_FORWARD_SYMBOL:GT_MAP_STRAND_REVERSE_SYMBOL);
  // Print position
  gt_gprintf(gprinter,GT_MAP_SEP_S"%"PRIu64 GT_MAP_SEP_S,gt_map_get_position(map));
  // Print mismatch string (compact it)
  register gt_map* map_it = map, *next_map=NULL;
  register bool cigar_pending = true;
  while (cigar_pending) {
    register const bool has_next_block = gt_map_has_next_block(map_it);
    gt_output_map_gprint_mismatch_string_(gprinter,map_it,next_map==NULL,!has_next_block);
    if (has_next_block) {
      next_map = gt_map_get_next_block(map_it);
      if ((cigar_pending=(gt_string_equals(map_it->seq_name,next_map->seq_name)))) {
        switch (gt_map_get_junction(map_it)) {
          case SPLICE:
            gt_gprintf(gprinter,">""%"PRIu64"*",gt_map_get_junction_distance(map_it));
            break;
          case POSITIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRIu64"+",gt_map_get_junction_distance(map_it));
            break;
          case NEGATIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRIu64"-",gt_map_get_junction_distance(map_it));
            break;
          case INSERT:
            cigar_pending=false;
            break;
          case NO_JUNCTION:
            break;
          default:
            gt_error(SELECTION_NOT_VALID);
            break;
        }
        map_it = next_map;
      }
    } else {
      cigar_pending = false;
    }
  }
  // Print attributes (scores)
  if (print_scores && gt_map_get_global_score(map)!=GT_MAP_NO_SCORE) {
    gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,gt_map_get_global_score(map));
  }
  // Print possible next blocks (out of the current sequence => split-maps across chromosomes)
  if (gt_map_has_next_block(map_it)) { // FIXME: trimmings, do this really occurs?
    gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
    gt_output_map_gprint_map_(gprinter,next_map,print_scores,false);
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map,gt_map* const map,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_map(gt_generic_printer* const gprinter,gt_map* const map,const bool print_scores) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  return gt_output_map_gprint_map_(gprinter,map,print_scores,true);
}

/*
 * Print maps unsorted
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_template_maps(
    gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores) {
  GT_NULL_CHECK(gprinter); GT_TEMPLATE_CHECK(template);
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sorting must be done beforehand
  register uint64_t i = 0;
  if (gt_expect_false(gt_template_get_num_mmaps(template)==0)) {
    gt_gprintf(gprinter,"-");
  } else {
    GT_TEMPLATE__ATTR_ITERATE(template,map_array,map_array_attr) {
      if (i>=max_printable_maps) break;
      if ((i++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
      GT_MULTIMAP_ITERATE(map_array,map,end_position) {
        if (end_position>0) gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
        gt_output_map_gprint_map_(gprinter,map,print_scores,true);
        if (print_scores && map_array_attr!=NULL && map_array_attr->score!=GT_MAP_NO_SCORE) {
          gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,map_array_attr->score);
        }
      }
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores) {
  GT_NULL_CHECK(gprinter); GT_ALIGNMENT_CHECK(alignment);
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sort beforehand
  register uint64_t i = 0;
  if (gt_expect_false(gt_alignment_get_num_maps(alignment)==0)) {
    gt_gprintf(gprinter,"-");
  } else {
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if (i>=max_printable_maps) break;
      if ((i++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
      gt_output_map_gprint_map_(gprinter,map,print_scores,true);
    }
  }
  return 0;
}
/*
 * Print map sorted by distance
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps_sorted,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_template_maps_sorted(
    gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores) {
  GT_NULL_CHECK(gprinter); GT_TEMPLATE_CHECK(template);
  register gt_status error_code = 0;
  if (gt_expect_false(gt_template_get_num_mmaps(template)==0 || max_printable_maps==0)) {
    gt_gprintf(gprinter,"-");
  } else {
    register const uint64_t num_maps = gt_template_get_num_mmaps(template);
    uint64_t strata = 0, pending_maps = 0, total_maps_printed = 0;
    while (gt_template_get_next_matching_strata(template,strata,&strata,&pending_maps)) {
      GT_TEMPLATE__ATTR_ITERATE(template,map_array,map_array_attr) {
        if (map_array_attr->distance!=strata) continue;
        // Print mmap
        --pending_maps;
        if ((total_maps_printed++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
        GT_MULTIMAP_ITERATE(map_array,map,end_position) {
          if (end_position>0) gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
          gt_output_map_gprint_map_(gprinter,map,false,true);
        }
        if (print_scores && map_array_attr!=NULL && map_array_attr->score!=GT_MAP_NO_SCORE) {
          gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,map_array_attr->score);
        }
        if (total_maps_printed>=max_printable_maps || total_maps_printed>=num_maps) return 0;
        if (pending_maps==0) break;
      }
      if (pending_maps>0) {
        error_code = GT_MOE_INCONSISTENT_COUNTERS;
        gt_error(TEMPLATE_INCONSISTENT_COUNTERS);
      }
      ++strata;
    }
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps_sorted,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps_sorted(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores) {
  GT_NULL_CHECK(gprinter); GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code = 0;
  if (gt_expect_false(gt_alignment_get_num_maps(alignment)==0 || max_printable_maps==0)) {
    gt_gprintf(gprinter,"-");
  } else {
    register const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
    uint64_t strata = 0, pending_maps = 0, total_maps_printed = 0;
    while (gt_alignment_get_next_matching_strata(alignment,strata,&strata,&pending_maps)) {
      GT_ALIGNMENT_ITERATE(alignment,map) {
        if (gt_map_get_global_distance(map)!=strata) continue;
        // Print map
        --pending_maps;
        if ((total_maps_printed++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
        gt_output_map_gprint_map_(gprinter,map,print_scores,true);
        if (total_maps_printed>=max_printable_maps || total_maps_printed>=num_maps) return 0;
        if (pending_maps==0) break;
      }
      if (pending_maps>0) {
        error_code = GT_MOE_INCONSISTENT_COUNTERS;
        gt_error(ALIGNMENT_INCONSISTENT_COUNTERS);
      }
      ++strata;
    }
  }
  return error_code;
}

/*
 * Specific High-level MAP Printers {Alignment/Template}
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_template(
    gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  // Print TAG
  gt_gprintf(gprinter,"%s",gt_template_get_tag(template));
  // Print READ(s)
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i = 0;
  gt_gprintf(gprinter,"\t%s",gt_alignment_get_read(gt_template_get_block(template,i)));
  while (++i<num_blocks) {
    gt_gprintf(gprinter," %s",gt_alignment_get_read(gt_template_get_block(template,i)));
  }
  // Print QUALITY
  i = 0;
  if (gt_alignment_has_qualities(gt_template_get_block(template,i))) {
    gt_gprintf(gprinter,"\t%s",gt_alignment_get_qualities(gt_template_get_block(template,i)));
    while (++i<num_blocks) {
      gt_gprintf(gprinter," %s",gt_alignment_get_qualities(gt_template_get_block(template,i)));
    }
  }
  // Print COUNTERS
  if (gt_expect_false(gt_template_get_not_unique_flag(template))) {
    gt_gprintf(gprinter,"\t" GT_MAP_COUNTS_NOT_UNIQUE_S);
  } else {
    gt_gprintf(gprinter,"\t");
    gt_output_map_gprint_counters(gprinter,gt_template_get_counters_vector(template),gt_template_get_mcs(template),false);
  }
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    error_code = gt_output_map_gprint_alignment_maps_sorted(gprinter,alignment,max_printable_maps,print_scores); // _sorted
    gt_gprintf(gprinter,"\n");
    return error_code;
  } GT_TEMPLATE_END_REDUCTION;
  error_code = gt_output_map_gprint_template_maps_sorted(gprinter,template,max_printable_maps,print_scores); // _sorted
  gt_gprintf(gprinter,"\n");
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_alignment(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,const uint64_t max_printable_maps,const bool print_scores) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  register gt_status error_code;
  // Print TAG
  gt_gprintf(gprinter,"%s",gt_alignment_get_tag(alignment));
  // Print READ(s)
  gt_gprintf(gprinter,"\t%s",gt_alignment_get_read(alignment));
  // Print QUALITY
  if (gt_alignment_has_qualities(alignment)) {
    gt_gprintf(gprinter,"\t%s",gt_alignment_get_qualities(alignment));
  }
  // Print COUNTERS
  if (gt_expect_false(gt_alignment_get_not_unique_flag(alignment))) {
    gt_gprintf(gprinter,"\t"GT_MAP_COUNTS_NOT_UNIQUE_S);
  } else {
    gt_gprintf(gprinter,"\t");
    gt_output_map_gprint_counters(gprinter,gt_alignment_get_counters_vector(alignment),gt_alignment_get_mcs(alignment),false);
  }
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  error_code = gt_output_map_gprint_alignment_maps_sorted(gprinter,alignment,max_printable_maps,print_scores); // _sorted
  gt_gprintf(gprinter,"\n");
  return error_code;
}


/*
 * GEM printer
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,max_printable_maps,print_scores
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_gem_template,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores);
GT_INLINE gt_status gt_output_map_gprint_gem_template(
    gt_generic_printer* const gprinter,gt_template* const template,const uint64_t max_printable_maps,const bool print_scores) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  if (gt_template_get_num_mmaps(template)>0) {
    return gt_output_map_gprint_template(gprinter,template,max_printable_maps,print_scores);
  } else {
    register gt_status error_code = 0;
    GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
      if ((error_code=gt_output_map_gprint_alignment(gprinter,alignment,max_printable_maps,print_scores))) return error_code;
    }
    return error_code;
  }
}


