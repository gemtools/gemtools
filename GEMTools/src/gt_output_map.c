/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_map.h"

#define GT_OUTPUT_MAP_COMPACT_COUNTERS_ZEROS_TH 5

GT_INLINE gt_output_map_attributes* gt_output_map_attributes_new() {
	gt_output_map_attributes* attr = malloc(sizeof(gt_output_map_attributes));
	gt_cond_fatal_error(!attr,MEM_HANDLER);
	gt_output_map_attributes_reset_defaults(attr);
	return attr;
}
GT_INLINE void gt_output_map_attributes_delete(gt_output_map_attributes* attributes) {
  GT_NULL_CHECK(attributes);
	free(attributes);
}
GT_INLINE void gt_output_map_attributes_reset_defaults(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  /* TAG */
	attributes->print_extra = true;
	attributes->print_casava = true;
	/* COUNTERS */
	attributes->compact = false;
	/* MAPS */
	attributes->print_scores = true;
	attributes->max_printable_maps = GT_ALL;
}

GT_INLINE bool gt_output_map_attributes_is_print_scores(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->print_scores;
}
GT_INLINE void gt_output_map_attributes_set_print_scores(gt_output_map_attributes* const attributes,const bool print_scores) {
  GT_NULL_CHECK(attributes);
	attributes->print_scores = print_scores;
}
GT_INLINE bool gt_output_map_attributes_is_print_extra(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->print_extra;
}
GT_INLINE void gt_output_map_attributes_set_print_extra(gt_output_map_attributes* const attributes,const bool print_extra) {
  GT_NULL_CHECK(attributes);
	attributes->print_extra = print_extra;
}
GT_INLINE bool gt_output_map_attributes_is_print_casava(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->print_casava;
}
GT_INLINE void gt_output_map_attributes_set_print_casava(gt_output_map_attributes* const attributes,const bool print_casava) {
  GT_NULL_CHECK(attributes);
	attributes->print_casava = print_casava;
}
GT_INLINE uint64_t gt_output_map_attributes_get_max_printable_maps(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->max_printable_maps;
}
GT_INLINE void gt_output_map_attributes_set_max_printable_maps(gt_output_map_attributes* const attributes,const uint64_t max_printable_maps) {
  GT_NULL_CHECK(attributes);
	attributes->max_printable_maps = max_printable_maps;
}
/*
 * TAG building block printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,attributes,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_tag,gt_string* const tag,gt_shash* const attributes,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_tag(
    gt_generic_printer* const gprinter,gt_string* const tag,gt_shash* const attributes,gt_output_map_attributes* const output_map_attributes) {
  // PRIgts needed as this calls gt_string_get_string downstream, which returns the full buffer not trimmed to length
  gt_gprintf(gprinter,PRIgts,PRIgts_content(tag));
  // Check if we have casava attributes
  if (gt_output_map_attributes_is_print_casava(output_map_attributes) && gt_shash_is_contained(attributes,GT_TAG_CASAVA)) {
    // Print casava
    gt_gprintf(gprinter," "PRIgts,PRIgts_content(gt_shash_get(attributes,GT_TAG_CASAVA,gt_string)));
  } else {
    // Append /1 /2 if paired
    if (gt_shash_is_contained(attributes,GT_TAG_PAIR)) {
        int64_t p = *gt_shash_get(attributes,GT_TAG_PAIR,int64_t);
      if (p > 0) {
        gt_gprintf(gprinter,"/%d",p);
      }
    }
  }
  if(gt_output_map_attributes_is_print_extra(output_map_attributes) && gt_shash_is_contained(attributes,GT_TAG_EXTRA)) {
    // Print additional
    gt_gprintf(gprinter," "PRIgts, PRIgts_content(gt_shash_get(attributes,GT_TAG_EXTRA,gt_string)));
  }
  return 0;
}
/*
 * Internal MAP printers (take parameters as to control flow/format options)
 */
GT_INLINE gt_status gt_output_map_gprint_mismatch_string_(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes,
    const bool begin_trim,const bool end_trim) {
  GT_NULL_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
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
GT_INLINE gt_status gt_output_map_gprint_map_block_(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes,
    const bool begin_trim,const bool end_trim) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  /*
   * FORMAT => chr11:-:51590050:(5)43T46A9>24*
   */
  // Print sequence name
  gt_gprintf(gprinter,PRIgts,PRIgts_content(map->seq_name));
  // Print strand
  gt_gprintf(gprinter,GT_MAP_SEP_S"%c",gt_map_get_strand(map)==FORWARD?GT_MAP_STRAND_FORWARD_SYMBOL:GT_MAP_STRAND_REVERSE_SYMBOL);
  // Print position
  gt_gprintf(gprinter,GT_MAP_SEP_S"%"PRIu64 GT_MAP_SEP_S,gt_map_get_position(map));
  // Print CIGAR
  gt_output_map_gprint_mismatch_string_(gprinter,map,output_map_attributes,begin_trim,end_trim);
  return 0;
}
GT_INLINE gt_status gt_output_map_gprint_map_(
    gt_generic_printer* const gprinter,gt_map* const map,
    gt_output_map_attributes* const output_map_attributes,const bool print_scores,const bool print_trims) {
  GT_NULL_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  /*
   * FORMAT => chr11:-:51590050:(5)43T46A9>24*
   */
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
    gt_output_map_gprint_mismatch_string_(gprinter,map_it,output_map_attributes,next_map==NULL,!has_next_block);
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
    gt_output_map_gprint_map_(gprinter,next_map,output_map_attributes,print_scores,false);
  }
  return 0;
}
GT_INLINE gt_status gt_output_map_gprint_counters_(
    gt_generic_printer* const gprinter,gt_vector* const counters,gt_output_map_attributes* const output_map_attributes,
    const uint64_t max_complete_strata,const bool not_unique_flag) {
  register const uint64_t num_counters = gt_vector_get_used(counters);
  register uint64_t i;
  // Not unique
  if (not_unique_flag) {
    gt_gprintf(gprinter,GT_MAP_COUNTS_NOT_UNIQUE_S);
    return 0;
  }
  // No counters
  if (num_counters==0) {
    gt_gprintf(gprinter,"0");
    return 0;
  }
  // Print all counters
  for (i=0;i<num_counters;) {
    if (i>0) gt_gprintf(gprinter,"%c",gt_expect_false(i==max_complete_strata)?GT_MAP_MCS:GT_MAP_COUNTS_SEP);
    register const uint64_t counter = *gt_vector_get_elm(counters,i,uint64_t);
    if (gt_expect_false(output_map_attributes->compact && counter==0)) {
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
  // MCS (zeros)
  if (max_complete_strata < UINT64_MAX) {
    for (;i<max_complete_strata;++i) {
      if (i>0) {
        gt_gprintf(gprinter,"%c0",GT_MAP_COUNTS_SEP);
      } else {
        gt_gprintf(gprinter,"0");
      }
    }
  }
  return 0;
}
/*
 * MAP building block printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS counters
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_counters,gt_vector* const counters);
GT_INLINE gt_status gt_output_map_gprint_counters(gt_generic_printer* const gprinter,gt_vector* const counters) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(counters);
  gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_counters_(gprinter,counters,&output_map_attributes,UINT64_MAX,false);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_string,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_mismatch_string(gt_generic_printer* const gprinter,gt_map* const map) {
  GT_NULL_CHECK(gprinter);
  GT_MAP_CHECK(map);
  gt_output_map_attributes map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_mismatch_string_g(gprinter,map,&map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_block,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_map_block(gt_generic_printer* const gprinter,gt_map* const map) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  gt_output_map_attributes map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_map_block_g(gprinter,map,&map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_map(gt_generic_printer* const gprinter,gt_map* const map) {
  GT_NULL_CHECK(gprinter); GT_MAP_CHECK(map);
  gt_output_map_attributes map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_map_g(gprinter,map,&map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS maps
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_list,gt_vector* const maps);
GT_INLINE gt_status gt_output_map_gprint_map_list(gt_generic_printer* const gprinter,gt_vector* const maps) {
  GT_NULL_CHECK(gprinter); GT_VECTOR_CHECK(maps);
  gt_output_map_attributes map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_map_list_g(gprinter,maps,&map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps,gt_template* const template);
GT_INLINE gt_status gt_output_map_gprint_template_maps(gt_generic_printer* const gprinter,gt_template* const template) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_template_maps_g(gprinter,template,&output_map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps,gt_alignment* const alignment);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps(gt_generic_printer* const gprinter,gt_alignment* const alignment) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_alignment_maps_g(gprinter,alignment,&output_map_attributes);
}
/*
 * General MAP printers (more flexible taking @gt_output_map_attributes)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_string_g,gt_map* const map,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_mismatch_string_g(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  return gt_output_map_gprint_mismatch_string_(gprinter,map,output_map_attributes,true,true);
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS counters,attributes,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_counters_g,
    gt_vector* const counters,gt_shash* const attributes,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_counters_g(gt_generic_printer* const gprinter,gt_vector* const counters,
    gt_shash* const attributes,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(counters);
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(output_map_attributes);
  register uint64_t* const mcs_ptr = (uint64_t*)gt_attribute_get(attributes,GT_ATTR_MAX_COMPLETE_STRATA);
  register bool* const not_unique_flag = (bool*)gt_attribute_get(attributes,GT_ATTR_NOT_UNIQUE);
  return gt_output_map_gprint_counters_(gprinter,counters,output_map_attributes,
      ((mcs_ptr!=NULL) ? *mcs_ptr : UINT64_MAX),
      ((not_unique_flag!=NULL) ? *not_unique_flag : false));
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_block_g,gt_map* const map,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_block_g(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  return gt_output_map_gprint_map_block_(gprinter,map,output_map_attributes,true,true);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_g,gt_map* const map,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_g(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  return gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS maps,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_list_g,gt_vector* const maps,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_list_g(
    gt_generic_printer* const gprinter,gt_vector* const maps,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(maps);
  GT_NULL_CHECK(output_map_attributes);
  register gt_status error_code = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    error_code |= gt_output_map_gprint_map_(gprinter,*map,output_map_attributes,output_map_attributes->print_scores,true);
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps_g,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_template_maps_g(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_NULL_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_map_attributes);
  register gt_status error_code = 0;
  if (gt_expect_false(gt_template_get_num_mmaps(template)==0 || output_map_attributes->max_printable_maps==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
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
          gt_output_map_gprint_map_(gprinter,map,output_map_attributes,false,true);
        }
        if (output_map_attributes->print_scores && map_array_attr!=NULL && map_array_attr->score!=GT_MAP_NO_SCORE) {
          gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,map_array_attr->score);
        }
        if (total_maps_printed>=output_map_attributes->max_printable_maps || total_maps_printed>=num_maps) return 0;
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
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps_g,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps_g(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_map_attributes);
  register gt_status error_code = 0;
  if (gt_expect_false(gt_alignment_get_num_maps(alignment)==0 || output_map_attributes->max_printable_maps==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
  } else {
    register const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
    uint64_t strata = 0, pending_maps = 0, total_maps_printed = 0;
    while (gt_alignment_get_next_matching_strata(alignment,strata,&strata,&pending_maps)) {
      GT_ALIGNMENT_ITERATE(alignment,map) {
        if (gt_map_get_global_distance(map)!=strata) continue;
        // Print map
        --pending_maps;
        if ((total_maps_printed++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
        gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true);
        if (total_maps_printed>=output_map_attributes->max_printable_maps || total_maps_printed>=num_maps) return 0;
        if (pending_maps==0) break;
      }
      if (pending_maps > 0) {
        error_code = GT_MOE_INCONSISTENT_COUNTERS;
        gt_error(ALIGNMENT_INCONSISTENT_COUNTERS);
      }
      ++strata;
    }
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps_unsorted,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_template_maps_unsorted(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_map_attributes);
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sorting must be done beforehand
  if (gt_expect_false(gt_template_get_num_mmaps(template)==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
  } else {
    register uint64_t i = 0;
    GT_TEMPLATE__ATTR_ITERATE(template,map_array,map_array_attr) {
      if (i>=output_map_attributes->max_printable_maps) break;
      if ((i++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
      GT_MULTIMAP_ITERATE(map_array,map,end_position) {
        if (end_position>0) gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
        gt_output_map_gprint_map_(gprinter,map,output_map_attributes,false,true);
        if (output_map_attributes->print_scores && map_array_attr!=NULL && map_array_attr->score!=GT_MAP_NO_SCORE) {
          gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,map_array_attr->score);
        }
      }
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps_unsorted,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps_unsorted(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_map_attributes);
  // NOTE: No sorting performed. Written as laid in the vector.
  //       Thus, if you want a particular sorting (by score, by distance, ...) sort beforehand
  if (gt_expect_false(gt_alignment_get_num_maps(alignment)==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
  } else {
    register uint64_t i = 0;
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if (i>=output_map_attributes->max_printable_maps) break;
      if ((i++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
      gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true);
    }
  }
  return 0;
}
/*
 * High-level MAP Printers {Alignment/Template}
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_template(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_map_attributes);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_output_map_gprint_alignment(gprinter,alignment,output_map_attributes);
  } GT_TEMPLATE_END_REDUCTION;
  register gt_status error_code;
  // Print TAG
  gt_output_map_gprint_tag(gprinter,template->tag,template->attributes,output_map_attributes);
  // Print READ(s)
  register const uint64_t num_blocks = gt_template_get_num_blocks(template);
  register uint64_t i = 0;
  gt_gprintf(gprinter,"\t%s",gt_alignment_get_read(gt_template_get_block(template,i)));
  while (++i<num_blocks) {
    gt_gprintf(gprinter," %s",gt_alignment_get_read(gt_template_get_block(template,i)));
  }
  // Print QUALITY
  gt_gprintf(gprinter,"\t"); i = 0;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    if (gt_alignment_has_qualities(alignment)) {
      if (i > 0) {
        gt_gprintf(gprinter," %s",gt_alignment_get_qualities(alignment));
      } else {
        gt_gprintf(gprinter,"%s",gt_alignment_get_qualities(alignment));
      }
    } else if (i > 0) {
      gt_gprintf(gprinter," ");
    }
    ++i;
  }
  // Print COUNTERS
  gt_gprintf(gprinter,"\t");
  gt_output_map_gprint_counters_(gprinter,gt_template_get_counters_vector(template),
      output_map_attributes,gt_template_get_mcs(template),gt_template_get_not_unique_flag(template));
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  error_code = gt_output_map_gprint_template_maps_g(gprinter,template,output_map_attributes);
  gt_gprintf(gprinter,"\n");
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment,gt_alignment* const alignment, gt_output_map_attributes*  const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment(
    gt_generic_printer* const gprinter,gt_alignment* const alignment, gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_map_attributes);
  register gt_status error_code;
  // Print TAG
  gt_output_map_gprint_tag(gprinter,alignment->tag,alignment->attributes,output_map_attributes);
  // Print READ(s)
  gt_gprintf(gprinter,"\t%s",gt_alignment_get_read(alignment));
  // Print QUALITY
  if (gt_alignment_has_qualities(alignment)) {
    gt_gprintf(gprinter,"\t%s",gt_alignment_get_qualities(alignment));
  } else {
    gt_gprintf(gprinter,"\t");
  }
  // Print COUNTERS
  gt_gprintf(gprinter,"\t");
  gt_output_map_gprint_counters_(gprinter,gt_alignment_get_counters_vector(alignment),
        output_map_attributes,gt_alignment_get_mcs(alignment),gt_alignment_get_not_unique_flag(alignment));
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  error_code = gt_output_map_gprint_alignment_maps_g(gprinter,alignment,output_map_attributes);
  gt_gprintf(gprinter,"\n");
  return error_code;
}
/*
 * GEM printer
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_gem_template,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_gem_template(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  if (gt_template_get_num_mmaps(template)>0) {
    return gt_output_map_gprint_template(gprinter,template,output_map_attributes);
  } else {
    register gt_status error_code = 0;
    GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
      if ((error_code=gt_output_map_gprint_alignment(gprinter,alignment,output_map_attributes))) return error_code;
    }
    return error_code;
  }
}
/*
 * Misc. Handy printers
 */
GT_INLINE gt_status gt_output_map_gprint_mismatch_summary_(gt_generic_printer* const gprinter,gt_map* const map,const uint64_t block_num) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  gt_gprintf(gprinter,"{Block %lu}\n",block_num);
  gt_gprintf(gprinter,"  --> GEM-Distance %lu\n",gt_map_get_distance(map));
  gt_gprintf(gprinter,"  --> Levenshtein-Distance %lu\n",gt_map_get_levenshtein_distance(map));
  gt_gprintf(gprinter,"  --> Mismatches %lu\n",gt_map_get_num_mismatch_bases(map));
  gt_gprintf(gprinter,"  --> Indels %lu\n",gt_map_get_num_indels(map));
  gt_gprintf(gprinter,"    --> Deletions %lu\n",gt_map_get_num_deletions(map));
  gt_gprintf(gprinter,"    --> Insertions %lu\n",gt_map_get_num_insertions(map));
  gt_gprintf(gprinter,"  --> Mismatch.list::",block_num);
  GT_MISMS_ITERATE(map,misms) {
    switch (misms->misms_type) {
      case MISMS:
        gt_gprintf(gprinter,"  MIS.%02lu.%c",misms->position,misms->base);
        break;
      case DEL:
        gt_gprintf(gprinter,"  DEL.%02lu<-%02lu>",misms->position,misms->size);
        break;
      case INS:
        gt_gprintf(gprinter,"  INS.%02lu<+%02lu>",misms->position,misms->size);
        break;
    }
  }
  gt_gprintf(gprinter,"\n");
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_summary,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_mismatch_summary(gt_generic_printer* const gprinter,gt_map* const map) {
  gt_gprintf(gprinter,"[Mismatch/Indel Summary]\n");
  gt_gprintf(gprinter,"  --> Total.Blocks %lu\n",gt_map_get_num_blocks(map));
  register uint64_t block_pos = 0;
  GT_MAP_ITERATE(map,map_block) {
    gt_output_map_gprint_mismatch_summary_(gprinter,map_block,block_pos);
    ++block_pos;
  }
  gt_gprintf(gprinter,"  --> Stringify \t ");
  gt_output_map_gprint_map(gprinter,map);
  gt_gprintf(gprinter,"\n");
  return 0;
}
#define GT_OUTPUT_MAP_RELOAD_MISMS(map,misms_offset,misms_ptr) { \
  if (misms_offset<gt_map_get_num_misms(map)) { \
    misms_ptr=gt_map_get_misms(map,misms_offset); \
  } else { \
    misms_ptr=NULL; \
  } \
}
GT_INLINE gt_status gt_output_map_gprint_pretty_alignment_(
    gt_generic_printer* const gprinter,gt_map* const map,
    char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,const uint64_t block_num) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  /*
   * Print the alignment (in short)
   */
  gt_gprintf(gprinter,"#%lu{",block_num);
  gt_output_map_gprint_map_block(gprinter,map);
  gt_gprintf(gprinter,"}\n");
  /*
   * Construct the schemes
   */
  register const uint64_t approx_bound_length = pattern_length+sequence_length;
  register gt_string* const pattern_scheme = gt_string_new(approx_bound_length);
  register gt_string* const sequence_scheme = gt_string_new(approx_bound_length);
  register gt_string* const relation_scheme = gt_string_new(approx_bound_length);
  register uint64_t pattern_centinel=0, sequence_centinel=0;
  // Init misms
  register uint64_t misms_offset=0, i;
  register gt_misms* misms;
  GT_OUTPUT_MAP_RELOAD_MISMS(map,misms_offset,misms);
  // Some printing artifacts
  for (i=0;i<4;++i) {
    gt_string_append_char(sequence_scheme,'-');
    gt_string_append_char(relation_scheme,' ');
    gt_string_append_char(pattern_scheme,' ');
  }
  // Traverse the sequence
  register bool exception = false, go_on = true;
  while (go_on && (pattern_centinel<pattern_length || sequence_centinel<sequence_length)) {
    if (misms!=NULL && misms->position==pattern_centinel) { // Misms
      switch (misms->misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
            exception=true; go_on=false; break;
          }
          gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
          gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) {
            gt_string_append_char(relation_scheme,'O'); exception=true;
          } else if (misms->base!=sequence[sequence_centinel]) {
            gt_string_append_char(relation_scheme,'#'); exception=true;
          } else {
            gt_string_append_char(relation_scheme,'X');
          }
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms->size>sequence_length) {
            exception=true; go_on=false; break;
          }
          for (i=0;i<misms->size;++i) {
            gt_string_append_char(sequence_scheme,sequence[sequence_centinel++]);
            gt_string_append_char(relation_scheme,'-');
            gt_string_append_char(pattern_scheme,' ');
          }
          break;
        case DEL:
          if (pattern_centinel+misms->size>pattern_length) {
            exception=true; go_on=false; break;
          }
          for (i=0;i<misms->size;++i) {
            gt_string_append_char(sequence_scheme,' ');
            gt_string_append_char(relation_scheme,'-');
            gt_string_append_char(pattern_scheme,pattern[pattern_centinel++]);
          }
          break;
      }
      ++misms_offset;
      GT_OUTPUT_MAP_RELOAD_MISMS(map,misms_offset,misms);
    } else { // Match
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
        exception=true; break;
      }
      gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
      gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) {
        gt_string_append_char(relation_scheme,'*');
      } else {
        gt_string_append_char(relation_scheme,'|');
      }
      ++pattern_centinel;
      ++sequence_centinel;
    }
  }
  // Fill the rest (in case of exceptions)
  while (pattern_centinel<pattern_length || sequence_centinel<sequence_length) {
    if (sequence_centinel<sequence_length) {
      gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
    } else {
      gt_string_append_char(sequence_scheme,'?');
    }
    gt_string_append_char(relation_scheme,'!');
    if (pattern_centinel<pattern_length) {
      gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
    } else {
      gt_string_append_char(pattern_scheme,'?');
    }
    ++pattern_centinel; ++sequence_centinel;
  }
  // Some printing artifacts
  for (i=0;i<4;++i) {
    gt_string_append_char(sequence_scheme,'-');
    gt_string_append_char(relation_scheme,' ');
    gt_string_append_char(pattern_scheme,' ');
  }
  // Append EOS
  gt_string_append_eos(sequence_scheme);
  gt_string_append_eos(relation_scheme);
  gt_string_append_eos(pattern_scheme);
  /*
   * Print the alignment (pretty)
   */
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(sequence_scheme));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(relation_scheme));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(pattern_scheme));
  // Free
  gt_string_delete(sequence_scheme);
  gt_string_delete(relation_scheme);
  gt_string_delete(pattern_scheme);
  return exception?1:0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,pattern,pattern_length,sequence,sequence_length
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_pretty_alignment,
    gt_map* const map,char* const pattern,const uint64_t pattern_length,char* const sequence,const uint64_t sequence_length);
GT_INLINE gt_status gt_output_map_gprint_pretty_alignment(
    gt_generic_printer* const gprinter,gt_map* const map,
    char* const pattern,const uint64_t pattern_length,char* const sequence,const uint64_t sequence_length) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  // Print map (short)
  gt_gprintf(gprinter,"[");
  gt_output_map_gprint_map(gprinter,map);
  gt_gprintf(gprinter,"] TotalBlocks=%lu\n",gt_map_get_num_blocks(map));
  // Print all the blocks (pretty)
  register uint64_t block_pos = 0;
  GT_MAP_ITERATE(map,map_block) {
    gt_output_map_gprint_pretty_alignment_(gprinter,map,
        pattern,pattern_length,sequence,sequence_length,block_pos++);
  }
  return 0;
}

