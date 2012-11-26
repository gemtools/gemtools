/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.h
 * DATE: 01/06/2012
 * DESCRIPTION: Data structure modeling sequences' templates.
 *   That is, set of alignments and relationships between their maps.
 */

#ifndef GT_TEMPLATE_H_
#define GT_TEMPLATE_H_

#include "gt_commons.h"
#include "gt_alignment.h"

// Codes gt_status
#define GT_TEMPLATE_OK 1
#define GT_TEMPLATE_FAIL 0

typedef struct {
  uint32_t template_id;
  uint32_t in_block_id;
  gt_string* tag;
  gt_vector* blocks; /* (gt_alignment*) */ /* paired::blocks->used==2 */
  gt_vector* counters; /* (uint64_t) */
  gt_vector* mmaps; /* (gt_map*) */
  gt_vector* mmaps_attributes; /* ( (gt_mmap_attributes) ) */
  char* maps_txt;
  gt_shash* attributes;
} gt_template;

typedef struct {
  uint64_t distance;
  int64_t score;
} gt_mmap_attributes;

// Iterators
typedef struct {
  gt_template* template;
  uint64_t next_pos;
} gt_template_alignment_iterator;
typedef struct {
  gt_template* template;
  uint64_t next_pos;
} gt_template_maps_iterator;

/*
 * Checkers & Errors
 */
#define GT_TEMPLATE_CHECK(template) \
  gt_fatal_check(template==NULL,NULL_HANDLER); \
  GT_STRING_CHECK(template->tag); \
  GT_VECTOR_CHECK(template->blocks); \
  GT_VECTOR_CHECK(template->counters); \
  GT_VECTOR_CHECK(template->mmaps); \
  GT_VECTOR_CHECK(template->mmaps_attributes); \
  GT_HASH_CHECK(template->attributes)

#define GT_TEMPLATE_CONSISTENCY_CHECK(template) \
  GT_TEMPLATE_CHECK(template); \
  gt_fatal_check(gt_vector_get_used(template->blocks)==0,TEMPLATE_ZERO_BLOCKS); \
  gt_fatal_check(gt_vector_get_used(template->mmaps)%gt_vector_get_used(template->blocks) != 0, \
           TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION); \
  gt_fatal_check((gt_vector_get_used(template->mmaps)/gt_vector_get_used(template->blocks)) != \
            gt_vector_get_used(template->mmaps_attributes),TEMPLATE_INCONSISTENT_MMAPS_ATTRB_RELATION)

#define GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_A,template_B) \
  gt_cond_fatal_error(gt_template_get_num_blocks(template_A) != \
    gt_template_get_num_blocks(template_B),TEMPLATE_INCONSISTENT_NUM_BLOCKS)


/*
 *  TODO: Scheduled for v2.0 (all lazy parsing)
 *  #define GT_TEMPLATE_EDITABLE_CHECK(template) \
 *  GT_TEMPLATE_CONSISTENCY_CHECK(template); \
 *  gt_fatal_check(template->maps_txt!=NULL,TEMPLATE_MAPS_NOT_PARSED)
 *  TODO: Scheduled for v2.0 (high level handling modules)
 */

/*
 * Reduction to single alignment
 */
#define GT_TEMPLATE_REDUCTION(template,reduced_alignment) \
  register gt_alignment* const reduced_alignment = gt_template_get_block((template),0)
#define GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,reduced_alignment) \
  if (gt_expect_false(gt_template_get_num_blocks((template))==1)) { \
    GT_TEMPLATE_REDUCTION(template,reduced_alignment);
#define GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT_(template) { \
  if (gt_expect_false(gt_template_get_num_blocks((template))==1))
#define GT_TEMPLATE_END_REDUCTION }
#define GT_TEMPLATE_END_REDUCTION__RETURN return;}

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new(void);
GT_INLINE void gt_template_clear(gt_template* const template,const bool delete_alignments);
GT_INLINE void gt_template_clear_alignments(gt_template* const template);
GT_INLINE void gt_template_delete(gt_template* const template);

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template);
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag,const uint64_t length);
GT_INLINE uint64_t gt_template_get_total_length(gt_template* const template);
/* Blocks (single alignments) */
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template);
GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position);
GT_INLINE gt_alignment* gt_template_get_block_dyn(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_clear_blocks(gt_template* const template);
GT_INLINE void gt_template_delete_blocks(gt_template* const template);
/* Counters */
GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template);
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters);
GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum);

/*
 * Attribute accessors
 */
GT_INLINE void* gt_template_get_attr(
    gt_template* const template,char* const attribute_id);
GT_INLINE void gt_template_set_attr(
    gt_template* const template,char* const attribute_id,
    void* const attribute,const size_t element_size);
#define gt_template_get_attribute(template,attribute_id,type) ((type*)gt_template_get_attr(template,attribute_id))
#define gt_template_set_attribute(template,attribute_id,attribute,type) gt_template_set_attr(template,attribute_id,attribute,sizeof(type))

GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template);
GT_INLINE void gt_template_set_mcs(gt_template* const template,uint64_t max_complete_strata);
GT_INLINE bool gt_template_has_qualities(gt_template* const template);
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template);
GT_INLINE void gt_template_set_not_unique_flag(gt_template* const template,bool is_not_unique);

/*
 * Multi-maps handlers (Map's relation)
 */
GT_INLINE void gt_template_mmap_attr_new(gt_mmap_attributes* const mmap_attr);
GT_INLINE gt_mmap_attributes* gt_template_get_mmap_attr(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_set_mmap_attr(gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr);
/* */
GT_INLINE uint64_t gt_template_get_num_mmaps(gt_template* const template);
GT_INLINE void gt_template_clear_mmaps(gt_template* const template);
/* */
GT_INLINE void gt_template_add_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE gt_map** gt_template_get_mmap(
    gt_template* const template,const uint64_t position,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_set_mmap(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attr);
/* */
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr);
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const mmap,gt_mmap_attributes* const mmap_attr);
/* */
GT_INLINE void gt_template_add_mmap_v(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,va_list v_args);
GT_INLINE void gt_template_add_mmap_va(
    gt_template* const template,gt_mmap_attributes* const mmap_attr,...);

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_copy(gt_template* const template,const bool copy_maps,const bool copy_mmaps);

/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE uint64_t gt_template_next_alignment_pos(gt_template_alignment_iterator* const template_alignment_iterator);
/*
 * Template's Multi-maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
 */
GT_INLINE void gt_template_new_mmap_iterator(gt_template* const template,gt_template_maps_iterator* const template_maps_iterator);
GT_INLINE gt_status gt_template_next_mmap(
    gt_template_maps_iterator* const template_maps_iterator,
    gt_map*** const mmap_ref,gt_mmap_attributes** const mmap_attr);
GT_INLINE uint64_t gt_template_next_mmap_pos(gt_template_maps_iterator* const template_maps_iterator);

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
#define GT_TEMPLATE__ATTR_ITERATE_(template,map_array,map_array_attr) \
  gt_map** map_array; \
  gt_mmap_attributes *map_array_attr = NULL; \
  gt_template_maps_iterator __##template##_maps_iterator; \
  gt_template_new_mmap_iterator(template,&(__##template##_maps_iterator)); \
  while (gt_template_next_mmap(&(__##template##_maps_iterator),&map_array,&map_array_attr))
#define GT_TEMPLATE__ATTR_ITERATE(template,map_array,map_array_attr) \
  register const uint64_t __map_array##_num_blocks = gt_template_get_num_blocks(template); \
  GT_TEMPLATE__ATTR_ITERATE_(template,map_array,map_array_attr)

/*
 * Iterate over array of maps provided by GT_TEMPLATE_ITERATE
 *   map_array = (end1.map1,end2.map1)
 *   GT_MULTIMAP_ITERATE(map_array) := {end1.map1,end2.map1}
 */
#define GT_MULTIMAP_ITERATE_BLOCKS(mmap_array,num_blocks,map,end_position) \
  register const uint64_t _num_blocks_##mmap_array = num_blocks; \
  register uint64_t end_position; \
  register gt_map* map; \
  for (end_position=0,map=*mmap_array; \
       end_position<(_num_blocks_##mmap_array);map=*(mmap_array+(++end_position)))
#define GT_MULTIMAP_ITERATE(mmap_array,map,end_position) \
  GT_MULTIMAP_ITERATE_BLOCKS(mmap_array,__map_array##_num_blocks,map,end_position)

/*
 * Iterate over the alignment of a template (individual blocks)
 *   Template = Alignment_end1 + Alignment_end2 + {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
 *   GT_TEMPLATE_ALIGNMENT_ITERATE(template) := {Alignment_end1,Alignment_end2}
 */
#define GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) \
  /* Template. Iterate over all alignments */ \
  gt_template_alignment_iterator alignment##_iterator; \
  register gt_alignment* alignment; \
  gt_template_new_alignment_iterator(template,&(alignment##_iterator)); \
  while ((alignment=gt_template_next_alignment(&(alignment##_iterator))))

#endif /* GT_TEMPLATE_H_ */
