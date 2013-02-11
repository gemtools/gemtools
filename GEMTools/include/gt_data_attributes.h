/*
 * PROJECT: GEM-Tools library
 * FILE: gt_data_attributes.h
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes (from general attributes to specific to certain file formats)
 */

#ifndef GT_DATA_ATTRIBUTES_H_
#define GT_DATA_ATTRIBUTES_H_

#include "gt_commons.h"
#include "gt_sequence_archive.h"

/*
 * Checkers
 */
#define GT_SAM_HEADERS_CHECK(sam_headers) \
  GT_VECTOR_CHECK(sam_headers->header); \
  GT_VECTOR_CHECK(sam_headers->read_group); \
  GT_VECTOR_CHECK(sam_headers->program); \
  GT_VECTOR_CHECK(sam_headers->comments); \
  GT_SEQUENCE_ARCHIVE_CHECK(sam_headers->sequence_archive)


/*
 * General Attributes
 */
GT_INLINE gt_shash* gt_attribute_new(void);
GT_INLINE void gt_attribute_clear(gt_shash* const attributes);
GT_INLINE void gt_attribute_delete(gt_shash* const attributes);

GT_INLINE void* gt_attribute_get(gt_shash* const attributes,char* const attribute_id);
GT_INLINE void gt_attribute_set_string(
    gt_shash* const attributes,char* const attribute_id,gt_string* const attribute_string);
GT_INLINE void gt_attribute_set_(
    gt_shash* const attributes,char* const attribute_id,void* const attribute,const size_t element_size);

#define gt_attribute_set(attributes,attribute_id,attribute,element_type) \
    gt_attribute_set_(attributes,attribute_id,(void*)attribute,sizeof(element_type))

/*
 * MAP File specifics Attribute
 */
typedef struct {
  bool contains_qualities;
} gt_map_file_format;

/*
 * FASTQ/FASTA/MULTIFASTA File specifics Attribute
 */
typedef enum { F_FASTA, F_FASTQ, F_MULTI_FASTA } gt_file_fasta_format;
typedef struct {
  gt_file_fasta_format fasta_format;
} gt_fasta_file_format;

/*
 * SAM File specifics Attribute (headers)
 */
typedef struct {
  char tag[2];
  gt_string* value;
} gt_sam_header_record; // SAM Header Record
typedef struct {
  gt_vector* header; // @HD /* (gt_sam_header_record) */
  gt_sequence_archive* sequence_archive; // @SQ
  gt_vector* read_group; // @RG  /* (gt_sam_header_record) */
  gt_vector* program; // @PG /* (gt_sam_header_record) */
  gt_vector* comments; // @ CO /* (gt_sam_header_record) */
} gt_sam_headers; // SAM Headers

// SAM Optional Fields
#define GT_ATTR_SAM "SAM_ATTR"

typedef enum { SAM_ATTR_INT, SAM_ATTR_FLOAT, SAM_ATTR_STRING } gt_sam_attribute_t;
typedef struct {
  char tag[2];
  gt_sam_attribute_t attribute_type;
  char type_id;
  union {
    int64_t i_value;
    double d_value;
    gt_string* s_value;
  };
} gt_sam_attribute;


GT_INLINE gt_sam_headers* gt_sam_header_new(void);
GT_INLINE void gt_sam_header_clear(gt_sam_headers* const sam_headers);
GT_INLINE void gt_sam_header_delete(gt_sam_headers* const sam_headers);

GT_INLINE bool gt_attribute_sam_has_attr_vector(gt_shash* const attributes);
GT_INLINE gt_vector* gt_attribute_sam_get_attr_vector(gt_shash* const attributes);

GT_INLINE void gt_attribute_sam_add_ivalue(gt_shash* const attributes,char* const tag,char type_id,const int64_t value);
GT_INLINE void gt_attribute_sam_add_fvalue(gt_shash* const attributes,char* const tag,char type_id,const double value);
GT_INLINE void gt_attribute_sam_add_svalue(gt_shash* const attributes,char* const tag,char type_id,char* const text,const int64_t length);

GT_INLINE gt_sam_attribute* gt_attribute_sam_get(gt_shash* const attributes,char* const tag);

#define GT_ATTRIBUTES_SAM_ITERATE(attributes,attr_element) \
  register const gt_vector* sam_attr_##attributes = gt_attribute_get(attributes,GT_ATTR_SAM);  \
  GT_VECTOR_ITERATE(sam_attr_##attributes,attr_element,attr_element_##counter,gt_sam_attribute)


#endif /* GT_DATA_ATTRIBUTES_H_ */
