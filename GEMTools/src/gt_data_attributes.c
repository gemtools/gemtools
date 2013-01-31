/*
 * PROJECT: GEM-Tools library
 * FILE: gt_data_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes (from general attributes to specific to certain file formats)
 */


#include "gt_data_attributes.h"

/*
 * General Attribute accessors
 */
GT_INLINE gt_shash* gt_attribute_new(void) {
  return gt_shash_new();
}

GT_INLINE void gt_attribute_clear(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  // Free string SAM attributes
  if (gt_attribute_sam_has_attr_vector(attributes)) {
    GT_ATTRIBUTES_SAM_ITERATE(attributes,attr) {
      if (attr->attribute_type==SAM_ATTR_STRING) gt_string_delete(attr->s_value);
    }
  }
  // Actual clear of the hash
  gt_shash_clear(attributes,true);
}

GT_INLINE void gt_attribute_delete(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  if (gt_attribute_sam_has_attr_vector(attributes)) { // Free string attributes
    GT_ATTRIBUTES_SAM_ITERATE(attributes,attr) {
      if (attr->attribute_type==SAM_ATTR_STRING) gt_string_delete(attr->s_value);
    }
  }
  // Actual free of the hash
  gt_shash_delete(attributes,true);
}

GT_INLINE void* gt_attribute_get(gt_shash* const attributes,char* const attribute_id) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_get(attributes,attribute_id,void);
}
GT_INLINE void gt_attribute_set(
    gt_shash* const attributes,char* const attribute_id,
    void* const attribute,const size_t element_size) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_ZERO_CHECK(element_size);
  // NOTE: We do a copy of the element as to handle it ourselves from here
  register void* attr = malloc(element_size); // Allocate attribute
  gt_cond_fatal_error(!attr,MEM_ALLOC);
  memcpy(attr,attribute,element_size); // Copy attribute
  // Test attribute
  register gt_shash_element* shash_element = gt_shash_get_shash_element(attributes,attribute_id);
  if (gt_expect_false(shash_element!=NULL)) {
    free(shash_element->element);
    shash_element->element = attr;
  } else {
    // Insert attribute
    gt_shash_insert_element(attributes,attribute_id,attr,element_size);
  }
}

/*
 * SAM Attributes
 */

#define GT_ATTR_SAM_INIT_ELEMENTS 5

GT_INLINE gt_sam_headers* gt_sam_header_new(void) {
  gt_sam_headers* sam_headers = malloc(sizeof(gt_sam_headers));
  gt_cond_fatal_error(!sam_headers,MEM_HANDLER);
  sam_headers->header = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_header_record)); // @HD
  sam_headers->read_group = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_header_record)); // @RG
  sam_headers->program = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_header_record)); // @PG
  sam_headers->comments = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_header_record)); // @ CO
  sam_headers->sequence_archive = gt_sequence_archive_new(); // @SQ
  return sam_headers;
}
GT_INLINE void gt_sam_header_clear(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Clear vectors + hash
  gt_vector_clear(sam_headers->header);
  gt_vector_clear(sam_headers->read_group);
  gt_vector_clear(sam_headers->program);
  gt_vector_clear(sam_headers->comments);
  gt_sequence_archive_clear(sam_headers->sequence_archive);
}

GT_INLINE void gt_sam_header_delete(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Delete vectors + hash
  gt_vector_delete(sam_headers->header);
  gt_vector_delete(sam_headers->read_group);
  gt_vector_delete(sam_headers->program);
  gt_vector_delete(sam_headers->comments);
  gt_sequence_archive_delete(sam_headers->sequence_archive);
}

#define GT_ATTR_SAM "SAM_ATTR"

#define GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag_src) \
  sam_attribute.tag[0]=tag_src[0]; \
  sam_attribute.tag[1]=tag_src[1]

#define GT_ATTRIBUTE_SAM_CMP_TAG(sam_attribute_ptr,tag_src) (sam_attribute_ptr->tag[0]==tag_src[0] && sam_attribute_ptr->tag[1]==tag_src[1])

GT_INLINE bool gt_attribute_sam_has_attr_vector(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  register gt_vector* vector = gt_attribute_get(attributes,GT_ATTR_SAM);
  return (vector!=NULL);
}
GT_INLINE gt_vector* gt_attribute_sam_get_attr_vector(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  register gt_vector* vector = gt_attribute_get(attributes,GT_ATTR_SAM);
  if (vector==NULL) {
    vector = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_header_record));
  }
  return vector;
}

GT_INLINE void gt_attribute_sam_add_ivalue(gt_shash* const attributes,char* const tag,char type_id,const int64_t value) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_INT;
  sam_attribute.type_id = type_id;
  sam_attribute.i_value = value;
  register gt_vector* vector = gt_attribute_sam_get_attr_vector(attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_fvalue(gt_shash* const attributes,char* const tag,char type_id,const double value) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_FLOAT;
  sam_attribute.type_id = type_id;
  sam_attribute.d_value = value;
  register gt_vector* vector = gt_attribute_sam_get_attr_vector(attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_svalue(gt_shash* const attributes,char* const tag,char type_id,char* const text,const int64_t length) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_STRING;
  sam_attribute.type_id = type_id;
  sam_attribute.s_value = gt_string_new(length+1);
  gt_string_set_nstring(sam_attribute.s_value,text,length);
  register gt_vector* vector = gt_attribute_sam_get_attr_vector(attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}

GT_INLINE gt_sam_attribute* gt_attribute_sam_get(gt_shash* const attributes,char* const tag) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(tag);
  // Iterate over all attributes till we find it !
  GT_ATTRIBUTES_SAM_ITERATE(attributes,attr) {
    if (GT_ATTRIBUTE_SAM_CMP_TAG(attr,tag)) return attr;
  }
  return NULL;
}

