/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sam_data_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to SAM format data structures
 */


#include "gt_sam_data_attributes.h"

/*
 * SAM File specifics Attribute (SAM Headers)
 */
#define GT_ATTR_SAM_INIT_ELEMENTS 2
GT_INLINE gt_sam_headers* gt_sam_header_new(void) {
  gt_sam_headers* sam_headers = gt_alloc(gt_sam_headers);
  sam_headers->header = gt_string_new(50); // @HD
  sam_headers->read_group = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @RG
  sam_headers->program = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @PG
  sam_headers->comments = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @ CO
  sam_headers->sequence_archive = NULL; // @SQ
  return sam_headers;
}
GT_INLINE void gt_sam_header_clear(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Header
  gt_string_clear(sam_headers->header);
  // Read group
  GT_VECTOR_ITERATE(sam_headers->read_group,record_h,nh,gt_string*) { gt_string_delete(*record_h); }
  gt_vector_clear(sam_headers->read_group);
  // Program
  GT_VECTOR_ITERATE(sam_headers->program,record_p,np,gt_string*) { gt_string_delete(*record_p); }
  gt_vector_clear(sam_headers->program);
  // Comments
  GT_VECTOR_ITERATE(sam_headers->comments,comment,nc,gt_string*) { gt_string_delete(*comment); }
  gt_vector_clear(sam_headers->comments);
  // Seq Archive
  sam_headers->sequence_archive=NULL;
}
GT_INLINE void gt_sam_header_delete(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Clear
  gt_sam_header_clear(sam_headers);
  // Delete
  gt_string_delete(sam_headers->header);
  gt_vector_delete(sam_headers->read_group);
  gt_vector_delete(sam_headers->program);
  gt_vector_delete(sam_headers->comments);
}

GT_INLINE gt_status gt_sam_header_set_header_record(gt_sam_headers* const sam_headers,gt_string* const header_line) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(header_line);
  gt_string_copy(sam_headers->header,header_line);
  return GT_SAM_HEADER_OK;
}
GT_INLINE gt_status gt_sam_header_add_read_group_record(gt_sam_headers* const sam_headers,gt_string* const read_group_record) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(read_group_record);
  gt_vector_insert(sam_headers->read_group,read_group_record,gt_string*);
  return GT_SAM_HEADER_OK;
}
GT_INLINE gt_status gt_sam_header_add_program_record(gt_sam_headers* const sam_headers,gt_string* const program_record) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(program_record);
  gt_vector_insert(sam_headers->program,program_record,gt_string*);
  return GT_SAM_HEADER_OK;
}
GT_INLINE gt_status gt_sam_header_add_comment(gt_sam_headers* const sam_headers,gt_string* const comment) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(comment);
  gt_vector_insert(sam_headers->comments,comment,gt_string*);
  return GT_SAM_HEADER_OK;
}

/*
 * SAM Optional Fields
 *   - SAM Attributes(optional fields) are just a vector of @gt_sam_attribute
 *     embedded into the general attributes(@gt_shash) of any object(@template,@alignment,@map,...)
 */
// Setup
GT_INLINE gt_vector* gt_attribute_sam_get_attributes(gt_shash* const general_attributes) {
  GT_HASH_CHECK(general_attributes);
  register gt_vector* vector = gt_attribute_get(general_attributes,GT_ATTR_ID_SAM);
  if (vector==NULL) vector = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_sam_attribute));
  return vector;
}
GT_INLINE gt_vector* gt_attribute_sam_fetch_attributes(gt_shash* const general_attributes) {
  GT_HASH_CHECK(general_attributes);
  return gt_attribute_get(general_attributes,GT_ATTR_ID_SAM);
}
GT_INLINE void gt_attribute_sam_delete_attributes(gt_shash* const general_attributes) {
  GT_HASH_CHECK(general_attributes);
  GT_ATTRIBUTES_SAM_ITERATE(general_attributes,attr) {
    if (attr->attribute_type==SAM_ATTR_STRING_VALUE) gt_string_delete(attr->s_value);
  }
}
GT_INLINE void gt_attribute_sam_clear_attributes(gt_shash* const general_attributes) {
  GT_HASH_CHECK(general_attributes);
  GT_ATTRIBUTES_SAM_ITERATE(general_attributes,attr) {
    if (attr->attribute_type==SAM_ATTR_STRING_VALUE) gt_string_delete(attr->s_value);
  }
}
GT_INLINE bool gt_attribute_has_sam_attributes(gt_shash* const general_attributes) {
  GT_HASH_CHECK(general_attributes);
  register gt_vector* vector = gt_attribute_get(general_attributes,GT_ATTR_ID_SAM);
  return (vector!=NULL);
}
/*
 * Accessors
 */
#define GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag_src) \
  sam_attribute.tag[0]=tag_src[0]; \
  sam_attribute.tag[1]=tag_src[1]
#define GT_ATTRIBUTE_SAM_CMP_TAG(sam_attribute_ptr,tag_src) (sam_attribute_ptr->tag[0]==tag_src[0] && sam_attribute_ptr->tag[1]==tag_src[1])
GT_INLINE gt_sam_attribute* gt_attribute_sam_get(gt_shash* const general_attributes,char* const tag) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  // Iterate over all attributes till we find it !
  GT_ATTRIBUTES_SAM_ITERATE(general_attributes,attr) {
    if (GT_ATTRIBUTE_SAM_CMP_TAG(attr,tag)) return attr;
  }
  return NULL;
}
/*
 * Add values (optional fields values)
 */
GT_INLINE void gt_attribute_sam_add_ivalue(gt_shash* const general_attributes,char* const tag,char type_id,const int64_t value) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_INT_VALUE;
  sam_attribute.type_id = type_id;
  sam_attribute.i_value = value;
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_fvalue(gt_shash* const general_attributes,char* const tag,char type_id,const double value) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_FLOAT_VALUE;
  sam_attribute.type_id = type_id;
  sam_attribute.d_value = value;
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_svalue(gt_shash* const general_attributes,char* const tag,char type_id,char* const text,const int64_t length) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_STRING_VALUE;
  sam_attribute.type_id = type_id;
  sam_attribute.s_value = gt_string_new(length+1);
  gt_string_set_nstring(sam_attribute.s_value,text,length);
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}

/*
 * Add function (optional fields values generated by the function)
 */
GT_INLINE void gt_attribute_sam_add_ifunc(gt_shash* const general_attributes,char* const tag,char type_id,int64_t (*i_func)(gt_template*,gt_alignment*,gt_map**)) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_INT_FUNC;
  sam_attribute.type_id = type_id;
  sam_attribute.i_func = i_func;
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_ffunc(gt_shash* const general_attributes,char* const tag,char type_id,double (*d_func)(gt_template*,gt_alignment*,gt_map**)) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_FLOAT_FUNC;
  sam_attribute.type_id = type_id;
  sam_attribute.d_func = d_func;
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}
GT_INLINE void gt_attribute_sam_add_sfunc(gt_shash* const general_attributes,char* const tag,char type_id,gt_string* (*s_func)(gt_template*,gt_alignment*,gt_map**)) {
  GT_HASH_CHECK(general_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute sam_attribute;
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute.attribute_type = SAM_ATTR_STRING_FUNC;
  sam_attribute.type_id = type_id;
  sam_attribute.s_func = s_func;
  register gt_vector* vector = gt_attribute_sam_get_attributes(general_attributes);
  gt_vector_insert(vector,sam_attribute,gt_sam_attribute);
}

