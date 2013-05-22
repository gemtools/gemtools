/*
 * PROJECT: GEM-Tools library
 * FILE: gt_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes (from general attributes to specific to certain file formats)
 */


#include "gt_attributes.h"
#include "gt_sam_attributes.h"

/*
 * General Attribute accessors
 */
GT_INLINE gt_attributes* gt_attributes_new(void) {
  return gt_shash_new();
}
GT_INLINE void gt_attributes_clear(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_shash_clear(attributes,true);
}
GT_INLINE void gt_attributes_delete(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_shash_delete(attributes,true);
}
GT_INLINE void* gt_attributes_get(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_get(attributes,attribute_id,void);
}
GT_INLINE bool gt_attributes_is_contained(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_is_contained(attributes,attribute_id);
}
GT_INLINE void gt_attributes_add_string(
    gt_attributes* const attributes,char* const attribute_id,gt_string* const attribute_string) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_STRING_CHECK(attribute_string);
  // Insert attribute
  gt_shash_insert_string(attributes,attribute_id,attribute_string);
}
GT_INLINE void gt_attributes_add_primitive(
    gt_attributes* const attributes,char* const attribute_id,void* const attribute,const size_t element_size) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_ZERO_CHECK(element_size);
  // We do a copy of the element as to handle it ourselves from here
  void* attribute_cp = gt_malloc(element_size); // Allocate attribute
  memcpy(attribute_cp,attribute,element_size); // Copy attribute
  // Insert attribute
  gt_shash_insert_primitive(attributes,attribute_id,attribute_cp,element_size);
}
GT_INLINE void gt_attributes_add_object(
    gt_attributes* const attributes,char* const attribute_id,
    void* const attribute,void* (*attribute_dup_fx)(),void (*attribute_free_fx)()) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_NULL_CHECK(attribute_dup_fx);
  GT_NULL_CHECK(attribute_free_fx);
  // Insert attribute
  gt_shash_insert_object(attributes,attribute_id,attribute,attribute_dup_fx,attribute_free_fx);
}
GT_INLINE void gt_attributes_remove(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_shash_remove(attributes,attribute_id,true);
}
GT_INLINE gt_attributes* gt_attributes_dup(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  return gt_shash_dup(attributes);
}
GT_INLINE void gt_attributes_copy(gt_attributes* const attributes_dst,gt_attributes* const attributes_src) {
  GT_ATTRIBUTES_CHECK(attributes_dst);
  GT_ATTRIBUTES_CHECK(attributes_src);
  gt_shash_copy(attributes_dst,attributes_src);
}

