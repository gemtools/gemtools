/*
 * PROJECT: GEM-Tools library
 * FILE: gt_data_attributes.h
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes (from general attributes to specific to certain file formats)
 */

#ifndef GT_DATA_ATTRIBUTES_H_
#define GT_DATA_ATTRIBUTES_H_

#include "gt_essentials.h"

/*
 * General Attributes
 */
GT_INLINE gt_shash* gt_attribute_new(void);
GT_INLINE void gt_attribute_clear(gt_shash* const attributes);
GT_INLINE void gt_attribute_delete(gt_shash* const attributes);

GT_INLINE void* gt_attribute_get(gt_shash* const attributes,char* const attribute_id);
GT_INLINE void gt_attribute_add_string(
    gt_shash* const attributes,char* const attribute_id,gt_string* const attribute_string);
GT_INLINE void gt_attribute_add_primitive(
    gt_shash* const attributes,char* const attribute_id,void* const attribute,const size_t element_size);
GT_INLINE void gt_attribute_add_object(
    gt_shash* const attributes,char* const attribute_id,
    void* const attribute,void* (*attribute_dup_fx)(),void (*attribute_free_fx)());

#define gt_attribute_add(attributes,attribute_id,attribute,element_type) \
    gt_attribute_add_primitive(attributes,attribute_id,(void*)attribute,sizeof(element_type))

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

#endif /* GT_DATA_ATTRIBUTES_H_ */
