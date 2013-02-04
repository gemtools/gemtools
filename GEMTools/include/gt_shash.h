/*
 * PROJECT: GEM-Tools library
 * FILE: gt_shash.h
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_SHASH_H_
#define GT_SHASH_H_

#include "gt_commons.h"
#include "uthash.h"

/*
 * String Key Hash
 */
typedef struct {
  char* key;
  void* element;
  size_t element_size;
  UT_hash_handle hh;
} gt_shash_element;
typedef struct {
  gt_shash_element* shash_head;
} gt_shash;

/*
 * Constructor
 */
GT_INLINE gt_shash* gt_shash_new(void);
GT_INLINE void gt_shash_clear(gt_shash* const shash,const bool free_element);
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element);

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_shash_element* gt_shash_get_shash_element(gt_shash* const shash,char* const key);
GT_INLINE char* gt_shash_insert_element(gt_shash* const shash,char* const key,void* const element,const size_t element_size);
GT_INLINE char* gt_shash_get_key(gt_shash* const shash,char* const key);
GT_INLINE void* gt_shash_get_element(gt_shash* const shash,char* const key);
GT_INLINE void gt_shash_remove_element(gt_shash* const shash,char* const key);

/*
 * Type-safe Accessors
 */
#define gt_shash_get(shash,string_key,type) ((type*)gt_shash_get_element(shash,string_key))
#define gt_shash_insert(shash,string_key,element,type) gt_shash_insert_element(shash,string_key,(void*)element,sizeof(type))
#define gt_shash_remove(shash,string_key) gt_shash_remove_element(shash,string_key)
GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key);
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash);

/*
 * Miscellaneous
 */
GT_INLINE gt_shash* gt_shash_copy(gt_shash* const shash);
GT_INLINE void gt_shash_deep_copy(gt_shash* const shash_dst,gt_shash* const shash_src);

/*
 * Iterator
 */
#define GT_SHASH_BEGIN_ITERATE(shash,it_skey,it_element,type) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    register type* const it_element = (type*)(shash_##sh_element->element); \
    register char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_BEGIN_ELEMENT_ITERATE(shash,it_element,type) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    register type* const it_element = (type*)(shash_##sh_element->element);

#define GT_SHASH_BEGIN_KEY_ITERATE(shash,it_skey,type) { \
  gt_shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,shash_##tmp) { \
    register char* const it_skey = shash_##sh_element->key;

#define GT_SHASH_END_ITERATE }}

#endif /* GT_SHASH_H_ */
