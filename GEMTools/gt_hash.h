/*
 * PROJECT: GEM-Tools library
 * FILE: gt_hash.h
 * DATE: 10/07/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_HASH_H_
#define GT_HASH_H_

#include "gt_commons.h"
#include "uthash.h"

/*
 * String Key Hash
 */
typedef struct {
  char* key;
  void* element;
  UT_hash_handle hh;
} gt_shash_element;
typedef struct {
  gt_shash_element* shash_head;
} gt_shash;

GT_INLINE gt_shash* gt_shash_new(void);
GT_INLINE void gt_shash_clean(gt_shash* const shash,const bool free_element,const bool free_key);
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element,const bool free_key);

GT_INLINE void gt_shash_insert_element(gt_shash* const shash,char* const key,void* element);
GT_INLINE void* gt_shash_get_key(gt_shash* const shash,char* const key);
GT_INLINE void* gt_shash_get_element(gt_shash* const shash,char* const key);
GT_INLINE void gt_shash_remove_element(gt_shash* const shash,char* const key);
/* Type-safe functions */
#define gt_shash_get(shash,key,type) ((type*)gt_shash_get_element(shash,key))
#define gt_shash_insert(shash,key,element) gt_shash_insert_element(shash,key,(void*)element)
#define gt_shash_remove(shash,key) gt_shash_remove_element(shash,key)

GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key);
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash);

#define GT_SHASH_BEGIN_ITERATE(shash,skey,element,type) { \
  gt_shash_element *shash_##sh_element, *tmp; \
  HASH_ITER(hh,shash->shash_head,shash_##sh_element,tmp) { \
    register type* const element = (type*)(shash_##sh_element->element); \
    register char* const skey = shash_##sh_element->key;
#define GT_SHASH_END_ITERATE }}

#endif /* GT_HASH_H_ */
