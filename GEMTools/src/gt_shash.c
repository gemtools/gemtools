/*
 * PROJECT: GEM-Tools library
 * FILE: gt_shash.c
 * DATE: 10/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_hash.h"

/*
 * Constructor
 */
GT_INLINE gt_shash* gt_shash_new(void) {
  gt_shash* shash=malloc(sizeof(gt_shash));
  gt_cond_fatal_error(!shash,MEM_HANDLER);
  shash->shash_head = NULL; // uthash initializer
  return shash;
}
GT_INLINE void gt_shash_clear(gt_shash* const shash,const bool free_element) {
  GT_HASH_CHECK(shash);
  gt_shash_element *shash_element, *tmp;
  HASH_ITER(hh,shash->shash_head,shash_element,tmp) {
    HASH_DEL(shash->shash_head,shash_element);
    free(shash_element->key);
    if (free_element) free(shash_element->element);
    free(shash_element);
  }
}
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element) {
  GT_HASH_CHECK(shash);
  gt_shash_clear(shash,free_element);
  free(shash);
}

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_shash_element* gt_shash_get_shash_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element;
  HASH_FIND_STR(shash->shash_head,key,shash_element);
  return shash_element;
}
GT_INLINE char* gt_shash_insert_element(gt_shash* const shash,char* const key,void* const element,const size_t element_size) {
  GT_HASH_CHECK(shash);
  register gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  if (gt_expect_true(shash_element==NULL)) {
    shash_element = malloc(sizeof(gt_shash_element));
    gt_cond_fatal_error(!shash_element,MEM_ALLOC);
    shash_element->key = gt_strndup(key,strlen(key));
    shash_element->element = element;
    shash_element->element_size = element_size;
    HASH_ADD_KEYPTR(hh,shash->shash_head,
        shash_element->key,strlen(shash_element->key),shash_element);
  } else {
    shash_element->element = element;
  }
  return shash_element->key;
}
GT_INLINE char* gt_shash_get_key(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  return gt_expect_true(shash_element!=NULL) ? shash_element->key : NULL;
}
GT_INLINE void* gt_shash_get_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  return gt_expect_true(shash_element!=NULL) ? shash_element->element : NULL;
}
GT_INLINE void gt_shash_remove_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  if (shash_element) {
    free(shash_element->key);
    HASH_DEL(shash->shash_head,shash_element);
  }
}

/*
 * Type-safe Accessors
 */
GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  return (gt_shash_get_shash_element(shash,key)!=NULL);
}
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash) {
  GT_HASH_CHECK(shash);
  return (uint64_t)HASH_COUNT(shash->shash_head);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_shash* gt_shash_copy(gt_shash* const shash) {
  register gt_shash* const shash_cpy =  gt_shash_new();
  GT_SHASH_BEGIN_ITERATE(shash,skey,shash_element,void) {
    // Insert element into the copy
    gt_shash_insert_element(shash_cpy,skey,shash_element,shash_sh_element->element_size);
  } GT_SHASH_END_ITERATE
  return shash_cpy;
}
GT_INLINE gt_shash* gt_shash_deep_copy(gt_shash* const shash) {
  register gt_shash* const shash_cp =  gt_shash_new();
  GT_SHASH_BEGIN_ITERATE(shash,skey,shash_element,void) {
    // Copy element
    register void* const shash_element_cp = malloc(shash_sh_element->element_size);
    gt_cond_fatal_error(!shash_element_cp,MEM_ALLOC);
    memcpy(shash_element_cp,shash_element,shash_sh_element->element_size);
    // Insert element into the copy
    gt_shash_insert_element(shash_cp,skey,shash_element_cp,shash_sh_element->element_size);
  } GT_SHASH_END_ITERATE
  return shash_cp;
}

