/*
 * PROJECT: GEM-Tools library
 * FILE: gt_hash.c
 * DATE: 10/07/2012
 * DESCRIPTION: // TODO
 */

#include "gt_hash.h"

/*
 * Checkers
 */
#define GT_HASH_CHECK(hash) gt_fatal_check(hash==NULL,NULL_HANDLER)

GT_INLINE gt_shash* gt_shash_new(void) {
  gt_shash* shash=malloc(sizeof(gt_shash));
  gt_cond_fatal_error(!shash,MEM_HANDLER);
  shash->shash_head = NULL; // uthash initializer
  return shash;
}
GT_INLINE void gt_shash_clean(gt_shash* const shash,const bool free_element,const bool free_key) {
  GT_HASH_CHECK(shash);
  gt_shash_element *shash_element, *tmp;
  HASH_ITER(hh,shash->shash_head,shash_element,tmp) {
    HASH_DEL(shash->shash_head,shash_element);
    if (free_element) free(shash_element->element);
    if (free_key) free(shash_element->key);
    free(shash_element);
  }
}
GT_INLINE void gt_shash_delete(gt_shash* const shash,const bool free_element,const bool free_key) {
  GT_HASH_CHECK(shash);
  gt_shash_clean(shash,free_element,free_key);
  free(shash);
}
GT_INLINE void gt_shash_insert_element(gt_shash* const shash,char* const key,void* element) {
  GT_HASH_CHECK(shash);
  gt_shash_element *shash_element = malloc(sizeof(gt_shash_element));
  shash_element->key = key;
  shash_element->element = element;
  HASH_ADD_KEYPTR(hh,shash->shash_head,
      shash_element->key,strlen(shash_element->key),shash_element);
}
GT_INLINE gt_shash_element* gt_shash_get_shash_element(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element;
  HASH_FIND_STR(shash->shash_head,key,shash_element);
  return shash_element;
}
GT_INLINE void* gt_shash_get_key(gt_shash* const shash,char* const key) {
  if (!shash) {
    printf("ss");
  }
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
    HASH_DEL(shash->shash_head,shash_element);
  }
}
GT_INLINE bool gt_shash_is_contained(gt_shash* const shash,char* const key) {
  GT_HASH_CHECK(shash);
  GT_NULL_CHECK(key);
  gt_shash_element *shash_element = gt_shash_get_shash_element(shash,key);
  return (shash_element!=NULL);
}
GT_INLINE uint64_t gt_shash_get_num_elements(gt_shash* const shash) {
  GT_HASH_CHECK(shash);
  return (uint64_t)HASH_COUNT(shash->shash_head);
}
