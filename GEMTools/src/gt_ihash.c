/*
 * PROJECT: GEM-Tools library
 * FILE: gt_ihash.c
 * DATE: 2/09/2012
 * DESCRIPTION: // TODO
 */

#include "gt_hash.h"

/*
 * Constructor
 */
GT_INLINE gt_ihash* gt_ihash_new(void) {
  gt_ihash* ihash=malloc(sizeof(gt_ihash));
  gt_cond_fatal_error(!ihash,MEM_HANDLER);
  ihash->ihash_head = NULL; // uthash initializer
  return ihash;
}
GT_INLINE void gt_ihash_clean(gt_ihash* const ihash,const bool free_element) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element, *tmp;
  HASH_ITER(hh,ihash->ihash_head,ihash_element,tmp) {
    HASH_DEL(ihash->ihash_head,ihash_element);
    if (free_element) free(ihash_element->element);
    free(ihash_element);
  }
}
GT_INLINE void gt_ihash_delete(gt_ihash* const ihash,const bool free_element) {
  GT_HASH_CHECK(ihash);
  gt_ihash_clean(ihash,free_element);
  free(ihash);
}

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE gt_ihash_element* gt_ihash_get_ihash_element(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element;
  HASH_FIND_INT(ihash->ihash_head,&key,ihash_element);
  return ihash_element;
}
GT_INLINE void gt_ihash_insert_element(gt_ihash* const ihash,const int64_t key,void* const element,const size_t element_size) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element* ihash_element = gt_ihash_get_ihash_element(ihash,key);
  if (gt_expect_true(ihash_element==NULL)) {
    ihash_element = malloc(sizeof(gt_ihash_element));
    gt_cond_fatal_error(!ihash_element,MEM_ALLOC);
    ihash_element->key = key;
    ihash_element->element = element;
    ihash_element->element_size = element_size;
    HASH_ADD_INT(ihash->ihash_head,key,ihash_element);
  } else {
    ihash_element->element = element;
  }
}
GT_INLINE void* gt_ihash_get_element(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  register gt_ihash_element* const ihash_element = gt_ihash_get_ihash_element(ihash,key);
  return gt_expect_true(ihash_element!=NULL) ? ihash_element->element : NULL;
}
GT_INLINE void gt_ihash_remove_element(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  gt_ihash_element *ihash_element = gt_ihash_get_ihash_element(ihash,key);
  if (ihash_element) {
    HASH_DEL(ihash->ihash_head,ihash_element);
  }
}

/*
 * Type-safe Accessors
 */
GT_INLINE bool gt_ihash_is_contained(gt_ihash* const ihash,const int64_t key) {
  GT_HASH_CHECK(ihash);
  return (gt_ihash_get_ihash_element(ihash,key)!=NULL);
}
GT_INLINE uint64_t gt_ihash_get_num_elements(gt_ihash* const ihash) {
  GT_HASH_CHECK(ihash);
  return (uint64_t)HASH_COUNT(ihash->ihash_head);
}

/*
 * Miscellaneous
 */
GT_INLINE gt_ihash* gt_ihash_copy(gt_ihash* const ihash) {
  register gt_ihash* const ihash_cp =  gt_ihash_new();
  GT_IHASH_BEGIN_ITERATE(ihash,ikey,ihash_element,void) {
    // Insert element into the copy
    gt_ihash_insert_element(ihash_cp,ikey,ihash_element,ihash_ih_element->element_size);
  } GT_IHASH_END_ITERATE
  return ihash_cp;
}
GT_INLINE gt_ihash* gt_ihash_deep_copy(gt_ihash* const ihash) {
  register gt_ihash* const ihash_cp =  gt_ihash_new();
  GT_IHASH_BEGIN_ITERATE(ihash,ikey,ihash_element,void) {
    // Copy element
    register void* const ihash_element_cp = malloc(ihash_ih_element->element_size);
    gt_cond_fatal_error(!ihash_element_cp,MEM_ALLOC);
    memcpy(ihash_element_cp,ihash_element,ihash_ih_element->element_size);
    // Insert element into the copy
    gt_ihash_insert_element(ihash_cp,ikey,ihash_element_cp,ihash_ih_element->element_size);
  } GT_IHASH_END_ITERATE
  return ihash_cp;
}
