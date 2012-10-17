/*
 * PROJECT: GEM-Tools library
 * FILE: gt_ihash.h
 * DATE: 2/09/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_IHASH_H_
#define GT_IHASH_H_

#include "gt_commons.h"
#include "uthash.h"

/*
 * Integer Key Hash
 */
typedef struct {
  int64_t key;
  void* element;
  size_t element_size;
  UT_hash_handle hh;
} gt_ihash_element;
typedef struct {
  gt_ihash_element* ihash_head;
} gt_ihash;

/*
 * Constructor
 */
GT_INLINE gt_ihash* gt_ihash_new(void);
GT_INLINE void gt_ihash_clean(gt_ihash* const ihash,const bool free_element);
GT_INLINE void gt_ihash_delete(gt_ihash* const ihash,const bool free_element);

/*
 * Basic (Type-unsafe) Accessors
 */
GT_INLINE void gt_ihash_insert_element(gt_ihash* const ihash,const int64_t key,void* const element,const size_t element_size);
GT_INLINE void* gt_ihash_get_element(gt_ihash* const ihash,const int64_t key);
GT_INLINE void gt_ihash_remove_element(gt_ihash* const ihash,const int64_t key);

/*
 * Type-safe Accessors
 */
#define gt_ihash_get(ihash,integer_key,type) ((type*)gt_ihash_get_element(ihash,integer_key))
#define gt_ihash_insert(ihash,integer_key,element,type) gt_ihash_insert_element(ihash,integer_key,(void*)element,sizeof(type))
#define gt_ihash_remove(ihash,integer_key) gt_ihash_remove_element(ihash,integer_key)
GT_INLINE bool gt_ihash_is_contained(gt_ihash* const ihash,const int64_t key);
GT_INLINE uint64_t gt_ihash_get_num_elements(gt_ihash* const ihash);

/*
 * Miscellaneous
 */
GT_INLINE gt_ihash* gt_ihash_copy(gt_ihash* const ihash);
GT_INLINE gt_ihash* gt_ihash_deep_copy(gt_ihash* const ihash);

/*
 * Iterator
 */
#define GT_IHASH_BEGIN_ITERATE(ihash,it_ikey,it_element,type) { \
  gt_ihash_element *ihash_##ih_element, *ihash_##tmp; \
  HASH_ITER(hh,ihash->ihash_head,ihash_##ih_element,ihash_##tmp) { \
    register type* const it_element = (type*)(ihash_##ih_element->element); \
    register int64_t const it_ikey = ihash_##ih_element->key;
#define GT_IHASH_END_ITERATE }}

#endif /* GT_IHASH_H_ */
