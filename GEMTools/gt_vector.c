/*
 * PROJECT: GEM-Tools library
 * FILE: gt_vector.c
 * DATE: 01/06/2012
 * DESCRIPTION: This file implements vectors based on raw memory buffers.
 * AUTHOR:
 *   (C) 2008-2011 P. Ribeca <paolo.ribeca@gmail.com>, all rights reserved
 *   (C) 2011      S. Marco Sola <santiagomsola@gmail.com>, all rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "gt_vector.h"

#define GT_VECTOR_EXPAND_FACTOR (3.0/2.0)

#define GT_VECTOR_CHECK(vector) \
  GT_NULL_CHECK(vector); \
  GT_NULL_CHECK(vector->memory); \
  GT_ZERO_CHECK(vector->element_size)

GT_INLINE gt_vector* gt_vector_new(size_t num_initial_elements,size_t element_size) {
  GT_ZERO_CHECK(element_size);
  gt_vector* vector=malloc(sizeof(gt_vector));
  gt_cond_fatal_error(!vector,MEM_HANDLER);
  vector->element_size=element_size;
  vector->elements_allocated=num_initial_elements;
  vector->memory=malloc(num_initial_elements*element_size);
  gt_cond_fatal_error(!vector->memory,MEM_ALLOC);
  vector->used=0;
  return vector;
}
GT_INLINE gt_status gt_vector_reserve(gt_vector* vector,size_t num_elements,bool zero_mem) {
  GT_VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    size_t proposed=(float)vector->elements_allocated*GT_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    vector->memory=realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) return GT_VECTOR_FAIL;
  }
  if (gt_expect_false(zero_mem)) {
    memset(vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
  return GT_VECTOR_OK;
}
GT_INLINE gt_status gt_vector_resize__clean(gt_vector* vector,size_t num_elements) {
  GT_VECTOR_CHECK(vector);
  if (vector->elements_allocated < num_elements) {
    size_t proposed=(float)vector->elements_allocated*GT_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    free(vector->memory);
    vector->memory=malloc(vector->elements_allocated*vector->element_size);
    if (!vector->memory) return GT_VECTOR_FAIL;
  }
  vector->used=0;
  return GT_VECTOR_OK;
}

GT_INLINE void gt_vector_cast__clean(gt_vector* vector,size_t element_size) {
  GT_VECTOR_CHECK(vector); GT_ZERO_CHECK(element_size);
  vector->elements_allocated=(vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size=element_size;
  vector->used=0;
}
GT_INLINE void gt_vector_delete(gt_vector* vector) {
  GT_VECTOR_CHECK(vector);
  free(vector->memory);
  free(vector);
}
GT_INLINE void gt_vector_copy(gt_vector* vector_to,gt_vector* vector_from) {
  GT_VECTOR_CHECK(vector_to); GT_VECTOR_CHECK(vector_from);
  gt_vector_cast__clean(vector_to,vector_from->element_size);
  gt_vector_reserve(vector_to,vector_from->used,false);
  memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
}

GT_INLINE void* gt_vector_get_mem_element(gt_vector* vector,size_t position,size_t element_size) {
  GT_VECTOR_CHECK(vector); GT_ZERO_CHECK(element_size);
  GT_VECTOR_RANGE_CHECK(vector,position);
  return vector->memory+(position*element_size);
}
