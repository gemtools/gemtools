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

inline gt_vector* gt_vector_new(size_t num_initial_elements,size_t element_size) {
  gt_vector* vector=calloc(1,sizeof(vector));
  gt_cond_fatal_error(!vector,MEM_HANDLER);
  vector->element_size=element_size;
  vector->elements_allocated=num_initial_elements;
  vector->memory=malloc(num_initial_elements*element_size);
  gt_cond_fatal_error(!vector->memory,MEM_ALLOC);
  return vector;
}
inline gt_status gt_vector_reserve(gt_vector* vector,size_t num_elements,bool zero_mem) {
  if (vector->elements_allocated < num_elements) {
    size_t proposed=(float)vector->elements_allocated*GT_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated=num_elements>proposed?num_elements:proposed;
    vector->memory=realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) return GT_VECTOR_FAIL;
  }
  if (__builtin_expect(zero_mem,0)) {
    memset(vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
  return GT_VECTOR_OK;
}
inline gt_status gt_vector_resize__clean(gt_vector* vector,size_t num_elements) {
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

inline void gt_vector_cast__clean(gt_vector* vector,size_t element_size) {
  vector->elements_allocated=(vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size=element_size;
  vector->used=0;
}
inline void gt_vector_delete(gt_vector* vector) {
  free(vector->memory);
  free(vector);
}
