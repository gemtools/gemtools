/*
 * PROJECT: GEM-Tools library
 * FILE: gt_vector.h
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

#ifndef _GT_VECTOR_H_GUARD_
#define _GT_VECTOR_H_GUARD_

#include "gt_commons.h"

// Codes gt_status
#define GT_VECTOR_OK 1
#define GT_VECTOR_FAIL 0

// Simple linear vector for generic type elements
typedef struct {
  void* memory;
  size_t used;
  size_t element_size;
  size_t elements_allocated;
} gt_vector;

// Get the content of the vector
#define gt_vector_get_mem(vector,type) ((type*)((vector)->memory))
#define gt_vector_get_elm(vector,position,type) (gt_vector_get_mem(vector,type)+position)
#define gt_vector_get_free_elm(vector,type) (gt_vector_get_mem(vector,type)+(vector)->used)
// Getters/Setter used elements in the vector
#define gt_vector_get_used(vector) ((vector)->used)
#define gt_vector_set_used(vector,total_used) (vector)->used=(total_used)
#define gt_vector_inc_used(vector) (++((vector)->used))
#define gt_vector_dec_used(vector) (--((vector)->used))
#define gt_vector_add_used(vector,additional) gt_vector_set_used(vector,gt_vector_get_used(vector)+additional)
#define gt_vector_clean(vector) (vector)->used=0
// Initialization and allocation
#define gt_vector_reserve_additional(vector,additional) gt_vector_reserve(vector,gt_vector_get_used(vector)+additional)
#define gt_vector_prepare(vector,data_type,num_elements) \
  ({ gt_vector_cast__clean(vector,sizeof(data_type)); \
     gt_vector_reserve(vector,num_elements); })
// Macro generic iterator
//  GT_VECTOR_ITERATE(vector_of_ints,elm_iterator,elm_counter,int) {
//    ..code..
//  }
#define GT_VECTOR_ITERATE(vector,element,counter,type) \
  register uint64_t counter; \
  register type* element = gt_vector_get_mem(vector,type); \
  for (counter=0;counter<vector_get_used(vector);++element,++counter)
// Add element to the vector (at the end)
#define gt_vector_insert(vector,element,type) \
  gt_vector_reserve_additional(vector,1); \
  *(gt_vector_get_elm(vector,gt_vector_get_used(vector),type))=element; \
  gt_vector_inc_used(vector)

inline gt_vector* gt_vector_new(size_t num_initial_elements,size_t element_size);
inline gt_status gt_vector_reserve(gt_vector* vector,size_t num_elements,bool zero_mem);
inline gt_status gt_vector_resize__clean(gt_vector* vector,size_t num_elements);

inline void gt_vector_cast__clean(gt_vector* vector,size_t element_size);
inline void gt_vector_delete(gt_vector* vector);

#endif /* _GT_VECTOR_H_GUARD_ */
