/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_buffer.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"

#define GT_OUTPUT_BUFFER_INITIAL_SIZE GT_BUFFER_SIZE_2M
#define GT_OUTPUT_BUFFER_FORCE_DUMP_SIZE GT_BUFFER_SIZE_8M

/*
 * Setup
 */
GT_INLINE gt_output_buffer* gt_output_buffer_new(void) {
  gt_output_buffer* output_buffer = malloc(sizeof(gt_output_buffer));
  gt_cond_fatal_error(!output_buffer,MEM_HANDLER);
  gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_FREE);
  output_buffer->buffer=gt_vector_new(GT_OUTPUT_BUFFER_INITIAL_SIZE,sizeof(char));
  output_buffer->buffered_output_file=NULL;
  return output_buffer;
}
GT_INLINE void gt_output_buffer_clear(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->mayor_block_id=UINT32_MAX;
  output_buffer->minor_block_id=0;
  output_buffer->is_final_block=true;
  gt_vector_clean(output_buffer->buffer);
}
GT_INLINE void gt_output_buffer_initiallize(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_output_buffer_clear(output_buffer);
  gt_output_buffer_set_state(output_buffer,buffer_state);
}
GT_INLINE void gt_output_buffer_delete(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_vector_delete(output_buffer->buffer);
  free(output_buffer);
}

/*
 * Accessors
 */
GT_INLINE void gt_output_buffer_set_state(gt_output_buffer* const output_buffer,const gt_output_buffer_state buffer_state) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->buffer_state=buffer_state;
}
GT_INLINE gt_output_buffer_state gt_output_buffer_get_state(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->buffer_state;
}
GT_INLINE void gt_output_buffer_set_mayor_block_id(gt_output_buffer* const output_buffer,const uint32_t mayor_block_id) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->mayor_block_id=mayor_block_id;
}
GT_INLINE uint32_t gt_output_buffer_get_mayor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->mayor_block_id;
}
GT_INLINE void gt_output_buffer_set_minor_block_id(gt_output_buffer* const output_buffer,const uint32_t minor_block_id) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->minor_block_id=minor_block_id;
}
GT_INLINE void gt_output_buffer_inc_minor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  ++output_buffer->minor_block_id;
}
GT_INLINE uint32_t gt_output_buffer_get_minor_block_id(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->minor_block_id;
}
GT_INLINE void gt_output_buffer_set_partial_block(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  output_buffer->is_final_block=false;
}

/*
 * Adaptors
 */
GT_INLINE char* gt_output_buffer_to_char(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  gt_vector_insert(output_buffer->buffer,EOS,char);
  return gt_vector_get_mem(output_buffer->buffer,char);
}
GT_INLINE gt_vector* gt_output_buffer_to_vchar(gt_output_buffer* const output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  return output_buffer->buffer;
}

/*
 * Buffer printer
 */
GT_INLINE void gt_output_buffer_check_safety_dump(gt_output_buffer** output_buffer) {
  GT_OUTPUT_BUFFER_CHECK(*output_buffer);
  if ((*output_buffer)->buffered_output_file!=NULL &&
      gt_vector_get_used((*output_buffer)->buffer)>GT_OUTPUT_BUFFER_FORCE_DUMP_SIZE) {
    gt_output_buffer_set_partial_block(*output_buffer);
    *output_buffer = gt_buffered_output_file_dump_buffer(
      (gt_buffered_output_file*)(*output_buffer)->buffered_output_file,*output_buffer);
    gt_cond_fatal_error(*output_buffer==NULL,BUFFER_SAFETY_DUMP);
    gt_output_buffer_clear(*output_buffer);
    gt_output_buffer_inc_minor_block_id(*output_buffer);
  }
}
GT_INLINE gt_status gt_vbprintf(gt_output_buffer** const output_buffer,const char *template,va_list v_args) {
  GT_OUTPUT_BUFFER_CHECK(*output_buffer);
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Calculate buffer size required to dump template{v_args}
  register int64_t mem_required = gt_calculate_memory_required_v(template,v_args);
  return gt_vbprintf_mem(output_buffer,mem_required,template,v_args);
}
GT_INLINE gt_status gt_bprintf(gt_output_buffer** const output_buffer,const char *template,...) {
  GT_OUTPUT_BUFFER_CHECK(*output_buffer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  return gt_vbprintf(output_buffer,template,v_args);
}
GT_INLINE gt_status gt_vbprintf_mem(
    gt_output_buffer** const output_buffer,const uint64_t expected_mem_usage,
    const char *template,va_list v_args) {
  GT_OUTPUT_BUFFER_CHECK(*output_buffer);
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Check buffer size. In case, we do an emergency dump
  gt_output_buffer_check_safety_dump(output_buffer);
  // Calculate buffer size required to dump template{v_args}
  gt_vector_reserve_additional((*output_buffer)->buffer,expected_mem_usage);
  register int64_t mem_used=vsprintf(gt_vector_get_free_elm((*output_buffer)->buffer,char),template,v_args);
  if (gt_expect_true(mem_used>0)) {
    gt_vector_add_used((*output_buffer)->buffer,mem_used);
    return mem_used;
  } else {
    return GT_OUT_BUFFER_FAIL;
  }
}
GT_INLINE gt_status gt_bprintf_mem(
    gt_output_buffer** const output_buffer,const uint64_t expected_mem_usage,const char *template,...) {
  GT_OUTPUT_BUFFER_CHECK(*output_buffer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  return gt_vbprintf_mem(output_buffer,expected_mem_usage,template,v_args);
}

/*
 * Memory usage helper functions
 */
GT_INLINE uint64_t gt_calculate_memory_required_v(const char *template,va_list v_args) {
  GT_NULL_CHECK(template); GT_NULL_CHECK(v_args);
  // Calculate memory required to print the template{v_args}
  register uint64_t mem_required = 0;
  register const char* centinel;
  for (centinel=template;*centinel!=EOS;++centinel,++mem_required) {
    if (*centinel==FORMAT) {
      ++centinel;
      // Read modifiers
      while (gt_is_number(*centinel)) ++centinel;
      if (*centinel==DOT){
        ++centinel;
        if (*centinel==STAR) {
          ++centinel;
        } else {
          while (gt_is_number(*centinel)) ++centinel;
        }
      }
      gt_check(centinel==EOS,PRINT_FORMAT);
      // Check format
      switch (*centinel) {
        case 's': { // String requires fetching the argument length
          register char* const string = va_arg(v_args,char*);
          mem_required+=strlen(string);
          break;
        }
        default:
          // As for the rest, we estimate the memory usage
          // Also we assume an upper bound over the possible formats (int, chars, floats, ...)
          va_arg(v_args,int);
          mem_required+=20;
          break;
      }
    }
  }
  return mem_required;
}
GT_INLINE uint64_t gt_calculate_memory_required_va(const char *template,...) {
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  return gt_calculate_memory_required_v(template,v_args);
}
