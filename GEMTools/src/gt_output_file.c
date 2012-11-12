/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_output_file.h"

/*
 * Checkers
 */
#define GT_OUTPUT_FILE_CHECK(output_file) \
  gt_fatal_check(output_file==NULL|| \
    output_file->file==NULL||output_file->file_name==NULL|| \
    output_file->buffer==NULL,NULL_HANDLER)

#define GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file) \
  GT_OUTPUT_FILE_CHECK(output_file); \
  gt_fatal_check( \
    output_file->buffer_busy>GT_MAX_OUTPUT_BUFFERS|| \
    output_file->buffer_write_pending>GT_MAX_OUTPUT_BUFFERS,OUTPUT_FILE_INCONSISTENCY)

/*
 * Setup
 */
GT_INLINE void gt_output_file_init_buffers(gt_output_file* const output_file) {
  /* Output Buffers */
  register uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS;++i) {
    output_file->buffer[i]=NULL;
  }
  output_file->buffer_busy=0;
  output_file->buffer_write_pending=0;
  /* Block ID (for synchronization purposes) */
  output_file->mayor_block_id=0;
  output_file->minor_block_id=0;
  /* Mutexes */
  gt_cond_fatal_error(pthread_mutex_init(&output_file->out_buffer_mutex, NULL),SYS_MUTEX_INIT);
  gt_cond_fatal_error(pthread_cond_init(&output_file->out_buffer_cond,NULL),SYS_COND_VAR_INIT);
  gt_cond_fatal_error(pthread_mutex_init(&output_file->out_file_mutex, NULL),SYS_MUTEX_INIT);
}

gt_output_file* gt_output_stream_new(FILE* const file,const gt_output_file_type output_file_type) {
  GT_NULL_CHECK(file);
  gt_output_file* output_file = malloc(sizeof(gt_output_file));
  gt_cond_fatal_error(!output_file,MEM_HANDLER);
  /* Output file */
  output_file->file_name=GT_STREAM_FILE_NAME;
  output_file->file=file;
  output_file->file_type=output_file_type;
  /* Setup buffers */
  gt_output_file_init_buffers(output_file);
  return output_file;
}
gt_output_file* gt_output_file_new(char* const file_name,const gt_output_file_type output_file_type) {
  GT_NULL_CHECK(file_name);
  gt_output_file* output_file = malloc(sizeof(gt_output_file));
  gt_cond_fatal_error(!output_file,MEM_HANDLER);
  /* Output file */
  output_file->file_name=file_name;
  gt_cond_fatal_error(!(output_file->file=fopen(file_name,"w")),FILE_OPEN,file_name);
  output_file->file_type=output_file_type;
  /* Setup buffers */
  gt_output_file_init_buffers(output_file);
  return output_file;
}
gt_status gt_output_file_close(gt_output_file* const output_file) {
  GT_BUFFERED_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  register gt_status error_code = 0;
  // Close file/stream
  gt_cond_error(error_code|=fclose(output_file->file),FILE_CLOSE,output_file->file_name);
  // Delete allocated buffers
  register uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&output_file->buffer[i]!=NULL;++i) {
    gt_output_buffer_delete(output_file->buffer[i]);
  }
  // Free mutex/CV
  gt_cond_error(error_code|=pthread_mutex_destroy(&output_file->out_buffer_mutex),SYS_MUTEX_DESTROY);
  gt_cond_error(error_code|=pthread_cond_destroy(&output_file->out_buffer_cond),SYS_COND_VAR_INIT);
  gt_cond_error(error_code|=pthread_mutex_destroy(&output_file->out_file_mutex),SYS_MUTEX_DESTROY);
  // Free handler
  free(output_file);
  return error_code;
}

/*
 * Internal Buffers Accessors
 */
GT_INLINE gt_output_buffer* __gt_buffered_output_file_request_buffer(gt_buffered_output_file* const buffered_output_file) {
  GT_BUFFERED_OUTPUT_FILE_CONSISTENCY_CHECK(buffered_output_file);
  // Conditional guard. Wait till there is any free buffer left
  while (buffered_output_file->buffer_busy==GT_MAX_OUTPUT_BUFFERS) {
    pthread_cond_wait(&buffered_output_file->out_buffer_cond,&buffered_output_file->out_buffer_mutex);
  }
  // There is at least one free buffer. Get it!
  register uint64_t i;
  for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&buffered_output_file->buffer[i]!=NULL;++i) {
    if (buffered_output_file->buffer[i]->buffer_state==GT_OUTPUT_BUFFER_FREE) {
      ++buffered_output_file->buffer_busy;
      gt_output_buffer_initiallize(buffered_output_file->buffer[i],GT_OUTPUT_BUFFER_BUSY);
      buffered_output_file->buffer[i]->buffered_output_file = buffered_output_file; // Attach output file
      return buffered_output_file->buffer[i];
    }
  }
  if (i<GT_MAX_OUTPUT_BUFFERS) { // There is at least one buffer available but not allocated yet
    ++buffered_output_file->buffer_busy;
    buffered_output_file->buffer[i] = gt_output_buffer_new();
    gt_output_buffer_initiallize(buffered_output_file->buffer[i],GT_OUTPUT_BUFFER_BUSY);
    buffered_output_file->buffer[i]->buffered_output_file = buffered_output_file; // Attach output file
    return buffered_output_file->buffer[i];
  }
  gt_fatal_error(ALG_INCONSISNTENCY);
  return NULL;
}
GT_INLINE gt_output_buffer* gt_buffered_output_file_request_buffer(gt_buffered_output_file* const buffered_output_file) {
  GT_BUFFERED_OUTPUT_FILE_CONSISTENCY_CHECK(buffered_output_file);
  register gt_output_buffer* fresh_buffer;
  GT_BEGIN_MUTEX_SECTION(buffered_output_file->out_buffer_mutex) {
    fresh_buffer = __gt_buffered_output_file_request_buffer(buffered_output_file);
  } GT_END_MUTEX_SECTION(buffered_output_file->out_buffer_mutex);
  return fresh_buffer;
}

GT_INLINE void __gt_buffered_output_file_release_buffer(
    gt_buffered_output_file* const buffered_output_file,gt_output_buffer* const output_buffer) {
  GT_BUFFERED_OUTPUT_FILE_CONSISTENCY_CHECK(buffered_output_file);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  // Broadcast. Wake up sleepy.
  if (buffered_output_file->buffer_busy==GT_MAX_OUTPUT_BUFFERS) {
    gt_cond_fatal_error(pthread_cond_broadcast(&buffered_output_file->out_buffer_cond),SYS_COND_VAR);
  }
  // Free buffer
  output_buffer->buffer_state=GT_OUTPUT_BUFFER_FREE;
  --buffered_output_file->buffer_busy;
  gt_output_buffer_clear(output_buffer);
}
GT_INLINE void gt_buffered_output_file_release_buffer(
    gt_buffered_output_file* const buffered_output_file,gt_output_buffer* const output_buffer) {
  GT_BUFFERED_OUTPUT_FILE_CONSISTENCY_CHECK(buffered_output_file);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  GT_BEGIN_MUTEX_SECTION(buffered_output_file->out_buffer_mutex) {
    __gt_buffered_output_file_release_buffer(buffered_output_file,output_buffer);
  } GT_END_MUTEX_SECTION(buffered_output_file->out_buffer_mutex);
}

GT_INLINE gt_output_buffer* gt_output_file_write_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  register int64_t bytes_written;
  gt_vector* vbuffer = gt_output_buffer_to_vchar(output_buffer);
  GT_BEGIN_MUTEX_SECTION(output_file->out_file_mutex)
  {
    bytes_written = fwrite(gt_vector_get_mem(vbuffer,char),1,
        gt_vector_get_used(vbuffer),output_file->file);
  }
  GT_END_MUTEX_SECTION(output_file->out_file_mutex);
  gt_cond_fatal_error(bytes_written!=gt_vector_get_used(vbuffer),OUTPUT_FILE_FAIL_WRITE);
  gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_BUSY);
  return output_buffer;
}
GT_INLINE gt_output_buffer* gt_output_file_sorted_write_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  // Set the block buffer as write pending and set the victim
  register bool victim;
  register uint32_t mayor_block_id, minor_block_id;
  GT_BEGIN_MUTEX_SECTION(output_file->out_buffer_mutex)
  {
    // Set the buffer as write pending
    ++output_file->buffer_write_pending;
    gt_output_buffer_set_state(output_buffer,GT_OUTPUT_BUFFER_WRITE_PENDING);
    victim = (output_file->mayor_block_id==gt_output_buffer_get_mayor_block_id(output_buffer) &&
              output_file->minor_block_id==gt_output_buffer_get_minor_block_id(output_buffer));
    if (!victim) {
      // Enqueue the buffer and continue (someone else will do the writing)
      output_buffer = __gt_buffered_output_file_request_buffer(buffered_output_file);
      GT_END_MUTEX_SECTION(buffered_output_file->out_buffer_mutex);
      return output_buffer;
    } else {
      mayor_block_id = gt_output_buffer_get_mayor_block_id(output_buffer);
      minor_block_id = gt_output_buffer_get_minor_block_id(output_buffer);
    }
  }
  GT_END_MUTEX_SECTION(output_file->out_buffer_mutex);
  // I'm the victim, I will output as much as I can
  do {
    // Write the current buffer
    gt_vector* vbuffer = gt_output_buffer_to_vchar(output_buffer);
    register const int64_t bytes_written =
        fwrite(gt_vector_get_mem(vbuffer,char),1,gt_vector_get_used(vbuffer),buffered_output_file->file);
    if (bytes_written!=gt_vector_get_used(vbuffer)) return NULL;
    // Update buffers' state
    GT_BEGIN_MUTEX_SECTION(buffered_output_file->out_buffer_mutex) {
      // Decrement write-pending blocks and update next block ID (mayorID,minorID)
      --buffered_output_file->buffer_write_pending;
      if (output_buffer->is_final_block) {
        ++mayor_block_id;
        minor_block_id = 0;
      } else {
        ++minor_block_id;
      }
      // Search for the next block buffer in order
      register bool next_found = false;
      if (buffered_output_file->buffer_write_pending>0) {
        register uint64_t i;
        for (i=0;i<GT_MAX_OUTPUT_BUFFERS&&buffered_output_file->buffer[i]!=NULL;++i) {
          if (mayor_block_id==gt_output_buffer_get_mayor_block_id(buffered_output_file->buffer[i]) &&
              minor_block_id==gt_output_buffer_get_minor_block_id(buffered_output_file->buffer[i])) {
            if (buffered_output_file->buffer[i]->buffer_state!=GT_OUTPUT_BUFFER_WRITE_PENDING) {
              // Cannot dump a busy buffer
              break;
            } else {
              // I'm still the victim, free the current buffer and output the new one
              __gt_buffered_output_file_release_buffer(buffered_output_file,output_buffer);
              next_found = true;
              output_buffer = buffered_output_file->buffer[i];
              break;
            }
          }
        }
      }
      if (!next_found) {
        // Fine, I'm done, let's get out of here ASAP
        buffered_output_file->mayor_block_id = mayor_block_id;
        buffered_output_file->minor_block_id = minor_block_id;
        gt_output_buffer_initiallize(output_buffer,GT_OUTPUT_BUFFER_BUSY);
        GT_END_MUTEX_SECTION(buffered_output_file->out_buffer_mutex);
        return output_buffer;
      }
    } GT_END_MUTEX_SECTION(buffered_output_file->out_buffer_mutex);
  } while (true);
}
GT_INLINE gt_output_buffer* gt_output_file_dump_buffer(
    gt_output_file* const output_file,gt_output_buffer* const output_buffer) {
  GT_OUTPUT_FILE_CONSISTENCY_CHECK(output_file);
  switch (output_file->file_type) {
    case SORTED_FILE:
      return gt_output_file_sorted_write_buffer(output_file,output_buffer);
      break;
    case UNSORTED_FILE:
      return gt_output_file_write_buffer(output_file,output_buffer);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_IMPLEMENTED);
      break;
  }
}





