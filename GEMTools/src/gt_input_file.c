/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_file.h"

// Internal constants
#define GT_INPUT_BUFFER_SIZE GT_BUFFER_SIZE_16M

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_input_file* input_file = malloc(sizeof(gt_input_file));
  gt_cond_fatal_error(!input_file,MEM_HANDLER);
  // Input file
  input_file->file_name = GT_STREAM_FILE_NAME;
  input_file->file_type = STREAM;
  input_file->file = stream;
  input_file->fildes = -1;
  input_file->eof = feof(stream);
  input_file->file_size = UINT64_MAX;
  input_file->file_format = UNKNOWN;
  gt_cond_fatal_error(pthread_mutex_init(&input_file->input_mutex, NULL),SYS_MUTEX_INIT);
  // Auxiliary Buffer (for synch purposes)
  input_file->file_buffer = malloc(GT_INPUT_BUFFER_SIZE);
  gt_cond_fatal_error(!input_file->file_buffer,MEM_ALLOC);
  input_file->buffer_size = 0;
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  // Detect file format
  gt_input_file_detect_file_format(input_file);
  return input_file;
}
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file) {
  GT_NULL_CHECK(file_name);
  // Allocate handler
  gt_input_file* input_file = malloc(sizeof(gt_input_file));
  gt_cond_fatal_error(!input_file,MEM_HANDLER);
  // Input file
  struct stat stat_info;
  gt_cond_fatal_error(stat(file_name,&stat_info)==-1,FILE_STAT,file_name);
  input_file->file_name = file_name;
  input_file->file_size = stat_info.st_size;
  input_file->eof = (input_file->file_size==0);
  input_file->file_format = UNKNOWN;
  gt_cond_fatal_error(pthread_mutex_init(&input_file->input_mutex,NULL),SYS_MUTEX_INIT);
  if (mmap_file) {
    input_file->file = NULL;
    input_file->fildes = open(file_name,O_RDONLY,0); // TODO: O_NOATIME condCompl (Thanks Jordi Camps)
    gt_cond_fatal_error(input_file->fildes==-1,FILE_OPEN,file_name);
    input_file->file_buffer =
      (uint8_t*) mmap(0,input_file->file_size,PROT_READ,MAP_PRIVATE,input_file->fildes,0);
    gt_cond_fatal_error(input_file->file_buffer==MAP_FAILED,SYS_MMAP,file_name);
    input_file->file_type = MAPPED_FILE;
    input_file->buffer_size = input_file->file_size;
  } else {
    input_file->fildes = -1;
    gt_cond_fatal_error(!(input_file->file=fopen(file_name,"r")),FILE_OPEN,file_name);
    input_file->file_type = REGULAR_FILE;
    input_file->file_buffer = malloc(GT_INPUT_BUFFER_SIZE);
    gt_cond_fatal_error(!input_file->file_buffer,MEM_ALLOC);
    input_file->buffer_size = 0;
  }
  // Auxiliary Buffer (for synch purposes)
  input_file->buffer_begin = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  // Detect file format
  gt_input_file_detect_file_format(input_file);
  return input_file;
}
/*
 * POST: Closes the gt_input_file
 * RETURN VALUE: Returns zero on success and error code
 */
gt_status gt_input_file_close(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_status status = GT_INPUT_FILE_OK;
  switch (input_file->file_type) {
    case REGULAR_FILE:
      free(input_file->file_buffer);
      if (fclose(input_file->file)) status = GT_INPUT_FILE_CLOSE_ERR;
      break;
    case MAPPED_FILE:
      gt_cond_error(munmap(input_file->file,input_file->file_size)==-1,SYS_UNMAP);
      if (close(input_file->fildes)) status = GT_INPUT_FILE_CLOSE_ERR;
      break;
    case STREAM:
      free(input_file->file_buffer);
      break;
  }
  free(input_file);
  return status;
}

/*
 * Internal buffer handlers
 */
GT_INLINE size_t gt_input_file_dump_to_buffer(
    gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  // Copy internal file buffer to buffer_dst
  register const uint64_t chunk_size = input_file->buffer_pos-input_file->buffer_begin;
  if (gt_expect_false(chunk_size==0)) return 0;
  gt_vector_reserve_additional(buffer_dst,chunk_size);
  memcpy(gt_vector_get_mem(buffer_dst,uint8_t)+gt_vector_get_used(buffer_dst),
      input_file->file_buffer+input_file->buffer_begin,chunk_size);
  gt_vector_add_used(buffer_dst,chunk_size);
  // Update position
  input_file->buffer_begin=input_file->buffer_pos;
  // Return number of written bytes
  return chunk_size;
}
GT_INLINE size_t gt_input_file_fill_buffer(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  input_file->global_pos += input_file->buffer_size;
  input_file->buffer_pos = 0;
  input_file->buffer_begin = 0;
  if (gt_expect_true(
      (input_file->file_type==STREAM && !feof(input_file->file)) ||
      (input_file->global_pos < input_file->file_size))) {
    input_file->buffer_size =
        fread(input_file->file_buffer,sizeof(uint8_t),GT_INPUT_BUFFER_SIZE,input_file->file);
    if (input_file->buffer_size==0) {
      input_file->eof = true;
    }
    return input_file->buffer_size;
  } else if (input_file->file_type==MAPPED_FILE && input_file->global_pos < input_file->file_size) {
    return input_file->file_size-input_file->global_pos;
  } else {
    input_file->eof = true;
    return 0;
  }
}
GT_INLINE size_t gt_input_file_next_line(
    gt_input_file* const input_file,gt_vector* const buffer_dst) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  GT_INPUT_FILE_CHECK__FILL_BUFFER(input_file,buffer_dst);
  if (input_file->eof) return GT_INPUT_FILE_EOF;
  while (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)!=EOL)) {
    GT_INPUT_FILE_NEXT_CHAR(input_file,buffer_dst);
  }
  if (!input_file->eof) {
    // GT_INPUT_FILE_CURRENT_CHAR(input_file) = EOL;
    GT_INPUT_FILE_NEXT_CHAR(input_file,buffer_dst); // Check DOS EOF
    if (gt_expect_true(!input_file->eof && GT_INPUT_FILE_CURRENT_CHAR(input_file)==DOS_EOL)) {
      // GT_INPUT_FILE_CURRENT_CHAR(input_file) = EOS;
      GT_INPUT_FILE_NEXT_CHAR(input_file,buffer_dst);
    }
  }
  return GT_INPUT_FILE_LINE_READ;
}

/*
 * Format detection
 */
// Forward declarations (gt_file_format_test_<FORMAT> in each logic module)
GT_INLINE bool gt_input_file_test_map(
    gt_input_file* const input_file,gt_map_file_format* const map_file_format,const bool show_errors);
// Format detector (cascade of checkers)
gt_file_format gt_input_file_detect_file_format(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  if (input_file->file_format != UNKNOWN) return input_file->file_format;
  // Try to determine the file format
  gt_input_file_fill_buffer(input_file);
  // MAP test
  if (gt_input_file_test_map(input_file,&(input_file->map_type),false)) {
    input_file->file_format = MAP;
    return MAP;
  }
  // TODO: Install SAM
  // TODO: Install FASTQ
  // Unknown format
  // gt_error(FILE_FORMAT);
  return UNKNOWN;
}

/*
 * Advanced I/O
 */
gt_input_file* gt_input_file_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t segment_number,const uint64_t total_segments) {
  // TODO
  return NULL;
}
gt_input_file* gt_input_file_reads_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t num_init_line,const uint64_t num_end_line) {
  // TODO
  return NULL;
}

/*
 * Mutex/ID functions
 */
GT_INLINE void gt_input_file_lock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_lock(&input_file->input_mutex),SYS_MUTEX);
}
GT_INLINE void gt_input_file_unlock(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  gt_cond_fatal_error(pthread_mutex_unlock(&input_file->input_mutex),SYS_MUTEX);
}
GT_INLINE uint64_t gt_input_file_next_id(gt_input_file* const input_file) {
  GT_INPUT_FILE_CHECK(input_file);
  return (input_file->processed_id)++;
}

/*
 * Reading from input (NOT thread safe, must call mutex functions before)
 */
GT_INLINE uint64_t gt_input_file_add_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  register uint64_t lines_read = 0;
  // Clear dst buffer
  gt_vector_clean(buffer_dst);
  // Read lines
  while (lines_read<num_lines && gt_input_file_next_line(input_file,buffer_dst)) {
    ++lines_read;
  }
  // Dump remaining content into the buffer
  gt_input_file_dump_to_buffer(input_file,buffer_dst);
  if (lines_read > 0 && *gt_vector_get_last_elm(buffer_dst,char) != EOL) {
    gt_vector_insert(buffer_dst,EOL,char);
  }
  input_file->processed_lines+=lines_read;
  return lines_read;
}
GT_INLINE uint64_t gt_input_file_get_lines(
    gt_input_file* const input_file,gt_vector* buffer_dst,const uint64_t num_lines) {
  GT_INPUT_FILE_CHECK(input_file);
  GT_VECTOR_CHECK(buffer_dst);
  // Clear dst buffer
  gt_vector_clean(buffer_dst);
  // Read lines
  return gt_input_file_add_lines(input_file,buffer_dst,num_lines);
}
