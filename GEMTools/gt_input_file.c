/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_file.h"

// Internal constants
#define GT_INPUT_STREAM_FILE_NAME "<<STREAM>>"
#define GT_INPUT_BUFFER_SIZE GT_BUFFER_SIZE_8K
#define GT_INPUT_BUFFER_NUM_LINES GT_NUM_LINES_10K

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream) {
  // Allocate handler
  gt_input_file* input_file = malloc(sizeof(gt_input_file));
  gt_cond_fatal_error(!input_file,MEM_HANDLER);
  // Input file
  input_file->file_name = GT_INPUT_STREAM_FILE_NAME;
  input_file->file_type = STREAM;
  input_file->file = stream;
  input_file->fildes = -1;
  input_file->eof = feof(stream);
  input_file->file_size = UINT64_MAX;
  input_file->file_format = UNKNOWN;
  pthread_mutex_init(&input_file->input_mutex, NULL);
  // Auxiliary Buffer (for synch purposes)
  input_file->file_buffer = malloc(GT_INPUT_BUFFER_SIZE);
  gt_cond_fatal_error(!input_file->file_buffer,MEM_ALLOC);
  input_file->buffer_size = 0;
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  return input_file;
}
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file) {
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
  pthread_mutex_init(&input_file->input_mutex,NULL);
  if (mmap_file) {
    input_file->file = NULL;
    input_file->fildes = open(file_name,O_RDONLY | O_NOATIME,0); // Thanks Jordi Camps
    gt_cond_fatal_error(input_file->fildes==-1,FILE_OPEN,file_name);
    input_file->file_buffer =
      (uint8_t*) mmap(0,input_file->file_size,PROT_READ,MAP_PRIVATE,input_file->fildes,0);
    gt_cond_fatal_error(input_file->mapped_file==MAP_FAILED,SYS_MMAP,file_name);
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
  input_file->buffer_pos = 0;
  input_file->global_pos = 0;
  input_file->processed_lines = 0;
  // ID generator
  input_file->processed_id = 0;
  return input_file;
}
/*
 * POST: Closes the gt_input_file
 * RETURN VALUE: Returns zero on success and error code
 */
void gt_input_file_close(gt_input_file* const input_file) {
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
 * Format detection
 */
gt_file_format gt_input_file_test_fastx(gt_input_file* const input_file) {
  // TODO
}
gt_file_format gt_input_file_test_map(gt_input_file* const input_file) {
  // TODO
}
gt_file_format gt_input_file_get_file_format(gt_input_file* const input_file) {
  if (input_file->file_format != UNKNOWN) return input_file->file_format;
  // Try to determine the file format
  register gt_file_format format;
  // MAP test
  format = gt_input_file_test_map(input_file);
  if (format!=UNKNOWN) {
    input_file->file_format = format;
    return format;
  }
  // FASTX test
  format = gt_input_file_test_fastx(input_file);
  if (format!=UNKNOWN) {
    input_file->file_format = format;
    return format;
  }
  // Unknown format
  gt_error(FILE_FORMAT);
  return UNKNOWN;
}

/*
 * Internal buffer handlers
 */
inline size_t gt_input_file_fill_buffer();
inline size_t gt_input_file_();

/*
 * Advanced I/O
 */
gt_input_file* gt_input_file_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t segment_number,const uint64_t total_segments) {
  // TODO
}
gt_input_file* gt_input_file_reads_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t num_init_line,const uint64_t num_end_line) {
  // TODO
}

/*
 * Mutex/ID functions
 */
inline gt_status gt_input_file_get_lock(gt_input_file* const input_file) {
  pthread_mutex_lock(&(input_file->input_mutex));
}
inline gt_status gt_input_file_release_lock(gt_input_file* const input_file) {
  pthread_mutex_unlock(&(input_file->input_mutex));
}
inline uint64_t gt_input_file_next_id(gt_input_file* const input_file) {
  register const uint64_t id = input_file->processed_id;
  ++input_file->processed_id;
  return id;
}

/*
 * Reading from input (NOT thread safe, must call mutex functions before)
 */
inline size_t gt_input_file_get_lines(
    gt_input_file* const input_file,
    gt_vector* buffer_dst,const uint64_t num_lines) {
  // TODO
}
