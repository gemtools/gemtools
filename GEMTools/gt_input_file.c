/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_file.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_input_file.h"

/*
 * Basic I/O functions
 */
gt_input_file* gt_input_stream_open(FILE* stream) {
  // TODO
}
gt_input_file* gt_input_file_open(char* const file_name,const bool mmap_file) {
  // TODO
}
void gt_input_file_close(gt_input_file* const input_file) {
  // TODO
}

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
 * Advanced I/O
 */
gt_input_file* gt_input_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t segment_number,const uint64_t total_segments) {
  // TODO
}
gt_input_file* gt_input_reads_segmented_file_open(
    char* const file_name,const bool mmap_file,
    const uint64_t num_init_line,const uint64_t num_end_line) {
  // TODO
}

/*
 * Mutex functions
 */
inline gt_status gt_input_get_lock(gt_input_file* const input_file) {
  // TODO
}
inline gt_status gt_input_release_lock(gt_input_file* const input_file) {
  // TODO
}

/*
 * Reading from input (NO thread safe, must call mutex functions before)
 */
inline size_t gt_input_get_lines(gt_input_file* const input_file,gt_vector* buffer_dst) {
  // TODO
}
