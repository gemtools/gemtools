/*
GEMTools python binding utilities
 */
#include "gemtools_binding.h"
#include <omp.h>

void gt_merge_files(char* const input_1, char* const input_2, char* const output_file_name, bool const mmap_input, bool const same_content, bool const paired_reads, uint64_t threads) {
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(input_1,mmap_input);
  gt_input_file* input_file_2 = gt_input_file_open(input_2,mmap_input);
  gt_output_file* output_file = gt_output_file_new(output_file_name, SORTED_FILE);
  // Mutex
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  // Parallel reading+process
  #pragma omp parallel num_threads(threads)
  {
    gt_merge_map_files(&input_mutex,input_file_1,input_file_2,paired_reads, same_content ,output_file);
  }
  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
