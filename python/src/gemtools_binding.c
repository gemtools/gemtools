/*
GEMTools python binding utilities
 */
#include "gemtools_binding.h"
#include <omp.h>

void gt_merge_files(char* const input_1, char* const input_2, char* const output_file_name, bool const same_content, uint64_t threads) {

  char* i1 = malloc(1+strlen(input_1) * sizeof(char));
  char* i2 = malloc(1+strlen(input_2) * sizeof(char));
  char* out = malloc(1+strlen(output_file_name) * sizeof(char));
  strcpy(i1, input_1);
  strcpy(i2, input_2);
  strcpy(out, output_file_name);
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(i1, false);
  gt_input_file* input_file_2 = gt_input_file_open(i2, false);
  gt_output_file* output_file = gt_output_file_new(out, SORTED_FILE);
  // Mutex
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  // Parallel reading+process
  #pragma omp parallel num_threads(threads)
  {
    gt_merge_map_files(&input_mutex,input_file_1,input_file_2,false, same_content ,output_file);
  }
  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
