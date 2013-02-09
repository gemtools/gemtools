
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

void gt_merge_files(char* const input_1, char* const input_2, char* const output_file,
    bool const mmap_input, bool const same_content, bool const paired_reads, uint64_t threads);

#endif /* GEMTOOLS_BINDING_H */