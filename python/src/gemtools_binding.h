
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

void gt_merge_files(char* const input_1, char* const input_2, char* const output_file,
    bool const same_content,  uint64_t threads);

#endif /* GEMTOOLS_BINDING_H */
