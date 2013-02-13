
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

void gt_merge_files_synch(gt_output_file* const output_file, uint64_t threads, const uint64_t num_files,  gt_input_file** files);

#endif /* GEMTOOLS_BINDING_H */
