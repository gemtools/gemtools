
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

void gt_merge_files_synch(gt_output_file* const output_file, uint64_t threads, const uint64_t num_files,  gt_input_file** files);
void gt_write_stream(gt_output_file* output, gt_input_file** inputs, uint64_t num_inputs, bool append_extra, bool clean_id, bool interleave, uint64_t threads, bool write_map);
void gt_stats_fill(gt_input_file* input_file, gt_stats* target_stats, uint64_t num_threads, bool paired_end, bool best_map);
bool gt_input_file_has_qualities(gt_input_file* file);
void gt_stats_print_stats(FILE* output, gt_stats* const stats, const bool paired_end);
#endif /* GEMTOOLS_BINDING_H */
