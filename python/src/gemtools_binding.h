
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

#define get_mapq(score) ((int)floor((sqrt(score)/256.0)*255))

void gt_write_stream(gt_output_file* output, gt_input_file** inputs, uint64_t num_inputs, bool append_extra, bool clean_id, bool interleave, uint64_t threads, bool write_map, bool remove_scores);
bool gt_input_file_has_qualities(gt_input_file* file);
#endif /* GEMTOOLS_BINDING_H */
