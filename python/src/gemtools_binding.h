
#ifndef GEMTOOLS_BINDING_H
#define GEMTOOLS_BINDING_H

#include "gem_tools.h"

#define get_mapq(score) ((int)floor((sqrt(score)/256.0)*255))

typedef struct {
  uint64_t max_matches;
  uint64_t min_event_distance;
  uint64_t max_event_distance;
  uint64_t min_levenshtein_distance;
  uint64_t max_levenshtein_distance;
  /* Filter-pairs */
  int64_t max_inss;
  int64_t min_inss;
  bool filter_by_strand;
  bool keep_unique;
  uint64_t min_score;
  bool filter_groups;
  bool group_1;
  bool group_2;
  bool group_3;
  bool group_4;
  bool paired;
  bool close_output;
  char* annotation;
} gt_filter_params;


void gt_merge_files_synch(gt_output_file* const output_file, uint64_t threads, const uint64_t num_files,  gt_input_file** files);
void gt_write_stream(gt_output_file* output, gt_input_file** inputs, uint64_t num_inputs, bool append_extra, bool clean_id, bool interleave, uint64_t threads, bool write_map, bool remove_scores);
void gt_stats_fill(gt_input_file* input_file, gt_stats* target_all_stats, gt_stats* target_best_stats, uint64_t num_threads, bool paired_end);
bool gt_input_file_has_qualities(gt_input_file* file);
void gt_stats_print_stats(FILE* output, gt_stats* const stats, const bool paired_end);
void gt_filter_stream(gt_input_file* input, gt_output_file* output, uint64_t threads, gt_filter_params* params);
#endif /* GEMTOOLS_BINDING_H */
