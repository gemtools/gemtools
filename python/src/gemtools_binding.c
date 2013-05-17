/*
   GEMTools python binding utilities
   */
#include "gemtools_binding.h"
#include <omp.h>


void gt_stats_fill(gt_input_file* input_file, gt_stats* target_all_stats, gt_stats* target_best_stats, uint64_t num_threads, bool paired_end){
  // Stats info
  gt_stats** all_stats = malloc(num_threads*sizeof(gt_stats*));
  all_stats[0] = target_all_stats;
  gt_stats** best_stats = malloc(num_threads*sizeof(gt_stats*));
  best_stats[0] = target_best_stats;

  gt_stats_analysis params_all = GT_STATS_ANALYSIS_DEFAULT();
  params_all.best_map = false;
  gt_stats_analysis params_best = GT_STATS_ANALYSIS_DEFAULT();
  params_best.best_map = true;

  //params_all.indel_profile = true
  // Parallel reading+process
#pragma omp parallel num_threads(num_threads)
  {
    uint64_t tid = omp_get_thread_num();
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

    gt_status error_code;
    gt_template *template = gt_template_new();
    if(tid > 0){
      all_stats[tid] = gt_stats_new();
      best_stats[tid] = gt_stats_new();
    }
    gt_generic_parser_attr* generic_parser_attr =  gt_input_generic_parser_attributes_new(paired_end);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template, generic_parser_attr))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file\n");
      }
      // Extract all_stats
      if(target_all_stats != NULL){
        gt_stats_calculate_template_stats(all_stats[tid],template,NULL, &params_all);
      }
      if(target_best_stats != NULL){
        gt_stats_calculate_template_stats(best_stats[tid],template,NULL, &params_best);
      }
    }
    // Clean
    gt_template_delete(template);
    // gt_template_delete(template_copy);
    gt_buffered_input_file_close(buffered_input);
  }

  // Merge all_stats
  if(target_all_stats != NULL){
    gt_stats_merge(all_stats, num_threads);
  }
  if(target_best_stats != NULL){
    gt_stats_merge(best_stats, num_threads);
  }

  // Clean
  free(all_stats);
  free(best_stats);
  gt_input_file_close(input_file);
}


bool gt_input_file_has_qualities(gt_input_file* file){
  return (file->file_format == FASTA && file->fasta_type.fasta_format == F_FASTQ) || (file->file_format == MAP && file->map_type.contains_qualities);
}

void gt_merge_files_synch(gt_output_file* const output_file, uint64_t threads, const uint64_t num_files,  gt_input_file** files) {
  // Mutex
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;
  // Parallel reading+process
#pragma omp parallel num_threads(threads)
  {
    //gt_merge_map_files(&input_mutex,input_file_1,input_file_2,false, same_content ,output_file);
    gt_merge_synch_map_files_a(&input_mutex, false, output_file, files, num_files);
  }
}

void gt_write_stream(gt_output_file* output, gt_input_file** inputs, uint64_t num_inputs, bool append_extra, bool clean_id, bool interleave, uint64_t threads, bool write_map, bool remove_scores){
  // prepare attributes

  gt_output_fasta_attributes* attributes = 0;
  gt_output_map_attributes* map_attributes = 0;
  if(!write_map){
    attributes = gt_output_fasta_attributes_new();
    gt_output_fasta_attributes_set_print_extra(attributes, append_extra);
    gt_output_fasta_attributes_set_print_casava(attributes, !clean_id);
    // check qualities
    if(!gt_input_file_has_qualities(inputs[0])){
      gt_output_fasta_attributes_set_format(attributes, F_FASTA);
    }
  }else{

    map_attributes = gt_output_map_attributes_new();
    gt_output_map_attributes_set_print_extra(map_attributes, append_extra);
    gt_output_map_attributes_set_print_casava(map_attributes, !clean_id);
    gt_output_map_attributes_set_print_scores(map_attributes, !remove_scores);
  }

  // generic parser attributes
  gt_generic_parser_attr* parser_attributes = gt_input_generic_parser_attributes_new(false); // do not force pairs
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  if(interleave){


    // main loop, interleave
#pragma omp parallel num_threads(threads)
    {
      register uint64_t i = 0;
      register uint64_t c = 0;
      gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output);
      gt_buffered_input_file** buffered_input = malloc(num_inputs * sizeof(gt_buffered_input_file*));

      for(i=0; i<num_inputs; i++){
        buffered_input[i] = gt_buffered_input_file_new(inputs[i]);
      }
      // attache first input to output
      gt_buffered_input_file_attach_buffered_output(buffered_input[0], buffered_output);

      gt_template* template = gt_template_new();
      gt_status status;
      i=0;
      while( gt_input_generic_parser_synch_blocks_a(&input_mutex, buffered_input, num_inputs, parser_attributes) == GT_STATUS_OK ){
        for(i=0; i<num_inputs; i++){
          if( (status = gt_input_generic_parser_get_template(buffered_input[i], template, parser_attributes)) == GT_STATUS_OK){
            if(write_map){
              gt_output_map_bofprint_template(buffered_output, template, map_attributes);
            }else{
              gt_output_fasta_bofprint_template(buffered_output, template, attributes);
            }
            c++;
          }
        }
      }
      gt_buffered_output_file_close(buffered_output);
      for(i=0; i<num_inputs; i++){
        gt_buffered_input_file_close(buffered_input[i]);
      }
      gt_template_delete(template);
      free(buffered_input);
    }
  }else{
    // main loop, cat
    #pragma omp parallel num_threads(threads)
    {
      register uint64_t i = 0;
      register uint64_t c = 0;
      register uint64_t last_id = 0;
      gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output);
      gt_buffered_input_file** buffered_input = malloc(1 * sizeof(gt_buffered_input_file*));
      gt_buffered_input_file* current_input = 0;
      gt_template* template = gt_template_new();
      gt_status status = 0;
      for(i=0; i<num_inputs;i++){
        // create input buffer
        if(i>0){
          gt_buffered_input_file_close(current_input);
          inputs[i]->processed_id = last_id;
        }
        current_input = gt_buffered_input_file_new(inputs[i]);
        buffered_input[0] = current_input;
        // attache the buffer
        gt_buffered_input_file_attach_buffered_output(current_input, buffered_output);

        // read
        while( gt_input_generic_parser_synch_blocks_a(&input_mutex, buffered_input, 1, parser_attributes) == GT_STATUS_OK ){
          if( (status = gt_input_generic_parser_get_template(current_input, template, parser_attributes)) == GT_STATUS_OK){
            if(write_map){
              gt_output_map_bofprint_template(buffered_output, template, map_attributes);
            }else{
              gt_output_fasta_bofprint_template(buffered_output, template, attributes);
            }
            c++;
          }
        }
        last_id = inputs[i]->processed_id;
      }
      gt_buffered_input_file_close(current_input);
      gt_buffered_output_file_close(buffered_output);
      gt_template_delete(template);
      free(buffered_input);
    }

  }
  if(attributes != NULL) gt_output_fasta_attributes_delete(attributes);
  if(map_attributes != NULL)gt_output_map_attributes_delete(map_attributes);
  gt_input_generic_parser_attributes_delete(parser_attributes);

  // register uint64_t i = 0;
  // for(i=0; i<num_inputs; i++){
  //     gt_input_file_close(inputs[i]);
  // }
  // gt_output_file_close(output);
}

void gt_stats_print_stats(FILE* output, gt_stats* const stats, const bool paired_end) {
  register uint64_t num_reads = stats->num_blocks;
  /*
   * General.Stats (Reads,Alignments,...)
   */
  fprintf(output,"[GENERAL.STATS]\n");
  gt_stats_print_general_stats(output,stats,num_reads,paired_end);
  /*
   * Maps
   */
  // if (parameters.maps_profile) {
  fprintf(output,"[MAPS.PROFILE]\n");
  gt_stats_print_maps_stats(output, stats,num_reads,paired_end);
  // }
  if (paired_end) {
    gt_stats_print_inss_distribution(output,stats->maps_profile->inss,stats->num_maps);
  }
  /*
   * Print Quality Scores vs Errors/Misms
   */
  // if (parameters.mismatch_quality) {
  {
    register const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(output,"[MISMATCH.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          output,maps_profile->qual_score_misms,maps_profile->total_mismatches);
    }
    if (maps_profile->total_errors_events > 0) {
      fprintf(output,"[ERRORS.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          output,maps_profile->qual_score_errors,maps_profile->total_errors_events);
    }
  }
  // }
  /*
   * Print Mismatch transition table
   */
  // if (parameters.mismatch_transitions) {
  {
    register const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(output,"[MISMATCH.TRANSITIONS]\n");
      fprintf(output,"MismsTransitions\n");
      gt_stats_print_misms_transition_table(
          output,maps_profile->misms_transition,maps_profile->total_mismatches);
      fprintf(output,"MismsTransitions.1-Nucleotide.Context");
      gt_stats_print_misms_transition_table_1context(
          output,maps_profile->misms_1context,maps_profile->total_mismatches);
    }
  }
  // }
  /*
   * Print Splitmaps profile
   */
  // if (parameters.splitmaps_profile) {
  fprintf(output,"[SPLITMAPS.PROFILE]\n");
  gt_stats_print_split_maps_stats(output,stats, paired_end);
  // }
}

void gt_score_filter(gt_template* template_dst,gt_template* template_src, gt_filter_params* params) {
  register bool is_4 = false;
  register bool best_printed = false;
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src, alignment_src) {
    GT_TEMPLATE_REDUCTION(template_dst, alignment_dst);
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      const int64_t score = get_mapq(map->gt_score);

      if(   (params->group_1 && 252 <= score && score <= 254)
          ||	(params->group_2 && 177 <= score && score <= 180)
          ||	(params->group_3 && 123 <= score && score <= 127)
          ||	(params->group_4 && (  (114 <= score && score <= 119)
              ||(95  <= score && score <= 110 && is_4)))

        ) { 
        if (!is_4 && 114 <= score && score <= 119){
          is_4= true;
        }
        if(!best_printed) gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
        if(score > 119){
          best_printed = true;
        }
      }
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;

  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  is_4 = false;
  best_printed = false;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attr) {
    register const int64_t score = get_mapq(mmap_attr->gt_score);

    if(   (params->group_1 && 252 <= score && score <= 254)
        ||	(params->group_2 && 177 <= score && score <= 180)
        ||	(params->group_3 && 123 <= score && score <= 127)
        ||	(params->group_4 && (  (114 <= score && score <= 119)
            ||(95  <= score && score <= 110 && is_4)))

      ) { 
      if (!is_4 && 114 <= score && score <= 119){
        is_4= true;
      }
      register gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
      if(!best_printed) gt_template_add_mmap(template_dst,mmap_copy);
      if(score > 119){
        best_printed = true;
      }
      free(mmap_copy);
    }
  }

}

void gt_annotation_filter(gt_template* template_dst,gt_template* template_src, gt_filter_params* params, gt_gtf* gtf) {
  printf("Doing the annotaiton filtering ...");
  gt_gtf_hits* hits = gt_gtf_hits_new();
  gt_gtf_search_template_for_exons(gtf, hits, template_src);

  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src, alignment_src) {
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;

  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attr) {
    GT_MMAP_ITERATE(mmap, map, end_position){
      printf("Annotation mapped :) \n");
    }
  }

}

void gt_template_filter(gt_template* template_dst,gt_template* template_src, gt_filter_params* params) {
  /*SE*/
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_REDUCTION(template_src,alignment_dst);
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      // Check SM contained
      register const uint64_t num_blocks = gt_map_get_num_blocks(map);
      register const int64_t total_distance = gt_map_get_global_distance(map);
      register const int64_t lev_distance = gt_map_get_global_levenshtein_distance(map);
      if (params->min_event_distance > total_distance || total_distance > params->max_event_distance) continue;
      if (params->min_levenshtein_distance > lev_distance || lev_distance > params->max_levenshtein_distance) continue;
      //if (0 > lev_distance || lev_distance > 2) continue;
      // Insert the map
      gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
      // Skip the rest if best
      //if (parameters.best_map) return;
    }
  } GT_TEMPLATE_END_REDUCTION__RETURN;

  /*
   * PE
   */
  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attr) {
    if(params->min_score > 0){

    }
    // Check strata
    register const int64_t total_distance = gt_map_get_global_distance(mmap[0])+gt_map_get_global_distance(mmap[1]);
    if (params->min_event_distance > total_distance || total_distance > params->max_event_distance) continue;

    // Check levenshtein distance
    register const int64_t lev_distance = gt_map_get_global_levenshtein_distance(mmap[0])+gt_map_get_global_levenshtein_distance(mmap[1]);
    if (params->min_levenshtein_distance > lev_distance || lev_distance > params->max_levenshtein_distance) continue;

    // Check inss
    if (params->min_inss > INT64_MIN || params->max_inss < INT64_MAX) {
      uint64_t gt_err;
      register const int64_t inss = gt_template_get_insert_size(mmap,&gt_err);
      if (params->min_inss > inss || inss > params->max_inss) continue;
    }
    // Check strandness
    if (params->filter_by_strand) {
      if (mmap[0]->strand==FORWARD && mmap[1]->strand==FORWARD) continue;
      if (mmap[0]->strand==REVERSE && mmap[1]->strand==REVERSE) continue;
    }
    // Add the mmap
    register gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
    gt_template_insert_mmap(template_dst,mmap_copy,mmap_attr);
    free(mmap_copy);
    // Skip the rest if best
    //      if (parameters.best_map) return;
  }
}

void gt_filter_stream(gt_input_file* input, gt_output_file* output, uint64_t threads, gt_filter_params* params){
  // prepare attributes
  gt_output_map_attributes* map_attributes = gt_output_map_attributes_new();
  gt_output_map_attributes_set_print_extra(map_attributes, true);
  gt_output_map_attributes_set_print_casava(map_attributes, true);

  // generic parser attributes
  gt_generic_parser_attr* parser_attributes = gt_input_generic_parser_attributes_new(false); // do not force pairs

  if(params->max_matches > 0){
    parser_attributes->map_parser_attr.max_parsed_maps = params->max_matches;
  }
  gt_gtf* gtf;

  if(params->annotation != NULL && strlen(params->annotation) > 0){
    FILE* of = fopen(params->annotation, "r");
    if(of == NULL){
      printf("ERROR opening annotation !\n");
      return;
    }
    printf("Read annotation \n");
    gtf = gt_gtf_read(of);
    fclose(of);
  }
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  // main loop, cat
  #pragma omp parallel num_threads(threads)
  {
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output);
    gt_buffered_input_file_attach_buffered_output(buffered_input, buffered_output);
    gt_template* template = gt_template_new();
    gt_status status = 0;
    register bool is_mapped = false;
    

    while( (status = gt_input_generic_parser_get_template(buffered_input,template, parser_attributes)) == GT_STATUS_OK ){
      gt_template_sort_by_distance__score(template);
      if(gtf != NULL){
        gt_template *template_filtered = gt_template_copy(template,false,false);
        gt_annotation_filter(template_filtered, template, params, gtf);
        gt_template_delete(template);
        template = template_filtered;
        gt_template_recalculate_counters(template);
      }


      if(params->filter_groups){
        gt_template *template_filtered = gt_template_copy(template,false,false);
        gt_score_filter(template_filtered, template, params);
        gt_template_delete(template);
        template = template_filtered;
        gt_template_recalculate_counters(template);
      }
      is_mapped = gt_template_is_mapped(template);
      register const bool is_unique = gt_template_get_not_unique_flag(template);

      if(is_mapped && (!is_unique || !params->keep_unique )){
        gt_template *template_filtered = gt_template_copy(template,false,false);
        gt_template_filter(template_filtered,template, params);
        gt_template_delete(template);
        template = template_filtered;
        gt_template_recalculate_counters(template);
      }
      gt_output_map_bofprint_template(buffered_output, template, map_attributes);
    }
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
    gt_template_delete(template);
  }
  gt_output_map_attributes_delete(map_attributes);
  gt_input_generic_parser_attributes_delete(parser_attributes);
  if(params->close_output){
    gt_output_file_close(output);
  }
}
