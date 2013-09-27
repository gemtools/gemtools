/*
   GEMTools python binding utilities
   */
#include "gemtools_binding.h"
#include <omp.h>


bool gt_input_file_has_qualities(gt_input_file* file){
  return (file->file_format == FASTA && file->fasta_type.fasta_format == F_FASTQ) || (file->file_format == MAP && file->map_type.contains_qualities);
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
  gt_generic_parser_attributes* parser_attributes = gt_input_generic_parser_attributes_new(false); // do not force pairs
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