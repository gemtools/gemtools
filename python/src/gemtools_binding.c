/*
GEMTools python binding utilities
 */
#include "gemtools_binding.h"
#include <omp.h>

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

void gt_write_stream(gt_output_file* output, gt_input_file** inputs, uint64_t num_inputs, bool append_extra, bool clean_id, bool interleave, uint64_t threads, bool write_map){
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
    }

    // generic parser attributes
    gt_generic_parser_attr* parser_attributes = gt_input_generic_parser_attributes_new(false); // do not force pairs


    if(true){
        pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

        // main loop, interleave
        #pragma omp parallel num_threads(threads)
        {
            register uint64_t i = 0;
            register uint64_t c = 0;
            gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output);
            gt_buffered_input_file** buffered_input = malloc(num_inputs * sizeof(gt_buffered_input_file*));

            for(i=0; i<num_inputs; i++){
                buffered_input[i] = gt_buffered_input_file_new(inputs[ num_inputs - i - 1 ]);
            }
            // attache first input to output
            gt_buffered_input_file_attach_buffered_output(buffered_input[num_inputs-1], buffered_output);

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
        // main loop cat
        register uint64_t i = 0;
        register uint64_t c = 0;
        gt_status status;
        gt_template* template = gt_template_new();
        gt_buffered_output_file* buffered_output;
        gt_buffered_input_file** buffered_input = malloc(num_inputs * sizeof(gt_buffered_input_file*));

        for(i=0; i<num_inputs; i++){
            buffered_input[i] = gt_buffered_input_file_new(inputs[i]);
        }
        for(i=0; i<num_inputs; i++){
            buffered_output = gt_buffered_output_file_new(output);
            gt_buffered_input_file_attach_buffered_output(buffered_input[i], buffered_output);
            while((status = gt_input_generic_parser_get_template(buffered_input[i], template, parser_attributes)) == GT_STATUS_OK){
                gt_output_fasta_bofprint_template(buffered_output, template, attributes);
                c++;
            }
            gt_buffered_input_file_close(buffered_input[i]);
            gt_buffered_output_file_close(buffered_output);
        }
        // printf("Written %d\n", c);

        gt_template_delete(template);
        // for(i=0; i<num_inputs; i++){
        //     gt_buffered_input_file_close(buffered_input[i]);
        // }
        //gt_buffered_output_file_close(buffered_output);
        free(buffered_input);

    }
    printf("COMPLETE WRITING DONE >>>>>>>>>>>>>>>>\n");
    gt_output_fasta_attributes_delete(attributes);
    gt_input_generic_parser_attributes_delete(parser_attributes);
    printf("CLOSING IN AND OUT>>>>>>>>\n");
    register uint64_t i = 0;
    for(i=0; i<num_inputs; i++){
        gt_input_file_close(inputs[i]);
    }
    gt_output_file_close(output);
}
