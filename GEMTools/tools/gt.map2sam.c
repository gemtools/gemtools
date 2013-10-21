/*
 * PROJECT: GEM-Tools library
 * FILE: gt.map2sam.c
 * DATE: 02/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Converter from MAP to SAM
 */

#define GT_MAP2SAM "gt.map2sam"
#define GT_MAP2SAM_VERSION "1.1"

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

typedef struct {
  /* I/O */
  char* name_input_file;
  char* name_output_file;
  char* name_reference_file;
  char* name_gem_index_file;
  char* sam_header_file;
  bool mmap_input;
  bool paired_end;
  /* Headers */
  char *read_group_id;
  /* SAM format */
  bool compact_format;
  /* Optional Fields */
  bool optional_field_NH;
  bool optional_field_NM;
  bool optional_field_XT;
  bool optional_field_XS;
  bool optional_field_md;
  bool calc_phred;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool load_index;
  bool load_index_sequences;
  gt_map_score_attributes map_score_attr;
} gt_stats_args;

gt_stats_args parameters = {
  /* I/O */
  .name_input_file=NULL,
  .name_output_file=NULL,
  .name_reference_file=NULL,
  .name_gem_index_file=NULL,
  .sam_header_file=NULL,
  .mmap_input=false,
  .paired_end=false,
  /* Headers */
  .read_group_id=NULL,
  /* SAM format */
  .compact_format=false,
  /* Optional Fields */
  .optional_field_NH=false,
  .optional_field_NM=false,
  .optional_field_XT=false,
  .optional_field_XS=false,
  .optional_field_md=false,
  .calc_phred=false,
  /* Misc */
  .num_threads=1,
  .verbose=false,
  /* Control flags */
  .load_index=false,
  .load_index_sequences=false,
  .map_score_attr.quality_format=GT_QUALS_OFFSET_33,
  .map_score_attr.max_strata_searched=0,
};

gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(parameters.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  return sequence_archive;
}

#define PHRED_KONST -0.23025850929940456840 // -log(10)/10;

void gt_map2sam_calc_phred(gt_template *template,gt_map_score_attributes *ms_attr)
{
	gt_sam_attribute *ys=NULL,*yq=NULL;
	ys=gt_attributes_get_sam_attribute(template->attributes,"ms");
	if(!ys || ys->type_id!='B') return;
	yq=gt_attributes_get_sam_attribute(template->attributes,"mx");
	ms_attr->mapping_cutoff=(yq && yq->type_id=='i')?yq->i_value:0;
	uint64_t max_complete_strata[2]={0,0};
	char *p=gt_string_get_string(ys->s_value);
	if(p && *p && p[1]==',') {
		char *p1;
		max_complete_strata[0]=(uint64_t)strtoul(p+2,&p1,10);
		if(*p1==',') max_complete_strata[1]=(uint64_t)strtoul(p1+1,&p,10);
	}
	uint64_t rd;
	for(rd=0;rd<2;rd++) {
		gt_alignment *al=gt_template_get_block(template,rd);
		gt_attributes_add(al->attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA,&max_complete_strata[rd],uint64_t);
	}
	gt_map_calculate_template_mapq_score(template,ms_attr);
}

void gt_map2sam_read__write() {
  // Open file IN/OUT
  gt_input_file* const input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_map_open(stdin) : gt_input_file_map_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* const output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
  gt_sam_headers* const sam_headers = gt_sam_header_new(); // SAM headers

  if(parameters.sam_header_file) {
  	gt_input_file* const sam_headers_input_file = gt_input_file_sam_open(parameters.sam_header_file,false);
  	uint64_t characters_read = 0, lines_read = 0;
  	gt_status error_code=gt_input_file_sam_read_headers((char *)sam_headers_input_file->file_buffer,sam_headers_input_file->buffer_size,sam_headers,&characters_read,&lines_read);
  	if(error_code) gt_error(PARSE_SAM_HEADER_NOT_SAM,sam_headers_input_file->file_name);
  	gt_input_file_close(sam_headers_input_file);
  	uint64_t xx=0;
  	gt_string *pg_id=gt_string_new(64);
		gt_sprintf(pg_id,"GToolsLib#%"PRIu64,++xx);
		gt_string *prev_id=NULL;
  	if(sam_headers->program_id_hash) {
  		gt_sam_header_record* hr=*(gt_sam_header_record **)gt_vector_get_last_elm(sam_headers->program,gt_sam_header_record*);
			prev_id=gt_sam_header_record_get_tag(hr,"ID");
  		do {
  			if(!gt_shash_get_element(sam_headers->program_id_hash,gt_string_get_string(pg_id))) break;
  			gt_sprintf(pg_id,"GToolsLib#%"PRIu64,++xx);
  		} while(xx<10000);
  	}
  	if(xx<10000) {
  		gt_sam_header_record *hr=gt_sam_header_record_new();
  		gt_sam_header_record_add_tag(hr,"ID",pg_id);
  		gt_string *pn_st=gt_string_set_new(GT_MAP2SAM);
  		gt_sam_header_record_add_tag(hr,"PN",pn_st);
  		gt_string *vn_st=gt_string_set_new(GT_MAP2SAM_VERSION);
  		if(prev_id) gt_sam_header_record_add_tag(hr,"PP",prev_id);
  		gt_sam_header_record_add_tag(hr,"VN",vn_st);
  		gt_sam_header_add_program_record(sam_headers,hr);
  	}
  }

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) {
    sequence_archive = gt_filter_open_sequence_archive(parameters.load_index_sequences);
    gt_sam_header_load_sequence_archive(sam_headers,sequence_archive);
  }

  // Print SAM headers
  gt_output_sam_ofprint_headers_sh(output_file,sam_headers);

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    // I/O attributes
    gt_map_parser_attributes* const input_map_attributes = gt_input_map_parser_attributes_new(parameters.paired_end);
    gt_output_sam_attributes* const output_sam_attributes = gt_output_sam_attributes_new();
    // Set out attributes
    gt_output_sam_attributes_set_compact_format(output_sam_attributes,parameters.compact_format);
    gt_output_sam_attributes_set_qualities_offset(output_sam_attributes,parameters.map_score_attr.quality_format);
  	gt_output_sam_attributes_set_print_mismatches(output_sam_attributes,false);

    if (parameters.optional_field_NH) gt_sam_attributes_add_tag_NH(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_NM) gt_sam_attributes_add_tag_NM(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XT) gt_sam_attributes_add_tag_XT(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_md) gt_sam_attributes_add_tag_md(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XS) gt_sam_attributes_add_tag_XS(output_sam_attributes->sam_attributes);
  	gt_sam_attributes_add_tag_SA(output_sam_attributes->sam_attributes);
    if (parameters.calc_phred) {
  		gt_sam_attributes_add_tag_MQ(output_sam_attributes->sam_attributes);
  		gt_sam_attributes_add_tag_XP(output_sam_attributes->sam_attributes);
  		gt_sam_attributes_add_tag_UQ(output_sam_attributes->sam_attributes);
  		gt_sam_attributes_add_tag_PQ(output_sam_attributes->sam_attributes);
  		gt_sam_attributes_add_tag_TQ(output_sam_attributes->sam_attributes);
  		gt_sam_attributes_add_tag_TP(output_sam_attributes->sam_attributes);
    }
  	if(sam_headers->read_group_id_hash) {
  		gt_sam_header_record *hr=NULL;
  		if (parameters.read_group_id) {
      		size_t* ix=gt_shash_get_element(sam_headers->read_group_id_hash,parameters.read_group_id);
      		if(ix) {
      			hr=*(gt_sam_header_record **)gt_vector_get_elm(sam_headers->read_group,*ix,gt_sam_header_record*);
      		} else gt_error(SAM_OUTPUT_UNKNOWN_RG_ID,parameters.read_group_id);
  		} else {
  			hr=*(gt_sam_header_record **)gt_vector_get_last_elm(sam_headers->read_group,gt_sam_header_record*);
  		}
  		if(hr) {
  			gt_string *id_tag=gt_sam_header_record_get_tag(hr,"ID");
  			if(!id_tag) gt_fatal_error(PARSE_SAM_HEADER_MISSING_TAG,"RG","ID"); // Should have been detected before, but we check again anyway
  			gt_sam_attributes_add_tag_RG(output_sam_attributes->sam_attributes,id_tag);
  			gt_string *lib_tag=gt_sam_header_record_get_tag(hr,"LB");
  			if(lib_tag) gt_sam_attributes_add_tag_LB(output_sam_attributes->sam_attributes,lib_tag);
  		}
    } else if(parameters.read_group_id) gt_error(SAM_OUTPUT_NO_HEADER_FOR_RG);
    gt_template* template = gt_template_new();
    while ((error_code=gt_input_map_parser_get_template(buffered_input,template,input_map_attributes))) {
      if (error_code!=GT_IMP_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
        continue;
      }
      if(parameters.calc_phred || parameters.optional_field_XT) gt_map2sam_calc_phred(template,&parameters.map_score_attr);
      // Print SAM template
      gt_output_sam_bofprint_template(buffered_output,template,output_sam_attributes);
    }

    // Clean
    gt_template_delete(template);
    gt_input_map_parser_attributes_delete(input_map_attributes);
    gt_output_sam_attributes_delete(output_sam_attributes);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_sam_header_delete(sam_headers);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

void usage(const gt_option* const options,char* groups[],const bool print_inactive) {
  fprintf(stderr, "USE: ./gt.map2sam [ARGS]...\n");
  gt_options_fprint_menu(stderr,options,groups,false,print_inactive);
}

void parse_arguments(int argc,char** argv) {
  struct option* gt_map2sam_getopt = gt_options_adaptor_getopt(gt_map2sam_options);
  gt_string* const gt_map2sam_short_getopt = gt_options_adaptor_getopt_short(gt_map2sam_options);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_map2sam_short_getopt),gt_map2sam_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      parameters.load_index = true;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      parameters.load_index = true;
      break;
    case 's':
    	parameters.sam_header_file = optarg;
    	break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 'Q':
    	parameters.calc_phred = true;
    	break;
    case 200:
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    /* Headers */
    case 300: // Read-group ID
    	parameters.read_group_id = optarg;
    	break;
      // TODO
    /* Alignments */
    case 'q':
      if (gt_streq(optarg,"offset-64")) {
        parameters.map_score_attr.quality_format=GT_QUALS_OFFSET_64;
      } else if (gt_streq(optarg,"offset-33")) {
        parameters.map_score_attr.quality_format=GT_QUALS_OFFSET_33;
      } else {
        gt_fatal_error_msg("Quality format not recognized: '%s'",optarg);
      }
      break;
    /* Optional Fields */
    case 500: // NH
      parameters.optional_field_NH = true;
      break;
    case 501: // NM
      parameters.optional_field_NM = true;
      break;
    case 502: // XT
      parameters.optional_field_XT = true;
      break;
    case 503: // XS
      parameters.optional_field_XS = true;
      break;
    case 504: // md
      parameters.optional_field_md = true;
      break;
    /* Format */
    case 'c':
      parameters.compact_format = true;
      break;
      /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage(gt_map2sam_options,gt_map2sam_groups,false);
      exit(1);
      break;
    case 'H':
      usage(gt_map2sam_options,gt_map2sam_groups,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_map2sam_options,gt_map2sam_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters check
   */
  if (parameters.load_index && parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
  if(!parameters.load_index && parameters.optional_field_XS){
    gt_fatal_error_msg("Reference file required to compute XS field in SAM");
  }
  // Free
  gt_string_delete(gt_map2sam_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // map2sam !!
  gt_map2sam_read__write();

  return 0;
}

