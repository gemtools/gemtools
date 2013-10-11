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
  gt_qualities_offset_t quality_format;
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
  int uniq_mapq_thresh;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool load_index;
  bool load_index_sequences;
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
  .quality_format=GT_QUALS_OFFSET_33,
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
  .uniq_mapq_thresh=30,
  /* Misc */
  .num_threads=1,
  .verbose=false,
  /* Control flags */
  .load_index=false,
  .load_index_sequences=false
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

static int cmp_int(const void *s1,const void *s2)
{
	return (*(int *)s1)-(*(int *)s2);
}

void gt_map2sam_calc_phred(gt_template *template)
{
	typedef struct {
		char *key;
		gt_map *map;
		double prob;
		double prob1;
		uint64_t score;
		uint8_t phred;
		uint8_t phred1;
		UT_hash_handle hh;
	} map_hash;

	map_hash *mhash[3]={0,0,0}, *mp_hash, *tmp;
	int rd;
	size_t buf_len=1024;
	gt_sam_attribute *ys=NULL,*yq=NULL;
	ys=gt_attributes_get_sam_attribute(template->attributes,"ms");
	if(!ys || ys->type_id!='B') return;
	yq=gt_attributes_get_sam_attribute(template->attributes,"mx");
	uint32_t map_cutoff=(yq && yq->type_id=='i')?yq->i_value:0;
	uint32_t max_complete_strata[2]={0,0};
	char *p=gt_string_get_string(ys->s_value);
	if(p && *p && p[1]==',') {
		char *p1;
		max_complete_strata[0]=(uint32_t)strtoul(p+2,&p1,10);
		if(*p1==',') max_complete_strata[1]=(uint32_t)strtoul(p1+1,&p,10);
	}
	char *buf=malloc(buf_len);
	gt_cond_fatal_error(!buf,MEM_HANDLER);
	bool unpaired=false;
	uint64_t min_score[4]={0xffff,0xffff,0xffff,0xffff};
	{
		GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,maps,maps_attr) {
			if(!maps_attr) gt_fatal_error(TEMPLATE_NOT_SCORED);
			uint64_t score=maps_attr->gt_score;
			if(score==GT_MAP_NO_GT_SCORE) gt_fatal_error(TEMPLATE_NOT_SCORED);
			uint64_t seq_like[2],interval_like;
			seq_like[0]=score&0xffff;
			seq_like[1]=(score>>16)&0xffff;
			interval_like=(score>>32)&0xff;
			// Build up list of single end alignments (need this so we can scale the MAPQ score)
			// Use hash to avoid counting a single end alignment twice if it occurs in two paired alignments
			for(rd=0;rd<2;rd++) if(maps[rd]) {
				size_t ssize=gt_string_get_length(maps[rd]->seq_name);
				size_t key_size=ssize+sizeof(maps[rd]->position);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				memcpy(buf,gt_string_get_string(maps[rd]->seq_name),ssize);
				memcpy(buf+ssize,&maps[rd]->position,sizeof(maps[rd]->position));
				HASH_FIND(hh,mhash[rd],buf,key_size,mp_hash);
				if(!mp_hash) {
					mp_hash=malloc(sizeof(map_hash));
					gt_cond_fatal_error(!mp_hash,MEM_HANDLER);
					mp_hash->key=malloc(key_size);
					gt_cond_fatal_error(!mp_hash->key,MEM_HANDLER);
					memcpy(mp_hash->key,buf,key_size);
					mp_hash->score=seq_like[rd];
					mp_hash->map=maps[rd];
					HASH_ADD_KEYPTR(hh,mhash[rd],mp_hash->key,key_size,mp_hash);
					if(seq_like[rd]<min_score[rd]) min_score[rd]=seq_like[rd];
				}
			}
			if(maps[0] && maps[1]) { // True paired alignments.  Shouldn't need to check for duplicates, but we will anyway
				size_t ssize=gt_string_get_length(maps[0]->seq_name);
				size_t key_size=ssize+2*sizeof(maps[0]->position);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				memcpy(buf,gt_string_get_string(maps[0]->seq_name),ssize);
				memcpy(buf+ssize,&maps[0]->position,sizeof(maps[0]->position));
				memcpy(buf+ssize+sizeof(maps[0]->position),&maps[1]->position,sizeof(maps[0]->position));
				HASH_FIND(hh,mhash[2],buf,key_size,mp_hash);
				if(!mp_hash) {
					mp_hash=malloc(sizeof(map_hash));
					gt_cond_fatal_error(!mp_hash,MEM_HANDLER);
					mp_hash->key=malloc(key_size);
					gt_cond_fatal_error(!mp_hash->key,MEM_HANDLER);
					memcpy(mp_hash->key,buf,key_size);
					uint64_t sc=seq_like[0]+seq_like[1]+interval_like;
					mp_hash->score=sc;
					HASH_ADD_KEYPTR(hh,mhash[2],mp_hash->key,key_size,mp_hash);
					if(sc<min_score[2]) min_score[2]=sc;
				}
			} else unpaired=true;
		}
	}
	// Get score of best possible aligning read that we didn't look for with the mapping parameters
	uint64_t fake_sc[3]={0,0,0};
	for(rd=0;rd<2;rd++) if(max_complete_strata[rd]) {
		gt_alignment *al=gt_template_get_block(template,rd);
		if(!al) continue;
	  gt_string* const quals=al->qualities;
	  uint64_t i;
	  int *qvs=malloc(sizeof(int)*quals->length);
	  size_t k=0;
	  int quals_offset=parameters.quality_format==GT_QUALS_OFFSET_33?33:64;
	  for(i=0;i<quals->length;i++) {
	  	int q=gt_string_get_string(quals)[i]-quals_offset;
	  	if(q>=map_cutoff) qvs[k++]=q;
	  }
	  if(k>=max_complete_strata[rd]) {
	  	qsort(qvs,k,sizeof(int),cmp_int);
	  	for(i=0;i<max_complete_strata[rd];i++) fake_sc[rd]+=qvs[i];
	  }
	  free(qvs);
	}
	// For paired reads, the best aligning read that we didn't find would be the best read found from one end combined
	// with the fake score from the other with the most likely interval size (i.e. 0)
	uint64_t f1=fake_sc[0]+min_score[1];
	uint64_t f2=fake_sc[1]+min_score[0];
	fake_sc[2]=f1<f2?f1:f2;
	// Insert fake scores into hash structures
	for(rd=0;rd<3;rd++) if(fake_sc[rd]) {
	  size_t key_size=6;
	  mp_hash=malloc(sizeof(map_hash));
	  gt_cond_fatal_error(!mp_hash,MEM_HANDLER);
	  mp_hash->key=malloc(key_size);
	  gt_cond_fatal_error(!mp_hash->key,MEM_HANDLER);
	  memcpy(mp_hash->key,"_FAKE_",key_size);
	  mp_hash->score=fake_sc[rd];
	  mp_hash->map=0;
	  HASH_ADD_KEYPTR(hh,mhash[rd],mp_hash->key,key_size,mp_hash);
		if(fake_sc[rd]<min_score[rd]) min_score[rd]=fake_sc[rd];
	}
	// Now we can calculate the single and paired end MAPQ values
	double z;
	for(rd=0;rd<3;rd++) if(mhash[rd]) {
		z=0.0;
		for(mp_hash=mhash[rd];mp_hash;mp_hash=mp_hash->hh.next) {
			mp_hash->prob=exp(PHRED_KONST*((double)(mp_hash->score)-(double)min_score[rd]));
			z+=mp_hash->prob;
		}
		for(mp_hash=mhash[rd];mp_hash;mp_hash=mp_hash->hh.next) {
			mp_hash->prob/=z;
			if(1.0-mp_hash->prob<1.0e-255) mp_hash->phred=254;
			else {
				int tp=(int)(0.5+log(1.0-mp_hash->prob)/PHRED_KONST);
				if(tp>254) tp=254;
				mp_hash->phred=tp;
			}
		}
	}
	// To see how well supported our paired reads are compared to all unpaired combinations,
	// if unpaired alignments exist we look at all combinations of single end reads
	if(unpaired) {
		map_hash *mp1,*mp2;
		for(mp1=mhash[0];mp1;mp1=mp1->hh.next) {
			for(mp2=mhash[1];mp2;mp2=mp2->hh.next) {
				uint64_t sc=mp1->score+mp2->score;
				if(sc<min_score[3]) min_score[3]=sc;
			}
		}
		z=0.0;
		for(mp1=mhash[0];mp1;mp1=mp1->hh.next) {
			for(mp2=mhash[1];mp2;mp2=mp2->hh.next) {
				uint64_t sc=mp1->score+mp2->score;
				double prb=exp(PHRED_KONST*((double)sc-(double)min_score[3]));
				z+=prb;
				if(mp1->map && mp2->map) {
					size_t ssize=gt_string_get_length(mp1->map->seq_name);
					size_t key_size=ssize+2*sizeof(mp1->map->position);
					memcpy(buf,gt_string_get_string(mp1->map->seq_name),ssize);
					memcpy(buf+ssize,&mp1->map->position,sizeof(mp1->map->position));
					memcpy(buf+ssize+sizeof(mp1->map->position),&mp2->map->position,sizeof(mp1->map->position));
					HASH_FIND(hh,mhash[2],buf,key_size,mp_hash);
					if(mp_hash) mp_hash->prob1=prb;
				}
			}
		}
		for(mp_hash=mhash[2];mp_hash;mp_hash=mp_hash->hh.next) {
			mp_hash->prob1/=z;
			if(1.0-mp_hash->prob1<1.0e-255) mp_hash->phred1=254;
			else {
				int tp=(int)(0.5+log(1.0-mp_hash->prob1)/PHRED_KONST);
				if(tp>254) tp=254;
				mp_hash->phred1=tp>mp_hash->phred?mp_hash->phred:tp;
			}
		}
	}
	// And now we have to enter the MAPQ values in the map structures
	{
		GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,maps,maps_attr) {
			for(rd=0;rd<2;rd++) if(maps[rd]) {
				size_t ssize=gt_string_get_length(maps[rd]->seq_name);
				size_t key_size=ssize+sizeof(maps[rd]->position);
				memcpy(buf,gt_string_get_string(maps[rd]->seq_name),ssize);
				memcpy(buf+ssize,&maps[rd]->position,sizeof(maps[rd]->position));
				HASH_FIND(hh,mhash[rd],buf,key_size,mp_hash);
				assert(mp_hash);
				maps[rd]->phred_score=fake_sc[rd]?mp_hash->phred:255;
			}
			if(maps[0] && maps[1]) { // True paired alignments.  Shouldn't need to check for duplicates, but we will anyway
				// seq_name should be the same for the two ends in a paired alignment, but we're not taking any chances
				size_t ssize=gt_string_get_length(maps[0]->seq_name);
				size_t key_size=ssize+2*sizeof(maps[0]->position);
				memcpy(buf,gt_string_get_string(maps[0]->seq_name),ssize);
				memcpy(buf+ssize,&maps[0]->position,sizeof(maps[0]->position));
				memcpy(buf+ssize+sizeof(maps[0]->position),&maps[1]->position,sizeof(maps[0]->position));
				HASH_FIND(hh,mhash[2],buf,key_size,mp_hash);
				assert(mp_hash);
				maps_attr->phred_score=mp_hash->phred;
				maps_attr->pair_score=unpaired?mp_hash->phred1:mp_hash->phred;
			}
		}
	}
	free(buf);
	for(rd=0;rd<3;rd++) if(mhash[rd]) {
		HASH_ITER(hh,mhash[rd],mp_hash,tmp) {
			HASH_DEL(mhash[rd],mp_hash);
			free(mp_hash->key);
			free(mp_hash);
		}
	}
}

void gt_map2sam_read__write() {
  // Open file IN/OUT
  gt_input_file* const input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* const output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
  gt_sam_headers* const sam_headers = gt_sam_header_new(); // SAM headers

  if(parameters.sam_header_file) {
  	gt_input_file* const sam_headers_input_file = gt_input_file_open(parameters.sam_header_file,false);
  	if(sam_headers_input_file->file_format!=SAM) gt_error(PARSE_SAM_HEADER_NOT_SAM,sam_headers_input_file->file_name);
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
    gt_output_sam_attributes_set_qualities_offset(output_sam_attributes,parameters.quality_format);
    if (parameters.optional_field_NH) gt_sam_attributes_add_tag_NH(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_NM) gt_sam_attributes_add_tag_NM(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XT) gt_sam_attributes_add_tag_XT(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_md) gt_sam_attributes_add_tag_md(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XS) gt_sam_attributes_add_tag_XS(output_sam_attributes->sam_attributes);
    if (parameters.calc_phred) {
    	gt_sam_attributes_add_tag_MQ(output_sam_attributes->sam_attributes);
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
      if(parameters.calc_phred || parameters.optional_field_XT) gt_map2sam_calc_phred(template);
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
        parameters.quality_format=GT_QUALS_OFFSET_64;
      } else if (gt_streq(optarg,"offset-33")) {
        parameters.quality_format=GT_QUALS_OFFSET_33;
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

