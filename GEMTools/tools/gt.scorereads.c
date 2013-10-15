/*
 * gt.scorereads.c
 *
 *  Created on: 8 Jul 2013
 *      Author: heath
 */


#include <getopt.h>
#include <ctype.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "gem_tools.h"

#define MAP_THRESHOLD 3
#define AP_BUF_SIZE 16384
#define QUAL_FASTQ 33
#define QUAL_SOLEXA 64
#define SOLEXA_BAD_QUAL 2
#define DEFAULT_QUAL (QUAL_FASTQ)
#define MISSING_QUAL 40 // Value to use in alignment score if qualities not available
#define INDEL_QUAL 40 // Value to use in alignment score for indels
#define MAX_GT_SCORE 0xFFFF
#define MAX_QUAL 42
#define ALIGN_NORM 0
#define ALIGN_BS_POS 1
#define ALIGN_BS_NEG 2
#define ALIGN_TYPES 3
#define AL_FORWARD 4
#define AL_REVERSE 8
#define AL_DIRECTIONS ((AL_FORWARD)|(AL_REVERSE))
#define AL_USED 128
#define GT_SCOREREADS_NUM_LINES GT_IMP_NUM_LINES

typedef struct {
  char *input_files[2];
  char *output_file;
  char *dist_file;
  bool mmap_input;
  bool verbose;
  gt_output_file_compression compress;
  gt_generic_parser_attributes *parser_attr;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output[2];
  int64_t min_insert;
  int64_t max_insert;
  uint8_t *ins_phred;
  int num_threads;
  int mapping_cutoff;
  int max_strata_searched;
  int indel_penalty;
  int split_penalty;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina)
} sr_param;

typedef struct {
  int64_t x;
  uint64_t ct;
  UT_hash_handle hh;
} sr_dist_element;

sr_param param = {
		.input_files={NULL,NULL},
		.output_file=NULL,
		.dist_file=NULL,
		.mmap_input=false,
		.compress=NONE,
		.parser_attr=NULL,
		.verbose=false,
		.num_threads=1,
		.indel_penalty=INDEL_QUAL,
		.split_penalty=INDEL_QUAL,
		.qual_offset=DEFAULT_QUAL,
		.mapping_cutoff=0,
		.max_strata_searched=0,
		.min_insert=0,
		.max_insert=0,
		.ins_phred=NULL
};

void usage(const gt_option* const options,char* groups[],const bool print_inactive) {
  fprintf(stderr, "USE: ./gt_scorereads [ARGS]...\n");
  gt_options_fprint_menu(stderr,options,groups,false,print_inactive);
}

static sr_dist_element *sr_find_insert_counter(sr_dist_element **de,int64_t x)
{
	sr_dist_element *new_de;
	HASH_FIND(hh,*de,&x,sizeof(int64_t),new_de);
	if(!new_de) {
		new_de=gt_alloc(sr_dist_element);
		new_de->ct=0;
		new_de->x=x;
		HASH_ADD(hh,*de,x,sizeof(int64_t),new_de);
	}
	return new_de;
}

static sr_dist_element *sr_increase_insert_count(sr_dist_element **de,int64_t x,uint64_t ct)
{
	sr_dist_element *new_de=sr_find_insert_counter(de,x);
	new_de->ct+=ct;
	return new_de;
}

static uint64_t calculate_dist_score(gt_alignment *al, gt_map *map, sr_param *param)
{
	int qual_offset=param->qual_offset;
	int indel_penalty=param->indel_penalty;
	int split_penalty=param->split_penalty;
	int map_cutoff=param->mapping_cutoff;
	register gt_string* const quals = al->qualities;
	register gt_string* const read = al->read;
	register const bool has_qualities = gt_alignment_has_qualities(al);
	uint64_t score=0;
	GT_MAP_ITERATE(map,map_block) {
		GT_MISMS_ITERATE(map_block,misms) {
			int quality_misms;
			if (has_qualities) {
				quality_misms = gt_string_get_string(quals)[misms->position]-qual_offset;
				if(quality_misms>MAX_QUAL) quality_misms=MAX_QUAL;
				else if(quality_misms<map_cutoff || gt_string_get_string(read)[misms->position]=='N') quality_misms=0;
			} else quality_misms=MISSING_QUAL;
			switch (misms->misms_type) {
			case MISMS:
				score+=quality_misms;
				break;
			case INS:
			case DEL:
				score+=indel_penalty;
				break;
			}
		}
		if(gt_map_get_next_block(map_block)) score+=split_penalty;
	}
	if(score>MAX_GT_SCORE) score=MAX_GT_SCORE;
	return score;
}

static gt_status gt_scorereads_get_insert_size(gt_template *template,gt_alignment *al[],int64_t *size)
{
	gt_status err=GT_STATUS_FAIL;
	int j;
	uint64_t nmap[2]={0,0};
	for(j=0;j<2;j++) {
		uint64_t i,k;
		k=gt_alignment_get_num_counters(al[j]);
		for(i=0;i<k;i++) {
			uint64_t x=gt_alignment_get_counter(al[j],i);
			if(x) {
				nmap[j]+=x;
			} else if(nmap[j]) break;
		}
	}
	if(nmap[0]==1 && nmap[1]==1) {
		gt_map *mmap[2];
		mmap[0]=gt_alignment_get_map(al[0],0);
		mmap[1]=gt_alignment_get_map(al[1],0);
		int64_t x=gt_template_get_insert_size(mmap,&err,0,0);
		if(err==GT_TEMPLATE_INSERT_SIZE_OK) {
			err=GT_STATUS_OK;
			*size=x;
		}
	}
	return err;
}

static void gt_scorereads_free_insert_hash(sr_dist_element *ins) {
	sr_dist_element *next_elem;
	while(ins) {
		next_elem=ins->hh.next;
		free(ins);
		ins=next_elem;
	}
}

static void gt_scorereads_pair_template(gt_template *template,gt_alignment *al[2],sr_param *param)
{
	gt_mmap_attributes attr;
	gt_map *mmap[2];
	uint64_t nmap[2],rd;
	for(rd=0;rd<2;rd++) {
		gt_alignment_recalculate_counters(al[rd]);
		nmap[rd]=gt_alignment_get_num_maps(al[rd]);
		if(nmap[rd]) {
			GT_ALIGNMENT_ITERATE(al[rd],map) {
				map->gt_score=calculate_dist_score(al[rd],map,param);
			}
		}
	}
	if(nmap[0]+nmap[1]) {
		uint64_t i=0;
		char *map_flag[2];
		map_flag[0]=gt_calloc((size_t)(nmap[0]+nmap[1]),char,true);
		map_flag[1]=map_flag[0]+nmap[0];
		GT_ALIGNMENT_ITERATE(al[0],map1) {
			mmap[0]=map1;
			uint64_t j=0;
			GT_ALIGNMENT_ITERATE(al[1],map2) {
				mmap[1]=map2;
				gt_status gt_err;
				int64_t x=gt_template_get_insert_size(mmap,&gt_err,0,0);
				if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK && x>=param->min_insert && x<=param->max_insert) {
					attr.distance=gt_map_get_global_distance(map1)+gt_map_get_global_distance(map2);
					attr.gt_score=map1->gt_score|(map2->gt_score<<16);
					if(param->ins_phred) attr.gt_score|=((uint64_t)param->ins_phred[x-param->min_insert]<<32);
					attr.phred_score=255;
					gt_template_inc_counter(template,attr.distance);
					gt_template_add_mmap_ends(template,map1,map2,&attr);
					map_flag[0][i]=map_flag[1][j]=1;
				}
				j++;
			}
			i++;
		}
		for(rd=0;rd<2;rd++) {
			for(i=0;i<nmap[rd];i++) {
				if(!map_flag[rd][i]) {
					gt_map *map=gt_alignment_get_map(al[rd],i);
					attr.distance=gt_map_get_global_distance(map);
					attr.gt_score=map->gt_score;
					attr.phred_score=255;
					gt_template_inc_counter(template,attr.distance);
					gt_template_add_mmap_ends(template,map,0,&attr);
				}
			}
		}
		free(map_flag[0]);
	}
	gt_attributes_remove(template->attributes,GT_ATTR_ID_TAG_PAIR);
}

void gt_add_extra_tags(gt_attributes *attributes,uint32_t mcs[],bool paired,int map_cutoff) {
	char *tag;
	if(paired) asprintf(&tag,"ms:B:I,%" PRIu32 ",%" PRIu32 " mx:i:%d",mcs[0],mcs[1],map_cutoff);
	else asprintf(&tag,"ms:B:I,%" PRIu32 " mx:i:%d",mcs[0],map_cutoff);
	gt_cond_fatal_error(!tag,MEM_HANDLER);
	gt_string *extra_string=gt_string_set_new(tag);
	free(tag);
	gt_string *old=gt_attributes_get(attributes,GT_ATTR_ID_TAG_EXTRA);
	if(!old) {
		gt_attributes_add_string(attributes,GT_ATTR_ID_TAG_EXTRA,extra_string);
	} else {
		gt_string_append_char(old,SPACE);
		gt_string_append_gt_string(old,extra_string);
		gt_string_delete(extra_string);
	}
}

void gt_template_add_mcs_tags(gt_template *template,gt_alignment *al[2],int map_cutoff)
{
	uint32_t rd,mcs[2];
	for(rd=0;rd<2;rd++) {
		mcs[rd]=al[rd]?gt_alignment_get_mcs(al[rd]):0;
		if(param.max_strata_searched) {
			uint32_t limit=param.max_strata_searched+1;
			if(mcs[rd]>limit) mcs[rd]=limit;
		}
	}
	gt_add_extra_tags(template->attributes,mcs,true,map_cutoff);
}

void gt_alignment_add_mcs_tags(gt_alignment *al,int map_cutoff)
{
	uint32_t mcs[2];
	mcs[0]=gt_alignment_get_mcs(al);
	if(param.max_strata_searched) {
		uint32_t limit=param.max_strata_searched+1;
		if(mcs[0]>limit) mcs[0]=limit;
	}
	gt_add_extra_tags(al->attributes,mcs,false,map_cutoff);
}

static int sort_dist_element(void *a,void *b)
{
	return ((sr_dist_element *)a)->x-((sr_dist_element *)b)->x;
}

int estimate_insert_histogram(sr_dist_element *ins,sr_param *param)
{
	uint64_t tot=0;
	sr_dist_element *elem;

	// Get interquartile range
	HASH_SORT(ins,sort_dist_element);
	for(elem=ins;elem;elem=elem->hh.next) tot+=elem->ct;
	uint64_t k=0,i=(tot+2)>>2;
	int64_t q1,q3;
	for(elem=ins;elem;elem=elem->hh.next) {
		k+=elem->ct;
		if(k>=i) break;
	}
	if(!elem) return GT_STATUS_FAIL;
	q1=elem->x;
	i=(tot*3+2)>>2;
	elem=elem->hh.next;
	for(;elem;elem=elem->hh.next) {
		k+=elem->ct;
		if(k>=i) break;
	}
	if(!elem) return GT_STATUS_FAIL;
	q3=elem->x;
	// Robust estimate of standard deviation from interquartile range
	double sd=(q3>q1)?(double)(q3-q1)/1.349:1.0;
	// Select bandwidth for kde (Silverman's rule of thumb)
	double h=exp(log(sd)+.2*log(4.0/(3.0*(double)tot)));
	size_t l=param->max_insert-param->min_insert+1;
	double *p=gt_malloc(sizeof(double)*2*l);
	double *cp=p+l;
	int *ph=gt_malloc(sizeof(int)*l);
	for(i=0;i<l;i++) p[i]=0.0;
	for(elem=ins;elem;elem=elem->hh.next) {
		double x=(double)elem->x;
		double n=(double)elem->ct;
		for(i=0;i<l;i++) {
			double d=(x-(double)(param->min_insert+i))/h;
			p[i]+=n*exp(-.5*d*d);
		}
	}
	double zmax=0.0,ztot=0.0;
	for(i=0;i<l;i++) {
		if(p[i]>zmax) {
			zmax=p[i];
			k=i;
		}
		ztot+=p[i];
	}
	cp[0]=p[0]/ztot;
	for(i=1;i<l;i++) cp[i]=cp[i-1]+p[i]/ztot;
	for(i=0;i<l;i++) {
		double z=p[i]/zmax;
		int phred=255;
		if(z>1.0e-26) {
			phred=(int)(log(z)*-10.0/log(10.0)+.5);
			if(phred>255) phred=255;
		}
		ph[i]=phred;
	}
	uint64_t i1=0,i2=l;
	for(i=k;i>0;i--) if(ph[i-1]==255) {
		i1=i;
		break;
	}
	for(i=k+1;i<l;i++) if(ph[i]==255) {
		i2=i-1;
		break;
	}
	while(i1 && i2<l && (cp[i2]-cp[i1])<.95) {
		for(i=i1;i>0;i--) if(ph[i-1]>255) break;
		if(i) {
			for(;i>0;i--) if(ph[i-1]==255) {
				i1=i;
				break;
			}
		}
		for(i=i2+1;i<l;i++) if(ph[i]>255) break;
		if(i<l) {
			for(;i<l;i++) if(ph[i]==255) {
				i2=i-1;
				break;
			}
		}
	}
	param->max_insert=param->min_insert+i2;
	param->min_insert+=i1;
	l=i2-i1+1;
	param->ins_phred=gt_malloc(l);
	for(i=i1,k=0;i<=i2;i++) param->ins_phred[k++]=ph[i];
//	for(i=0;i<l;i++) printf("%"PRId64"\t%d\n",param->min_insert+i,param->ins_phred[i]);
	free(p);
	free(ph);
	return GT_STATUS_OK;
}

void read_dist_file(sr_param *param)
{
	gt_input_file* file=gt_input_file_open(param->dist_file,false);
	gt_buffered_input_file* bfile=gt_buffered_input_file_new(file);
	gt_status nl=0;
	typedef enum {V1,V2,UNKNOWN} dist_file_type;
	dist_file_type ftype=UNKNOWN;
	sr_dist_element *ins_dist=NULL;
	do {
		nl=gt_buffered_input_file_get_block(bfile,0);
		int i;
		char *cur=(char *)gt_vector_get_mem(bfile->block_buffer,char);
		for(i=0;i<(int)nl;i++) {
			char *p=cur;
			cur=strchr(p,'\n');
			assert(cur);
			while(*p!='\n' && isspace(*p)) p++;
			if(*p && *p!='n') {
				char *p1=p+1;
				while(!isspace(*p1)) p1++;
				if(*p1!='\n') {
					*(p1++)=0;
					while(*p1!='\n' && isspace(*p1)) p1++;
					if(*p1!='\n') {
						char *p2=p1+1;
						while(!isspace(*p2)) p2++;
						*p2=0;
						if(ftype==UNKNOWN) {
							if(!strcmp(p,"Size")) {
								if(!strcmp(p1,"DS_count")) ftype=V1;
								else if(!strcmp(p1,"Paired")) ftype=V2;
							}
						} else {
							if(p[1]=='=' && (p[0]=='<' || p[0]=='>')) p+=2;
							int ef=0;
							int64_t x=strtoul(p,&p2,10);
							if(*p2) ef=1;
							else {
								int64_t y=strtoul(p1,&p2,10);
								if(*p2 || y<0) ef=2;
								else {
									(void)sr_increase_insert_count(&ins_dist,x,y);
								}
							}
							if(ef) {
								fprintf(stderr,"Invalid entry in insert distribution file: %s\t%s\n",p,p1);
								break;
							}
						}
					}
				}
			}
			cur++;
		}
		if(i<nl) break;
	} while(nl);
	if(ftype==UNKNOWN) fprintf(stderr,"Insert distribution file format not recognized\n");
	else estimate_insert_histogram(ins_dist,param);
	gt_scorereads_free_insert_hash(ins_dist);
	gt_buffered_input_file_close(bfile);
	gt_input_file_close(file);

}


gt_status parse_arguments(int argc,char** argv) {
	gt_status err=GT_STATUS_OK;
	param.parser_attr=gt_input_generic_parser_attributes_new(false);
  struct option* gt_scorereads_getopt = gt_options_adaptor_getopt(gt_scorereads_options);
  gt_string* const gt_scorereads_short_getopt = gt_options_adaptor_getopt_short(gt_scorereads_options);
  int option, option_index;
  char *p;
  int insert_set[2]={0,0};
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_scorereads_short_getopt),gt_scorereads_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i':
    	param.dist_file = optarg;
    	break;
    case 'o':
    	param.output_file = optarg;
    	break;
    case 300:
      param.input_files[0] = optarg;
      break;
    case 301:
      param.input_files[1] = optarg;
      break;
  	case 'p':
  		gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
  		break;
  	case 'z':
#ifdef HAVE_ZLIB
  		param.compress=GZIP;
#endif
  		break;
  	case 'j':
#ifdef HAVE_BZLIB
  		param.compress=BZIP2;
#endif
  		break;
  	case 'Z':
  		param.compress=NONE;
  		break;
 		/* Score function */
  	case 'q':
      if (gt_streq(optarg,"offset-64")) {
        param.qual_offset=64;
      } else if (gt_streq(optarg,"offset-33")) {
        param.qual_offset=33;
      } else {
        gt_fatal_error_msg("Quality format not recognized: '%s'",optarg);
      }
      break;
  	case 401:
  		param.min_insert=(int)strtol(optarg,&p,10);
  		if(*p || param.min_insert<0) {
  			fprintf(stderr,"Illegal minimum insert size: '%s'\n",optarg);
  			err=-7;
  		} else insert_set[0]=1;
  		break;
  	case 402:
  		param.max_insert=(int)strtol(optarg,&p,10);
  		if(*p || param.max_insert<0) {
  			fprintf(stderr,"Illegal maximum insert size: '%s'\n",optarg);
  			err=-7;
  		} else insert_set[1]=1;
  		break;
  	case 403:
  		param.indel_penalty=(int)strtol(optarg,&p,10);
  		if(*p || param.indel_penalty<0) {
  			fprintf(stderr,"Illegal indel penalty: '%s'\n",optarg);
  			err=-7;
  		}
  		break;
  	case 's':
  		param.split_penalty=(int)strtol(optarg,&p,10);
  		if(*p || param.split_penalty<0) {
  			fprintf(stderr,"Illegal split penalty score: '%s'\n",optarg);
  			err=-7;
  		}
  		break;
  	case 'M':
  		param.mapping_cutoff=(int)strtol(optarg,&p,10);
  		if(*p || param.mapping_cutoff<0) {
  			fprintf(stderr,"Illegal mapping cutoff: '%s'\n",optarg);
  			err=-7;
  		}
  		break;
  	case 'm':
  		param.max_strata_searched=(int)strtol(optarg,&p,10);
  		if(*p || param.max_strata_searched<0) {
  			fprintf(stderr,"Illegal mismatch limit: '%s'\n",optarg);
  			err=-7;
  		}
  		break;
    /* Misc */
    case 'v':
      param.verbose = true;
      break;
    case 't':
#ifdef HAVE_OPENMP
      param.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage(gt_scorereads_options,gt_scorereads_groups,false);
      exit(1);
      break;
    case 'H':
      usage(gt_scorereads_options,gt_scorereads_groups,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_map2sam_options,gt_map2sam_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      usage(gt_scorereads_options,gt_scorereads_groups,false);
      gt_fatal_error_msg("Option '%c' %d not recognized",option,option);
    }
  }
  /*
   * Parameters check
   */
	if(param.input_files[1]) gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
	if(err==GT_STATUS_OK && insert_set[0] && insert_set[1] && param.min_insert>param.max_insert) {
		fputs("Minimum insert size > maximum insert size\n",stderr);
    usage(gt_scorereads_options,gt_scorereads_groups,false);
		err=-15;
	}
	if(err==GT_STATUS_OK) {
		if(param.output_file && param.compress!=NONE) {
			size_t l=strlen(param.output_file);
			switch(param.compress) {
			case GZIP:
				if(l<3 || strcmp(param.output_file+l-3,".gz")) {
					char *s;
					asprintf(&s,"%s.gz",param.output_file);
					param.output_file=s;
				}
				break;
			case BZIP2:
				if(l<4 || strcmp(param.output_file+l-4,".bz2")) {
					char *s;
					asprintf(&s,"%s.bz2",param.output_file);
					param.output_file=s;
				}
				break;
			default:
				break;
			}
		}
		if(!insert_set[1]) {
				if(param.min_insert<=1000) param.max_insert=1000;
				else param.max_insert=param.min_insert+1000;
		}
		if(gt_input_generic_parser_attributes_is_paired(param.parser_attr) && param.dist_file) read_dist_file(&param);
	}
  // Free
  gt_string_delete(gt_scorereads_short_getopt);
  return err;
}

gt_status gt_input_map_parser_synch_get_template(gt_buffered_input_file *buffered_input_file1,gt_buffered_input_file *buffered_input_file2,gt_template *template,pthread_mutex_t *mutex) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file1);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input_file2);
  GT_TEMPLATE_CHECK(template);
	gt_status code;
	if((code=gt_input_map_parser_synch_blocks(buffered_input_file1,buffered_input_file2,mutex))) {
		code=gt_input_map_parser_get_template(buffered_input_file1,template,NULL);
		if(code!=GT_IMP_OK) {
			gt_input_map_parser_get_template(buffered_input_file2,template,NULL);
			return code;
		}
		code=gt_input_map_parser_get_alignment(buffered_input_file2,gt_template_get_block_dyn(template,1),NULL);
	}
	return code;
}

gt_status gt_scorereads_process(sr_param *param)
{
	gt_status err=GT_STATUS_OK;

	// Open out file
	gt_output_file *output_file;
	if(param->output_file) {
		output_file=gt_output_file_new_compress(param->output_file,UNSORTED_FILE,param->compress);
	} else {
		output_file=gt_output_stream_new_compress(stdout,UNSORTED_FILE,param->compress);
	}
	gt_cond_fatal_error(!output_file,FILE_OPEN,param->output_file);
	param->printer_attr=gt_generic_printer_attributes_new(MAP);
	param->printer_attr->output_map_attributes->print_casava=true;
	param->printer_attr->output_map_attributes->print_extra=true;
	param->printer_attr->output_map_attributes->print_scores=true;
	param->printer_attr->output_map_attributes->hex_print_scores=true;
	// Do we have two map files as input (one for each read)?
	if(param->input_files[1]) {
		pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
		gt_input_file* input_file1=gt_input_file_open(param->input_files[0],param->mmap_input);
		gt_input_file* input_file2=gt_input_file_open(param->input_files[1],param->mmap_input);
		if(input_file1->file_format!=MAP || input_file2->file_format!=MAP) {
			gt_fatal_error_msg("Fatal error: paired files '%s','%s' are not in MAP format\n",param->input_files[0],param->input_files[1]);
		}
		gt_buffered_input_file* buffered_input_a=gt_buffered_input_file_new(input_file1);
		gt_buffered_input_file* buffered_input_b=gt_buffered_input_file_new(input_file2);
		// If no insert size file supplied, estimate from the first GT_SCOREREADS_NUM_LINES records
		if(!param->ins_phred) {
			uint64_t nread=0;
			gt_template *template=gt_template_new();
			sr_dist_element *ins_dist=NULL;
			gt_status error_code;
			while((error_code=gt_input_map_parser_synch_get_template(buffered_input_a,buffered_input_b,template,&mutex))) {
				if(error_code!=GT_IMP_OK) continue;
				gt_alignment *al[2];
				al[0]=gt_template_get_block(template,0);
				al[1]=gt_template_get_block(template,1);
				int64_t size;
				error_code=gt_scorereads_get_insert_size(template,al,&size);
				if(error_code==GT_STATUS_OK) (void)sr_increase_insert_count(&ins_dist,size,1);
				if(++nread==GT_SCOREREADS_NUM_LINES) break;
			}
			gt_template_delete(template);
			// Reset input buffers to the beginning of the file
			buffered_input_a->cursor=gt_vector_get_mem(buffered_input_a->block_buffer,char);
			buffered_input_a->current_line_num=0;
			buffered_input_b->cursor=gt_vector_get_mem(buffered_input_b->block_buffer,char);
			buffered_input_b->current_line_num=0;
			(void)estimate_insert_histogram(ins_dist,param);
			gt_scorereads_free_insert_hash(ins_dist);
		}

#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(param->num_threads)
#endif
		{
#ifdef HAVE_OPENMP
			uint64_t tid=omp_get_thread_num();
#else
			uint64_t tid=0;
#endif

			gt_buffered_input_file* buffered_input1=tid?gt_buffered_input_file_new(input_file1):buffered_input_a;
			gt_buffered_input_file* buffered_input2=tid?gt_buffered_input_file_new(input_file2):buffered_input_b;
			gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
			gt_buffered_input_file_attach_buffered_output(buffered_input1,buffered_output);
			gt_status error_code;
			gt_template *template=gt_template_new();
			while((error_code=gt_input_map_parser_synch_get_template(buffered_input1,buffered_input2,template,&mutex))) {
				if(error_code!=GT_IMP_OK) continue;
				gt_alignment *al[2];
				al[0]=gt_template_get_block(template,0);
				al[1]=gt_template_get_block(template,1);
				if(!(gt_string_nequals(template->tag,al[1]->tag,gt_string_get_length(template->tag)))) {
					gt_error_msg("Fatal ID mismatch ('%*s','%*s') parsing files '%s','%s'\n",PRIgts_content(template->tag),PRIgts_content(al[1]->tag),param->input_files[0],param->input_files[1]);
					break;
				}
				gt_scorereads_pair_template(template,al,param);
				gt_template_add_mcs_tags(template,al,param->mapping_cutoff);
				if (gt_output_generic_bofprint_template(buffered_output,template,param->printer_attr)) {
					gt_error_msg("Fatal error outputting read '"PRIgts"'\n",PRIgts_content(gt_template_get_string_tag(template)));
				}
			}
			gt_template_delete(template);
			gt_buffered_input_file_close(buffered_input1);
			gt_buffered_input_file_close(buffered_input2);
			gt_buffered_output_file_close(buffered_output);
		}
		gt_input_file_close(input_file1);
		gt_input_file_close(input_file2);
	} else { // Single input file (could be single end or interleaved paired end
		gt_input_file* input_file=param->input_files[0]?gt_input_file_open(param->input_files[0],param->mmap_input):gt_input_stream_open(stdin);
		gt_buffered_input_file* buffered_input_a=gt_buffered_input_file_new(input_file);
		if(gt_input_generic_parser_attributes_is_paired(param->parser_attr)) {
			// If no insert size file supplied, estimate from the first GT_SCOREREADS_NUM_LINES records
			if(!param->ins_phred) {
				uint64_t nread=0;
				gt_template *template=gt_template_new();
				sr_dist_element *ins_dist=NULL;
				gt_status error_code;
				while ((error_code=gt_input_generic_parser_get_template(buffered_input_a,template,param->parser_attr))) {
					if(error_code!=GT_IMP_OK) continue;
					gt_alignment *al[2];
					al[0]=gt_template_get_block(template,0);
					al[1]=gt_template_get_block(template,1);
					int64_t size;
					error_code=gt_scorereads_get_insert_size(template,al,&size);
					if(error_code==GT_STATUS_OK) (void)sr_increase_insert_count(&ins_dist,size,1);
					if(++nread==GT_SCOREREADS_NUM_LINES) break;
				}
				gt_template_delete(template);
				// Reset input buffers to the beginning of the file
				buffered_input_a->cursor=gt_vector_get_mem(buffered_input_a->block_buffer,char);
				buffered_input_a->current_line_num=0;
				(void)estimate_insert_histogram(ins_dist,param);
				gt_scorereads_free_insert_hash(ins_dist);
			}
#ifdef OPENMP
#pragma omp parallel num_threads(param->num_threads)
#endif
			{
#ifdef HAVE_OPENMP
				uint64_t tid=omp_get_thread_num();
#else
				uint64_t tid=0;
#endif
				gt_buffered_input_file* buffered_input=tid?gt_buffered_input_file_new(input_file):buffered_input_a;
				gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
				gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
				gt_status error_code;
				gt_template *template=gt_template_new();
				while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param->parser_attr))) {
					if (error_code!=GT_IMP_OK) continue;
					gt_alignment *al[2];
					al[0]=gt_template_get_block(template,0);
					al[1]=gt_template_get_block(template,1);
					gt_scorereads_pair_template(template,al,param);
					gt_template_add_mcs_tags(template,al,param->mapping_cutoff);
					if (gt_output_generic_bofprint_template(buffered_output,template,param->printer_attr)) {
						gt_error_msg("Fatal error outputting read '"PRIgts"'\n",PRIgts_content(gt_template_get_string_tag(template)));
					}
				}
				// Clean
				gt_template_delete(template);
				gt_buffered_input_file_close(buffered_input);
				gt_buffered_output_file_close(buffered_output);
			}
		} else {
#ifdef OPENMP
#pragma omp parallel num_threads(param->num_threads)
#endif
			{
				gt_buffered_input_file* buffered_input=gt_buffered_input_file_new(input_file);
				gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
				gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
				gt_status error_code;
				gt_template *template=gt_template_new();
				while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param->parser_attr))) {
					if (error_code!=GT_IMP_OK) {
						gt_error_msg("Error parsing file '%s'\n",param->input_files[0]);
						continue;
					}
					gt_alignment *alignment=gt_template_get_block(template,0);
					gt_alignment_recalculate_counters(alignment);
					GT_ALIGNMENT_ITERATE(alignment,map) {
						if(map->gt_score==GT_MAP_NO_GT_SCORE) map->gt_score=calculate_dist_score(alignment,map,param);
						map->phred_score=255;
					}
					gt_alignment_add_mcs_tags(alignment,param->mapping_cutoff);
					if (gt_output_generic_bofprint_alignment(buffered_output,alignment,param->printer_attr)) {
						gt_error_msg("Fatal error outputting read '"PRIgts"'\n",PRIgts_content(gt_template_get_string_tag(template)));
					}
				}
				// Clean
				gt_template_delete(template);
				gt_buffered_input_file_close(buffered_input);
				gt_buffered_output_file_close(buffered_output);
			}
		}
		gt_input_file_close(input_file);
	}
	gt_output_file_close(output_file);
	gt_generic_printer_attributes_delete(param->printer_attr);
	if(param->ins_phred) free(param->ins_phred);
	return err;
}

int main(int argc,char *argv[])
{
	gt_status err=0;
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  err=parse_arguments(argc,argv);
  if(err==GT_STATUS_OK) err=gt_scorereads_process(&param);
	return err=GT_STATUS_OK?0:-1;
}

