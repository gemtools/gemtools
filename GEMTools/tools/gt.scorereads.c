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

#define DEFAULT_INS_CUTOFF 0.01 /* Insert sizes in the upper or lower cutoff percentiles will not be used */
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

typedef struct {
  char *input_files[2];
  char *output_file;
  char *dist_file;
  double ins_cutoff;
  bool mmap_input;
  bool verbose;
  gt_output_file_compression compress;
  gt_generic_parser_attributes *parser_attr;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output[2];
  int64_t min_insert;
  int64_t max_insert;
  double *ins_dist;
  uint8_t *ins_phred;
  int num_threads;
  int mapping_cutoff;
  int indel_quality;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina)
} sr_param;

sr_param param = {
		.input_files={NULL,NULL},
		.output_file=NULL,
		.dist_file=NULL,
		.ins_cutoff=DEFAULT_INS_CUTOFF,
		.mmap_input=false,
		.compress=NONE,
		.parser_attr=NULL,
		.verbose=false,
		.num_threads=1,
		.indel_quality=INDEL_QUAL,
		.qual_offset=DEFAULT_QUAL,
		.mapping_cutoff=0,
		.min_insert=0,
		.max_insert=0,
		.ins_dist=NULL,
		.ins_phred=NULL
};

void usage(const gt_option* const options,char* groups[],const bool print_inactive) {
  fprintf(stderr, "USE: ./gt_scorereads [ARGS]...\n");
  gt_options_fprint_menu(stderr,options,groups,false,print_inactive);
}

static void *sr_malloc(size_t s)
{
	void *p;

	p=malloc(s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
}

static void *sr_calloc(size_t n,size_t s)
{
	void *p;

	p=calloc(n,s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
}

static void *sr_realloc(void *ptr,size_t s)
{
	void *p;

	p=realloc(ptr,s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
}

typedef struct {
	int64_t x;
	u_int64_t y;
} hist_entry;

void read_dist_file(sr_param *param,int iset[2])
{
	gt_input_file* file=gt_input_file_open(param->dist_file,false);
	gt_buffered_input_file* bfile=gt_buffered_input_file_new(file);
	gt_status nl=0;
	typedef enum {V1,V2,UNKNOWN} dist_file_type;
	dist_file_type ftype=UNKNOWN;
	bool first=true;
	size_t sz=1024;
	size_t ct=0;
	u_int64_t total=0;
	hist_entry *hist=sr_malloc(sizeof(hist_entry)*sz);
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
									if(first==false && x<=param->max_insert) {
										ef=3;
									} else if(x>=0 && y>0) {
										param->max_insert=x;
										if(first==true) {
											first=false;
											param->min_insert=x;
										}
										if(ct==sz) {
											sz*=1.5;
											hist=sr_realloc(hist,sizeof(hist_entry)*sz);
										}
										hist[ct].x=x;
										hist[ct++].y=y;
										total+=y;
									}
								}
							}
							if(ef) {
								if(ef==3) fprintf(stderr,"Insert distribution file not sorted\n");
								else fprintf(stderr,"Invalid entry in insert distribution file: %s\t%s\n",p,p1);
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
	if(first==false && ct>2) {
		double z1=param->ins_cutoff*(double)total;
		double z2=(1.0-param->ins_cutoff)*(double)total;
		double z=0.0;
		int i;
		for(i=0;i<ct;i++) {
			z+=hist[i].y;
			if(z>=z1) break;
		}
		int i1=i;
		if(!iset[0] || hist[i].x>param->min_insert) param->min_insert=hist[i].x;
		for(i++;i<ct;i++) {
			z+=hist[i].y;
			if(z>=z2) break;
		}
		int i2=i-1;
		if(!iset[1] || hist[i2].x<param->max_insert) param->max_insert=hist[i2].x;
		fprintf(stderr,"Insert distribution %"PRId64" - %"PRId64"\n",param->min_insert,param->max_insert);
		int k=param->max_insert-param->min_insert+1;
		param->ins_dist=sr_malloc(sizeof(double)*k);
		param->ins_phred=sr_malloc((size_t)k);
		for(i=0;i<k;i++) param->ins_dist[i]=0.0;
		for(i=i1;i<=i2;i++) {
			double zt=(double)hist[i].y/(double)total;
			param->ins_dist[hist[i].x-param->min_insert]=zt;
			int phred=255;
			if(zt>0.0) {
				phred=(int)(log(zt)*-10.0/log(10.0)+.5);
				if(phred>255) phred=255;
			}
			param->ins_phred[hist[i].x-param->min_insert]=phred;
		}
	} else {
		if(ftype==UNKNOWN) fprintf(stderr,"Insert distribution file format not recognized\n");
		else fprintf(stderr,"No valid lines read in from insert distribution file\n");
	}
	free(hist);
	gt_buffered_input_file_close(bfile);
	gt_input_file_close(file);
}

static uint64_t calculate_dist_score(gt_alignment *al, gt_map *map, int qual_offset,int qual_penalty)
{
	register gt_string* const quals = al->qualities;
	register const bool has_qualities = gt_alignment_has_qualities(al);
	uint64_t score=0;
	GT_MAP_ITERATE(map,map_block) {
		GT_MISMS_ITERATE(map_block,misms) {
			int quality_misms;
			if (has_qualities) {
				quality_misms = gt_string_get_string(quals)[misms->position]-qual_offset;
				if(quality_misms>MAX_QUAL) quality_misms=MAX_QUAL;
				else if(quality_misms<0) quality_misms=0;
			} else quality_misms=MISSING_QUAL;
			switch (misms->misms_type) {
			case MISMS:
				score+=quality_misms;
				break;
			case INS:
			case DEL:
				score+=qual_penalty;
				break;
			}
		}
	}
	if(score>MAX_GT_SCORE) score=MAX_GT_SCORE;
	return score;
}

static void pair_read(gt_template *template,gt_alignment *alignment1,gt_alignment *alignment2,sr_param *param)
{
	gt_alignment_recalculate_counters(alignment1);
	gt_alignment_recalculate_counters(alignment2);
	gt_mmap_attributes attr;
	gt_map *mmap[2];
	uint64_t nmap[2];
	nmap[0]=gt_alignment_get_num_maps(alignment1);
	nmap[1]=gt_alignment_get_num_maps(alignment2);
	if(nmap[0]+nmap[1]) {
		char *map_flag[2];
		map_flag[0]=sr_calloc((size_t)(nmap[0]+nmap[1]),sizeof(char));
		map_flag[1]=map_flag[0]+nmap[0];
		uint64_t i=0;
		GT_ALIGNMENT_ITERATE(alignment1,map1) {
			if(map1->gt_score==GT_MAP_NO_GT_SCORE) map1->gt_score=calculate_dist_score(alignment1,map1,param->qual_offset,param->indel_quality);
			mmap[0]=map1;
			uint64_t j=0;
			GT_ALIGNMENT_ITERATE(alignment2,map2) {
				if(map2->gt_score==GT_MAP_NO_GT_SCORE) map2->gt_score=calculate_dist_score(alignment2,map2,param->qual_offset,param->indel_quality);
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
		for(i=0;i<nmap[0];i++) {
			if(!map_flag[0][i]) {
				gt_map *map=gt_alignment_get_map(alignment1,i);
				attr.distance=gt_map_get_global_distance(map);
				attr.gt_score=map->gt_score;
				attr.phred_score=255;
				gt_template_inc_counter(template,attr.distance);
				gt_template_add_mmap_ends(template,map,0,&attr);
			}
		}
		for(i=0;i<nmap[1];i++) {
			if(!map_flag[1][i]) {
				gt_map *map=gt_alignment_get_map(alignment2,i);
				attr.distance=gt_map_get_global_distance(map);
				attr.gt_score=(map->gt_score<<16);
				attr.phred_score=255;
				gt_template_inc_counter(template,attr.distance);
				gt_template_add_mmap_ends(template,0,map,&attr);
			}
		}
		free(map_flag[0]);
	}
	gt_attributes_remove(template->attributes,GT_ATTR_ID_TAG_PAIR);
}

int parse_arguments(int argc,char** argv) {
	int err=0;
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
  		param.indel_quality=(int)strtol(optarg,&p,10);
  		if(*p || param.indel_quality<0) {
  			fprintf(stderr,"Illegal indel score: '%s'\n",optarg);
  			err=-7;
  		}
  		break;
  	case 'x':
  		param.ins_cutoff=strtod(optarg,&p);
  		if(*p || param.ins_cutoff>0.5 || param.ins_cutoff<0.0) {
  			fprintf(stderr,"Illegal insert distribution cutoff percentile: '%s'\n",optarg);
  			err=-6;
  		}
  		break;
  	case 'm':
  		param.mapping_cutoff=(int)strtol(optarg,&p,10);
  		if(*p || param.mapping_cutoff<0) {
  			fprintf(stderr,"Illegal mapping cutoff: '%s'\n",optarg);
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
	if(!err && insert_set[0] && insert_set[1] && param.min_insert>param.max_insert) {
		fputs("Minimum insert size > maximum insert size\n",stderr);
    usage(gt_scorereads_options,gt_scorereads_groups,false);
		err=-15;
	}
	if(!err) {
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
		if(gt_input_generic_parser_attributes_is_paired(param.parser_attr) && param.dist_file) read_dist_file(&param,insert_set);
		else if(!insert_set[1]) {
				if(param.min_insert<=1000) param.max_insert=1000;
				else param.max_insert=param.min_insert+1000;
		}
	}
  // Free
  gt_string_delete(gt_scorereads_short_getopt);
  return err;
}

int main(int argc,char *argv[])
{
	int err=0;
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  err=parse_arguments(argc,argv);
  if(!err) {
		// Open out file
		gt_output_file *output_file;
		if(param.output_file) {
			output_file=gt_output_file_new_compress(param.output_file,UNSORTED_FILE,param.compress);
		} else {
			output_file=gt_output_stream_new_compress(stdout,UNSORTED_FILE,param.compress);
		}
		gt_cond_fatal_error(!output_file,FILE_OPEN,param.output_file);
		param.printer_attr=gt_generic_printer_attributes_new(MAP);
		param.printer_attr->output_map_attributes->print_casava=true;
		param.printer_attr->output_map_attributes->print_extra=true;
		param.printer_attr->output_map_attributes->print_scores=true;
		param.printer_attr->output_map_attributes->hex_print_scores=true;
		// Do we have two map files as input (one for each read)?
		if(param.input_files[1]) {
			pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
			gt_input_file* input_file1=gt_input_file_open(param.input_files[0],param.mmap_input);
			gt_input_file* input_file2=gt_input_file_open(param.input_files[1],param.mmap_input);
			if(input_file1->file_format!=MAP || input_file2->file_format!=MAP) {
				gt_fatal_error_msg("Fatal error: paired files '%s','%s' are not in MAP format\n",param.input_files[0],param.input_files[1]);
			}
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(param.num_threads)
#endif
			{
				gt_buffered_input_file* buffered_input1=gt_buffered_input_file_new(input_file1);
				gt_buffered_input_file* buffered_input2=gt_buffered_input_file_new(input_file2);
				gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
				gt_buffered_input_file_attach_buffered_output(buffered_input1,buffered_output);
				gt_status error_code;
				gt_template *template=gt_template_new();
				while(gt_input_map_parser_synch_blocks(buffered_input1,buffered_input2,&mutex)) {
					error_code=gt_input_map_parser_get_template(buffered_input1,template,NULL);
					if(error_code!=GT_IMP_OK) {
						gt_input_map_parser_get_template(buffered_input2,template,NULL);
						gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
						continue;
					}
					if(gt_template_get_num_blocks(template)!=1) {
						gt_error_msg("Error parsing files '%s','%s': wrong number of blocks\n",param.input_files[0],param.input_files[1]);
						continue;
					}
					gt_alignment *alignment1=gt_template_get_block(template,0);
					gt_alignment *alignment2=gt_template_get_block_dyn(template,1);
					error_code=gt_input_map_parser_get_alignment(buffered_input2,alignment2,NULL);
					if (error_code!=GT_IMP_OK) {
						gt_error_msg("Error parsing file '%s'\n",param.input_files[1]);
						continue;
					}
					if(!(gt_string_nequals(template->tag,alignment2->tag,gt_string_get_length(template->tag)))) {
						gt_error_msg("Fatal ID mismatch ('%*s','%*s') parsing files '%s','%s'\n",PRIgts_content(template->tag),PRIgts_content(alignment2->tag),param.input_files[0],param.input_files[1]);
						break;
					}
					pair_read(template,alignment1,alignment2,&param);
					if (gt_output_generic_bofprint_template(buffered_output,template,param.printer_attr)) {
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
			gt_input_file* input_file=param.input_files[0]?gt_input_file_open(param.input_files[0],param.mmap_input):gt_input_stream_open(stdin);
#ifdef OPENMP
#pragma omp parallel num_threads(param.num_threads)
#endif
			{
				gt_buffered_input_file* buffered_input=gt_buffered_input_file_new(input_file);
				gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
				gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
				gt_status error_code;
				gt_template *template=gt_template_new();
				if(gt_input_generic_parser_attributes_is_paired(param.parser_attr)) {
					while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param.parser_attr))) {
						if (error_code!=GT_IMP_OK) {
							gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
							continue;
						}
						gt_alignment *alignment1=gt_template_get_block(template,0);
						gt_alignment *alignment2=gt_template_get_block(template,1);
						pair_read(template,alignment1,alignment2,&param);
						if (gt_output_generic_bofprint_template(buffered_output,template,param.printer_attr)) {
							gt_error_msg("Fatal error outputting read '"PRIgts"'\n",PRIgts_content(gt_template_get_string_tag(template)));
						}
					}
				} else {
					while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param.parser_attr))) {
						if (error_code!=GT_IMP_OK) {
							gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
							continue;
						}
						gt_alignment *alignment=gt_template_get_block(template,0);
						gt_alignment_recalculate_counters(alignment);
						GT_ALIGNMENT_ITERATE(alignment,map) {
							if(map->gt_score==GT_MAP_NO_GT_SCORE) map->gt_score=calculate_dist_score(alignment,map,param.qual_offset,param.indel_quality);
							map->phred_score=255;
						}
						if (gt_output_generic_bofprint_alignment(buffered_output,alignment,param.printer_attr)) {
							gt_error_msg("Fatal error outputting read '"PRIgts"'\n",PRIgts_content(gt_template_get_string_tag(template)));
						}
					}
				}
				// Clean
				gt_template_delete(template);
				gt_buffered_input_file_close(buffered_input);
				gt_buffered_output_file_close(buffered_output);
			}
			gt_input_file_close(input_file);
		}
		gt_output_file_close(output_file);
		gt_generic_printer_attributes_delete(param.printer_attr);
		if(param.ins_dist) {
			free(param.ins_dist);
			free(param.ins_phred);
		}
	}
	return err;
}

