/*
 * sel_reads.c
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
#include "sel_reads.h"

static void usage(FILE *f)
{
	fputs("usage:\n align_stats\n",f);
	fputs("  -i|--insert_dist <insert size distribution file>   (mandatory)\n",f);
	fputs("  -o|--output <output file>                          (default=stdout)\n",f);
#ifdef HAVE_ZLIB
	fputs("  -z|--gzip output with gzip\n",f);
	fputs("  -Z|--no-compress                                  (default)\n",f);
#endif
#ifdef HAVE_BZLIB
	fputs("  -j|--bzip2 output with bzip2\n",f);
#endif
		fprintf(f,"  -x|--insert_dist_cutoff <cutoff>                   (default=%g)\n",DEFAULT_INS_CUTOFF);
	fputs("  -r|--reads <reads file or file pair>\n",f);
#ifdef HAVE_OPENMP
	fputs("  -t|--threads <number of threads>\n",f);
#endif
	fputs("  -m|--mmap      Mmap input files\n",f);
	fputs("  -p|--paired    Paired mapping input file\n",f);
	fprintf(f,"  -F|--fastq     select fastq quality coding         %s\n",DEFAULT_QUAL==QUAL_FASTQ?"(default)":"");
	fprintf(f,"  -S|--solexa    select ilumina quality coding       %s\n",DEFAULT_QUAL==QUAL_SOLEXA?"(default)":"");
	fprintf(f,"  -C|--combined  combined output file for paired reads\n");
	fprintf(f,"  -c|--separate  separate output files for paired reads (default)\n");
	fprintf(f,"  -q|--qual_val  select quality value adjustment     (default=%d)\n",DEFAULT_QUAL);
	fprintf(f,"  -L|--left <n> left trim reads from base 1-n, 0 = do not trim (default=%d)\n",DEFAULT_LTRIM);
	fprintf(f,"  -R|--right <n> right trim reads from base n-end, 0 = do not trim (default=%d)\n",DEFAULT_RTRIM);
	fputs("  -B|--bad_qual <qual> quality value considered as especially bad\n",f);
	fputs("  -h|help|usage                                      (print this file\n\n",f);
}

static void set_opt(char *opt,char **opt_p,char *val)
{
	if(*opt_p) {
		fprintf(stderr,"multiple %s options: '%s' overwriting previous definition '%s'\n",opt,val,*opt_p);
		free(*opt_p);
	}
	*opt_p=strdup(val);
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

void read_dist_file(sr_param *param)
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
		param->min_insert=hist[i].x;
		for(i++;i<ct;i++) {
			z+=hist[i].y;
			if(z>=z2) break;
		}
		int i2=i-1;
		param->max_insert=hist[i2].x;
		fprintf(stderr,"Insert distribution %"PRId64" - %"PRId64"\n",param->min_insert,param->max_insert);
		int k=param->max_insert-param->min_insert+1;
		param->ins_dist=sr_malloc(sizeof(double)*k);
		for(i=0;i<k;i++) param->ins_dist[i]=0.0;
		for(i=i1;i<=i2;i++) {
			param->ins_dist[hist[i].x-param->min_insert]=(double)hist[i].y/(double)total;
		}
	} else {
		if(ftype==UNKNOWN) fprintf(stderr,"Insert distribution file format not recognized\n");
		else fprintf(stderr,"No valid lines read in from insert distribution file\n");
	}
	free(hist);
	gt_buffered_input_file_close(bfile);
	gt_input_file_close(file);
}

void as_select_reads(gt_template *template,sr_param *param)
{
	uint64_t nrd;
	bool paired_file=false; // Was the input file from a paired mapping
	if(gt_input_generic_parser_attributes_is_paired(param->parser_attr)) {
		if(gt_template_get_num_blocks(template)!=2) {
			gt_fatal_error_msg("Fatal error: Expecting paired reads\n");
		}
		nrd=2;
		if(!param->input_files[1]) paired_file=true;
	} else {
		if(gt_template_get_num_blocks(template)!=1) {
			gt_fatal_error_msg("Fatal error: Expecting unpaired reads\n");
		}
		nrd=1;
	}
	printf("nrd=%"PRIu64"\n",nrd);
}

int main(int argc,char *argv[])
{
	int err=0;

	static struct option longopts[]={
			{"reads",required_argument,0,'r'},
			{"insert_dist",required_argument,0,'i'},
			{"insert_dist_cutoff",required_argument,0,'x'},
			{"fastq",no_argument,0,'F'},
			{"solexa",no_argument,0,'S'},
			{"bad_qual",required_argument,0,'B'},
			{"qual_val",required_argument,0,'q'},
			{"output",required_argument,0,'o'},
			{"paired",no_argument,0,'p'},
			{"help",no_argument,0,'h'},
			{"usage",no_argument,0,'h'},
			{"combined",no_argument,0,'C'},
			{"separate",no_argument,0,'c'},
			{"gzip",no_argument,0,'z'},
			{"bzip2",no_argument,0,'j'},
			{"no-compress",no_argument,0,'Z'},
			{"left",required_argument,0,'L'},
			{"right",required_argument,0,'R'},
			{"threads",required_argument,0,'t'},
			{"mmap",no_argument,0,'m'},
			{0,0,0,0}
	};

	sr_param param = {
			.input_files={NULL,NULL},
			.output_file=NULL,
			.dist_file=NULL,
			.ins_cutoff=DEFAULT_INS_CUTOFF,
			.mmap_input=false,
			.combined=false,
			.compress=NONE,
			.left=DEFAULT_LTRIM,
			.right=DEFAULT_RTRIM,
			.parser_attr=gt_input_generic_parser_attributes_new(false),
			.num_threads=1,
			.qual_offset=DEFAULT_QUAL,
			.bad_qual=0,
			.min_insert=0,
			.max_insert=10000,
	};
	char c,*p;
	while(!err && (c=getopt_long(argc,argv,"i:t:r:o:q:B:L:R:x:wpzjCcZFSh?",longopts,0))!=-1) {
		switch(c) {
		case 'i':
			set_opt("insert_dist",&param.dist_file,optarg);
			break;
		case 'p':
			gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
			break;
		case 'o':
			set_opt("output",&param.output_file,optarg);
			break;
		case 'q':
			param.qual_offset=(int)strtol(optarg,&p,10);
			if(*p || param.qual_offset<0 || param.qual_offset>255) {
				fprintf(stderr,"Illegal quality value adjustment: '%s'\n",optarg);
				err=-7;
			}
			break;
		case 'B':
			param.bad_qual=(int)strtol(optarg,&p,10);
			if(*p || param.bad_qual<0 || param.bad_qual>255) {
				fprintf(stderr,"Illegal 'bad' quality value: '%s'\n",optarg);
				err=-7;
			}
			break;
    case 'L':
      param.left=(int)strtol(optarg,&p,10);
      if(*p || param.left<0) {
        fprintf(stderr,"Illegal left trim value: '%s'\n",optarg);
        err=-8;
      }
      break;
    case 'R':
      param.right=(int)strtol(optarg,&p,10);
      if(*p || param.right<0) {
        fprintf(stderr,"Illegal right trim value: '%s'\n",optarg);
        err=-8;
      }
      break;
    case 'x':
      param.ins_cutoff=strtod(optarg,&p);
      if(*p || param.ins_cutoff>0.5 || param.ins_cutoff<0.0) {
        fprintf(stderr,"Illegal insert distribution cutoff percentile: '%s'\n",optarg);
        err=-6;
      }
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
		case 'C':
			param.combined=true;
			break;
		case 'c':
			param.combined=false;
			break;
		case 'F':
			param.qual_offset=QUAL_FASTQ;
			break;
		case 'S':
			param.qual_offset=QUAL_SOLEXA;
			break;
		case 'w':
			param.mmap_input=true;
			break;
		case 't':
#ifdef HAVE_OPENMP
			param.num_threads=atoi(optarg);
#endif
			break;
		case 'r':
			if(param.input_files[0]) {
				fprintf(stderr,"multiple reads options: '%s' overwriting previous definition\n",optarg);
				param.input_files[0]=0;
				param.input_files[1]=0;
			}
			p=strchr(optarg,',');
			if(p) {
				*p++=0;
				if(strchr(p,',')) {
					fprintf(stderr,"Alignment files should be specified either in comma separated pairs (paired end) or individually (single end or paired alignment)\n");
					err=-10;
				} else {
					param.input_files[0]=strdup(optarg);
					param.input_files[1]=strdup(p);
				}
			} else {
				param.input_files[0]=strdup(optarg);
			}
			break;
			fprintf(stderr,"Alignment files should be specified either in comma separated pairs (paired end) or individually (single end or paired alignment)\n");
			break;
		case 'h':
		case '?':
			usage(stdout);
			exit(0);
		}
	}
	if(param.input_files[1]) {
		gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
	}
  if(!err && !param.dist_file && gt_input_generic_parser_attributes_is_paired(param.parser_attr)) {
    fputs("No insert distribution file specified; use -i or --insert_dist options\n",stderr);
    err=-15;
  }
  if(!err && !param.input_files[0]) {
    fputs("No alignment files specified; use -r or --reads options\n",stderr);
    err=-20;
  }
  if(!err && param.left && param.right && param.right<=param.left) {
    fprintf(stderr,"Incompatible trim parameteres (%d,%d)\nNote that both are measured counting from the left of the read\n",param.left,param.right);
    err=-21;
  }
  if(!err) {
  	if(gt_input_generic_parser_attributes_is_paired(param.parser_attr)) read_dist_file(&param);
  	char *qs="";
  	if(param.qual_offset==QUAL_FASTQ) qs="<FASTQ>";
  	else if(param.qual_offset==QUAL_SOLEXA) qs="<SOLEXA>";
  	if(gt_input_generic_parser_attributes_is_paired(param.parser_attr)) {
  		fprintf(stderr,"sel_reads:\ndist_file = %s\ninsert distribution cutoff percentile %g\nquality value adjustment: %d %s\n",param.dist_file,param.ins_cutoff,param.qual_offset,qs);
  	} else {
  		printf("sel_reads: \nquality value adjustment: %d %s\n",param.qual_offset,qs);
  	}
  	// Do we have two map files as input (one for each read)?
  	if(param.input_files[1]) {
  		gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
  		pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
  		gt_input_file* input_file1=gt_input_file_open(param.input_files[0],param.mmap_input);
  		gt_input_file* input_file2=gt_input_file_open(param.input_files[1],param.mmap_input);
  		if(input_file1->file_format!=MAP || input_file2->file_format!=MAP) {
  			gt_fatal_error_msg("Fatal error: paired files '%s','%s' are not in MAP format\n",param.input_files[0],param.input_files[1]);
  		}
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(param.num_threads)
  		{
  			uint64_t tid=omp_get_thread_num();
#else
  		{
  			uint64_t tid=0;
#endif
  			gt_buffered_input_file* buffered_input1=gt_buffered_input_file_new(input_file1);
  			gt_buffered_input_file* buffered_input2=gt_buffered_input_file_new(input_file2);
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
  				gt_alignment *alignment1=gt_template_get_block(template,0);
  				gt_mmap_attributes attr;
  				gt_map *mmap[2];
  				GT_ALIGNMENT_ITERATE(alignment1,map1) {
  					mmap[0]=map1;
  					GT_ALIGNMENT_ITERATE(alignment2,map2) {
  						mmap[1]=map2;
  						gt_status gt_err;
  						int64_t x=gt_template_get_insert_size(mmap,&gt_err,0,0);
  						if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK && x>=param.min_insert && x<=param.max_insert) {
  							attr.distance=gt_map_get_global_distance(map1)+gt_map_get_global_distance(map2);
  							attr.gt_score=GT_MAP_NO_GT_SCORE;
  							gt_template_inc_counter(template,attr.distance+1);
  							gt_template_add_mmap_ends(template,map1,map2,&attr);
  						}
  					}
  				}
  				as_select_reads(template,&param);
  			}
  			gt_template_delete(template);
  			gt_buffered_input_file_close(buffered_input1);
  			gt_buffered_input_file_close(buffered_input2);
  		}
  		gt_input_file_close(input_file1);
  		gt_input_file_close(input_file2);
  	} else { // Single input file (could be singled end or paired end)
  		gt_input_file* input_file=param.input_files[0]?gt_input_file_open(param.input_files[0],param.mmap_input):gt_input_stream_open(stdin);
#ifdef OPENMP
  #pragma omp parallel num_threads(param.num_threads)
  		{
  			uint64_t tid=omp_get_thread_num();
#else
  		{
    		uint64_t tid=0;
#endif
  			gt_buffered_input_file* buffered_input=gt_buffered_input_file_new(input_file);
  			gt_status error_code;
  			gt_template *template=gt_template_new();
  			while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param.parser_attr))) {
  				if (error_code!=GT_IMP_OK) {
  					gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
  					continue;
  				}
  				as_select_reads(template,&param);
  			}
  			// Clean
  			gt_template_delete(template);
  			gt_buffered_input_file_close(buffered_input);
  		}
  		gt_input_file_close(input_file);
  	}
  }
	return err;
}
