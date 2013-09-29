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
#include "score_reads.h"

static void usage(FILE *f)
{
	fputs("usage:\n align_stats <READ1_FILE> [<READ2_FILE>]\n",f);
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
#ifdef HAVE_OPENMP
	fputs(    "  -t|--threads <number of threads>\n",f);
#endif
	fputs(    "  -m|--mmap          Mmap input files\n",f);
	fprintf(f,"  -F|--fastq         select fastq quality coding         %s\n",DEFAULT_QUAL==QUAL_FASTQ?"(default)":"");
	fprintf(f,"  -S|--solexa        select ilumina quality coding       %s\n",DEFAULT_QUAL==QUAL_SOLEXA?"(default)":"");
	fprintf(f,"  -p|--paired        paired read input\n");
	fprintf(f,"  -q|--qual_val      select quality value adjustment     (default=%d)\n",DEFAULT_QUAL);
	fprintf(f,"  -I|--indel_score   indel score\n");
	fprintf(f,"  -r|--min_insert    minimum insert size\n");
	fprintf(f,"  -s|--max_insert    maximum insert size\n");
	fputs(    "  -h|help|usage      (print this file\n\n",f);
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
					attr.gt_score=map1->gt_score+map2->gt_score;
					if(param->ins_phred) attr.gt_score+=param->ins_phred[x-param->min_insert];
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
				attr.gt_score=map->gt_score;
				attr.phred_score=255;
				gt_template_inc_counter(template,attr.distance);
				gt_template_add_mmap_ends(template,0,map,&attr);
			}
		}
		free(map_flag[0]);
	}
	gt_attributes_remove(template->attributes,GT_ATTR_ID_TAG_PAIR);
}

int main(int argc,char *argv[])
{
	int err=0;

	static struct option longopts[]={
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
			{"separate",no_argument,0,'c'},
			{"gzip",no_argument,0,'z'},
			{"bzip2",no_argument,0,'j'},
			{"no-compress",no_argument,0,'Z'},
			{"threads",required_argument,0,'t'},
			{"mmap",no_argument,0,'m'},
			{"indel_score",required_argument,0,'I'},
			{"min_insert",required_argument,0,'r'},
			{"max_insert",required_argument,0,'s'},
			{0,0,0,0}
	};

	sr_param param = {
			.input_files={NULL,NULL},
			.output_file=NULL,
			.dist_file=NULL,
			.ins_cutoff=DEFAULT_INS_CUTOFF,
			.mmap_input=false,
			.compress=NONE,
			.parser_attr=gt_input_generic_parser_attributes_new(false),
			.num_threads=1,
			.indel_quality=INDEL_QUAL,
			.qual_offset=DEFAULT_QUAL,
			.min_insert=0,
			.max_insert=0,
			.ins_dist=NULL,
			.ins_phred=NULL
	};
	char c,*p;
	int insert_set[2]={0,0};
	while(!err && (c=getopt_long(argc,argv,"i:t:o:q:B:x:r:s:wzjpcIZFSh?",longopts,0))!=-1) {
		switch(c) {
		case 'i':
			set_opt("insert_dist",&param.dist_file,optarg);
			break;
		case 'o':
			set_opt("output",&param.output_file,optarg);
			break;
		case 'p':
			gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
			break;
		case 'q':
			param.qual_offset=(int)strtol(optarg,&p,10);
			if(*p || param.qual_offset<0 || param.qual_offset>255) {
				fprintf(stderr,"Illegal quality value adjustment: '%s'\n",optarg);
				err=-7;
			}
			break;
		case 'r':
			param.min_insert=(int)strtol(optarg,&p,10);
			if(*p || param.min_insert<0) {
				fprintf(stderr,"Illegal minimum insert size: '%s'\n",optarg);
				err=-7;
			} else insert_set[0]=1;
			break;
		case 'I':
			param.indel_quality=(int)strtol(optarg,&p,10);
			if(*p || param.indel_quality<0) {
				fprintf(stderr,"Illegal indel score: '%s'\n",optarg);
				err=-7;
			}
			break;
		case 's':
			param.max_insert=(int)strtol(optarg,&p,10);
			if(*p || param.max_insert<0) {
				fprintf(stderr,"Illegal maximum insert size: '%s'\n",optarg);
				err=-7;
			} else insert_set[1]=1;
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
		case 'h':
		case '?':
			usage(stdout);
			exit(0);
		}
	}
	int i;
	for(i=optind;i<argc;i++) {
		if(i-optind>1) {
			fputs("More than two input files specified on command line; extra arguments ignored\n",stderr);
			break;
		}
		param.input_files[i-optind]=strdup(argv[i]);
	}
	if(param.input_files[1]) gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
	if(!err && insert_set[0] && insert_set[1] && param.min_insert>param.max_insert) {
		fputs("Minimum insert size > maximum insert size\n",stderr);
		usage(stderr);
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
					free(param.output_file);
					param.output_file=s;
				}
				break;
			case BZIP2:
				if(l<4 || strcmp(param.output_file+l-4,".bz2")) {
					char *s;
					asprintf(&s,"%s.bz2",param.output_file);
					free(param.output_file);
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
		char *qs="";
		if(param.qual_offset==QUAL_FASTQ) qs="<FASTQ>";
		else if(param.qual_offset==QUAL_SOLEXA) qs="<SOLEXA>";
		if(gt_input_generic_parser_attributes_is_paired(param.parser_attr)) {
			fputs("sel_reads:\n",stderr);
			if(param.dist_file) fprintf(stderr,"dist_file = %s\n",param.dist_file);
			fprintf(stderr,"insert distribution cutoff percentile %g\nquality value adjustment: %d %s\n",param.ins_cutoff,param.qual_offset,qs);
		} else {
			printf("sel_reads: \nquality value adjustment: %d %s\n",param.qual_offset,qs);
		}
		// Open out file
		gt_output_file *output_file;
		if(param.output_file) {
			output_file=gt_output_file_new_compress(param.output_file,UNSORTED_FILE,param.compress);
		} else {
			output_file=gt_output_stream_new_compress(stdout,UNSORTED_FILE,param.compress);
		}
		gt_cond_fatal_error(!output_file,FILE_OPEN,param.output_file);
		gt_buffered_output_file *buffered_output=gt_buffered_output_file_new(output_file);
		param.printer_attr=gt_generic_printer_attributes_new(MAP);
		param.printer_attr->output_map_attributes->print_casava=true;
		param.printer_attr->output_map_attributes->print_extra=true;
		param.printer_attr->output_map_attributes->print_scores=true;
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
			}
			gt_input_file_close(input_file);
		}
		gt_buffered_output_file_close(buffered_output);
		gt_output_file_close(output_file);
		gt_generic_printer_attributes_delete(param.printer_attr);
		if(param.ins_dist) {
			free(param.ins_dist);
			free(param.ins_phred);
		}
	}
	return err;
}
