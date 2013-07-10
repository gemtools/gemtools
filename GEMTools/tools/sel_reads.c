/*
 * sel_reads.c
 *
 *  Created on: 8 Jul 2013
 *      Author: heath
 */


#include <getopt.h>
#include <omp.h>
#include <pthread.h>
#include "gem_tools.h"
#include "sel_reads.h"

static void usage(FILE *f)
{
	fputs("usage:\n align_stats\n",f);
	fputs("  -i|--insert_dist <insert size distribution file>   (mandatory)\n",f);
	fputs("  -o|--output <output file>                          (default=stdout)\n",f);
	fputs("  -z|--compress\n",f);
	fputs("  -Z|--no-compress                                  (default)\n",f);
	fprintf(f,"  -x|--insert_dist_cutoff <cutoff>                   (default=%g)\n",DEFAULT_INS_CUTOFF);
	fputs("  -r|--reads <reads file or file pair>\n",f);
	fputs("  -t|--threads <number of threads>\n",f);
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

static void *as_malloc(size_t s)
{
	void *p;

	p=malloc(s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
}

static void *as_calloc(size_t n,size_t s)
{
	void *p;

	p=calloc(n,s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
}

static void *as_realloc(void *ptr,size_t s)
{
	void *p;

	p=realloc(ptr,s);
	gt_cond_fatal_error(!p,MEM_HANDLER);
	return p;
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
			{"compress",no_argument,0,'z'},
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
			.compress=false,
			.left=DEFAULT_LTRIM,
			.right=DEFAULT_RTRIM,
			.parser_attr=gt_input_generic_parser_attributes_new(false),
			.num_threads=1,
			.qual_offset=DEFAULT_QUAL,
			.bad_qual=0,
	};
	char c,*p;
	while(!err && (c=getopt_long(argc,argv,"i:t:r:o:q:B:L:R:x:wpzCcZFSh?",longopts,0))!=-1) {
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
			param.compress=true;
			break;
		case 'Z':
			param.compress=false;
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
			param.num_threads=atoi(optarg);
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

	return err;
}
