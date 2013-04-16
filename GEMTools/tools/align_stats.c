#include <getopt.h>
#include <omp.h>
#include <pthread.h>
#include "gem_tools.h"
#include "align_stats.h"

// Every entry will be set to zero be default
static const char base_tab[256]= {
		['A']=2, ['C']=3, ['G']=4, ['T']=5, ['N']=1,
		['a']=2, ['c']=3, ['g']=4, ['t']=5, ['n']=1};

static void usage(FILE *f)
{
	fputs("usage:\n align_stats\n",f);
	fputs("  -r|--reads <reads file or file pair>\n",f);
	fputs("  -o|--output <output stats file>\n",f);
	fputs("  -d|--insert_dist <output insert size distribution file>\n",f);
	fputs("  -M|--min_insert <minimum insert size> (for pairing of single end files: default=0)\n",f);
	fprintf(f,"  -m|--max_insert <maximum insert size> (for pairing of single end files: default=%d)\n",DEFAULT_MAX_INSERT);
	fputs("  -t|--threads <number of threads>\n",f);
	fputs("  -l|--read_length <untrimmed read length>\n",f);
	fputs("  -V|--variable  Variable length reads\n",f);
	fputs("  -p|--paired    Paired mapping input file\n",f);
	fputs("  -i|--ignore_id Do not attempt to parse read IDs\n",f);
	fputs("  -m|--mmap      Mmap input files",f);
	fprintf(f,"  -L|--max_read_length <maximum valid read length>   (default=%u)\n",MAX_READ_LENGTH);
	fprintf(f,"  -F|--fastq     select fastq quality coding         %s\n",DEFAULT_QUAL_OFFSET==QUAL_FASTQ?"(default)":"");
	fprintf(f,"  -S|--solexa    select ilumina quality coding       %s\n",DEFAULT_QUAL_OFFSET==QUAL_SOLEXA?"(default)":"");
	fprintf(f,"  -q|--qual_off  select quality value offset         (default=%d)\n",DEFAULT_QUAL_OFFSET);
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

static as_stats* as_stats_new(bool paired)
{
	as_stats* stats=as_calloc((size_t)1,sizeof(as_stats));
	stats->max_indel_length=50; // We can expand if necessary
	stats->paired=paired;
	int i,j=(paired==true?2:1);
	for(i=0;i<j;i++) {
		int k;
		for(k=0;k<2;k++) stats->indel_length[i*2+k]=as_calloc(sizeof(uint64_t),stats->max_indel_length);
	}
	stats->insert_size=0;
	stats->loc_hash=0;
	return stats;
}

static void as_stats_free(as_stats *stats)
{
	uint64_t i,j=(stats->paired==true?2:1);
	for(i=0;i<j;i++) {
		free(stats->indel_length[i]);
		if(stats->curr_read_store[i]) {
			free(stats->read_length_stats[i]);
			uint64_t k;
			for(k=0;k<stats->curr_read_store[i];k++) free(stats->base_counts_by_cycle[i][k]);
			free(stats->base_counts_by_cycle[i]);
		}
	}
	free(stats);
}

static void as_stats_resize(as_stats *stats,uint64_t rd,uint64_t l)
{
	stats->max_read_length[rd]=l;
	if(l>=stats->curr_read_store[rd]) {
		uint64_t i,nlen=l*1.5; // Allocate a bit more space than we need now to avoid un-necessary re-sizing in future
		if(stats->curr_read_store[rd]) {
			stats->read_length_stats[rd]=as_realloc(stats->read_length_stats[rd],nlen*sizeof(uint64_t));
			stats->base_counts_by_cycle[rd]=as_realloc(stats->base_counts_by_cycle,sizeof(void *)*nlen);
			for(i=stats->curr_read_store[rd];i<nlen;i++) {
				stats->read_length_stats[rd][i]=0;
				stats->base_counts_by_cycle[rd][i]=as_calloc((size_t)(MAX_QUAL+1)*5,sizeof(uint64_t));
			}
		} else {
			stats->read_length_stats[rd]=as_calloc((size_t)nlen,sizeof(uint64_t));
			stats->base_counts_by_cycle[rd]=as_malloc(sizeof(void *)*nlen);
			for(i=0;i<nlen;i++) {
				stats->base_counts_by_cycle[rd][i]=as_calloc((size_t)(MAX_QUAL+1)*5,sizeof(uint64_t));
			}
		}
		stats->curr_read_store[rd]=nlen;
	}
}

static void as_collect_stats(gt_template* template,as_stats* stats,as_param *param) 
{
	stats->nreads++;
	uint64_t nrd;
	bool paired_file=false; // Was the input file from a paired mapping
	if(param->paired_read) {
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
	uint64_t j;
	gt_alignment *al[2];
	char *rd[2],*ql[2];
	uint64_t len[2];
	register int qual_offset=param->qual_offset;
	// Get alignments, qualities, reads and lengths for both ends
	for(j=0;j<nrd;j++) {
		al[j]=gt_template_get_block(template,j);
		rd[j]=gt_alignment_get_read(al[j]);
		ql[j]=gt_alignment_get_qualities(al[j]);
		len[j]=strlen(rd[j]);
		// Update yield max_read_length and resize stats arrays if necessary
		if(stats->max_read_length[j]<len[j]) as_stats_resize(stats,j,len[j]);
		stats->read_length_stats[j][len[j]]++;
		uint64_t i,yld=0;
		char *p=rd[j];
		char *q=ql[j];
		uint64_t **bc=stats->base_counts_by_cycle[j];
		for(i=0;i<len[j];i++) {
			int base=base_tab[(int)p[i]]-1;
			int qual=q[i]-qual_offset;
			if(qual<0 || qual>MAX_QUAL || base<0) {
				gt_fatal_error_msg("Illegal base or quality character '%c %c' in read\n",p[i],q[i]);
			}
			if(base) yld++;
			bc[i][qual*5+base]++;
		}
		stats->yield[j]+=yld;
	}
	// Filter maps (both single and paired end) to remove maps after first zero strata after the first hit
	uint64_t nmaps[3]={0,0,0};
	uint64_t max_dist[3]={0,0,0};
	bool ambig[2]={false,false};
	// Single end alignments
	for(j=0;j<nrd;j++) {
		uint64_t i,k;
		k=gt_alignment_get_num_counters(al[j]);
		for(i=1;i<=k;i++) {
			uint64_t x=gt_alignment_get_counter(al[j],i);
			if(x) {
				nmaps[j]+=x;
				max_dist[j]=i-1;
			} else if(nmaps[j]) break;
		}
		if(i==k && paired_file==false) ambig[j]=true;
	}
	if(nrd==2) {
		// Paired end alignments
		uint64_t i,k;
		k=gt_template_get_num_counters(template);
		for(i=1;i<=k;i++) {
			uint64_t x=gt_template_get_counter(template,i);
			if(x) {
				nmaps[2]+=x;
				max_dist[2]=i-1;
			} else if(nmaps[2]) break;
		}
		// Now count number of split maps for each end
		uint64_t nsplit[3]={0,0,0};
		for(j=0;j<2;j++) {
			bool flg=false;
			GT_ALIGNMENT_ITERATE(al[j],map) {
				if(gt_map_get_distance(map)<=max_dist[j]) {
					if(gt_map_get_num_blocks(map)>1) nsplit[j]++;
					else flg=true;
				}
			}
			if(nsplit[j]) {
				stats->reads_with_splitmaps[j]++;
				if(flg==false) stats->reads_only_with_splitmaps[j]++;
			}
		}
		// And for paired alignments
		bool flg=false;
		GT_TEMPLATE__ATTR_ITERATE_(template,maps,maps_attr) {
			if(maps_attr->distance<=max_dist[2]) {
				if(gt_map_get_num_blocks(maps[0])>1 || gt_map_get_num_blocks(maps[1])>1) nsplit[2]++;
				else flg=true;
			}
		}
		if(nsplit[2]) {
			stats->reads_with_splitmaps[2]++;
			if(flg==false) stats->reads_only_with_splitmaps[2]++;
		}
		for(j=0;j<2;j++) {
			if(nmaps[j]) {
				stats->mapped[j]++;
				if(nmaps[j]==1) {
					if(ambig[j]==false) stats->unique[j]++;
					else stats->ambiguous[j]++;
				}
			}
		}
		if(nmaps[2]) {
			stats->paired_mapped++;
			if(nmaps[2]==1 && (nmaps[0]<=gt_alignment_get_num_maps(al[0])) && (nmaps[1]<=gt_alignment_get_num_maps(al[1]))) stats->paired_unique++;
		}
	}
}

static void as_merge_stats(as_stats **st,uint64_t nt,bool paired)
{
	uint64_t i,j,k,k1,nr,len[2];

	nr=paired?2:1;
	for(j=0;j<nr;j++) {
		len[j]=st[0]->max_read_length[j];
		for(i=1;i<nt;i++) {
			if(st[i]->max_read_length[j]>len[j]) len[j]=st[i]->max_read_length[j];
		}
		if(len[j]>st[0]->max_read_length[j]) as_stats_resize(st[0],j,len[j]);
	}
	for(i=1;i<nt;i++) {
		st[0]->nreads+=st[i]->nreads;
		for(j=0;j<nr;j++) {
			st[0]->yield[j]+=st[i]->yield[j];
			st[0]->mapped[j]+=st[i]->mapped[j];
			st[0]->unique[j]+=st[i]->unique[j];
			st[0]->ambiguous[j]+=st[i]->ambiguous[j];
			len[j]=st[i]->max_read_length[j];
			uint64_t **bc0=st[0]->base_counts_by_cycle[j];
			uint64_t **bc1=st[i]->base_counts_by_cycle[j];
			for(k=0;k<=len[j];k++) {
				st[0]->read_length_stats[j][k]+=st[i]->read_length_stats[j][k];
				for(k1=0;k1<5*(MAX_QUAL+1);k1++) bc0[k][k1]+=bc1[k][k1];
			}
		}
		for(j=0;j<3;j++) {
			st[0]->paired_type[j]+=st[i]->paired_type[j];
			st[0]->bis_stats[j]+=st[i]->bis_stats[j];
			st[0]->reads_with_splitmaps[j]+=st[i]->reads_with_splitmaps[j];
			st[0]->reads_only_with_splitmaps[j]+=st[i]->reads_only_with_splitmaps[j];
		}
		st[0]->paired_mapped+=st[i]->paired_mapped;
		st[0]->paired_unique+=st[i]->paired_unique;
		as_stats_free(st[i]);
	}
}

static void as_print_yield_summary(FILE *f,as_stats *st,as_param *param)
{
	bool paired=param->paired_read;
	uint64_t trimmed[2]={0,0},yield[2]={0,0},min_rl[2],i,j,k;
	fputs("Yield summary\n\n",f);
	j=paired?2:1;
	for(i=0;i<j;i++) {
		uint64_t l=st->max_read_length[i];
		for(k=0;k<=l;k++) if(st->read_length_stats[i][k]) break;
		min_rl[i]=k;
		for(;k<=l;k++) {
			uint64_t x=st->read_length_stats[i][k];
			trimmed[i]+=(l-k)*x;
			yield[i]+=k*x;
		}
	}
	if(paired) {
		fprintf(f,"Paired end reads.  No. pairs =\t%" PRIu64 "\n",st->nreads);
		if(!param->variable_read_length) {
			fprintf(f,"Read lengths:\tRead 1 =\t%" PRIu64 "\tRead 2 =\t%" PRIu64 "\n",st->max_read_length[0],st->max_read_length[1]);
			fprintf(f,"Yield PF:\tRead 1 = \t%" PRIu64 "\tRead 2 = \t%" PRIu64 "\tTotal = \t%" PRIu64 "\n",yield[0]+trimmed[0],yield[1]+trimmed[1],yield[0]+yield[1]+trimmed[0]+trimmed[1]);
			fprintf(f,"Bases trimmed:\tRead 1 = \t%" PRIu64 "\t(%.2f%%)\tRead 2 = \t%" PRIu64 "\t(%.2f%%)\tTotal = \t%" PRIu64 "\t(%.2f%%)\n",
					trimmed[0],100.0*(double)trimmed[0]/(double)(yield[0]+trimmed[0]),
					trimmed[1],100.0*(double)trimmed[1]/(double)(yield[1]+trimmed[1]),
					trimmed[0]+trimmed[1],100.0*(double)(trimmed[0]+trimmed[1])/(double)(yield[0]+yield[1]+trimmed[0]+trimmed[1]));
		} else {
			fprintf(f,"Read lengths:\tRead 1 =\t%" PRIu64 " - %" PRIu64 "\tRead 2 =\t%" PRIu64 " - %" PRIu64 "\n",min_rl[0],st->max_read_length[0],min_rl[1],st->max_read_length[1]);
			fprintf(f,"Yield PF:\tRead 1 =\t%" PRIu64 "\tRead 2 =\t%" PRIu64 "\tTotal =\t%" PRIu64 "\n",yield[0],yield[1],yield[0]+yield[1]);
		}
		fprintf(f,"No calls:\tRead 1 =\t%" PRIu64 "\t(%.2f%%)\tRead 2 =\t%" PRIu64 "\t(%.2f%%)\tTotal =\t%" PRIu64 "\t(%.2f%%)\n",
				yield[0]-st->yield[0],100.0*(double)(yield[0]-st->yield[0])/(double)yield[0],
				yield[1]-st->yield[1],100.0*(double)(yield[1]-st->yield[1])/(double)yield[1],
				yield[0]+yield[1]-st->yield[0]-st->yield[1],100.0*(double)(yield[0]+yield[1]-st->yield[0]-st->yield[1])/(double)(yield[0]+yield[1]));
		fprintf(f,"%s yield:\tRead 1 =\t%" PRIu64 "\t(%.2f%%)\tRead 2 =\t%" PRIu64 "\t(%.2f%%)\tTotal =\t%" PRIu64 "\t(%.2f%%)\n",param->variable_read_length?"Clean":"Trimmed",
				st->yield[0],100.0*(double)st->yield[0]/(double)(yield[0]+trimmed[0]),
				st->yield[1],100.0*(double)st->yield[1]/(double)(yield[1]+trimmed[1]),
				st->yield[0]+st->yield[1],100.0*(double)(st->yield[0]+st->yield[1])/(double)(yield[0]+yield[1]+trimmed[0]+trimmed[1]));
	} else {
		fprintf(f,"Single end reads.  No. reads =\t%" PRIu64 "\n",st->nreads);
		if(!param->variable_read_length) {
			fprintf(f,"Read length:\t%" PRIu64 "\n",st->max_read_length[0]);
			fprintf(f,"Yield PF:\t%" PRIu64 "\n",yield[0]+trimmed[0]);
			fprintf(f,"Bases trimmed:\t=\t%" PRIu64 "\t(%.2f%%)\n",
					trimmed[0],100.0*(double)trimmed[0]/(double)(yield[0]+trimmed[0]));
		} else {
			fprintf(f,"Read length:\t%" PRIu64 " - %" PRIu64 "\n",min_rl[0],st->max_read_length[0]);
			fprintf(f,"Yield PF:\t%" PRIu64 "\n",yield[0]);
		}
		fprintf(f,"No calls =\t%" PRIu64 "\t(%.2f%%)\n",yield[0]-st->yield[0],100.0*(double)(yield[0]-st->yield[0])/(double)yield[0]);
		fprintf(f,"%s yield:\t%" PRIu64 "\t(%.2f%%)\n",param->variable_read_length?"Clean":"Trimmed",st->yield[0],100.0*(double)st->yield[0]/(double)(yield[0]+trimmed[0]));
	}
}

static void as_print_mapping_summary(FILE *f,as_stats *st,as_param *param)
{
	bool paired=param->paired_read;
	bool paired_file=false; // Was the input file from a paired mapping
	if(paired==true && !param->input_files[1]) paired_file=true;
	fputs("\nSingle end mapping summary\n\n",f);
	double z=(double)st->nreads;
	if(paired==true) {
		fprintf(f,"Uniquely mapping reads:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n",
				st->unique[0],100.0*(double)st->unique[0]/z,
				st->unique[1],100.0*(double)st->unique[1]/z,
				st->unique[0]+st->unique[1],100.0*(double)(st->unique[0]+st->unique[1])/(z+z));
		uint64_t mult[3];
		mult[0]=st->mapped[0]-st->unique[0]-st->ambiguous[0];
		mult[1]=st->mapped[1]-st->unique[1]-st->ambiguous[1];
		fprintf(f,"Multiply mapping reads:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n",
				mult[0],100.0*(double)mult[0]/z,
				mult[1],100.0*(double)mult[1]/z,
				mult[0]+mult[1],100.0*(double)(mult[0]+mult[1])/(z+z));
		uint64_t unmap[3];
		unmap[0]=st->nreads-st->mapped[0];
		unmap[1]=st->nreads-st->mapped[1];
		fprintf(f,"Unmapped reads:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n",
				unmap[0],100.0*(double)unmap[0]/z,
				unmap[1],100.0*(double)unmap[1]/z,
				unmap[0]+unmap[1],100.0*(double)(unmap[0]+unmap[1])/(z+z));
		if(paired_file==false) {
			fprintf(f,"Ambiguous reads:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n\n",
					st->ambiguous[0],100.0*(double)st->ambiguous[0]/z,
					st->ambiguous[1],100.0*(double)st->ambiguous[1]/z,
					st->ambiguous[0]+st->ambiguous[1],100.0*(double)(st->ambiguous[0]+st->ambiguous[1])/(z+z));
		}
		fprintf(f,"Reads with splitmaps:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n",
				st->reads_with_splitmaps[0],100.0*(double)st->reads_with_splitmaps[0]/z,
				st->reads_with_splitmaps[1],100.0*(double)st->reads_with_splitmaps[1]/z,
				st->reads_with_splitmaps[0]+st->reads_with_splitmaps[1],100.0*(double)(st->reads_with_splitmaps[0]+st->reads_with_splitmaps[1])/(z+z));
		fprintf(f,"Reads with only splitmaps:\tRead 1 =\t%" PRIu64 "\t(%g%%)\tRead 2 =\t%" PRIu64 "\t(%g%%)\tTotal =\t%" PRIu64 "\t(%g%%)\n",
				st->reads_only_with_splitmaps[0],100.0*(double)st->reads_only_with_splitmaps[0]/z,
				st->reads_only_with_splitmaps[1],100.0*(double)st->reads_only_with_splitmaps[1]/z,
				st->reads_only_with_splitmaps[0]+st->reads_only_with_splitmaps[1],100.0*(double)(st->reads_only_with_splitmaps[0]+st->reads_only_with_splitmaps[1])/(z+z));
		fputs("\nPaired end mapping summary\n\n",f);
		fprintf(f,"Uniquely mapping read pairs:\t%" PRIu64 "\t(%g%%)\n",st->paired_unique,100.0*(double)st->paired_unique/z);
		mult[2]=st->paired_mapped-st->paired_unique;
		fprintf(f,"Multiply mapping reads:\t%" PRIu64 "\t(%g%%)\n",mult[2],100.0*(double)mult[2]/z);
		unmap[2]=st->nreads-st->paired_mapped;
		fprintf(f,"Unmapped read pairs:\t%" PRIu64 "\t(%g%%)\n",unmap[2],100.0*(double)unmap[2]/z);
		fprintf(f,"Read pairs with splitmaps:\t%" PRIu64 "\t(%g%%)\n",st->reads_with_splitmaps[2],100.0*(double)st->reads_with_splitmaps[2]/z);
		fprintf(f,"Read pairs with only splitmaps:\t%" PRIu64 "\t(%g%%)\n",st->reads_only_with_splitmaps[2],100.0*(double)st->reads_only_with_splitmaps[2]/z);
	} else {
		fprintf(f,"Uniquely mapping reads:\t%" PRIu64 "\t(%g%%)\n",st->unique[0],100.0*(double)st->unique[0]/z);
		uint64_t mult=st->mapped[0]-st->unique[0]-st->ambiguous[0];
		fprintf(f,"Multiply mapping reads:\t%" PRIu64 "\t(%g%%)\n",mult,100.0*(double)mult/z);
		uint64_t unmap=st->nreads-st->mapped[0];
		fprintf(f,"Unmapped reads:\t%" PRIu64 "\t(%g%%)\n",unmap,100.0*(double)unmap/z);
		fprintf(f,"Ambiguous mapping reads:\t%" PRIu64 "\t(%g%%)\n\n",st->ambiguous[0],100.0*(double)st->ambiguous[0]/z);
		fprintf(f,"Reads with splitmaps:\t%" PRIu64 "\t(%g%%)\n",st->reads_with_splitmaps[0],100.0*(double)st->reads_with_splitmaps[0]/z);
		fprintf(f,"Reads with only splitmaps:\t%" PRIu64 "\t(%g%%)\n",st->reads_only_with_splitmaps[0],100.0*(double)st->reads_only_with_splitmaps[0]/z);
	}
}

static void as_print_read_lengths(FILE *f,as_stats *st,bool paired)
{
	if(!st->max_read_length[0]) return;
	fprintf(f,"\n\nDistribution of reads lengths after trimming\n");
	uint64_t i,l,j,k;
	if(paired) {
		fputs("Read length\tR1:n_reads\tR1:p\tR2:nreads\tR2:p\n",f);
		j=2;
		l=st->max_read_length[0]>st->max_read_length[1]?st->max_read_length[0]:st->max_read_length[1];
	} else {
		fputs("Read length\tn_reads\tp\n",f);
		j=1;
		l=st->max_read_length[0];
	}
	double tot[2]={0.0,0.0};
	for(i=0;i<=l;i++) {
		for(k=0;k<j;k++) {
			if(i<=st->max_read_length[k]) tot[k]+=(double)st->read_length_stats[k][i];
		}
	}
	uint64_t x[2];
	for(i=0;i<=l;i++) {
		for(k=0;k<j;k++) x[k]=(i<=st->max_read_length[k]?st->read_length_stats[k][i]:0);
		if(paired) {
			if(x[0]||x[1]) {
				fprintf(f,"%" PRIu64 "\t%" PRIu64 "\t%.4f\t%" PRIu64 "\t%.4f\n",i,x[0],(double)x[0]/tot[0],x[1],(double)x[1]/tot[1]);
			}
		} else if(x[0]) {
			fprintf(f,"%" PRIu64 "\t%" PRIu64 "\t%.4f\n",i,x[0],(double)x[0]/tot[0]);
		}
	}
}

static void as_print_stats(as_stats *st,as_param *param)
{
	FILE *fout;
	if(param->output_file) {
		fout=fopen(param->output_file,"w");
	} else {
		fout=stdout;
	}
	as_print_yield_summary(fout,st,param);
	as_print_mapping_summary(fout,st,param);
	as_print_read_lengths(fout,st,param->paired_read);
}

int main(int argc,char *argv[])
{
	int err=0,c;
	char *p,*p1;

	static struct option longopts[]={
			{"reads",required_argument,0,'r'},
			{"insert_dist",required_argument,0,'d'},
			{"max_insert",required_argument,0,'m'},
			{"min_insert",required_argument,0,'M'},
			{"insert_dist",required_argument,0,'d'},
			{"paired",no_argument,0,'p'},
			{"variable",no_argument,0,'V'},
			{"ignore_id",no_argument,0,'i'},
			{"mmap",no_argument,0,'m'},
			{"fastq",no_argument,0,'F'},
			{"solexa",no_argument,0,'S'},
			{"threads",required_argument,0,'t'},
			{"qual_off",required_argument,0,'q'},
			{"output",required_argument,0,'o'},
			{"read_length",required_argument,0,'l'},
			{"max_read_length",required_argument,0,'L'},
			{"help",no_argument,0,'h'},
			{"usage",no_argument,0,'h'},
			{0,0,0,0}
	};

	as_param param = {
			.input_files={NULL,NULL},
			.output_file=NULL,
			.dist_file=NULL,
			.mmap_input=false,
			.parser_attr=GENERIC_PARSER_ATTR_DEFAULT(false),
			.paired_read=false,
			.ignore_id=false,
			.min_insert=0,
			.max_insert=DEFAULT_MAX_INSERT,
			.variable_read_length=false,
			.read_length={0,0},
			.max_read_length=MAX_READ_LENGTH,
			.num_threads=1,
			.qual_offset=DEFAULT_QUAL_OFFSET,
	};

	while(!err && (c=getopt_long(argc,argv,"d:t:r:o:q:m:M:l:L:x:FSVpi?",longopts,0))!=-1) {
		switch(c) {
		case 'd':
			set_opt("insert_dist",&param.dist_file,optarg);
			break;
		case 'p':
			param.paired_read=true;
			break;
		case 'o':
			set_opt("output",&param.output_file,optarg);
			break;
		case 'L':
			param.max_read_length=(uint64_t)strtoul(optarg,&p,10);
			break;
		case 'l':
			param.read_length[0]=(uint64_t)strtoul(optarg,&p,10);
			if(*p==',') param.read_length[1]=(uint64_t)strtoul(p+1,&p1,10);
			else param.read_length[1]=param.read_length[0];
			break;
		case 'q':
			param.qual_offset=(int)strtol(optarg,&p,10);
			if(*p || param.qual_offset<0 || param.qual_offset>255) {
				fprintf(stderr,"Illegal quality value adjustment: '%s'\n",optarg);
				err=-7;
			}
			break;
		case 'm':
			param.max_insert=(uint64_t)strtoul(optarg,&p,10);
			break;
		case 'M':
			param.min_insert=(uint64_t)strtoul(optarg,&p,10);
			break;
		case 'F':
			param.qual_offset=QUAL_FASTQ;
			break;
		case 'S':
			param.qual_offset=QUAL_SOLEXA;
			break;
		case 'V':
			param.variable_read_length=true;
			break;
		case 'i':
			param.ignore_id=true;
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
	as_stats** stats=as_malloc(param.num_threads*sizeof(void *));
	// Do we have two map files as input (one for each read)?
	if(param.input_files[1]) {
		param.paired_read=true;
		pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
		gt_input_file* input_file1=gt_input_file_open(param.input_files[0],param.mmap_input);
		gt_input_file* input_file2=gt_input_file_open(param.input_files[1],param.mmap_input);
		if(input_file1->file_format!=MAP || input_file2->file_format!=MAP) {
			gt_fatal_error_msg("Fatal error: paired files '%s','%s' are not in MAP format\n",param.input_files[0],param.input_files[1]);
		}
#pragma omp parallel num_threads(param.num_threads)
		{
			uint64_t tid=omp_get_thread_num();
			gt_buffered_input_file* buffered_input1=gt_buffered_input_file_new(input_file1);
			gt_buffered_input_file* buffered_input2=gt_buffered_input_file_new(input_file2);
			gt_status error_code;
			gt_template *template=gt_template_new();
			stats[tid]=as_stats_new(param.paired_read);
			while(gt_input_map_parser_synch_blocks(buffered_input1,buffered_input2,&mutex)) {
				error_code=gt_input_map_parser_get_template(buffered_input1,template);
				if(error_code!=GT_IMP_OK) {
					gt_input_map_parser_get_template(buffered_input2,template);
					gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
					continue;
				}
				if(gt_template_get_num_blocks(template)!=1) {
					gt_error_msg("Error parsing files '%s','%s': wrong number of blocks\n",param.input_files[0],param.input_files[1]);
					continue;
				}
				gt_alignment *alignment2=gt_template_get_block_dyn(template,1);
				error_code=gt_input_map_parser_get_alignment(buffered_input2,alignment2);
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
						uint64_t gt_err;
						int64_t x=gt_template_get_insert_size(mmap,&gt_err);
						if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK && x>=param.min_insert && x<=param.max_insert) {
							attr.distance=gt_map_get_global_distance(map1)+gt_map_get_global_distance(map2);
							attr.gt_score=GT_MAP_NO_GT_SCORE;
							gt_template_inc_counter(template,attr.distance+1);
							gt_template_add_mmap_va(template,&attr,map1,map2);
						}
					}
				}
				as_collect_stats(template,stats[tid],&param);
			}
			gt_template_delete(template);
			gt_buffered_input_file_close(buffered_input1);
			gt_buffered_input_file_close(buffered_input2);
		}
		gt_input_file_close(input_file1);
		gt_input_file_close(input_file2);
	} else { // Single file (could be single or paired end)
		gt_input_file* input_file=param.input_files[0]?gt_input_file_open(param.input_files[0],param.mmap_input):gt_input_stream_open(stdin);
#pragma omp parallel num_threads(param.num_threads)
		{
			uint64_t tid=omp_get_thread_num();
			gt_buffered_input_file* buffered_input=gt_buffered_input_file_new(input_file);
			gt_status error_code;
			gt_template *template=gt_template_new();
			stats[tid]=as_stats_new(param.paired_read);
			while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,&param.parser_attr))) {
				if (error_code!=GT_IMP_OK) {
					gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
					continue;
				}
				// For paired reads, insert single end mappings into alignments
				if(param.paired_read) {
					if(gt_template_get_num_blocks(template)!=2) {
						gt_fatal_error_msg("Fatal error: Expecting paired reads\n");
					}
					gt_alignment *al[2];
					al[0]=gt_template_get_block(template,0);
					al[1]=gt_template_get_block(template,1);
					gt_alignment_recalculate_counters(al[0]);
					gt_alignment_recalculate_counters(al[1]);
				}
				as_collect_stats(template,stats[tid],&param);
			}
			// Clean
			gt_template_delete(template);
			gt_buffered_input_file_close(buffered_input);
		}
		gt_input_file_close(input_file);
	}
	as_merge_stats(stats,param.num_threads,param.paired_read);
	as_print_stats(stats[0],&param);
	as_stats_free(stats[0]);
	free(stats);
	return err;
}
