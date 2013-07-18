#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "gem_tools.h"
#include "align_stats.h"

// Every entry will be set to zero be default
static const char base_tab[256]= {
		['A']=2, ['C']=3, ['G']=4, ['T']=5, ['N']=1,
		['a']=2, ['c']=3, ['g']=4, ['t']=5, ['n']=1};

static const char id_brk_char[256]= {
		[':']=ID_COLON_CHAR, ['#']=ID_HASH_CHAR, [' ']=ID_SPACE_CHAR, ['/']=ID_SLASH_CHAR, [0]=ID_END_CHAR
};

static void usage(FILE *f)
{
	fputs("usage:\n align_stats\n",f);
	fputs("  -r|--reads <reads file or file pair>\n",f);
	fputs("  -o|--output <output stats file>\n",f);
	fputs("  -d|--insert_dist <output insert size distribution file>\n",f);
	fputs("  -M|--min_insert <minimum insert size> (for pairing of single end files: default=0)\n",f);
	fprintf(f,"  -m|--max_insert <maximum insert size> (for pairing of single end files: default=%d)\n",DEFAULT_MAX_INSERT);
	fputs("  -t|--threads <number of threads>\n",f);
	fputs("  -z|--gzip (compress output files with gzip\n",f);
	fputs("  -j|--bzip2 (compress output files with bzip2\n",f);
	fputs("  -Z|--no=compress (default)\n",f);
	fputs("  -l|--read_length <untrimmed read length>\n",f);
	fputs("  -V|--variable  Variable length reads\n",f);
	fputs("  -p|--paired    Paired mapping input file\n",f);
	fputs("  -i|--ignore_id Do not attempt to parse read IDs\n",f);
	fputs("  -w|--mmap      mmap input files\n",f);
	fprintf(f,"  -P|--phage_lambda <identifier for phage lambda>    (default='%s')\n",PHAGE_LAMBDA);
	fprintf(f,"  -X|--phix174 <identifier for phiX174>    (default='%s')\n",PHIX174);
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

static void as_set_output_files(as_param *param)
{
	char *csuff[3]={"",".gz",".bz2"};
	char *cs;

	if(!(param->output_file && param->dist_file)) {
		switch(param->compress) {
		case GZIP:
			cs=csuff[1];
			break;
		case BZIP2:
			cs=csuff[2];
			break;
		default:
			cs=csuff[0];
			break;
		}
		if(param->input_files[0]) {
			char *root=strdup(param->input_files[0]);
			char *p=strrchr(root,'/');
			if(p) root=p+1;
			p=strchr(root,'.');
			if(p) *p=0;
			if(param->input_files[1]) {
				if(p-root>2 && (p[-1]=='1' || p[-1]=='2')) {
					if(p[-2]=='.' || p[-2]=='_') p[-2]=0;
					else p[-1]=0;
				}
				else if(p-root>6 && !strncmp(p-4,"_map",4) && (p[-5]=='1' || p[-5]=='2') && p[-6]=='_') p[-6]=0;
			} else {
				if(p-root>4 && !strncmp(p-4,"_map",4)) p[-4]=0;
			}
			if(!param->output_file) asprintf(&param->output_file,"%s_report.txt%s",root,cs);
			if(!param->dist_file) asprintf(&param->dist_file,"%s_frag_dist.txt%s",root,cs);
		} else {
			if(!param->output_file) asprintf(&param->output_file,"align_stats_report.txt%s",cs);
			if(!param->dist_file) asprintf(&param->dist_file,"align_stats_frag_dist.txt%s",cs);
		}
	}
	if(!strcmp(param->output_file,"-")) param->output_file=0;
	if(!strcmp(param->dist_file,"-")) param->dist_file=0;
}

static id_tag *new_id_tag(void)
{
	id_tag *idt;

	idt=as_malloc(sizeof(id_tag));
	idt->instrument_name=gt_string_new(128);
	idt->run=gt_string_new(128);
	idt->flowcell=gt_string_new(128);
	idt->index=gt_string_new(128);
	return idt;
}

static void clear_id_tag(id_tag *idt)
{
	gt_string_clear(idt->instrument_name);
	gt_string_clear(idt->run);
	gt_string_clear(idt->flowcell);
	gt_string_clear(idt->index);
}

static void free_id_tag(id_tag *idt)
{
	gt_string_delete(idt->instrument_name);
	gt_string_delete(idt->run);
	gt_string_delete(idt->flowcell);
	gt_string_delete(idt->index);
	free(idt);
}

static uint64_t parse_id_tag(gt_string *tag,id_tag *idt)
{

	// @HWUSI-EAS100R:6:73:941:1973#0/1   Old style Casava tag
  // Machine:lane:tile:x:y#multiplex tag/read
	//
	// @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG  New style (v1.8) tag
	// Machine:run:flowcell:lane:tile:x:y read:filter:flags:multiplex tag
  //
	// The end of the tag (/read for old style or the part after the space for the new style) may have been trimmed
	//
	char type[11],*tg,*p;
	int fields[11],len[11];
	int i=0,ix=0,err=ID_TAG_OK;

	clear_id_tag(idt);
	tg=gt_string_get_string(tag);
	ix=0;
	fields[0]=0;
	register char c;
	while(tg[ix] && i<11) {
		while(!(c=id_brk_char[(int)tg[ix]])) ix++;
		len[i]=ix+1-fields[i];
		type[i++]=c;
		ix++;
		if(c==ID_END_CHAR || i==11) break;
		if(c==ID_SPACE_CHAR) while((c=id_brk_char[(int)tg[ix]])==ID_SPACE_CHAR) ix++; // Skip multiple white space characters
		fields[i]=ix;
	}
	int j;
	for(j=0;j<i;j++) if(type[j]!=ID_COLON_CHAR) break;
	int sidx;
	if((i==6 || i==7) && j==4 && type[j]==ID_HASH_CHAR) {
		gt_string_copy_substr(idt->instrument_name,tag,fields[0],len[0]);
		sidx=1;
	} else {
		gt_string_copy_substr(idt->instrument_name,tag,fields[0],len[0]);
		gt_string_copy_substr(idt->run,tag,fields[1],len[1]);
		gt_string_copy_substr(idt->flowcell,tag,fields[2],len[2]);
		sidx=3;
	}
	idt->lane=(uint32_t)strtoul(tg+fields[sidx],&p,10);
	if(idt->lane<1 || idt->lane>MAX_LANE_ID) {
		err=ID_TAG_ERROR_BAD_LANE;
	}
	if(err==ID_TAG_OK) {
		idt->tile=(uint32_t)strtoul(tg+fields[sidx+1],&p,10);
		if(idt->tile<1) {
			err=ID_TAG_ERROR_BAD_TILE;
		}
	}
	if(err==ID_TAG_OK) {
		idt->x=(uint32_t)strtoul(tg+fields[sidx+2],&p,10);
		if(idt->x<1) {
			err=ID_TAG_ERROR_BAD_COORD;
		}
	}
	if(err==ID_TAG_OK) {
		idt->y=(uint32_t)strtoul(tg+fields[sidx+3],&p,10);
		if(idt->y<1) {
			err=ID_TAG_ERROR_BAD_COORD;
		}
	}
// if(err==ID_TAG_OK) {
// printf("Machine: " PRIgts "\tRun: " PRIgts "\tFC: " PRIgts "\tLane: %d\tTile: %" PRIu32 "\tX,Y: %" PRIu32 ",%" PRIu32 "\n",PRIgts_content(idt->instrument_name),PRIgts_content(idt->run),PRIgts_content(idt->flowcell),idt->lane,idt->tile,idt->x,idt->y);
//	}
	return err;
}

static as_stats* as_stats_new(bool paired)
{
	as_stats* stats=as_calloc((size_t)1,sizeof(as_stats));
	stats->max_indel_length=50; // We can expand if necessary
	stats->paired=paired;
	int i,j=(paired==true?2:1);
	for(i=0;i<j;i++) {
		int k;
		for(k=0;k<2;k++) stats->indel_length[i*2+k]=as_calloc(sizeof(uint64_t),stats->max_indel_length+1);
	}
	stats->insert_size=0;
	stats->loc_hash=0;
	return stats;
}

static void as_stats_free(as_stats *stats)
{
	uint64_t i,j=(stats->paired==true?2:1);
	for(i=0;i<j;i++) {
		if(stats->curr_read_store[i]) {
			free(stats->read_length_stats[i]);
			uint64_t k;
			for(k=0;k<2;k++) {
				free(stats->indel_length[i*2+k]);
				free(stats->indel_stats[i*2+k]);
			}
			for(k=0;k<stats->curr_read_store[i];k++) free(stats->base_counts_by_cycle[i][k]);
			for(k=i*(MAX_QUAL+1);k<(i+1)*(MAX_QUAL+1);k++) {
				free(stats->mm_stats[k]);
				free(stats->qual_stats[k]);
			}
			free(stats->base_counts_by_cycle[i]);
		}
	}
	if(stats->paired==true) {
		HASH_CLEAR(hh,stats->insert_size);
	}
	free(stats);
}

static void add_indel_stats(as_stats *stats,uint64_t size,uint64_t cycle,uint64_t ix)
{
	stats->indel_stats[ix][cycle]++;
	if(size<=MAX_INDEL_SIZE) {
		if(size>stats->max_indel_length) {
			uint64_t nsize=size*1.2;
			if(nsize>MAX_INDEL_SIZE) nsize=MAX_INDEL_SIZE;
			uint64_t i,j=(stats->paired==true?2:1);
			for(i=0;i<j;i++) {
				int k;
				for(k=0;k<2;k++) {
					stats->indel_length[i*2+k]=as_realloc(stats->indel_length[i*2+k],sizeof(uint64_t)*(nsize+1));
					uint64_t sz;
					for(sz=stats->max_indel_length+1;sz<=nsize;sz++) stats->indel_length[i*2+k][sz]=0;
				}
			}
			stats->max_indel_length=nsize;
		}
		stats->indel_length[ix][size]++;
	}
}

static dist_element *as_find_insert_counter(dist_element **de,int64_t x)
{
	dist_element *new_de;
	HASH_FIND(hh,*de,&x,sizeof(int64_t),new_de);
	if(!new_de) {
		new_de=as_malloc(sizeof(dist_element));
		int i;
		for(i=0;i<4;i++) new_de->ct[i]=0;
		new_de->x=x;
		HASH_ADD(hh,*de,x,sizeof(int64_t),new_de);
	}
	return new_de;
}

static dist_element *as_increase_insert_count(dist_element **de,int ix,int64_t x)
{
	dist_element *new_de=as_find_insert_counter(de,x);
	new_de->ct[ix]++;
	return new_de;
}

#define LH_BIN_SIZE 1024
static pthread_rwlock_t loc_hash_rwlock;

static void insert_loc(as_stats *stats,uint64_t x,int64_t ins_size,uint32_t tile,gt_string *ctg)
{
	loc_hash *lh;
	unsigned int k=x/LH_BIN_SIZE;
	uint16_t loc=x%LH_BIN_SIZE;
	pthread_rwlock_rdlock(&loc_hash_rwlock);
	HASH_FIND_STR(*stats->loc_hash,gt_string_get_string(ctg),lh);
	pthread_rwlock_unlock(&loc_hash_rwlock);
	if(!lh) {
		lh=as_malloc(sizeof(loc_hash));
		lh->ctg=strdup(gt_string_get_string(ctg));
		lh->lblock=0;
		pthread_rwlock_init(&lh->rwlock,NULL);
		pthread_rwlock_wrlock(&loc_hash_rwlock);
		HASH_ADD_KEYPTR(hh,*stats->loc_hash,lh->ctg,(int)strlen(lh->ctg),lh);
		pthread_rwlock_unlock(&loc_hash_rwlock);
	}
	loc_block *lb;
	pthread_rwlock_rdlock(&lh->rwlock);
	HASH_FIND_INT(lh->lblock,&k,lb);
	pthread_rwlock_unlock(&lh->rwlock);
	if(!lb) {
		lb=as_malloc(sizeof(loc_block));
		lb->n_elem=0;
		lb->size=INIT_LB_SIZE;
		lb->x=k;
		lb->elem=as_malloc(lb->size*sizeof(loc_elem));
		pthread_mutex_init(&lb->mutex,NULL);
		pthread_rwlock_wrlock(&lh->rwlock);
		HASH_ADD_INT(lh->lblock,x,lb);
		pthread_rwlock_unlock(&lh->rwlock);
	}
	pthread_mutex_lock(&lb->mutex);
	if(lb->n_elem==lb->size) {
		lb->size*=1.5;
		lb->elem=as_realloc(lb->elem,lb->size*sizeof(loc_elem));
	}
	loc_elem* le=lb->elem+(lb->n_elem++);
	pthread_mutex_unlock(&lb->mutex);
	le->loc=loc;
	le->tile=tile;
	le->dist=ins_size;
}

static void as_stats_resize(as_stats *stats,uint64_t rd,uint64_t l)
{
	stats->max_read_length[rd]=l;
	if(l>=stats->curr_read_store[rd]) {
		uint64_t i,nlen=l*1.5; // Allocate a bit more space than we need now to avoid un-necessary re-sizing in future
		if(stats->curr_read_store[rd]) {
			stats->read_length_stats[rd]=as_realloc(stats->read_length_stats[rd],nlen*sizeof(uint64_t));
			stats->base_counts_by_cycle[rd]=as_realloc(stats->base_counts_by_cycle[rd],sizeof(void *)*nlen);
			uint64_t j;
			for(j=rd*2;j<rd*2+2;j++) stats->indel_stats[j]=as_realloc(stats->indel_stats[j],sizeof(uint64_t)*nlen);
			for(j=rd*(MAX_QUAL+1);j<(rd+1)*(MAX_QUAL+1);j++) {
				stats->mm_stats[j]=as_realloc(stats->mm_stats[j],sizeof(uint64_t)*nlen);
				stats->qual_stats[j]=as_realloc(stats->qual_stats[j],sizeof(uint64_t)*nlen);
			}
			for(i=stats->curr_read_store[rd];i<nlen;i++) {
				stats->read_length_stats[rd][i]=0;
				stats->base_counts_by_cycle[rd][i]=as_calloc((size_t)(MAX_QUAL+1)*5,sizeof(uint64_t));
				for(j=rd*2;j<rd*2+2;j++) stats->indel_stats[j][i]=0;
				for(j=rd*(MAX_QUAL+1);j<(rd+1)*(MAX_QUAL+1);j++) stats->mm_stats[j][i]=0;
			}
		} else {
			stats->read_length_stats[rd]=as_calloc((size_t)nlen,sizeof(uint64_t));
			stats->base_counts_by_cycle[rd]=as_malloc(sizeof(void *)*nlen);
			uint64_t j;
			for(j=rd*2;j<rd*2+2;j++) stats->indel_stats[j]=as_calloc((size_t)nlen,sizeof(uint64_t));
			for(j=rd*(MAX_QUAL+1);j<(rd+1)*(MAX_QUAL+1);j++) {
				stats->mm_stats[j]=as_calloc((size_t)nlen,sizeof(uint64_t));
				stats->qual_stats[j]=as_calloc((size_t)nlen,sizeof(uint64_t));
			}
			for(i=0;i<nlen;i++) {
				stats->base_counts_by_cycle[rd][i]=as_calloc((size_t)(MAX_QUAL+1)*5,sizeof(uint64_t));
			}
		}
		stats->curr_read_store[rd]=nlen;
	}
}

static void get_error_profile(as_stats *stats,gt_alignment *al,uint64_t rd,int qual_offset)
{
	static int mis_type[]={0,2,1,2,2,0,2,1,1,2,0,2,2,1,2,0};
	if(!al->maps->used) return;
	// Get first map only from alignment
  register gt_string* const read = al->read;
  register gt_string* const quals = al->qualities;
  register const bool has_qualities = gt_alignment_has_qualities(al);
	gt_map *map=gt_alignment_get_map(al,0);
  register int quality_misms = 0;
  uint64_t i;
  for(i=0;i<read->length;i++) {
  	if(has_qualities) {
  		quality_misms = gt_string_get_string(quals)[i]-qual_offset;
  		if(quality_misms>MAX_QUAL) quality_misms=MAX_QUAL;
  	} else quality_misms=0;
  	stats->qual_stats[rd*(MAX_QUAL+1)+quality_misms][i]++;
  }
  GT_MAP_ITERATE(map,map_block) {
    GT_MISMS_ITERATE(map_block,misms) {
      if (has_qualities) {
        quality_misms = gt_string_get_string(quals)[misms->position]-qual_offset;
        if(quality_misms>MAX_QUAL) quality_misms=MAX_QUAL;
        else if(quality_misms<0) quality_misms=0;
      }
      switch (misms->misms_type) {
      case MISMS:
       	stats->mm_stats[rd*(MAX_QUAL+1)+quality_misms][misms->position]++;
      	int base=base_tab[(int)gt_string_get_string(read)[misms->position]]-2;
      	int rbase=base_tab[(int)misms->base]-2;
      	if(base>=0 && rbase>=0) {
      		// Get transition/transversion counts
      		int type=mis_type[(base<<2)|rbase];
      		if(type==1) stats->ts_stats[rd][quality_misms]++;
      		else if(type==2) stats->tv_stats[rd][quality_misms]++;
      	}
      	if(base>=0 && misms->position) {
      		int prev_base=base_tab[(int)gt_string_get_string(read)[misms->position-1]]-2;
      		if(prev_base>=0) {
      			if(base==prev_base) stats->pbc_stats[rd][quality_misms]++;
      			else stats->pbn_stats[rd][quality_misms]++;
      		}
      	}
      	break;
      case INS:
      	add_indel_stats(stats,misms->size,misms->position,rd);
      	break;
      case DEL:
      	add_indel_stats(stats,misms->size,misms->position,2+rd);
      	break;
      }
    }

  }
}

static void as_collect_stats(gt_template* template,as_stats* stats,as_param *param,id_tag *idt)
{
	stats->nreads++;
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
		for(i=0;i<k;i++) {
			uint64_t x=gt_alignment_get_counter(al[j],i);
			if(x) {
				nmaps[j]+=x;
				max_dist[j]=i;
			} else if(nmaps[j]) break;
		}
		if(i==k-1 && paired_file==false) ambig[j]=true;
		// Collect error states from first alignment only
		if(nmaps[j]) {
			get_error_profile(stats,al[j],j,qual_offset);
		}
	}
	if(nrd==2) {
		// Paired end alignments
		uint64_t i,k;
		k=gt_template_get_num_counters(template);
		for(i=0;i<k;i++) {
			uint64_t x=gt_template_get_counter(template,i);
			if(x) {
				nmaps[2]+=x;
				max_dist[2]=i;
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
		GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,maps,maps_attr) {
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
			if(nmaps[2]==1 && (nmaps[0]<=gt_alignment_get_num_maps(al[0])) && (nmaps[1]<=gt_alignment_get_num_maps(al[1]))) {
				stats->paired_unique++;
				maps=gt_template_get_mmap_array(template,0,NULL);
				gt_status gt_err;
				int64_t ins_size=gt_template_get_insert_size(maps,&gt_err,0,0);
				if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK) {
					dist_element* de=as_increase_insert_count(&stats->insert_size,AS_INSERT_TYPE_PAIRED,ins_size);
					if(nmaps[0]>1 || nmaps[1]>1) de->ct[AS_INSERT_TYPE_RECOVERED]++;
					if(nsplit[2]) de->ct[AS_INSERT_TYPE_SPLIT]++;
				}
			}
		}
		// Track insert sizes for all pairs where single end reads are uniquely mapping
		if(nmaps[0]==1 && nmaps[1]==1) {
			gt_status gt_err;
			gt_map *tmaps[2];
			tmaps[0]=gt_alignment_get_map(al[0],0);
			tmaps[1]=gt_alignment_get_map(al[1],0);
			uint64_t xx;
			gt_string *contig;
			int64_t ins_size=gt_template_get_insert_size(tmaps,&gt_err,&xx,&contig);
			if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK) {
				(void)as_increase_insert_count(&stats->insert_size,AS_INSERT_TYPE_ALL_UNIQUE,ins_size);
				stats->paired_type[PAIR_TYPE_DS]++;
				insert_loc(stats,xx,ins_size,idt->tile,contig);
			} else if(gt_err==GT_TEMPLATE_INSERT_SIZE_SAME_STRAND) stats->paired_type[PAIR_TYPE_SS]++;
			else if(gt_err==GT_TEMPLATE_INSERT_SIZE_DIFFERENT_CONTIGS) stats->paired_type[PAIR_TYPE_MM]++;
		}
	} else {
		/* We can still track duplicates for single end reads, although we will find too many */
		gt_map* tmap=gt_alignment_get_map(al[0],0);
		insert_loc(stats,tmap->position,0,idt->tile,tmap->seq_name);
	}
}

static int cmp_dist_elem(dist_element *a,dist_element *b)
{
	return a->x-b->x;
}

static int cmp_loc_elem(const void *s1,const void *s2)
{
  const loc_elem *le1,*le2;
  int x;

  le1=s1;
  le2=s2;
  if(!(x=le1->loc-le2->loc)) {
    if(!(x=abs(le1->dist)-abs(le2->dist))) {
      if(!(x=le1->tile-le2->tile)) {
        x=le1->dist-le2->dist;
      }
    }
  }
  return x;
}

static void *as_calc_duplicate_rate(void *ss)
{
	as_param* param=ss;
	as_stats* stats=param->stats[0];
	loc_hash* lh=*(stats->loc_hash);
	uint64_t (*dup_cnt)[DUP_LIMIT+1]=stats->duplicate_counts;
	uint64_t dcounts[2][DUP_LIMIT+1];
	uint64_t tot=0;
	int i,j;
	for(i=0;i<5;i++) for(j=0;j<=DUP_LIMIT;j++) dup_cnt[i][j]=0;
	for(;lh;lh=lh->hh.next) {
		loc_block *lb;
		for(lb=lh->lblock;lb;lb=lb->hh.next) {
			qsort(lb->elem,lb->n_elem,sizeof(loc_elem),cmp_loc_elem);
			u_int16_t tile,loc;
			int16_t dst=0;
			tile=loc=0;
			int k,k1,xx;
			k=k1=xx=0;
			uint64_t kk[4]={0,0,0,0};
			tot+=lb->n_elem;
			loc_elem* le=lb->elem;
			for(i=0;i<(int)lb->n_elem;i++,le++) {
				if(le->loc!=loc || abs(le->dist)!=abs(dst)) {
					if(k) {
						if(k>DUP_LIMIT) k=DUP_LIMIT+1;
						else if(k>1) {
							assert(xx<=DUP_LIMIT);
							dcounts[0][xx]=kk[0];
							dcounts[1][xx++]=kk[1];
							for(k1=0;k1<4;k1++) kk[k1]=0;
							for(k1=0;k1<xx;k1++) {
								int k2;
								for(k2=0;k2<2;k2++) {
									int k3=dcounts[k2][k1];
									if(k3>1) kk[0]+=k3*(k3-1);
								}
								kk[1]+=dcounts[0][k1]*dcounts[1][k1];
								for(k2=0;k2<k1;k2++) {
									kk[2]+=dcounts[0][k1]*dcounts[0][k2]+dcounts[1][k1]*dcounts[1][k2];
									kk[3]+=dcounts[0][k1]*dcounts[1][k2]+dcounts[1][k1]*dcounts[0][k2];
								}
							}
							kk[0]>>=1;
							for(k1=0;k1<4;k1++) dup_cnt[k1+1][k-1]+=kk[k1];
							xx=0;
						}
						dup_cnt[0][k-1]++;
					}
					k=1;
					xx=0;
					tile=le->tile;loc=le->loc;
					dst=le->dist;
					k1=(dst<0)?0:1;
					kk[k1]=1;
					kk[k1^1]=0;
				} else {
					k++;
					if(le->tile!=tile) {
						if(xx<DUP_LIMIT) {
							dcounts[0][xx]=kk[0];
							dcounts[1][xx++]=kk[1];;
							tile=le->tile;
							dst=le->dist;
							k1=(dst<0)?0:1;
							kk[k1]=1;
							kk[k1^1]=0;
						}
					} else {
						if(le->dist!=dst) {
							k1=1;
							dst=le->dist;
						}
						kk[k1]++;
					}
				}
			}
			if(k) {
				if(k>DUP_LIMIT) k=DUP_LIMIT+1;
				else if(k>1) {
					assert(xx<=DUP_LIMIT);
					dcounts[0][xx]=kk[0];
					dcounts[1][xx++]=kk[1];
					for(k1=0;k1<4;k1++) kk[k1]=0;
					for(k1=0;k1<xx;k1++) {
						int k2;
						for(k2=0;k2<2;k2++) {
							int k3=dcounts[k2][k1];
							if(k3>1) kk[0]+=k3*(k3-1);
						}
						kk[1]+=dcounts[0][k1]*dcounts[1][k1];
						for(k2=0;k2<k1;k2++) {
							kk[2]+=dcounts[0][k1]*dcounts[0][k2]+dcounts[1][k1]*dcounts[1][k2];
							kk[3]+=dcounts[0][k1]*dcounts[1][k2]+dcounts[1][k1]*dcounts[0][k2];
						}
					}
					kk[0]>>=1;
					for(k1=0;k1<4;k1++) dup_cnt[k1+1][k-1]+=kk[k1];
					xx=0;
				}
				dup_cnt[0][k-1]++;
			}
		}
	}
	double z1,z2,z3,z4,z5,z6;
	z1=z2=z3=z4=z5=z6=0.0;
  int k=0;
	for(i=0;i<=DUP_LIMIT;i++) {
		z1+=(double)dup_cnt[0][i];
		z2+=(double)dup_cnt[0][i]*(i+1);
		z3+=(double)dup_cnt[1][i];
		z4+=(double)dup_cnt[2][i];
		z5+=(double)dup_cnt[3][i];
		z6+=(double)dup_cnt[4][i];
		if(dup_cnt[0][i]) k=i;
	}
	if(z2 && (z3+z4+z5+z6)) {
		double z,z7,z8,z9;
    z=1.0-z1/z2;
    z1=(z3*z6>0.0)?(z3*(z5+z6)-z5*(z3+z4))/(z3*z6):1.0;
    z2=z3*z1/(z3+z4+z5+z6);
    z7=1.0-(1.0-z2)*z;
    z8=tot*1000.0;
    z9=z8;
    for(i=0;i<10000;i++) {
      z8=tot*z7/(1.0-exp(log(1.0-1.0/z8)*tot));
      if(fabs(z8-z9)<1.0e-2) break;
      z9=z8;
    }
    stats->duplicate_rate[0]=z;
    stats->duplicate_rate[1]=z2;
    stats->duplicate_reads_used=tot;
    stats->unique_fragment_estimate=(uint64_t)(z8+.5);
	} else {
		stats->duplicate_rate[0]=stats->duplicate_rate[1]=0.0;
    stats->duplicate_reads_used=stats->unique_fragment_estimate=0;
	}
	return 0;
}

static void *as_merge_stats(void *ss)
{
	uint64_t i,j,k,k1,nr,len[2],ins_size;

	as_param* param=ss;
	as_stats** st=param->stats;
	uint64_t nt=param->num_threads;
	bool paired=gt_input_generic_parser_attributes_is_paired(param->parser_attr);
	nr=paired?2:1;
	for(j=0;j<nr;j++) {
		len[j]=st[0]->max_read_length[j];
		for(i=1;i<nt;i++) {
			if(st[i]->max_read_length[j]>len[j]) len[j]=st[i]->max_read_length[j];
		}
		if(len[j]>st[0]->max_read_length[j]) as_stats_resize(st[0],j,len[j]);
	}
	ins_size=st[0]->max_indel_length;
	for(i=1;i<nt;i++) if(st[i]->max_indel_length>ins_size) ins_size=st[i]->max_indel_length;
	if(ins_size>st[0]->max_indel_length) {
		for(j=0;j<nr;j++) {
			int k;
			for(k=0;k<2;k++) {
				st[0]->indel_length[j*2+k]=as_realloc(st[0]->indel_length[j*2+k],sizeof(uint64_t)*(ins_size+1));
				uint64_t sz;
				for(sz=st[0]->max_indel_length+1;sz<=ins_size;sz++) st[0]->indel_length[j*2+k][sz]=0;
			}
		}
	}
	for(i=1;i<nt;i++) {
		st[0]->nreads+=st[i]->nreads;
		for(j=0;j<nr;j++) {
			st[0]->yield[j]+=st[i]->yield[j];
			st[0]->mapped[j]+=st[i]->mapped[j];
			st[0]->unique[j]+=st[i]->unique[j];
			st[0]->ambiguous[j]+=st[i]->ambiguous[j];
			for(k=0;k<=MAX_QUAL;k++) {
				st[0]->ts_stats[j][k]+=st[i]->ts_stats[j][k];
				st[0]->tv_stats[j][k]+=st[i]->tv_stats[j][k];
				st[0]->pbc_stats[j][k]+=st[i]->pbc_stats[j][k];
				st[0]->pbn_stats[j][k]+=st[i]->pbn_stats[j][k];
			}
			len[j]=st[i]->max_read_length[j];
			uint64_t **bc0=st[0]->base_counts_by_cycle[j];
			uint64_t **bc1=st[i]->base_counts_by_cycle[j];
			if(st[i]->curr_read_store[j]) {
				if(st[i]->curr_read_store[j]>st[0]->curr_read_store[j])
					as_stats_resize(st[0],j,st[i]->curr_read_store[j]);
				for(k=0;k<=len[j];k++) {
					st[0]->read_length_stats[j][k]+=st[i]->read_length_stats[j][k];
					for(k1=0;k1<5*(MAX_QUAL+1);k1++) bc0[k][k1]+=bc1[k][k1];
				}
				for(k1=j*(MAX_QUAL+1);k1<(j+1)*(MAX_QUAL+1);k1++) {
					for(k=0;k<len[j];k++) {
						st[0]->mm_stats[k1][k]+=st[i]->mm_stats[k1][k];
						st[0]->qual_stats[k1][k]+=st[i]->qual_stats[k1][k];
					}
				}
				for(k1=j*2;k1<j*2+2;k1++) {
					for(k=0;k<len[j];k++) st[0]->indel_stats[k1][k]+=st[i]->indel_stats[k1][k];
					for(k=0;k<=st[i]->max_indel_length;k++) st[0]->indel_length[k1][k]+=st[i]->indel_length[k1][k];
				}
			}
		}
		for(j=0;j<3;j++) {
			st[0]->paired_type[j]+=st[i]->paired_type[j];
			st[0]->bis_stats[j]+=st[i]->bis_stats[j];
			st[0]->reads_with_splitmaps[j]+=st[i]->reads_with_splitmaps[j];
			st[0]->reads_only_with_splitmaps[j]+=st[i]->reads_only_with_splitmaps[j];
		}
		if(paired) {
			st[0]->paired_mapped+=st[i]->paired_mapped;
			st[0]->paired_unique+=st[i]->paired_unique;
			dist_element *de,*de1;
			for(de=st[i]->insert_size;de!=NULL;de=de->hh.next) {
				de1=as_find_insert_counter(&st[0]->insert_size,de->x);
				for(j=0;j<4;j++) de1->ct[j]+=de->ct[j];
			}
		}
		as_stats_free(st[i]);
	}
	HASH_SORT(st[0]->insert_size,cmp_dist_elem);
	return 0;
}

static void as_print_yield_summary(FILE *f,as_param *param)
{
	bool paired=gt_input_generic_parser_attributes_is_paired(param->parser_attr);
	as_stats* st=param->stats[0];
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

static void as_print_mapping_summary(FILE *f,as_param *param)
{
	bool paired=gt_input_generic_parser_attributes_is_paired(param->parser_attr);
	as_stats *st=param->stats[0];
	bool paired_file=false; // Was the input file from a paired mapping
	if(paired==true && !param->input_files[1]) paired_file=true;
	uint64_t counts[4]={0,0,0,0};
	dist_element *de;
	int j;
	for(de=st->insert_size;de;de=de->hh.next) {
		for(j=0;j<4;j++) counts[j]+=de->ct[j];
	}
	double zcounts[4];
	for(j=0;j<4;j++) zcounts[j]=(double)counts[j];
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
		fputs("\nPair statistics (uniquely mapping read pairs only)\n\n",f);
		uint64_t cnt[4]={0,0,0,0};
		double lim[3]={.25,.5,.75};
		int state[4]={0,0,0,0};
		int64_t Q[3][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0}};
		dist_element *de;
		int ct=0;
		for(de=st->insert_size;de && ct<4;de=de->hh.next) {
			int i;
			for(i=0;i<4;i++) {
				if(state[i]<3) {
					if((double)cnt[i]/zcounts[i]>=lim[state[i]]) {
						Q[state[i]++][i]=de->x;
						if(state[i]==3) ct++;
					}
 				}
				cnt[i]+=de->ct[i];
			}
		}
		double ztot=(double)(st->paired_type[0]+st->paired_type[1]+st->paired_type[2]);

		fprintf(f,"Read pairs on different strand (DS):\t%" PRIu64 "\t(%g%%)\n",
					st->paired_type[PAIR_TYPE_DS],100.0*(double)st->paired_type[PAIR_TYPE_DS]/ztot);
		fprintf(f,"Read pairs on same strand (SS):\t%" PRIu64 "\t(%g%%)\n",
					st->paired_type[PAIR_TYPE_SS],100.0*(double)st->paired_type[PAIR_TYPE_SS]/ztot);
		fprintf(f,"Read pairs on different contigs:\t%" PRIu64 "\t(%g%%)\n",
					st->paired_type[PAIR_TYPE_MM],100.0*(double)st->paired_type[PAIR_TYPE_MM]/ztot);
		fputs("\nInsert size summary\n\n",f);
		fprintf(f,"Selected unique read pairs:\t(%g)\tQ1: %" PRId64 "\tMedian: %" PRId64 "\tQ3: %" PRId64 "\n",zcounts[0],Q[0][0],Q[1][0],Q[2][0]);
		fprintf(f,"All unique read pairs:\t(%g)\tQ1: %" PRId64 "\tMedian: %" PRId64 "\tQ3: %" PRId64 "\n",zcounts[1],Q[0][1],Q[1][1],Q[2][1]);
		fprintf(f,"Selected unique read pairs with recovered read:\t(%g)\tQ1: %" PRId64 "\tMedian: %" PRId64 "\tQ3: %" PRId64 "\n",zcounts[2],Q[0][2],Q[1][2],Q[2][2]);
		fprintf(f,"Selected unique read pairs with split reads:\t(%g)\tQ1: %" PRId64 "\tMedian: %" PRId64 "\tQ3: %" PRId64 "\n",zcounts[3],Q[0][3],Q[1][3],Q[2][3]);
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

static void as_print_read_lengths(FILE *f,as_param *param)
{
	as_stats *st=param->stats[0];
	bool paired=gt_input_generic_parser_attributes_is_paired(param->parser_attr);
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

static void as_print_distance_file(as_param *param)
{
	as_stats* st=param->stats[0];
	gt_output_file *file;

	if(param->dist_file) {
		file=gt_output_file_new_compress(param->dist_file,UNSORTED_FILE,param->compress);
	} else {
		file=gt_output_stream_new_compress(stdout,UNSORTED_FILE,param->compress);
	}
	gt_cond_fatal_error(!file,FILE_OPEN,param->dist_file);
	FILE *fp=file->file;

	uint64_t counts[4]={0,0,0,0};
	dist_element *de;
	int j;
	for(de=st->insert_size;de;de=de->hh.next) {
		for(j=0;j<4;j++) counts[j]+=de->ct[j];
	}
	double zcounts[4];
	for(j=0;j<4;j++) zcounts[j]=(double)counts[j];
	fputs("Fragment size distribution (uniquely mapping reads):\n\n",fp);
	fputs("Size\tPaired\tAll\tRecovered\tSplit\tPaired_freq\tAll_freq\tRecovered_freq\tSplit_freq\n",fp);
	for(de=st->insert_size;de;de=de->hh.next) {
		fprintf(fp,"%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%g\t%g\t%g\t%g\n",
				de->x,de->ct[0],de->ct[1],de->ct[2],de->ct[3],
				(double)de->ct[0]/zcounts[0],(double)de->ct[1]/zcounts[1],(double)de->ct[2]/zcounts[2],(double)de->ct[3]/zcounts[3]);
	}
	gt_output_file_close(file);
}

static void as_print_duplicate_summary(FILE *fp,as_param *param)
{
	as_stats* st=param->stats[0];
	fprintf(fp,"\nOverall duplicate percentage = %g%%\nOptical duplicate fraction = %g\nNo. read pairs used = %"PRIu64"\n",100.0*st->duplicate_rate[0],st->duplicate_rate[1],st->duplicate_reads_used);
	fprintf(fp,"Estimated number of unique fragments in library = %"PRIu64"\n",st->unique_fragment_estimate);
}

static void as_print_detailed_duplicate_report(FILE *fp,as_param *param)
{
	as_stats* st=param->stats[0];
	int k=0;
	int i;
	double z=0.0;
	for(i=0;i<=DUP_LIMIT;i++) {
		uint64_t c=st->duplicate_counts[0][i];
		if(c) {
			z+=(double)c;
			k=i;
		}
	}
	if(k) {
		fputs("\nDetailed duplicate report\nN_copies\tfreq\tprob\tW++\tW+-\tB++\tB+-\n",fp);
		for(i=0;i<=k;i++) {
			if(i==DUP_LIMIT) fputs(">=",fp);
			fprintf(fp,"%d\t%"PRIu64"\t%g\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",i+1,st->duplicate_counts[0][i],(double)st->duplicate_counts[0][i]/z,
					st->duplicate_counts[1][i],st->duplicate_counts[2][i],st->duplicate_counts[3][i],st->duplicate_counts[4][i]);
		}
	}
}

static void as_print_mismatch_report(FILE *fp,as_param *param)
{
  const double ln10=log(10.0);
	as_stats* st=param->stats[0];
	bool paired=gt_input_generic_parser_attributes_is_paired(param->parser_attr);
	// Collect overall error stats
	uint64_t mm_stats[2][MAX_QUAL+1],qual_stats[2][MAX_QUAL+1];
	uint64_t mm_total[2]={0,0};
	uint64_t qtotal[2]={0,0},total[2]={0,0},tv_stats[2]={0,0},ts_stats[2]={0,0},pbc_stats[2]={0,0},pbn_stats[2]={0,0};
	uint64_t base_ct[2][5],base_ct1[2][5*(MAX_QUAL+1)],base_ctt[2][MAX_QUAL+1];
	uint64_t i,j,k,nr;
	nr=paired?2:1;
	for(i=0;i<nr;i++) {
		for(j=0;j<5;j++) base_ct[i][j]=0;
		for(j=0;j<=MAX_QUAL;j++) {
			uint64_t *tp=st->mm_stats[i*(MAX_QUAL+1)+j];
			uint64_t *tq=st->qual_stats[i*(MAX_QUAL+1)+j];
			uint64_t tmp_mm=0,tmp=0;
			uint64_t tt[]={0,0,0,0,0};
			for(k=0;k<st->max_read_length[i];k++) {
				uint64_t *tb=st->base_counts_by_cycle[i][k]+j*5;
				int k1;
				for(k1=0;k1<5;k1++) tt[k1]+=tb[k1];
				tmp_mm+=tp[k];
				tmp+=tq[k];
			}
			uint64_t tt1=0;
			for(k=0;k<5;k++) {
				base_ct[i][k]+=tt[k];
				base_ct1[i][j*5+k]=tt[k];
				tt1+=tt[k];
			}
			base_ctt[i][j]=tt1;
			tv_stats[i]+=st->tv_stats[i][j];
			ts_stats[i]+=st->ts_stats[i][j];
			pbc_stats[i]+=st->pbc_stats[i][j];
			pbn_stats[i]+=st->pbn_stats[i][j];
			mm_stats[i][j]=tmp_mm;
			qual_stats[i][j]=tmp;
			mm_total[i]+=tmp_mm;
			qtotal[i]+=tmp;
		}
		for(j=0;j<5;j++) total[i]+=base_ct[i][j];
	}
	fputs("\nMismatch report (based on first alignment only)\n\n",fp);
	if(paired) {
		fprintf(fp,"Overall Mismatch percentage (Read 1 Read 2):\t%g%%\t%g%%\n",100.0*mm_total[0]/qtotal[0],100.0*mm_total[1]/qtotal[1]);
		fprintf(fp,"Overall base composition:\t(A:%.3f,C:%.3f,G:%.3f,T:%.3f,N:%.3f)\t(A:%.3f,C:%.3f,G:%.3f,T:%.3f,N:%.3f)\n",
				(double)base_ct[0][1]/total[0],(double)base_ct[0][2]/total[0],(double)base_ct[0][3]/total[0],(double)base_ct[0][4]/total[0],(double)base_ct[0][0]/total[0],
				(double)base_ct[1][1]/total[1],(double)base_ct[1][2]/total[1],(double)base_ct[1][3]/total[1],(double)base_ct[1][4]/total[1],(double)base_ct[1][0]/total[1]);
		fprintf(fp,"Overall transition:transversion ratio:\t%g\t%g\n",(double)ts_stats[0]/tv_stats[0],(double)ts_stats[1]/tv_stats[1]);
		fprintf(fp,"Overall probability of mismatch being a copy of previous base:\t%g\t%g\n",(double)pbc_stats[0]/(pbc_stats[0]+pbn_stats[0]),(double)pbc_stats[1]/(pbc_stats[1]+pbn_stats[1]));
	} else {
		fprintf(fp,"Overall mismatch percentage:\t%g%%\n",100.0*mm_total[0]/qtotal[0]);
		fprintf(fp,"Overall base composition:\t(A:%.3f,C:%.3f,G:%.3f,T:%.3f,N:%.3f)\n",
				(double)base_ct[0][1]/total[0],(double)base_ct[0][2]/total[0],(double)base_ct[0][3]/total[0],(double)base_ct[0][4]/total[0],(double)base_ct[0][0]/total[0]);
		fprintf(fp,"Overall transition:transversion ratio:\t%g\n",(double)ts_stats[0]/tv_stats[0]);
		fprintf(fp,"Overall probability of mismatch being a copy of previous base:\t%g\n",(double)pbc_stats[0]/(pbc_stats[0]+pbn_stats[0]));
	}
	for(i=0;i<nr;i++) {
		fprintf(fp,"\nMismatch quality profile - Read %"PRIu64"\n\n",i+1);
		fputs("Qual\tn_bases\tp(bases)\tcp(bases)\tn_mm\tp(mm)\t-log10_p(mm)\tts:tv\tp(pbc)\tp(A)\tp(C)\tp(G)\tp(T)\n",fp);
		uint64_t ttot=0;
		for(j=0;j<=MAX_QUAL;j++) {
			if(qual_stats[i][j]) {
				double z=(double)mm_stats[i][j]/qual_stats[i][j];
				ttot+=qual_stats[i][j];
				uint64_t tt=st->pbc_stats[i][j]+st->pbn_stats[i][j];
				fprintf(fp,"%"PRIu64"\t%"PRIu64"\t%g\t%g\t%"PRIu64"\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",j,qual_stats[i][j],(double)qual_stats[i][j]/qtotal[i],
						(double)ttot/qtotal[i],mm_stats[i][j],z,-log(z)/ln10,st->tv_stats[i][j]?(double)st->ts_stats[i][j]/st->tv_stats[i][j]:0.0,
								tt?(double)st->pbc_stats[i][j]/tt:0.0,(double)base_ct1[i][j*5+1]/base_ctt[i][j],(double)base_ct1[i][j*5+2]/base_ctt[i][j],
										(double)base_ct1[i][j*5+3]/base_ctt[i][j],(double)base_ct1[i][j*5+4]/base_ctt[i][j]);

			}
		}
		fprintf(fp,"\nMismatch read position profile - Read %"PRIu64"\n\n",i+1);
		fputs("Pos\tn_bases\tp(bases)\tavg_qual\tn_mm\tp(mm)\tp(A)\tp(C)\tp(G)\tp(T)\tp(N)\n",fp);
		uint64_t len=st->max_read_length[i];
		for(j=0;j<len;j++) {
			uint64_t mms[5]={0,0,0,0,0};
			uint64_t *bc=st->base_counts_by_cycle[i][j];
			uint64_t qs=0,mm=0;
			double qsmn=0.0;
			for(k=0;k<=MAX_QUAL;k++) {
				mm+=st->mm_stats[i*(MAX_QUAL+1)+k][j];
				qs+=st->qual_stats[i*(MAX_QUAL+1)+k][j];
				qsmn+=(double)k*st->qual_stats[i*(MAX_QUAL+1)+k][j];
				int k1;
				for(k1=0;k1<5;k1++) mms[k1]+=bc[k*5+k1];
			}
			if(qs) {
				qsmn/=(double)qs;
				uint64_t tt=0;
				for(k=0;k<5;k++) tt+=mms[k];
				fprintf(fp,"%"PRIu64"\t%"PRIu64"\t%g\t%.1f\t%"PRIu64"\t%g\t%g\t%g\t%g\t%g\t%g\n",j,qs,(double)qs/qtotal[i],qsmn,mm,(double)mm/qs,
						(double)mms[1]/tt,(double)mms[2]/tt,(double)mms[3]/tt,(double)mms[4]/tt,(double)mms[0]/tt);
			}
		}
	}
}

static void as_print_stats(as_param *param)
{
	gt_output_file *file;
	if(param->output_file) {
		file=gt_output_file_new_compress(param->output_file,UNSORTED_FILE,param->compress);
	} else {
		file=gt_output_stream_new_compress(stdout,UNSORTED_FILE,param->compress);
	}
	gt_cond_fatal_error(!file,FILE_OPEN,param->output_file);
	FILE *fp=file->file;
	if(gt_input_generic_parser_attributes_is_paired(param->parser_attr)) as_print_distance_file(param);
	as_print_yield_summary(fp,param);
	as_print_mapping_summary(fp,param);
	as_print_duplicate_summary(fp,param);
	as_print_mismatch_report(fp,param);
	as_print_read_lengths(fp,param);
	as_print_detailed_duplicate_report(fp,param);
	gt_output_file_close(file);
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
			{"phage_lambda",required_argument,0,'P'},
			{"phix174",required_argument,0,'X'},
			{"paired",no_argument,0,'p'},
			{"variable",no_argument,0,'V'},
			{"ignore_id",no_argument,0,'i'},
			{"mmap",no_argument,0,'w'},
			{"fastq",no_argument,0,'F'},
			{"solexa",no_argument,0,'S'},
			{"gzip",no_argument,0,'z'},
			{"bzip2",no_argument,0,'j'},
			{"no-compress",no_argument,0,'Z'},
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
			.phage_lambda=NULL,
			.phix174=NULL,
			.mmap_input=false,
			.parser_attr=gt_input_generic_parser_attributes_new(false),
			.ignore_id=false,
			.compress=NONE,
			.min_insert=0,
			.max_insert=DEFAULT_MAX_INSERT,
			.variable_read_length=false,
			.read_length={0,0},
			.max_read_length=MAX_READ_LENGTH,
			.num_threads=1,
			.qual_offset=DEFAULT_QUAL_OFFSET,
	};

	while(!err && (c=getopt_long(argc,argv,"d:t:r:o:q:m:M:l:L:x:P:X:FSVzjZwpi?",longopts,0))!=-1) {
		switch(c) {
		case 'd':
			set_opt("insert_dist",&param.dist_file,optarg);
			break;
		case 'p':
			gt_input_generic_parser_attributes_set_paired(param.parser_attr,true);
			break;
		case 'o':
			set_opt("output",&param.output_file,optarg);
			break;
		case 'P':
			set_opt("phage_lambda",&param.phage_lambda,optarg);
			break;
		case 'X':
			set_opt("phix174",&param.phix174,optarg);
			break;
		case 'L':
			param.max_read_length=(uint64_t)strtoul(optarg,&p,10);
			break;
		case 'l':
			param.read_length[0]=(uint64_t)strtoul(optarg,&p,10);
			if(*p==',') param.read_length[1]=(uint64_t)strtoul(p+1,&p1,10);
			else param.read_length[1]=param.read_length[0];
			break;
		case 'z':
			param.compress=GZIP;
			break;
		case 'j':
			param.compress=BZIP2;
			break;
		case 'Z':
			param.compress=NONE;
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
		case 'w':
			param.mmap_input=true;
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
	if(!param.phage_lambda) param.phage_lambda=strdup(PHAGE_LAMBDA);
	if(!param.phix174) param.phix174=strdup(PHIX174);
	as_set_output_files(&param);
	as_stats** stats=as_malloc(param.num_threads*sizeof(void *));
	param.stats=stats;
	loc_hash *lh=0;
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
			id_tag *idt=new_id_tag();
			stats[tid]=as_stats_new(gt_input_generic_parser_attributes_is_paired(param.parser_attr));
			stats[tid]->loc_hash=&lh;
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
				if(!param.ignore_id) {
					uint64_t idt_err=parse_id_tag(template->tag,idt);
					if(idt_err!=ID_TAG_OK) {
						gt_error_msg("Fatal error parsing ID '"PRIgts"'\n",PRIgts_content(template->tag));
						break;
					}
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
				as_collect_stats(template,stats[tid],&param,idt);
			}
			gt_template_delete(template);
			gt_buffered_input_file_close(buffered_input1);
			gt_buffered_input_file_close(buffered_input2);
			free_id_tag(idt);
		}
		gt_input_file_close(input_file1);
		gt_input_file_close(input_file2);
	} else { // Single file (could be single or paired end)
		gt_input_file* input_file=param.input_files[0]?gt_input_file_open(param.input_files[0],param.mmap_input):gt_input_stream_open(stdin);
#ifdef HAVE_OPENMP
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
			stats[tid]=as_stats_new(gt_input_generic_parser_attributes_is_paired(param.parser_attr));
			stats[tid]->loc_hash=&lh;
			id_tag *idt=new_id_tag();
			while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,param.parser_attr))) {
				if (error_code!=GT_IMP_OK) {
					gt_error_msg("Error parsing file '%s'\n",param.input_files[0]);
					continue;
				}
				// For paired reads, insert single end mappings into alignments
				if(gt_input_generic_parser_attributes_is_paired(param.parser_attr)) {
					if(gt_template_get_num_blocks(template)!=2) {
						gt_fatal_error_msg("Fatal error: Expecting paired reads\n");
					}
					gt_alignment *al[2];
					al[0]=gt_template_get_block(template,0);
					al[1]=gt_template_get_block(template,1);
					gt_alignment_recalculate_counters(al[0]);
					gt_alignment_recalculate_counters(al[1]);
				}
				if(!param.ignore_id) {
					uint64_t idt_err=parse_id_tag(template->tag,idt);
					if(idt_err!=ID_TAG_OK) {
						gt_error_msg("Fatal error parsing ID '"PRIgts"'\n",PRIgts_content(template->tag));
						break;
					}
				}
				as_collect_stats(template,stats[tid],&param,idt);
			}
			// Clean
			gt_template_delete(template);
			gt_buffered_input_file_close(buffered_input);
			free_id_tag(idt);
		}
		gt_input_file_close(input_file);
	}
	pthread_t stats_merge;
	if(pthread_create(&stats_merge,NULL,as_merge_stats,&param)) {
		gt_error_msg("Fatal error - could not create new thread\n");
		exit(-1);
	}
	pthread_t calc_dup;
	if(pthread_create(&calc_dup,NULL,as_calc_duplicate_rate,&param)) {
		gt_error_msg("Fatal error - could not create new thread\n");
		exit(-1);
	}
	pthread_join(calc_dup,NULL);
	pthread_join(stats_merge,NULL);
	as_print_stats(&param);
	as_stats_free(stats[0]);
	free(stats);
	return err;
}
