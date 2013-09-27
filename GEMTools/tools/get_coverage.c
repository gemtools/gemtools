/*
 * get_coverage.c
 *
 *  Created on: 25 Sep 2013
 *      Author: heath
 */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "gem_tools.h"
#include "uthash.h"
#include "get_coverage.h"

static char **input_files;
static int input_idx,n_input_files;
static pthread_mutex_t input_file_mutex;

static char *get_input_file(void)
{
  char *p;

  pthread_mutex_lock(&input_file_mutex);
  if(input_idx<n_input_files) {
    p=input_files[input_idx++];
  } else p=0;
  pthread_mutex_unlock(&input_file_mutex);
  return p;
}

static void set_opt(char *opt,char **opt_p,char *val) 
{
  if(*opt_p) {
    fprintf(stderr,"multiple %s options: '%s' overwriting previous definition '%s'\n",opt,val,*opt_p);
    free(*opt_p);
  }
  *opt_p=strdup(val);
}

static int read_ranges_file(char *fname,struct contig **contig_ptr,struct lk_compress *lkc)
{
  int i,err;
  u_int32_t x1,x1a,x2,x3,x;
  struct contig *ctg,*c;
  struct range_blk *r,*r1;
  char *estr;
  FILE *fptr;
  string *s;
  tokens *tok;
  void *tbuf;

  s=0;
  tok=0;
  tbuf=0;
  err=0;
  fptr=open_readfile_and_check(fname,&i,lkc);
  printf("Reading ranges file '%s'\n",fname);
  ctg=0;
  while(!err) {
    s=fget_string(fptr,s,&tbuf);
    if(!s->len) break;
    tok=tokenize(get_cstring(s),'\t',tok);
    if(tok->n_tok>=3) {
      HASH_FIND_STR(ctg,tok->toks[0],c);
      if(!c) {
	c=lk_malloc(sizeof(struct contig));
	c->name=strdup(tok->toks[0]);
	c->size=0;
	c->ranges=lk_malloc(sizeof(struct range_blk));
	c->ranges->idx=0;
	c->ranges->next=0;
	c->tranges=0;
	c->counts=0;
	c->tcounts=0;
	c->tot_count=0;
	c->tot_tcount=0;
	c->tot_bases=0;
	c->tot_tbases=0;
	HASH_ADD_KEYPTR(hh,ctg,c->name,strlen(c->name),c);
      }
      r=c->ranges;
      if(r->idx==RANGE_BLK_SIZE) {
	r=lk_malloc(sizeof(struct range_blk));
	r->idx=0;
	r->next=c->ranges;
	c->ranges=r;
      }
      x1=(u_int32_t)getlnumcolumn(tok,1,0,&estr);
      if(!estr) x2=(u_int32_t)getlnumcolumn(tok,2,0,&estr);
      if(estr || x1<1 || x2<x1) {
	err=1;
	fprintf(stderr,"Error in coordinate format: %s\n",estr);
      } else {
	r->name[r->idx]=getstrcolumn(tok,3,0,0);
	r->x1[r->idx]=x1;
	r->x2[r->idx++]=x2;
      }
    }
  }
  fclose(fptr);
  if(s) free_string(s);
  if(tok) free_tokens(tok);
  if(tbuf) free_fget_buffer(&tbuf);
  signal(SIGCHLD,SIG_DFL);
  x2=x3=0;
  for(c=ctg;c && !err;c=c->hh.next) {
    x1=0;
    r=c->ranges;
    c->ranges=0;
    while(r) {
      for(i=0;i<r->idx;i++) {
	if(r->x2[i]>x1) x1=r->x2[i];
      }
      r1=r->next;
      r->next=c->ranges;
      c->ranges=r;
      r=r1;
    }
    c->size=x1;
    x1a=0;
    if(!x1) err=2;
    else {
      x2+=x1;
      c->counts=lk_malloc(sizeof(count)*c->size);
      for(i=0;i<(int)c->size;i++) c->counts[i]=NO_COUNT;
      for(r=c->ranges;r;r=r->next) {
	for(i=0;i<r->idx;i++) {
	  for(x=r->x1[i]-1;x<r->x2[i];x++) c->counts[x]=0;
	}
      }
      for(i=0;i<(int)c->size;i++) if(!c->counts[i]) x1a++;
    }
    c->tot_bases=x1a;
    x3+=x1a;
  }
  printf("total=%u, in ranges=%u\n",x2,x3);
  while(waitpid(-1,&i,WNOHANG)>0);
  *contig_ptr=err?0:ctg;
  return err;
}

static int read_target_file(char *fname,struct contig *ctg,int extend,struct lk_compress *lkc)
{
  int i,err;
  u_int32_t x1,x2,x3,x,cs,x1a;
  struct contig *c;
  struct range_blk *r;
  char *estr;
  FILE *fptr;
  string *s;
  tokens *tok;
  void *tbuf;

  s=0;
  tok=0;
  tbuf=0;
  err=0;
  fptr=open_readfile_and_check(fname,&i,lkc);
  printf("Reading target file '%s'\n",fname);
  while(!err) {
    s=fget_string(fptr,s,&tbuf);
    if(!s->len) break;
    tok=tokenize(get_cstring(s),'\t',tok);
    if(tok->n_tok>=3) {
      HASH_FIND_STR(ctg,tok->toks[0],c);
      if(c) {
	r=c->tranges;
	if(!r) {
	  r=c->tranges=lk_malloc(sizeof(struct range_blk));
	  r->idx=0;
	} else {
	  if(r->idx==RANGE_BLK_SIZE) {
	    r=lk_malloc(sizeof(struct range_blk));
	    r->idx=0;
	    r->next=c->tranges;
	    c->tranges=r;
	  }
	}
	x1=(u_int32_t)getlnumcolumn(tok,1,0,&estr);
	if(!estr) x2=(u_int32_t)getlnumcolumn(tok,2,0,&estr);
	if(estr || x1<1 || x2<x1) {
	  err=1;
	  fprintf(stderr,"Error in coordinate format: %s\n",estr);
	} else {
	  if(x1<=(unsigned int)extend) x1=1;
	  else x1-=extend;
	  x2+=extend;
	  if(x1<1) x1=1;
	  r->name[r->idx]=getstrcolumn(tok,3,0,0);
	  r->x1[r->idx]=x1;
	  r->x2[r->idx++]=x2;
	}
      } else {
	printf("Not found %s\n",tok->toks[0]);
      }
    }
  }
  fclose(fptr);
  if(s) free_string(s);
  if(tok) free_tokens(tok);
  if(tbuf) free_fget_buffer(&tbuf);
  signal(SIGCHLD,SIG_DFL);
  x2=x3=0;
  for(c=ctg;c && !err;c=c->hh.next) {
    if(!c->tranges) continue;
    x1=0;
    for(r=c->tranges;r;r=r->next) {
      for(i=0;i<r->idx;i++) {
	if(r->x2[i]>x1) x1=r->x2[i];
      }
    }
    if(x1>c->size) x1=c->size;
    c->tsize=x1;
    x1a=0;
    if(!x1) {
      fprintf(stderr,"Illegal target range for %s\n",c->name);
      err=2;
      continue;
    } else {
      x2+=x1;
      /*      printf("%s\t%d\n",c->name,c->tsize);*/
      c->tcounts=lk_malloc(sizeof(count)*c->tsize);
      cs=c->tsize;
      for(i=0;i<(int)cs;i++) c->tcounts[i]=NO_COUNT;
      for(r=c->tranges;r;r=r->next) {
	for(i=0;i<r->idx;i++) {
	  for(x=r->x1[i]-1;x<r->x2[i];x++) {
	    if(x>=cs) break;
	    if(!c->counts[x]) c->tcounts[x]=0;
	  }
	}
      }
      for(i=0;i<(int)cs;i++) if(!c->tcounts[i]) x1a++;
    }
    c->tot_tbases=x1a;
    x3+=x1a;
  }
  printf("total=%u, on target=%u\n",x2,x3);
  while(waitpid(-1,&i,WNOHANG)>0);
  return err;
}

static char *parse_match(char *p,int fmt,struct match *m)
{
  char *p1;
	
  m->orientation=-1;
  if(fmt==GEM_FMT) {
    if(*p=='[') {
      while(*p && *p!='=') p++;
      if(*p) {
	p++;
      }
    }
    p1=p;
    while(*p1 && *p1!=':') {
      if(*p1=='#') {
	if(!strncmp(p1,"#C2T:",4) || !strncmp(p1,"#G2A:",4)) {
	  *p1=0;
	  p1+=4;
	  break;
	}
      } else if(*p1=='_') {
	if(!strncmp(p1,"_C2T:",4) || !strncmp(p1,"_G2A:",4)) {
	  *p1=0;
	  p1+=4;
	  break;
	}
      }
      p1++;
    }
    if(*p1) {
      *p1++=0;
      m->ctg=p;
      if(isdigit((int)*p1)) {
	m->pos=strtoul(p1,&p,10);
	if(p==p1) p=0;
	else if(*p=='F') m->orientation=0;
	else if(*p=='R') m->orientation=1;
	else p=0;
	m->cigar=0;
      } else {
	if(*p1=='F') m->orientation=0;
	else if(*p1=='R') m->orientation=1;
	else p=0;
	if(p) {
	  p1++;
	  if(*p1=='[') p1++;
	  m->pos=strtoul(p1,&p,10);
	  if(p==p1) p=0;
	  else {
	    m->cigar=p;
	    while(*p && *p!='@') p++;
	    if(*p) *p=0;
	  }
	}
      }
    } else p=0;
  } else {
    p1=strchr(p,'/');
    if(p1) p=p1+1;
    m->ctg=p;
    while(*p && *p!=':') p++;
    if(*p) {
      *p++=0;
      p1=p;
      m->pos=strtoul(p1,&p,10);
      if(p==p1) p=0;
      else if(*p=='F') m->orientation=0;
      else if(*p=='R') m->orientation=1;
      else p=0;
      m->cigar=0;
    } else p=0;
  }
  return p;
}

static int process_file(char *fname,struct contig *ctg,int block_size,int fmt,u_int64_t *number,struct lk_compress *lkc)
{
  int i,k,err,sz,nn;
  FILE *fptr;
  string *s;
  tokens *tok;
  void *tbuf;
  char *p,*p1;
  struct match m;
  struct contig *c;
  count *ct;
  u_int32_t x,x1;

  s=0;
  tok=0;
  tbuf=0;
  err=0;
  nn=(fmt==GEM_FMT?1:0);
  fptr=open_readfile(fname,&i,lkc);
  if(!fptr) return -1;
  printf("Reading matches file '%s'\n",fname);
  while(!err) {
    s=fget_string(fptr,s,&tbuf);
    if(!s->len) break;
    tok=tokenize(get_cstring(s),'\t',tok);
    if(tok->n_tok>=4+nn) {
      p=tok->toks[2+nn];
      i=0;
      if(!isdigit(*p)) continue;
      while(*p && !err) {
	k=(int)strtol(p,&p1,10);
	i+=k;
	if(i>1)  break;
	if(*p1) {
	  if(*p1!=':' && *p1!='+') err=5;
	  else p1++;
	}
	p=p1;
	if(!k && i) break;
      }
      if(err) {
	fprintf(stderr,"Bad format for match number column '%s'\n",tok->toks[3]);
	err=0;
	continue;
      } else if(i==1) {
	p=parse_match(tok->toks[3+nn],fmt,&m);
	if(!p) {
	  /*	  err=10;*/
	  fprintf(stderr,"Bad format for match\n");
	  continue;
	} else {
	  if(*number) {
	    (*number)--;
	    if(!(*number)) {
	      err=1;
	      break;
	    }
	  }
	  HASH_FIND_STR(ctg,m.ctg,c);
	  if(c) {
	    sz=strlen(tok->toks[1]);
	    ct=c->counts;
	    x=m.pos-1;
	    if(block_size>1) {
	      if(m.orientation) {
		for(i=0;i<sz;i++) {
		  if(x<c->size) {
		    x1=x/block_size;
		    x1*=block_size;
		    if(ct[x1]<MAX_COUNT) ct[x1]++;
		  }
		  if(!x) break;
		  x--;
		}
	      } else {
		for(i=0;i<sz && x<c->size;i++) {
		  x1=x/block_size;
		  x1*=block_size;
		  if(ct[x1]<MAX_COUNT) ct[x1]++;
		  x++;
		}
	      }
	    } else {
	      if(m.orientation) {
		for(i=0;i<sz;i++) {
		  if(x<c->size) {
		    if(ct[x]<MAX_COUNT) ct[x]++;
		  }
		  if(!x) break;
		  x--;
		}
	      } else {
		for(i=0;i<sz && x<c->size;i++) {
		  if(ct[x]<MAX_COUNT) ct[x]++;
		  x++;
		}
	      }
	    }
	    if((ct=c->tcounts)) {
	      x=m.pos-1;
	      if(block_size>1) {
		if(m.orientation) {
		  for(i=0;i<sz;i++) {
		    if(x<c->tsize) {
		      x1=x/block_size;
		      x1*=block_size;
		      if(ct[x1]<MAX_COUNT) ct[x1]++;
		    }
		    if(!x) break;
		    x--;
		  }
		} else {
		  for(i=0;i<sz && x<c->tsize;i++) {
		    x1=x/block_size;
		    x1*=block_size;
		    if(ct[x1]<MAX_COUNT) ct[x1]++;
		    x++;
		  }
		}
	      } else {
		if(m.orientation) {
		  for(i=0;i<sz;i++) {
		    if(x<c->tsize) {
		      if(ct[x]<MAX_COUNT) ct[x]++;
		    }
		    if(!x) break;
		    x--;
		  }
		} else {
		  for(i=0;i<sz && x<c->tsize;i++) {
		    if(ct[x]<MAX_COUNT) ct[x]++;
		    x++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  fclose(fptr);
  if(s) free_string(s);
  if(tok) free_tokens(tok);
  if(tbuf) free_fget_buffer(&tbuf);
  signal(SIGCHLD,SIG_DFL);
  while(waitpid(-1,&i,WNOHANG)>0);
  return err;
}

struct tdc_par {
  struct contig *ctg;
  struct lk_compress *lkc;
};

static void *read_det_cov(void *vd)
{
  int i,err;
  FILE *fptr;
  string *s;
  tokens *tok;
  void *tbuf;
  char *p,*p1,cc,*fname;
  struct contig *c,*ctg;
  struct tdc_par *tp;
  count *ct,ct1,ct2;
  struct lk_compress *lkc;
  u_int32_t x;

  tp=vd;
  ctg=tp->ctg;
  lkc=tp->lkc;
  while((fname=get_input_file())) {
    s=0;
    tok=0;
    tbuf=0;
    err=0;
    fptr=open_readfile(fname,&i,lkc);
    if(!fptr) return 0;
    printf("Reading detailed coverage file '%s'\n",fname);
    c=0;
    ct=0;
    x=0;
    while(!err) {
      s=fget_string(fptr,s,&tbuf);
      if(!s->len) break;
      tok=tokenize(get_cstring(s),'\t',tok);
      if(tok->n_tok==2) {
	p=tok->toks[0];
	if(*p=='>') {
	  if(!c || strcmp(c->name,p+1)) {
	    HASH_FIND_STR(ctg,p+1,c);
	  }
	  if(c) {
	    p=tok->toks[1];
	    x=strtoul(p,&p1,10)-1;
	    if(p==p1 || *p1) {
	      err=10;
	      fprintf(stderr,"(a) Bad file format\n");
	    }
	    ct=c->counts;
	  }
	}
      } else if(c) {
	if(tok->n_tok!=1) {
	  fprintf(stderr,"(b) Bad file format\n");
	  err=-1;
	}
	p=tok->toks[0];
	while(x<c->size && (cc=*p++)) {
	  if(cc>=33 && cc<=125) {
	    cc-=33;
	    ct1=(count)cc;
	  } else if(cc==126) {
	    ct1=(count)strtoul(p,&p1,16);
	    p=p1+1;
	  } else {
	    fprintf(stderr,"(c) Bad file format (%d) '%s'\n",cc,tok->toks[0]);
	    err=-1;
	  }
	  if(!err) {
	    pthread_mutex_lock(&c->mut);
	    ct2=ct[x];
	    if(ct2<MAX_COUNT) {
	      if(MAX_COUNT-ct2<ct1) ct[x]=MAX_COUNT;
	      else ct[x]=ct2+ct1;
	    }
	    pthread_mutex_unlock(&c->mut);
	    x++;
	  }
	}
      }
    }
    fclose(fptr);
    if(s) free_string(s);
    if(tok) free_tokens(tok);
    if(tbuf) free_fget_buffer(&tbuf);
    signal(SIGCHLD,SIG_DFL);
    while(waitpid(-1,&i,WNOHANG)>0);
  }
  return 0;
}

#define LN_92 4.52178857704904030964
static int write_cbuf(count *cbuf,int out_cnt,FILE *fptr)
{
  int i,j,k,k1;
  count c,c1,x,xx[3]={1,92,92*92};
  char buf[OUT_WIDTH+3];

  for(i=j=0;i<out_cnt;i++) if(cbuf[i]>j) j=cbuf[i];
  if(j) j=(int)(log((double)j)/LN_92+.999999);
  if(!j) j=1;
  assert(j>0 && j<4);
  buf[0]=32+j;
  for(i=k=1;i<out_cnt;i++) {
    c=cbuf[i];
    for(k1=j-1;k1>=0;k1--) {
      x=xx[k1];
      c1=c/x;
      buf[k++]=c1+33;
      c-=c1*x;
    }
    if(k+j>OUT_WIDTH) break;
  }
  buf[k++]='\n';
  buf[k]=0;
  fputs(buf,fptr);
  for(j=0;i<out_cnt;i++) cbuf[j++]=cbuf[i];
  return j;
}

int main(int argc,char *argv[])
{
  int i,j,k,c,err,nthr;
  int block_size,extend,detail,combine,out_cnt,format;
  struct lk_compress *lkc;
  struct contig *contigs,*ctg;
  struct range_blk *r;
  struct tdc_par tp;
  char *ranges_file,*target_file,*output_file,*compare_ref;
  char *filter,*suffix,*tn,*prefix,*pp;
  count *ct,cc,cbuf[OUT_WIDTH];
  u_int64_t *hist,*thist,nn,kk,tnn,tkk,number;
  u_int32_t x;
  pthread_t *read_threads;
  FILE *ofptr,*ctg_fptr;
  static struct option longopts[]={
    {"ranges_file",required_argument,0,'r'},
    {"target_file",required_argument,0,'t'},
    {"block_size",required_argument,0,'b'},
    {"number",required_argument,0,'n'},
    {"output",required_argument,0,'o'},
    {"extend_regions",required_argument,0,'x'},
    {"detailed_output",no_argument,0,'d'},
    {"combine",no_argument,0,'c'},
    {"compare",required_argument,0,'C'},
    {"eland",no_argument,0,'E'},
    {"gem",no_argument,0,'G'},
    {"prefix",no_argument,0,'p'},
    {0,0,0,0}
  };
  
  err=0;
  detail=combine=0;
  ranges_file=target_file=output_file=compare_ref=prefix=0;
  block_size=1;
  extend=0;
  contigs=0;
  number=0;
  format=GEM_FMT;
  while((c=getopt_long(argc,argv,"p:r:t:b:o:n:x:C:dcEG",longopts,0))!=-1) {
    switch(c) {
    case 'p':
      set_opt("prefix",&prefix,optarg);
      break;
    case 'r':
      set_opt("ranges_file",&ranges_file,optarg);
      break;
    case 't':
      set_opt("target_file",&target_file,optarg);
      break;
    case 'C':
      set_opt("compare",&compare_ref,optarg);
      break;
    case 'o':
      set_opt("output",&output_file,optarg);
      break;
    case 'b':
      block_size=atoi(optarg);
      break;
    case 'E':
      format=ELAND_FMT;
      break;
    case 'G':
      format=GEM_FMT;
      break;
    case 'x':
      extend=atoi(optarg);
      break;
    case 'n':
      number=strtoul(optarg,&pp,10);
      break;
    case 'd':
      detail=1;
      break;
    case 'c':
      combine=1;
      break;
    }
  }
  if(!ranges_file) {
    fprintf(stderr,"No ranges file specified.  Use -r or --ranges_file option\n");
    exit(-1);
  }
  if(!output_file) output_file=strdup("coverage_hist.txt");
  if(block_size<1) {
    fprintf(stderr,"Invalid block size specified - must be >0\n");
    exit(-1);
  }
  printf("Block size = %d, extend range = %d\n",block_size,extend);
  hist=lk_malloc(sizeof(u_int64_t)*(MAX_COUNT+1));
  if(target_file) thist=lk_malloc(sizeof(u_int64_t)*(MAX_COUNT+1));
  else thist=0;
  lkc=init_compress();
  err=read_ranges_file(ranges_file,&contigs,lkc);
  if(!err && target_file) err=read_target_file(target_file,contigs,extend,lkc);
  nthr=(int)sysconf(_SC_NPROCESSORS_ONLN);
  input_files=argv;
  input_idx=optind;
  n_input_files=argc;
  if(!nthr) nthr=1;
  read_threads=malloc(sizeof(pthread_t)*nthr);
  if(combine) {
    tp.ctg=contigs;
    tp.lkc=lkc;
    for(i=0;i<nthr;i++) {
      if((j=pthread_create(read_threads+i,NULL,read_det_cov,&tp))) abt(__FILE__,__LINE__,"Thread creation %d failed: %d\n",i+1,j);
    }
    for(i=0;i<nthr;i++) pthread_join(read_threads[i],NULL);
  } else {
    for(i=optind;!err && i<argc;i++) {
      err=process_file(argv[i],contigs,block_size,format,&number,lkc);
    }
  }
  printf("Generating histogram\n");
  filter=suffix=0;
  if(lkc->default_compress<COMPRESS_NONE) {
    i=lkc->default_compress;
    filter=lkc->comp_path[i][i==COMPRESS_ZIP?1:0];
    suffix=lkc->compress_suffix[i];
  }
  tn=0;
  if(prefix) {
    asprintf(&tn,"%s_contig_summ.txt",prefix);
    assert(tn);
    ctg_fptr=fopen(tn,"w");
    free(tn);
  } else ctg_fptr=fopen("contig_summ.txt","w");
  if(target_file) {
    for(i=0;i<=MAX_COUNT;i++) hist[i]=thist[i]=0;
    for(ctg=contigs;ctg;ctg=ctg->hh.next) {
      ct=ctg->counts;
      for(i=0;i<(int)ctg->size;i+=block_size) {
	cc=ct[i];
	if(cc<=MAX_COUNT) {
	  hist[cc]++;
	  ctg->tot_count+=cc;
	}
      }
      ct=ctg->tcounts;
      for(i=0;i<(int)ctg->tsize;i+=block_size) {
	cc=ct[i];
	if(cc<=MAX_COUNT) {
	  thist[cc]++;
	  ctg->tot_tcount+=cc;
	}
      }
      if(ctg_fptr) {
	fprintf(ctg_fptr,"%s\t%u\t%"PRIu64"\t%g",ctg->name,ctg->tot_bases,ctg->tot_count,(double)ctg->tot_count/(double)ctg->tot_bases);
	fprintf(ctg_fptr,"\t%u\t%"PRIu64"\t%g\n",ctg->tot_tbases,ctg->tot_tcount,ctg->tot_tbases?(double)ctg->tot_tcount/(double)ctg->tot_tbases:0.0);
      }
    }
    if(ctg_fptr) fclose(ctg_fptr);
    nn=tnn=0;
    for(i=0;i<=MAX_COUNT;i++) {
      nn+=hist[i];
      tnn+=thist[i];
    }
    kk=tkk=0;
    ofptr=fopen(output_file,"w");
    if(!ofptr) ofptr=stdout;
    for(i=0;i<=MAX_COUNT;i++) {
      kk+=hist[i];
      tkk+=thist[i];
      fprintf(ofptr,"%d\t%"PRIu64"\t%g\t%g\t",i,hist[i],(double)hist[i]/(double)nn,(double)kk/(double)nn);
      fprintf(ofptr,"%"PRIu64"\t%g\t%g\n",thist[i],(double)thist[i]/(double)tnn,(double)tkk/(double)tnn);
    }
    if(ofptr!=stdout) fclose(ofptr);
  } else {
    for(i=0;i<=MAX_COUNT;i++) hist[i]=0;
    for(ctg=contigs;ctg;ctg=ctg->hh.next) {
      ct=ctg->counts;
      for(i=0;i<(int)ctg->size;i+=block_size) {
	cc=ct[i];
	if(cc<=MAX_COUNT) {
	  hist[cc]++;
	  ctg->tot_count+=cc;
	}
      }
      if(ctg_fptr) {
	fprintf(ctg_fptr,"%s\t%u\t%"PRIu64"\t%g\n",ctg->name,ctg->tot_bases,ctg->tot_count,(double)ctg->tot_count/(double)ctg->tot_bases);
      }
    }
    if(ctg_fptr) fclose(ctg_fptr);
    nn=0;
    for(i=0;i<=MAX_COUNT;i++) nn+=hist[i];
    kk=0;
    ofptr=fopen(output_file,"w");
    if(!ofptr) ofptr=stdout;
    for(i=0;i<=MAX_COUNT;i++) {
      kk+=hist[i];
      fprintf(ofptr,"%d\t%"PRIu64"\t%g\t%g\n",i,hist[i],(double)hist[i]/(double)nn,(double)kk/(double)nn);
    }
    if(ofptr!=stdout) fclose(ofptr);
  }
  if(detail) {
    tn=0;
    ofptr=0;
    if(filter && suffix) {
      asprintf(&tn,"detailed_coverage.txt.%s",suffix);
      if(tn) {
	i=child_open(WRITE,tn,filter);
	ofptr=fdopen(i,"w");
      }
    }
    if(!ofptr) ofptr=fopen("detailed_coverage.txt","w");
    if(!ofptr) ofptr=stdout;
    printf("Writing detailed coverage information\n");
    for(ctg=contigs;ctg;ctg=ctg->hh.next) {
      k=0;
      ct=ctg->counts;
      /*      ct1=ctg->tcounts;*/
      for(r=ctg->ranges;r;r=r->next) {
	for(i=0;i<r->idx;i++) {
	  if(!k++) fprintf(ofptr,"*%s\t%x\n",ctg->name,r->x1[i]);
	  else fprintf(ofptr,"*\t%x\n",r->x1[i]);
	  out_cnt=0;
	  for(x=r->x1[i];x<=r->x2[i];x++) {
	    cbuf[out_cnt++]=ct[x-1];
	    if(out_cnt==OUT_WIDTH) {
 	      out_cnt=write_cbuf(cbuf,out_cnt,ofptr);
	    }
	  }
	  while(out_cnt) out_cnt=write_cbuf(cbuf,out_cnt,ofptr);
	}
      }
    }
    if(ofptr!=stdout) {
      fclose(ofptr);
      if(tn) {
	free(tn);
	while(waitpid(-1,&i,WNOHANG)>0);
      }
    }
  }
  return err;
}


