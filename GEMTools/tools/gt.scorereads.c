/*
 * gt.scorereads.c
 *
 *  Created on: 8 Jul 2013
 *      Author: heath
 */

#define GT_SCOREREADS "gt.scorereads"
#define GT_SCOREREADS_VERSION "1.0"

#include <getopt.h>
#include <ctype.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "gem_tools.h"
#include "gt_pipe_io.h"

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
#define GT_SCOREREADS_NUM_LINES (GT_IMP_NUM_LINES>>1)
#define GT_SCOREREADS_DEFAULT_OUTPUT_FORMAT MAP

typedef enum {all, unmapped, poorly_mapped, discordant, trim_maps, one_end_anchored} gt_scorereads_filter;

typedef struct {
  gt_string *name;
  gt_file_format format;
  gt_scorereads_filter filter;
  gt_output_file *output_file;
  gt_generic_printer_attributes *printer_attr;
} output_def;


typedef struct {
  gt_buffered_output_file *buffered_output;
  gt_output_sam_attributes *sam_attributes;
} gt_scorereads_buffered_output;

typedef struct {
  char *input_files[2];
  gt_vector *outputs; // output_def*
  char *dist_file;
  char* name_reference_file;
  char* name_gem_index_file;
  char* sam_header_file;
  char* read_group_id;
  bool mmap_input;
  bool verbose;
  bool compact_format; // for SAM output
  /* Control flags */
  bool load_index;
  bool load_index_sequences;
  bool is_paired;
  bool realign;
  bool bisulfite;
  gt_file_format output_format;
  gt_output_file_compression compress;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output[2];
  int num_threads;
  gt_map_score_attributes map_score_attr;
  gt_sequence_archive *sequence_archive;
} sr_param;

typedef struct {
  int64_t x;
  uint64_t ct;
  UT_hash_handle hh;
} sr_dist_element;

sr_param param = {
  .input_files={NULL,NULL},
  .outputs=NULL,
  .dist_file=NULL,
  .name_reference_file=NULL,
  .name_gem_index_file=NULL,
  .sam_header_file=NULL,
  .read_group_id=NULL,
  .mmap_input=false,
  .compress=NONE,
  .verbose=false,
  .compact_format=false,
  .load_index=false,
  .load_index_sequences=false,
  .is_paired=false,
  .realign=false,
  .bisulfite=false,
  .num_threads=1,
  .map_score_attr.indel_penalty=INDEL_QUAL,
  .map_score_attr.split_penalty=INDEL_QUAL,
  .map_score_attr.mapping_cutoff=0,
  .map_score_attr.max_strata_searched=0,
  .map_score_attr.max_pair_maps=GT_MAP_SCORE_MAX_PAIR_MAPS,
  .map_score_attr.max_orphan_maps=GT_MAP_SCORE_MAX_ORPHAN_MAPS,
  .map_score_attr.quality_format=GT_QUALS_OFFSET_33,
  .map_score_attr.minimum_insert=0,
  .map_score_attr.maximum_insert=0,
  .map_score_attr.insert_phred=NULL,
  .sequence_archive=NULL,
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

gt_sequence_archive* gt_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  if (param.name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(param.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(param.name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",param.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  return sequence_archive;
}

static gt_status gt_scorereads_get_insert_size(gt_template *template,int64_t *size)
{
  gt_status err=GT_STATUS_FAIL;
  int j;
  gt_alignment *al[2];
  uint64_t nmap[2]={0,0};
  for(j=0;j<2;j++) {
    al[j]=gt_template_get_block(template,j);
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

void gt_add_extra_tags(gt_attributes *attributes,uint32_t mcs[],bool paired,gt_map_score_attributes *ms_attr) {
  char *tag;
  if(paired) asprintf(&tag,"ms:B:I,%" PRIu32 ",%" PRIu32 " mx:i:%d",mcs[0],mcs[1],ms_attr->mapping_cutoff);
  else asprintf(&tag,"ms:B:I,%" PRIu32 " mx:i:%d",mcs[0],ms_attr->mapping_cutoff);
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

void gt_template_add_mcs_tags(gt_template *template,gt_map_score_attributes *ms_attr)
{
  uint32_t rd,mcs[2];
  for(rd=0;rd<2;rd++) {
    gt_alignment *al=gt_template_get_block(template,rd);
    mcs[rd]=al?gt_alignment_get_mcs(al):0;
    if(ms_attr->max_strata_searched) {
      uint32_t limit=ms_attr->max_strata_searched+1;
      if(mcs[rd]>limit) mcs[rd]=limit;
    }
  }
  gt_add_extra_tags(template->attributes,mcs,true,ms_attr);
}

void gt_alignment_add_mcs_tags(gt_alignment *al,gt_map_score_attributes *ms_attr)
{
  uint32_t mcs[2];
  mcs[0]=gt_alignment_get_mcs(al);
  if(ms_attr->max_strata_searched) {
    uint32_t limit=ms_attr->max_strata_searched+1;
    if(mcs[0]>limit) mcs[0]=limit;
  }
  gt_add_extra_tags(al->attributes,mcs,false,ms_attr);
}

static int sort_dist_element(void *a,void *b)
{
  return ((sr_dist_element *)a)->x-((sr_dist_element *)b)->x;
}

int estimate_insert_histogram(sr_dist_element *ins,gt_map_score_attributes *ms_attr)
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
  size_t l=ms_attr->maximum_insert-ms_attr->minimum_insert+1;
  double *p=gt_malloc(sizeof(double)*2*l);
  double *cp=p+l;
  int *ph=gt_malloc(sizeof(int)*l);
  for(i=0;i<l;i++) p[i]=0.0;
  for(elem=ins;elem;elem=elem->hh.next) {
    double x=(double)elem->x;
    double n=(double)elem->ct;
    for(i=0;i<l;i++) {
      double d=(x-(double)(ms_attr->minimum_insert+i))/h;
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
  ms_attr->maximum_insert=ms_attr->minimum_insert+i2;
  ms_attr->minimum_insert+=i1;
  l=i2-i1+1;
  ms_attr->insert_phred=gt_malloc(l);
  for(i=i1,k=0;i<=i2;i++) ms_attr->insert_phred[k++]=ph[i];
  //	for(i=0;i<l;i++) printf("%"PRId64"\t%d\n",ms_attr->minimum_insert+i,ms_attr->ins_phred[i]);
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
  else estimate_insert_histogram(ins_dist,&param->map_score_attr);
  gt_scorereads_free_insert_hash(ins_dist);
  gt_buffered_input_file_close(bfile);
  gt_input_file_close(file);
  
}

gt_output_sam_attributes *gt_scorereads_setup_sam_tags(gt_sam_headers *sam_headers,sr_param *param)
{
  // Set up optional record TAGs
  // Set out attributes
  gt_output_sam_attributes *output_sam_attributes = gt_output_sam_attributes_new();
  gt_output_sam_attributes_set_compact_format(output_sam_attributes,param->compact_format);
  gt_output_sam_attributes_set_qualities_offset(output_sam_attributes,param->map_score_attr.quality_format);
  gt_output_sam_attributes_set_print_mismatches(output_sam_attributes,false);
  gt_sam_attributes_add_tag_options(gt_scorereads_attribute_option_list,output_sam_attributes->sam_attributes);
  if(sam_headers->read_group_id_hash) {
    gt_sam_header_record *hr=NULL;
    if (param->read_group_id) {
      size_t* ix=gt_shash_get_element(sam_headers->read_group_id_hash,param->read_group_id);
      if(ix) {
	hr=*(gt_sam_header_record **)gt_vector_get_elm(sam_headers->read_group,*ix,gt_sam_header_record*);
      } else gt_error(SAM_OUTPUT_UNKNOWN_RG_ID,param->read_group_id);
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
  } else if(param->read_group_id) gt_error(SAM_OUTPUT_NO_HEADER_FOR_RG);
  if(param->bisulfite) gt_sam_attributes_add_tag_XB(output_sam_attributes->sam_attributes);
  return output_sam_attributes;
}

void gt_scorereads_print_template(gt_vector *outputs,gt_vector *buffered_outputs,gt_template *template,gt_map_score_attributes *map_score_attr)
{
  gt_scorereads_buffered_output *buf_out=gt_vector_get_mem(buffered_outputs,gt_scorereads_buffered_output);
  bool added_tags=false,calculated_mapq=false;
  GT_VECTOR_ITERATE(outputs,odef_p,idx,output_def*) {
    output_def* odef=*odef_p;
    gt_status print_code=GT_STATUS_OK;
    switch(odef->format) {
    case MAP:
      if(!added_tags) {
	gt_template_add_mcs_tags(template,map_score_attr);
	added_tags=true;
      }
      print_code=gt_output_generic_bofprint_template(buf_out->buffered_output,template,odef->printer_attr);
      break;
    case SAM:
      if(!calculated_mapq) {
	gt_map_calculate_template_mapq_score(template,map_score_attr);
	calculated_mapq=true;
      }
      print_code=gt_output_sam_bofprint_template(buf_out->buffered_output,template,buf_out->sam_attributes);
      break;
    case FASTA:
    case FASTQ:
      print_code=gt_output_fasta_bofprint_template(buf_out->buffered_output,template,odef->printer_attr->output_fasta_attributes);
      break;
    default:
      gt_fatal_error_msg("Fatal error - unsupported output format");
      break;
    }
    if(print_code) {
      gt_error_msg("Error outputting read '"PRIgts"' in %s format\n",PRIgts_content(gt_template_get_string_tag(template)),odef->format==MAP?"MAP":"SAM");
    }
    buf_out++;
  }
}

void gt_scorereads_print_alignment(gt_vector *outputs,gt_vector *buffered_outputs,gt_alignment *alignment,gt_map_score_attributes *map_score_attr)
{
  gt_scorereads_buffered_output *buf_out=gt_vector_get_mem(buffered_outputs,gt_scorereads_buffered_output);
  GT_VECTOR_ITERATE(outputs,odef_p,idx,output_def*) {
    output_def* odef=*odef_p;
    gt_status print_code=GT_STATUS_OK;
    switch(odef->format) {
    case MAP:
      gt_alignment_add_mcs_tags(alignment,map_score_attr);
      print_code=gt_output_generic_bofprint_alignment(buf_out->buffered_output,alignment,odef->printer_attr);
      break;
    case SAM:
      gt_map_calculate_alignment_mapq_score(alignment,map_score_attr);
      print_code=gt_output_sam_bofprint_alignment(buf_out->buffered_output,alignment,buf_out->sam_attributes);
      break;
    case FASTQ:
    case FASTA:
      print_code=gt_output_fasta_bofprint_alignment(buf_out->buffered_output,alignment,odef->printer_attr->output_fasta_attributes);
      break;
    default:
      gt_fatal_error_msg("Fatal error - unsupported output format");
      break;
    }
    if(print_code) {
      gt_error_msg("Error outputting read '"PRIgts"' in %s format\n",PRIgts_content(gt_alignment_get_string_tag(alignment)),odef->format==MAP?"MAP":"SAM");
    }
    buf_out++;
  }
}

gt_vector* gt_scorereads_open_and_attach_buffered_output_files(gt_vector *outputs,gt_buffered_input_file *buffered_input,gt_sam_headers *sam_headers,sr_param *param)
{
  GT_NULL_CHECK(outputs);
  GT_NULL_CHECK(buffered_input);
  gt_vector *output_buffers=gt_vector_new(gt_vector_get_used(outputs),sizeof(gt_scorereads_buffered_output));
  GT_VECTOR_ITERATE(outputs,odef_p,idx,output_def*) {
    output_def* odef=*odef_p;
    gt_scorereads_buffered_output *buf_out=gt_vector_get_free_elm(output_buffers,gt_scorereads_buffered_output);
    gt_vector_inc_used(output_buffers);
    buf_out->buffered_output=gt_buffered_output_file_new(odef->output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buf_out->buffered_output);
    buf_out->sam_attributes = odef->format==SAM?gt_scorereads_setup_sam_tags(sam_headers,param):NULL;
  }
  return output_buffers;
}

void gt_scorereads_close_buffered_output_files(gt_vector *buffered_outputs)
{
  GT_NULL_CHECK(buffered_outputs);
  GT_VECTOR_ITERATE(buffered_outputs,buf_out,idx,gt_scorereads_buffered_output) {
    if(buf_out->sam_attributes) gt_output_sam_attributes_delete(buf_out->sam_attributes);
    gt_buffered_output_file_close(buf_out->buffered_output);
  }
}

gt_sam_headers *gt_scorereads_setup_sam_headers(sr_param *param)
{
  gt_sam_headers* sam_headers = gt_sam_header_new();
  if(param->sam_header_file) {
    gt_input_file* const sam_headers_input_file = gt_input_file_sam_open(param->sam_header_file,false);
    uint64_t characters_read = 0, lines_read = 0;
    gt_status error_code=gt_input_file_sam_read_headers((char *)sam_headers_input_file->file_buffer,sam_headers_input_file->buffer_size,sam_headers,&characters_read,&lines_read);
    if(error_code) gt_error(PARSE_SAM_HEADER_NOT_SAM,sam_headers_input_file->file_name);
    gt_input_file_close(sam_headers_input_file);
  }
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
    gt_string *pn_st=gt_string_set_new(GT_SCOREREADS);
    gt_sam_header_record_add_tag(hr,"PN",pn_st);
    gt_string *vn_st=gt_string_set_new(GT_SCOREREADS_VERSION);
    if(prev_id) gt_sam_header_record_add_tag(hr,"PP",prev_id);
    gt_sam_header_record_add_tag(hr,"VN",vn_st);
    gt_sam_header_add_program_record(sam_headers,hr);
  }
  // Open reference file
  if (param->load_index) {
    param->sequence_archive = gt_open_sequence_archive(param->load_index_sequences);
    gt_sam_header_load_sequence_archive(sam_headers,param->sequence_archive);
  }
  return sam_headers;
}

static void alignment_strip_bisulfite_contig_names(gt_alignment *al)
{
  GT_ALIGNMENT_CHECK(al);
  GT_ALIGNMENT_ITERATE(al,map) {
    uint64_t bis_type=GT_BIS_TYPE_NORMAL;
    bool bt_flag=false;
    bool strand_misms=false;
    GT_MAP_ITERATE(map,mapb) {
      gt_string *seq=mapb->seq_name;
      uint64_t l=gt_string_get_length(seq);
      char *p=gt_string_get_string(seq);
      if(l>4) {
	if(p[l-4]=='#' || p[l-4]=='_') {
	  if(!strncmp(p+l-3,"C2T",3)) {
	    gt_string_trim_right(seq,4);
	    if(bt_flag==true) {
	      if(bis_type!=GT_BIS_TYPE_C2T) {
		strand_misms=true;
		break;
	      }
	    } else {
	      bis_type=GT_BIS_TYPE_C2T;
	      bt_flag=true;
	    }
	  } else if(!strncmp(p+l-3,"G2A",3)) {
	    gt_string_trim_right(seq,4);	  
	    if(bt_flag==true) {
	      if(bis_type!=GT_BIS_TYPE_G2A) {
		strand_misms=true;
		break;
	      }
	    } else {
	      bis_type=GT_BIS_TYPE_G2A;
	      bt_flag=true;
	    }
	  } else {
	    if(bis_type!=GT_BIS_TYPE_NORMAL) {
		strand_misms=true;
		break;
	    }
	    bt_flag=true;
	  }
	}
      }
    }
    if(strand_misms==true) bis_type=GT_BIS_TYPE_MISMATCH;
    gt_attributes *attr=map->attributes;
    if(attr==NULL) map->attributes=attr=gt_attributes_new();
    gt_attributes_add(attr,GT_ATTR_ID_BIS_TYPE,&bis_type,uint64_t);
  }
}

static void template_strip_bisulfite_contig_names(gt_template *tp)
{
  for(int j=0;j<2;j++) {
    gt_alignment *al=gt_template_get_block(tp,j);
    if(al) alignment_strip_bisulfite_contig_names(al);
  }
}

gt_status gt_scorereads_process(sr_param *param)
{
  gt_status err=GT_STATUS_OK;
  gt_sam_headers* sam_headers = NULL; // SAM headers
  // Open output file(s)
  GT_VECTOR_ITERATE(param->outputs,odef_p,idx,output_def*) {
    output_def* odef = *odef_p;
    odef->output_file=odef->name?gt_output_file_new_compress(gt_string_get_string(odef->name),UNSORTED_FILE,param->compress):
      gt_output_stream_new_compress(stdout,UNSORTED_FILE,param->compress);
    gt_cond_fatal_error(!odef->output_file,FILE_OPEN,odef->name?gt_string_get_string(odef->name):"<STDOUT>");
    switch(odef->format) {
    case MAP:
      odef->printer_attr=gt_generic_printer_attributes_new(MAP);
      odef->printer_attr->output_map_attributes->print_casava=true;
      odef->printer_attr->output_map_attributes->print_extra=true;
      odef->printer_attr->output_map_attributes->print_scores=true;
      odef->printer_attr->output_map_attributes->hex_print_scores=true;
      break;
    case FASTA:
    case FASTQ:
      odef->printer_attr=gt_generic_printer_attributes_new(odef->format);
      break;
    case SAM:
      if(sam_headers == NULL) sam_headers=gt_scorereads_setup_sam_headers(param);
      // Print SAM headers
      gt_output_sam_ofprint_headers_sh(odef->output_file,sam_headers);
      break;
    default:
      gt_fatal_error_msg("Fatal error - unsupported output format");
      break;
    }
  }
  // Do we have two map files as input (one for each read)?
  if(param->input_files[1]) {
    pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
    gt_input_file* input_file1=gt_input_file_map_open(param->input_files[0],param->mmap_input);
    gt_input_file* input_file2=gt_input_file_map_open(param->input_files[1],param->mmap_input);
    if(input_file1->file_format!=MAP || input_file2->file_format!=MAP) {
      gt_fatal_error_msg("Fatal error: paired files '%s','%s' are not in MAP format\n",param->input_files[0],param->input_files[1]);
    }
    gt_buffered_input_file* buffered_input_a=gt_buffered_input_file_new(input_file1);
    gt_buffered_input_file* buffered_input_b=gt_buffered_input_file_new(input_file2);
    // If no insert size file supplied, estimate from the first GT_SCOREREADS_NUM_LINES records
    if(!param->map_score_attr.insert_phred) {
      uint64_t nread=0;
      gt_template *template=gt_template_new();
      sr_dist_element *ins_dist=NULL;
      gt_status error_code;
      while((error_code=gt_input_map_parser_synch_get_template(buffered_input_a,buffered_input_b,template,&mutex))) {
	if(error_code!=GT_IMP_OK) continue;
	int64_t size;
	error_code=gt_scorereads_get_insert_size(template,&size);
	if(error_code==GT_STATUS_OK) (void)sr_increase_insert_count(&ins_dist,size,1);
	if(++nread==GT_SCOREREADS_NUM_LINES) break;
      }
      gt_template_delete(template);
      // Reset input buffers to the beginning of the file
      buffered_input_a->cursor=gt_vector_get_mem(buffered_input_a->block_buffer,char);
      buffered_input_a->current_line_num=0;
      buffered_input_b->cursor=gt_vector_get_mem(buffered_input_b->block_buffer,char);
      buffered_input_b->current_line_num=0;
      (void)estimate_insert_histogram(ins_dist,&param->map_score_attr);
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
      gt_vector *buffered_outputs=gt_scorereads_open_and_attach_buffered_output_files(param->outputs,buffered_input1,sam_headers,param);
      gt_map_score_tmap_buf *ms_tmap_buf=gt_map_score_new_map_score_tmap_buf(GT_MAX(param->map_score_attr.max_pair_maps,param->map_score_attr.max_orphan_maps));
      gt_status error_code;
      gt_template *template=gt_template_new();
      while((error_code=gt_input_map_parser_synch_get_template(buffered_input1,buffered_input2,template,&mutex))) {
	if(error_code!=GT_IMP_OK) continue;
	gt_map_pair_template(template,&param->map_score_attr,ms_tmap_buf);
	if(param->bisulfite) template_strip_bisulfite_contig_names(template);
	if(param->realign) gt_template_realign_levenshtein(template,param->sequence_archive);
	gt_scorereads_print_template(param->outputs,buffered_outputs,template,&param->map_score_attr);
      }
      gt_template_delete(template);
      gt_map_score_delete_map_score_tmap_buf(ms_tmap_buf);
      gt_buffered_input_file_close(buffered_input1);
      gt_buffered_input_file_close(buffered_input2);
      gt_scorereads_close_buffered_output_files(buffered_outputs);
    }
    gt_input_file_close(input_file1);
    gt_input_file_close(input_file2);
  } else { // Single input file
    gt_input_file* input_file=param->input_files[0]?gt_input_file_map_open(param->input_files[0],param->mmap_input):gt_input_stream_map_open(stdin);
    gt_buffered_input_file* buffered_input_a=gt_buffered_input_file_new(input_file);
    gt_map_parser_attributes* input_map_attributes_a = gt_input_map_parser_attributes_new(param->is_paired);
    if(param->is_paired) {
      /*
       * Paired input
       *
       * If no insert size file supplied, estimate from the first GT_SCOREREADS_NUM_LINES records
       *
       */
      if(!param->map_score_attr.insert_phred) {
	uint64_t nread=0;
	gt_template *template=gt_template_new();
	sr_dist_element *ins_dist=NULL;
	gt_status error_code;
	while ((error_code=gt_input_map_parser_get_template(buffered_input_a,template,input_map_attributes_a))) {
	  if(error_code!=GT_IMP_OK) continue;
	  int64_t size;
	  error_code=gt_scorereads_get_insert_size(template,&size);
	  if(error_code==GT_STATUS_OK) (void)sr_increase_insert_count(&ins_dist,size,1);
	  if(++nread==GT_SCOREREADS_NUM_LINES) break;
	}
	gt_template_delete(template);
	// Reset input buffers to the beginning of the file
	buffered_input_a->cursor=gt_vector_get_mem(buffered_input_a->block_buffer,char);
	buffered_input_a->current_line_num=0;
	(void)estimate_insert_histogram(ins_dist,&param->map_score_attr);
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
	gt_buffered_input_file* buffered_input=tid?gt_buffered_input_file_new(input_file):buffered_input_a;
	gt_vector *buffered_outputs=gt_scorereads_open_and_attach_buffered_output_files(param->outputs,buffered_input,sam_headers,param);
	gt_map_parser_attributes* input_map_attributes=tid?gt_input_map_parser_attributes_new(param->is_paired):input_map_attributes_a;
	gt_output_sam_attributes* const output_sam_attributes = param->output_format==SAM?gt_scorereads_setup_sam_tags(sam_headers,param):NULL;
	gt_map_score_tmap_buf *ms_tmap_buf=gt_map_score_new_map_score_tmap_buf(GT_MAX(param->map_score_attr.max_pair_maps,param->map_score_attr.max_orphan_maps));
	gt_status error_code;
	gt_template *template=gt_template_new();
	while ((error_code=gt_input_map_parser_get_template(buffered_input,template,input_map_attributes))) {
	  if (error_code!=GT_IMP_OK) continue;
	  gt_map_pair_template(template,&param->map_score_attr,ms_tmap_buf);
	  if(param->bisulfite) template_strip_bisulfite_contig_names(template);
	  if(param->realign) gt_template_realign_levenshtein(template,param->sequence_archive);
	  gt_scorereads_print_template(param->outputs,buffered_outputs,template,&param->map_score_attr);
	}
	// Clean
	gt_template_delete(template);
	gt_map_score_delete_map_score_tmap_buf(ms_tmap_buf);
	gt_input_map_parser_attributes_delete(input_map_attributes);
	if(output_sam_attributes) gt_output_sam_attributes_delete(output_sam_attributes);
	gt_buffered_input_file_close(buffered_input);
	gt_scorereads_close_buffered_output_files(buffered_outputs);
      }
    } else {
      /*
       * Single end reads
       */
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(param->num_threads)
#endif
      {
#ifdef HAVE_OPENMP
	uint64_t tid=omp_get_thread_num();
#else
	uint64_t tid=0;
#endif
	gt_buffered_input_file* buffered_input=tid?gt_buffered_input_file_new(input_file):buffered_input_a;
	gt_vector *buffered_outputs=gt_scorereads_open_and_attach_buffered_output_files(param->outputs,buffered_input,sam_headers,param);
	gt_map_parser_attributes* input_map_attributes=tid?gt_input_map_parser_attributes_new(false):input_map_attributes_a;
	gt_status error_code;
	gt_template *template=gt_template_new();
	while ((error_code=gt_input_map_parser_get_template(buffered_input,template,input_map_attributes))) {
	  if (error_code!=GT_IMP_OK) {
	    gt_error_msg("Error parsing file '%s'\n",param->input_files[0]);
	    continue;
	  }
	  gt_alignment *alignment=gt_template_get_block(template,0);
	  GT_ALIGNMENT_ITERATE(alignment,map) {
	    if(map->gt_score==GT_MAP_NO_GT_SCORE) map->gt_score=gt_map_calculate_gt_score(alignment,map,&param->map_score_attr);
	    map->phred_score=255;
	  }
	  if(param->bisulfite) alignment_strip_bisulfite_contig_names(alignment);
	  if(param->realign) gt_alignment_realign_levenshtein(alignment,param->sequence_archive);
	  gt_scorereads_print_alignment(param->outputs,buffered_outputs,alignment,&param->map_score_attr);
	}
	// Clean
	gt_template_delete(template);
	gt_input_map_parser_attributes_delete(input_map_attributes);
	gt_buffered_input_file_close(buffered_input);
	gt_scorereads_close_buffered_output_files(buffered_outputs);
      }
    }
    gt_input_file_close(input_file);
  }
  GT_VECTOR_ITERATE(param->outputs,odef_p1,idx1,output_def*) {
    output_def* odef = *odef_p1;
    gt_output_file_close(odef->output_file);
    if(odef->format==MAP) gt_generic_printer_attributes_delete(odef->printer_attr);
  }
  if(param->map_score_attr.insert_phred) free(param->map_score_attr.insert_phred);
  return err;
}

output_def *gt_scorereads_new_output_def(char *name,gt_file_format format,gt_scorereads_filter filter)
{
  output_def *od=gt_alloc(output_def);
  if(name==NULL || !strcmp(name,"-")) od->name=NULL;
  else od->name=gt_string_set_new(name);
  od->format=format;
  od->filter=filter;
  od->output_file=NULL;
  od->printer_attr=NULL;
  return od;
}

char *gt_scorereads_parse_output_option(char *opt)
{
  char *format_st[]={"MAP","FASTA","SAM","FASTQ",0};
  gt_file_format formats[]={MAP,FASTA,SAM,FASTQ};
  char *filter_st[]={"ALL","UNMAPPED","POORLY-MAPPED","DISCORDANT","TRIM-MAPS","ONE-END-ANCHORED",0};
  gt_scorereads_filter filters[]={all,unmapped,poorly_mapped,discordant,trim_maps,one_end_anchored};
  
  char *p=strchr(opt,',');
  while(p) {
    p++;
    int idx=0;
    while(format_st[idx] && strncasecmp(p,format_st[idx],strlen(format_st[idx]))) idx++;
    if(format_st[idx]) {
      char *p1=p+strlen(format_st[idx]);
      gt_file_format format=formats[idx];
      if(*p1==',') {
        idx=0;
        p1++;
        while(filter_st[idx] && strncasecmp(p1,filter_st[idx],strlen(filter_st[idx]))) idx++;
        if(filter_st[idx]) {
          p1+=strlen(format_st[idx]);
          if(!*p1) {
            size_t sz=p-opt;
            char *nm=gt_malloc(sz);
            strncpy(nm,opt,sz-1);
            nm[sz-1]=0;
            gt_vector_insert(param.outputs,gt_scorereads_new_output_def(nm,format,filters[idx]),output_def*);
            return NULL;
          }
        }
      } else if(!*p1) {
        size_t sz=p-opt;
        char *nm=NULL;
        if(sz!=2 || *opt!='-') {
	  nm=gt_malloc(sz);
	  strncpy(nm,opt,sz-1);
	  nm[sz-1]=0;
        }
        gt_vector_insert(param.outputs,gt_scorereads_new_output_def(nm,format,all),output_def*);
        return NULL;
      } else p=strchr(p,',');
    }
  }
  return opt;
}

gt_status parse_arguments(int argc,char** argv) {
  gt_status err=GT_STATUS_OK;
  struct option* gt_scorereads_getopt = gt_options_adaptor_getopt(gt_scorereads_options);
  gt_string* const gt_scorereads_short_getopt = gt_options_adaptor_getopt_short(gt_scorereads_options);
  int option, option_index;
  char *p,*output_file=NULL;
  gt_file_format output_format=GT_SCOREREADS_DEFAULT_OUTPUT_FORMAT;
  param.outputs=gt_vector_new(2,sizeof(output_def *));

  int insert_set[2]={0,0};
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
                            gt_string_get_string(gt_scorereads_short_getopt),gt_scorereads_getopt,&option_index))==-1) break;
    switch (option) {
      /* I/O */
    case 303:
      param.dist_file = optarg;
      break;
    case 'o':
      output_file = gt_scorereads_parse_output_option(optarg);
      break;
    case 'b':
      param.bisulfite=true;
      // Intentional fall through
    case 'R':
      param.realign=true;
      break;
    case 'r':
      param.name_reference_file = optarg;
      param.load_index = true;
      break;
    case 'I':
      param.name_gem_index_file = optarg;
      param.load_index = true;
      break;
    case 's':
      param.sam_header_file = optarg;
      break;
    case 300:
      param.input_files[0] = optarg;
      break;
    case 301:
      param.input_files[1] = optarg;
      break;
    case 'p':
      param.is_paired=true;
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
	param.map_score_attr.quality_format=GT_QUALS_OFFSET_64;
      } else if (gt_streq(optarg,"offset-33")) {
	param.map_score_attr.quality_format=GT_QUALS_OFFSET_33;
      } else {
	gt_fatal_error_msg("Quality format not recognized: '%s'",optarg);
      }
      break;
    case 704: // output-format
      if (gt_streq(optarg,"MAP")) {
	output_format = MAP;
      } else if (gt_streq(optarg,"SAM")) {
	output_format = SAM;
      } else if (gt_streq(optarg,"FASTA")) {
	output_format = FASTA;
      } else if (gt_streq(optarg,"FASTQ")) {
	output_format = FASTQ;
      } else {
	gt_fatal_error_msg("Output format '%s' not recognized",optarg);
      }
      break;
    case 401:
      param.map_score_attr.minimum_insert=(int)strtol(optarg,&p,10);
      if(*p || param.map_score_attr.minimum_insert<0) {
	fprintf(stderr,"Illegal minimum insert size: '%s'\n",optarg);
	err=-7;
      } else insert_set[0]=1;
      break;
    case 402:
      param.map_score_attr.maximum_insert=(int)strtol(optarg,&p,10);
      if(*p || param.map_score_attr.maximum_insert<0) {
	fprintf(stderr,"Illegal maximum insert size: '%s'\n",optarg);
	err=-7;
      } else insert_set[1]=1;
      break;
    case 403:
      param.map_score_attr.indel_penalty=(int)strtol(optarg,&p,10);
      if(*p || param.map_score_attr.indel_penalty<0) {
	fprintf(stderr,"Illegal indel penalty: '%s'\n",optarg);
	err=-7;
      }
      break;
    case 404:
      param.map_score_attr.max_pair_maps=strtoul(optarg,&p,10);
      if(*p) {
	fprintf(stderr,"Illegal max pair maps value: '%s'\n",optarg);
	err=-7;
      }
      break;
    case 405:
      param.map_score_attr.max_orphan_maps=strtoul(optarg,&p,10);
      if(*p) {
	fprintf(stderr,"Illegal max orphan maps value: '%s'\n",optarg);
	err=-7;
      }
      break;
    case 500:
      if(gt_sam_attributes_parse_tag_option_string(gt_scorereads_attribute_option_list,optarg)!=GT_STATUS_OK) {
	fprintf(stderr,"Unable to parse --tag option '%s'\n",optarg);
	err=-8;
      }
      break;
    case 'S':
      param.map_score_attr.split_penalty=(int)strtol(optarg,&p,10);
      if(*p || param.map_score_attr.split_penalty<0) {
	fprintf(stderr,"Illegal split penalty score: '%s'\n",optarg);
	err=-7;
      }
      break;
    case 'M':
      param.map_score_attr.mapping_cutoff=strtoul(optarg,&p,10);
      if(*p) {
	fprintf(stderr,"Illegal mapping cutoff: '%s'\n",optarg);
	err=-7;
      }
      break;
    case 'm':
      param.map_score_attr.max_strata_searched=strtoul(optarg,&p,10);
      if(*p) {
	fprintf(stderr,"Illegal mismatch limit: '%s'\n",optarg);
	err=-7;
      }
      break;
      /* Headers */
    case 600: // Read-group ID
      param.read_group_id = optarg;
      break;
      /* Format */
    case 'c':
      param.compact_format = true;
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
      gt_options_fprint_json_menu(stderr,gt_scorereads_options,gt_scorereads_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      usage(gt_scorereads_options,gt_scorereads_groups,false);
      gt_fatal_error_msg("Option '%c' %d not recognized",option,option);
      break;
    }
  }
  /*
   * Parameters check
   */
  if(param.input_files[1] && param.input_files[0]==NULL) {
    param.input_files[0]=param.input_files[1];
    param.input_files[1]=NULL;
  }
  if(param.input_files[1]) param.is_paired=true;
  if(err==GT_STATUS_OK && insert_set[0] && insert_set[1] && param.map_score_attr.minimum_insert>param.map_score_attr.maximum_insert) {
    fputs("Minimum insert size > maximum insert size\n",stderr);
    usage(gt_scorereads_options,gt_scorereads_groups,false);
    err=-15;
  } else if(param.realign) {
    if(!param.load_index) {
      fputs("Index or reference must be supplied for realign or bisulfite option\n",stderr);
      usage(gt_scorereads_options,gt_scorereads_groups,false);
      err=-16;
    } else {
      param.load_index_sequences=true;
    }
  }
  if(err==GT_STATUS_OK) {
    // Create output_def for normal (non specified) output using specified format.  If no output command at all specified, create a stdout stream
    if(output_file!=NULL || !gt_vector_get_used(param.outputs))
      gt_vector_insert(param.outputs,gt_scorereads_new_output_def(output_file,output_format,all),output_def*);
    bool use_stdout=false;
    GT_VECTOR_ITERATE(param.outputs,odef_p,idx,output_def*) {
      gt_string *com;
      if((com=(*odef_p)->name)) {
	uint64_t l=gt_string_get_length(com);
	switch(param.compress) {
	case GZIP:
	  if(l<3 || strcmp(gt_string_get_string(com)+l-3,".gz")) gt_string_append_string(com,".gz",3);
	  break;
	case BZIP2:
	  if(l<4 || strcmp(gt_string_get_string(com)+l-4,".bz2")) gt_string_append_string(com,".bz2",4);
	  break;
	default:
	  break;
	}
      } else {
	if(use_stdout) {
	  gt_error_msg("Error in output definitions: multiple outputs to stdout stream defined\n");
	  err=GT_STATUS_FAIL;
	} else use_stdout=true;
      }
    }
    if(!insert_set[1]) {
      if(param.map_score_attr.minimum_insert<=1000) param.map_score_attr.maximum_insert=1000;
      else param.map_score_attr.maximum_insert=param.map_score_attr.minimum_insert+1000;
    }
    if(param.is_paired && param.dist_file) read_dist_file(&param);
  }
  // Free
  gt_string_delete(gt_scorereads_short_getopt);
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

