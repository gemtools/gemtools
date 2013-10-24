/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_score.c
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_score.h"
#include "gt_sam_attributes.h"
/*
 * Map Score Accessors
 */
GT_INLINE uint64_t gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->gt_score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score) {
  GT_MAP_CHECK(map);
  map->gt_score = score;
}
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->phred_score;
}
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score) {
  GT_MAP_CHECK(map);
  map->phred_score = phred_score;
}

uint64_t gt_map_calculate_gt_score(gt_alignment *al, gt_map *map, gt_map_score_attributes *ms_attr)
{
	int qual_offset=ms_attr->quality_format==GT_QUALS_OFFSET_33?33:64;
	int indel_penalty=ms_attr->indel_penalty;
	int split_penalty=ms_attr->split_penalty;
	int map_cutoff=ms_attr->mapping_cutoff;
	register gt_string* const quals = al->qualities;
	register gt_string* const read = al->read;
	register const bool has_qualities = gt_alignment_has_qualities(al);
	uint64_t score=0;
	GT_MAP_ITERATE(map,map_block) {
		GT_MISMS_ITERATE(map_block,misms) {
			int quality_misms;
			if (has_qualities) {
				quality_misms = gt_string_get_string(quals)[misms->position]-qual_offset;
				if(quality_misms>GT_MAP_SCORE_MAX_QUALITY) quality_misms=GT_MAP_SCORE_MAX_QUALITY;
				else if(quality_misms<map_cutoff || gt_string_get_string(read)[misms->position]=='N') quality_misms=0;
			} else quality_misms=GT_MAP_SCORE_MISSING_QUALITY;
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
	if(score>GT_MAP_SCORE_MAX_GT_SCORE) score=GT_MAP_SCORE_MAX_GT_SCORE;
	return score;
}

void gt_map_pair_template(gt_template *template,gt_map_score_attributes *ms_attr)
{
	GT_TEMPLATE_CHECK(template);
	GT_MAP_SCORE_ATTRIBUTES_CHECK(ms_attr);

	gt_mmap_attributes attr;
	gt_alignment *al[2];
	gt_map *mmap[2];
	uint64_t rd,nmap[3]={0,0,0};
	for(rd=0;rd<2;rd++) {
		al[rd]=gt_template_get_block(template,rd);
		if(al[rd]) {
			uint64_t max_complete_strata=gt_alignment_get_mcs(al[rd]);
			if(ms_attr->max_strata_searched) {
				uint32_t limit=ms_attr->max_strata_searched+1;
				if(max_complete_strata>limit) max_complete_strata=limit;
			}
			nmap[rd]=gt_alignment_get_num_maps(al[rd]);
			if(nmap[rd]) {
				GT_ALIGNMENT_ITERATE(al[rd],map) {
					map->gt_score=gt_map_calculate_gt_score(al[rd],map,ms_attr);
				}
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
				if(gt_err==GT_TEMPLATE_INSERT_SIZE_OK && x>=ms_attr->minimum_insert && x<=ms_attr->maximum_insert) {
					attr.distance=gt_map_get_global_distance(map1)+gt_map_get_global_distance(map2);
					attr.gt_score=map1->gt_score|(map2->gt_score<<16);
					if(ms_attr->insert_phred) attr.gt_score|=((uint64_t)ms_attr->insert_phred[x-ms_attr->minimum_insert]<<32);
					attr.phred_score=255;
					gt_template_inc_counter(template,attr.distance);
					gt_template_add_mmap_ends(template,map1,map2,&attr);
					map_flag[0][i]=map_flag[1][j]=1;
					nmap[2]++;
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
					attr.gt_score=rd?(map->gt_score<<16):map->gt_score;
					attr.phred_score=255;
					gt_template_inc_counter(template,attr.distance);
					if(!rd)	gt_template_add_mmap_ends(template,map,0,&attr);
					else gt_template_add_mmap_ends(template,0,map,&attr);
				}
			}
		}
		free(map_flag[0]);
	} else {
		attr.distance=0;
		attr.gt_score=0;
		attr.phred_score=255;
		gt_template_add_mmap_ends(template,0,0,&attr);
	}
//	gt_attributes_remove(template->attributes,GT_ATTR_ID_TAG_PAIR);
}

/*
 * To check for duplicate alignments we make a key out of a concatenation of sequence name, position and strand.  If we have
 * multiple segments, keys for all segments are concatenated.  Keys for paired alignments are formed from the two keys
 * for each end joined with "::"
 *
 */

static size_t gt_map_score_get_key_size(gt_map *map)
{
	size_t ksize=0;
	uint64_t nbl=0;
	GT_MAP_SEGMENT_ITERATOR(map,map_segment_iterator) {
		nbl++;
		ksize+=gt_string_get_length(gt_map_segment_iterator_get_map(&map_segment_iterator)->seq_name);
	}
	return ksize+nbl*(sizeof(map->position)+sizeof(map->strand));
}

static size_t gt_map_score_get_pair_key_size(gt_map *map1,gt_map *map2)
{
	return gt_map_score_get_key_size(map1)+gt_map_score_get_key_size(map2)+2;
}

static size_t gt_map_score_make_key(void *buf,gt_map *map)
{
	size_t bpos=0;
  GT_MAP_SEGMENT_ITERATOR(map,map_segment_iterator) {
  	gt_map *mp=gt_map_segment_iterator_get_map(&map_segment_iterator);
  	gt_string *seq=mp->seq_name;
  	memcpy(buf+bpos,gt_string_get_string(seq),gt_string_get_length(seq));
  	bpos+=gt_string_get_length(seq);
  	memcpy(buf+bpos,&mp->position,sizeof(mp->position));
  	bpos+=sizeof(mp->position);
  	memcpy(buf+bpos,&mp->strand,sizeof(mp->strand));
  	bpos+=sizeof(mp->strand);
  }
  return bpos;
}

static size_t gt_map_score_make_pair_key(void *buf,gt_map *map1,gt_map *map2)
{
	size_t bpos=gt_map_score_make_key(buf,map1);
	memcpy(buf+bpos,"::",2);
	return bpos+2+gt_map_score_make_key(buf+bpos+2,map2);
}

static int gt_map_cmp_int(const void *s1,const void *s2)
{
	return (*(int *)s1)-(*(int *)s2);
}

#define PHRED_KONST -0.23025850929940456840 // -log(10)/10;

void gt_map_calculate_template_mapq_score(gt_template *template,gt_map_score_attributes *ms_attr)
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
	uint32_t map_cutoff=ms_attr->mapping_cutoff;
	uint32_t max_complete_strata[2]={0,0};
	for(rd=0;rd<2;rd++) {
		gt_alignment *al=gt_template_get_block(template,rd);
		max_complete_strata[rd]=al?gt_alignment_get_mcs(al):0;
		if(ms_attr->max_strata_searched) {
			uint32_t limit=ms_attr->max_strata_searched+1;
			if(max_complete_strata[rd]>limit) max_complete_strata[rd]=limit;
		}
	}
	size_t buf_len=1024;
	char *buf=gt_malloc(buf_len);
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
				size_t key_size=gt_map_score_get_key_size(maps[rd]);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				(void)gt_map_score_make_key(buf,maps[rd]);
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
				size_t key_size=gt_map_score_get_pair_key_size(maps[0],maps[1]);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				(void)gt_map_score_make_pair_key(buf,maps[0],maps[1]);
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
	  int quals_offset=ms_attr->quality_format==GT_QUALS_OFFSET_33?33:64;
	  for(i=0;i<quals->length;i++) {
	  	int q=gt_string_get_string(quals)[i]-quals_offset;
	  	if(q>=map_cutoff) qvs[k++]=q;
	  }
	  if(k>=max_complete_strata[rd]) {
	  	qsort(qvs,k,sizeof(int),gt_map_cmp_int);
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
				size_t key_size=gt_map_score_make_key(buf,maps[rd]);
				HASH_FIND(hh,mhash[rd],buf,key_size,mp_hash);
				assert(mp_hash);
	      GT_MAP_ITERATE(maps[rd],mapb) { mapb->phred_score=mp_hash->phred; }
			}
			if(maps[0] && maps[1]) {
				size_t key_size=gt_map_score_make_pair_key(buf,maps[0],maps[1]);
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

void gt_map_calculate_alignment_mapq_score(gt_alignment *alignment,gt_map_score_attributes *ms_attr)
{
	uint32_t map_cutoff=ms_attr->mapping_cutoff;
	uint32_t max_complete_strata=0;

	typedef struct {
		gt_map *map;
		double prob;
		uint64_t score;
	} pmap;

	pmap *maplist=NULL;
	max_complete_strata=gt_alignment_get_mcs(alignment);
	if(ms_attr->max_strata_searched) {
		uint32_t limit=ms_attr->max_strata_searched+1;
		if(max_complete_strata>limit) max_complete_strata=limit;
	}
	uint64_t min_score=0xffff;
	// Build up list of single end alignments (need this so we can scale the MAPQ score)
	uint64_t nmaps=gt_alignment_get_num_maps(alignment);
	if(nmaps) {
		maplist=gt_malloc(sizeof(pmap)*(nmaps+1));
		uint64_t k=0;
		GT_ALIGNMENT_ITERATE(alignment,map) {
			uint64_t score=map->gt_score;
			if(score==GT_MAP_NO_GT_SCORE) gt_fatal_error(ALIGNMENT_NOT_SCORED);
			maplist[k].score=score&0xffff;
			maplist[k].map=map;
			if(maplist[k].score<min_score) min_score=maplist[k].score;
			k++;
		}
		// Get score of best possible aligning read that we didn't look for with the mapping parameters
		uint64_t fake_sc=0;
		if(max_complete_strata) {
			gt_string* const quals=alignment->qualities;
			uint64_t i;
			int *qvs=malloc(sizeof(int)*quals->length);
			k=0;
			int quals_offset=ms_attr->quality_format==GT_QUALS_OFFSET_33?33:64;
			for(i=0;i<quals->length;i++) {
				int q=gt_string_get_string(quals)[i]-quals_offset;
				if(q>=map_cutoff) qvs[k++]=q;
			}
			if(k>=max_complete_strata) {
				qsort(qvs,k,sizeof(int),gt_map_cmp_int);
				for(i=0;i<max_complete_strata;i++) fake_sc+=qvs[i];
			}
			free(qvs);
			maplist[nmaps].score=fake_sc;
			maplist[nmaps].map=NULL;
			if(fake_sc<min_score) min_score=fake_sc;
		}
		double z=0.0;
		for(k=0;k<=nmaps;k++) {
			maplist[k].prob=exp(PHRED_KONST*((double)(maplist[k].score)-(double)min_score));
			z+=maplist[k].prob;
		}
		for(k=0;k<nmaps;k++) {
			maplist[k].prob/=z;
			int tp;
			if(1.0-maplist[k].prob<1.0e-255) tp=254;
			else {
				tp=(int)(0.5+log(1.0-maplist[k].prob)/PHRED_KONST);
				if(tp>254) tp=254;
			}
			GT_MAP_ITERATE(maplist[k].map,mapb) { mapb->phred_score=tp;	}
		}
		free(maplist);
	}
}

