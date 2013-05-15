#include "gt_gtf.h"


GT_INLINE gt_gtf_entry* gt_gtf_entry_new(const uint64_t start, const uint64_t end, const gt_gtf_strand strand, gt_string* type){
  gt_gtf_entry* entry = malloc(sizeof(gt_gtf_entry));
  entry->start = start;
  entry->end = end;
  entry->type = type;
  entry->strand = strand;
  return entry;
}

GT_INLINE void gt_gtf_entry_delete(gt_gtf_entry* const entry){
  free(entry);
}

GT_INLINE gt_gtf_ref* gt_gtf_ref_new(){
  gt_gtf_ref* ref = malloc(sizeof(gt_gtf_ref));
  ref->entries = gt_vector_new(10000, sizeof(gt_gtf_entry*));
  return ref;
}

GT_INLINE void gt_gtf_ref_delete(gt_gtf_ref* const ref){
  register uint64_t s = gt_vector_get_used(ref->entries);
  register uint64_t i = 0;
  for(i=0; i<s; i++){
    gt_gtf_entry_delete( (gt_vector_get_elm(ref->entries, i, gt_gtf_entry)));
  }
  gt_vector_delete(ref->entries);
  free(ref);
}

GT_INLINE gt_gtf* gt_gtf_new(){
  gt_gtf* gtf = malloc(sizeof(gt_gtf));
  gtf->refs = gt_shash_new();
  gtf->types = gt_shash_new();
  return gtf;
}

GT_INLINE void gt_gtf_delete(gt_gtf* const gtf){
  gt_shash_delete(gtf->refs, false);
  gt_shash_delete(gtf->types, false);
  free(gtf);
}

GT_INLINE gt_string* gt_gtf_get_type(gt_gtf* gtf, char* const type){
  if(!gt_shash_is_contained(gtf->types, type)){
    gt_string* s = gt_string_new(strlen(type) + 1);
    gt_string_set_string(s, type);
    gt_shash_insert_string(gtf->types, type, s);
  }
  return *gt_shash_get(gtf->types, type, gt_string*);
}

GT_INLINE gt_gtf_ref* gt_gtf_get_ref(gt_gtf* gtf, char* const name){
  if(!gt_shash_is_contained(gtf->refs, name)){
    gt_gtf_ref* rr = gt_gtf_ref_new();
    gt_shash_insert(gtf->refs, name, rr, gt_gtf_ref*);
  }
  return *gt_shash_get(gtf->refs, name, gt_gtf_ref*);
}

GT_INLINE void gt_gtf_read_line(char* line, gt_gtf* const gtf){
  char* ref;
  char* type;
  uint64_t start = 0;
  uint64_t end = 0;
  gt_gtf_strand strand = GTF_STRAND_UNKNOWN;


  char * pch;
  // name
  pch = strtok(line, "\t");
  if(pch == NULL){
    return;
  }
  ref = pch;
  // source
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  // type
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  type = pch;
  // start
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  start = atol(pch);
  // end
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  end = atol(pch);
  // score
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  // strand
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  if(pch[0] == '+'){
    strand = GTF_STRAND_FORWARD;
  }else if(pch[0] == '-'){
    strand = GTF_STRAND_REVERSE;
  }

  gt_string* tp = gt_gtf_get_type(gtf, type);
  gt_gtf_entry* e = gt_gtf_entry_new(start, end, strand, tp);	
  gt_gtf_ref* gtref = gt_gtf_get_ref(gtf, ref);
  gt_vector_insert(gtref->entries, e, gt_gtf_entry*);
}

GT_INLINE int gt_gtf_cmp_entries(const gt_gtf_entry** a, const gt_gtf_entry** b){
  uint64_t s_1 = (*a)->start;
  uint64_t s_2 = (*b)->start;
  if(s_1 < s_2){
    return -1;
  }else if (s_1 > s_2){
    return 1;
  }
  return 0;
}

GT_INLINE gt_gtf* gt_gtf_read(FILE* input){
  gt_gtf* gtf = gt_gtf_new();
  char line[2048];
  while ( fgets(line, 2048, input) != NULL ){
    gt_gtf_read_line(line, gtf);
  }

  // sort the refs

  GT_SHASH_BEGIN_ELEMENT_ITERATE(gtf->refs,shash_element,gt_gtf_ref) {
    qsort(gt_vector_get_mem(shash_element->entries, gt_gtf_entry*), gt_vector_get_used(shash_element->entries), sizeof(gt_gtf_entry**), (int (*)(const void *,const void *))gt_gtf_cmp_entries);
  } GT_SHASH_END_ITERATE
  return gtf;
}

GT_INLINE void gt_gtf_search(gt_gtf* gtf, gt_vector* target, char* const  ref, const uint64_t start, const uint64_t end, const gt_gtf_strand strand){
  gt_vector_clear(target);
  // make sure to no create the ref entries in serch
  if (! gt_shash_is_contained(gtf->refs, ref)){
    return;
  }
  gt_gtf_ref* source_ref = gt_gtf_get_ref(gtf, ref);
  gt_vector* entries = source_ref->entries;
  uint64_t len = gt_vector_get_used(entries);
  uint64_t hit = gt_gtf_bin_search(entries, start);
  register gt_gtf_entry* e;
  for(; hit<len; hit++){
    e =  *gt_vector_get_elm(entries, hit, gt_gtf_entry*);
    if(e->start > end){
      break;
    }
    if(e->strand == strand && e->start <= end && e->end >= start){
      gt_vector_insert(target, e, gt_gtf_entry*);
    }
  }
  return;
}

GT_INLINE void gt_gtf_search_template(gt_gtf* gtf, gt_vector* target, gt_template* template){
  GT_TEMPLATE_ITERATE(template, mmap)
  {
    GT_MULTIMAP_ITERATE(mmap, map, end_position)
    {
      GT_MAP_ITERATE(map, map_block)
      {

      }
    }
  }


}

GT_INLINE uint64_t gt_gtf_bin_search_first(gt_vector* entries, uint64_t i){
  register gt_gtf_entry* e =  *gt_vector_get_elm(entries, i, gt_gtf_entry*);
  register uint64_t c = e->end;
  while(c==e->end){
    if(i == 0){
      return 0;
    }
    i = i - 1;
    e =  *gt_vector_get_elm(entries, i, gt_gtf_entry*);
  }	
  return i+1;
}

GT_INLINE uint64_t gt_gtf_bin_search(gt_vector* entries, uint64_t t){	
  uint64_t used = gt_vector_get_used(entries);
  uint64_t l = 0;
  uint64_t h = used - 1;
  uint64_t m = 0;

  register gt_gtf_entry* e =  *gt_vector_get_elm(entries, h, gt_gtf_entry*);
  while(l < h ){
    m = (l + h) / 2;
    e = *gt_vector_get_elm(entries, m, gt_gtf_entry*);
    if(e->end < t){
      l = m + 1;
    }else{
      h = m;
    }
  }
  e = *gt_vector_get_elm(entries, l, gt_gtf_entry*);
  if (h == l){
    return gt_gtf_bin_search_first(entries, l);
  }else{
    return gt_gtf_bin_search_first(entries, m);
  }
}


