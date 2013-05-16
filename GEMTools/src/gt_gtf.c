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
  return gt_shash_get(gtf->types, type, gt_string);
}

GT_INLINE gt_gtf_ref* gt_gtf_get_ref(gt_gtf* gtf, char* const name){
  if(!gt_shash_is_contained(gtf->refs, name)){
    gt_gtf_ref* rr = gt_gtf_ref_new();
    gt_shash_insert(gtf->refs, name, rr, gt_gtf_ref*);
  }
  return gt_shash_get(gtf->refs, name, gt_gtf_ref);
}

GT_INLINE void gt_gtf_read_line(char* line, gt_gtf* const gtf){
  char* ref;
  char* type;
  char* gene_id;
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
  // typee
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
  // last thing where i can not remember what it was
  pch = strtok(NULL, "\t");
  if(pch == NULL){
    return;
  }
  // additional fields
  // search for gene_id
  register bool gid = false;
  while((pch = strtok(NULL, " ")) != NULL){
    if(strcmp("gene_id", pch) == 0){
      gid = true;
    }else{
      if(gid){
        gene_id = pch;
        register uint64_t l = strlen(gene_id);
        register uint64_t off = 1;
        if(gene_id[l-off] == ';'){
          gene_id[l-off] = '\0';
          off = 2;
        }
        if(gene_id[0] == '"'){
          gene_id++;
          gene_id[l-(off+1)] = '\0';
        }
        break;
      }
    }
  }

  gt_string* tp = gt_gtf_get_type(gtf, type);
  gt_gtf_entry* e = gt_gtf_entry_new(start, end, strand, tp);
  if(gene_id != NULL){
    gt_string* gids= gt_string_new(32);
    gt_string_set_string(gids, gene_id);
    e->gene_id = gids;
  }
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
  }else{
    uint64_t e_1 = (*a)->end;
    uint64_t e_2 = (*b)->end;
    if(e_1 < e_2){
      return -1;
    }else if (e_1 > e_2){
      return 1;
    }else{
      return gt_string_cmp((*a)->type, (*b)->type);
    }
  }
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

  uint64_t last_start = 0;
  uint64_t last_end = 0;
  gt_string* last_type = NULL;
  for(; hit<len; hit++){
    e =  *gt_vector_get_elm(entries, hit, gt_gtf_entry*);
    if(e->start > end){
      break;
    }

    // disabled strand check here
    // e->strand == strand
    if(e->start <= end && e->end >= start){
      // check for duplicates
      if(last_type != NULL){
        if(last_type == e->type && last_start == e->start && last_end == e->end){
          continue;
        }
      }
      last_type = e->type;
      last_start = e->start;
      last_end = e->end;
      gt_vector_insert(target, e, gt_gtf_entry*);
    }
  }
  return;
}


GT_INLINE bool gt_gtf_search_exact(gt_gtf* gtf, char* const  ref, const uint64_t start, const uint64_t end, char* const type, const gt_gtf_strand strand){
  if (! gt_shash_is_contained(gtf->refs, ref)){
    return false;
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
    // disabled strand check here
    if(e->start == start && 
       e->end == end && 
       e->strand == strand && 
       strcmp(type, gt_string_get_string(e->type)) == 0){
      return true;
    }
  }
  return false;
}

GT_INLINE void gt_gtf_search_template(gt_gtf* gtf, gt_vector* result, gt_template* template_src){
  // paird end alignments
  register const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  gt_string* exon_type = gt_gtf_get_type(gtf, "exon");
  gt_vector* target = gt_vector_new(32, sizeof(gt_gtf_entry*));
  gt_vector_clear(result);

  // process single alignments
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_REDUCTION(template_src,alignment_dst);
    gt_shash* gene_counts = gt_shash_new();
    GT_ALIGNMENT_ITERATE(alignment_src,map) {
      // search the map
      GT_BEGIN_MAP_BLOCKS_ITERATOR(map, map_it){
        gt_vector_clear(target);
        register uint64_t start = gt_map_get_begin_position(map_it);
        register uint64_t end   = gt_map_get_end_position(map_it); 
        register char* ref = gt_map_get_seq_name(map_it);
        register gt_gtf_strand strand = GTF_STRAND_UNKNOWN;
        if(gt_map_get_strand(map_it) == FORWARD) strand = GTF_STRAND_FORWARD;
        else strand = GTF_STRAND_REVERSE;
        gt_gtf_search(gtf, target, ref, start, end, strand);
        // get all the exons with same gene id
        GT_VECTOR_ITERATE(target, element, counter, gt_gtf_entry*){          
          if((*element)->gene_id != NULL && gt_string_equals(exon_type, (*element)->type)){
            char* gene_id = gt_string_get_string((*element)->gene_id);
            if(! gt_shash_is_contained(gene_counts, gene_id)){
              uint64_t c = 1;
              gt_shash_insert(gene_counts, gene_id, &c, uint64_t*);
            }else{
              uint64_t* c = gt_shash_get(gene_counts, gene_id, uint64_t*);
              *c = *c + 1;
            }
          }
        } 
        GT_SHASH_BEGIN_KEY_ITERATE(gene_counts,c,char*) { 
          gt_vector_insert(result, c, char*); 
        } GT_SHASH_END_ITERATE;
      } GT_END_MAP_BLOCKS_ITERATOR;
    }
    gt_shash_delete(gene_counts, false);
  } GT_TEMPLATE_END_REDUCTION__RETURN;


  
  GT_TEMPLATE__ATTR_ITERATE(template_src,mmap,mmap_attr) {
    gt_shash* gene_counts = gt_shash_new();
    GT_MULTIMAP_ITERATE(mmap, map, end_position){
      GT_BEGIN_MAP_BLOCKS_ITERATOR(map, map_it){
        gt_vector_clear(target);
        register uint64_t start = gt_map_get_begin_position(map_it);
        register uint64_t end   = gt_map_get_end_position(map_it); 
        register char* ref = gt_map_get_seq_name(map_it);
        register gt_gtf_strand strand = GTF_STRAND_UNKNOWN;
        if(gt_map_get_strand(map_it) == FORWARD) strand = GTF_STRAND_FORWARD;
        else strand = GTF_STRAND_REVERSE;
        gt_gtf_search(gtf, target, ref, start, end, strand);
        // get all the exons with same gene id
        GT_VECTOR_ITERATE(target, element, counter, gt_gtf_entry*){          
          if((*element)->gene_id != NULL && gt_string_equals(exon_type, (*element)->type)){
            char* gene_id = gt_string_get_string((*element)->gene_id);
            if(! gt_shash_is_contained(gene_counts, gene_id)){
              uint64_t c = 1;
              gt_shash_insert(gene_counts, gene_id, &c, uint64_t*);
            }else{
              uint64_t* c = gt_shash_get(gene_counts, gene_id, uint64_t*);
              *c = *c + 1;
            }
          }
          
          /*if((*element)->gene_id != NULL && gt_string_equals(exon_type, (*element)->type)){*/
            /*printf("%s : %d-%d : %s\n", */
                /*gt_string_get_string((*element)->type), */
                /*(*element)->start, */
                /*(*element)->end,*/
                /*gt_string_get_string((*element)->gene_id)*/
            /*);*/
          /*}*/
        } 

        GT_SHASH_BEGIN_KEY_ITERATE(gene_counts,c,char*) { 
          gt_vector_insert(result, c, char*); 
        } GT_SHASH_END_ITERATE;
      } GT_END_MAP_BLOCKS_ITERATOR;
    }
    gt_shash_delete(gene_counts, false);
  }
  gt_vector_delete(target);
}

GT_INLINE uint64_t gt_gtf_bin_search_first(gt_vector* entries, uint64_t i){
  if(i == 0){
    return 0;
  }
  register gt_gtf_entry* e =  *gt_vector_get_elm(entries, i, gt_gtf_entry*);
  register gt_gtf_entry* f =  *gt_vector_get_elm(entries, i-1, gt_gtf_entry*);
  while(f->end >= e->start){
    i = i - 1;
    if(i == 0){
      return 0;
    }
    f =  *gt_vector_get_elm(entries, i, gt_gtf_entry*);
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
    if(e->start < t){
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


