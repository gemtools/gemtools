#include "gt_gtf.h"

GT_INLINE gt_gtf_entry* gt_gtf_entry_new(const uint64_t start, const uint64_t end, const gt_strand strand, gt_string* const type){
  gt_gtf_entry* entry = malloc(sizeof(gt_gtf_entry));
  entry->uid = 0;
  entry->start = start;
  entry->end = end;
  entry->type = type;
  entry->strand = strand;
  return entry;
}
GT_INLINE void gt_gtf_entry_delete(gt_gtf_entry* const entry){
  free(entry);
}

GT_INLINE gt_gtf_ref* gt_gtf_ref_new(void){
  gt_gtf_ref* ref = malloc(sizeof(gt_gtf_ref));
  ref->entries = gt_vector_new(GTF_DEFAULT_ENTRIES, sizeof(gt_gtf_entry*));
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

GT_INLINE gt_gtf* gt_gtf_new(void){
  gt_gtf* gtf = malloc(sizeof(gt_gtf));
  gtf->refs = gt_shash_new();
  gtf->types = gt_shash_new();
  gtf->gene_ids = gt_shash_new();
  gtf->transcript_ids = gt_shash_new();
  gtf->gene_types = gt_shash_new();
  return gtf;
}

GT_INLINE void gt_gtf_delete(gt_gtf* const gtf){
  gt_shash_delete(gtf->refs, true);
  gt_shash_delete(gtf->types, true);
  gt_shash_delete(gtf->gene_ids, true);
  gt_shash_delete(gtf->transcript_ids, true);
  gt_shash_delete(gtf->gene_types, true);
  free(gtf);
}

GT_INLINE gt_gtf_hits* gt_gtf_hits_new(void){
  gt_gtf_hits* hits = malloc(sizeof(gt_gtf_hits));
  hits->exon_hits = gt_vector_new(16, sizeof(gt_gtf_hit*));
  hits->num_genes = 0;
  hits->num_protein_coding =0;
  hits->num_paired_genes =0;
  return hits;
}
GT_INLINE void gt_gtf_hits_delete(gt_gtf_hits* const hits){
  gt_gtf_hits_clear(hits);
  gt_vector_delete(hits->exon_hits);
  free(hits);
}
GT_INLINE void gt_gtf_hits_clear(gt_gtf_hits* const hits){
  uint64_t i = 0;
  for(i=0; i<gt_vector_get_used(hits->exon_hits); i++){
    gt_gtf_hit* hit = *gt_vector_get_elm(hits->exon_hits, i, gt_gtf_hit*);
    gt_gtf_hit_delete(hit);
  }
  hits->num_genes = 0;
  hits->num_protein_coding =0;
  hits->num_paired_genes =0;
  gt_vector_clear(hits->exon_hits);
}


GT_INLINE gt_string* gt_gtf_get_type(const gt_gtf* const gtf, char* const type){
  if(!gt_gtf_contains_type(gtf, type)){
    gt_string* s = gt_string_set_new(type);
    gt_shash_insert_string(gtf->types, type, s);
  }
  return gt_shash_get(gtf->types, type, gt_string);
}
GT_INLINE bool gt_gtf_contains_type(const gt_gtf* const gtf, char* const name){
	return gt_shash_is_contained(gtf->types, name);
}

GT_INLINE gt_gtf_ref* gt_gtf_get_ref(const gt_gtf* const gtf, char* const name){
  if(!gt_gtf_contains_ref(gtf, name)){
    gt_gtf_ref* rr = gt_gtf_ref_new();
    gt_shash_insert(gtf->refs, name, rr, gt_gtf_ref*);
  }
  return gt_shash_get(gtf->refs, name, gt_gtf_ref);
}
GT_INLINE bool gt_gtf_contains_ref(const gt_gtf* const gtf, char* const name){
	return gt_shash_is_contained(gtf->refs, name);
}

GT_INLINE gt_string* gt_gtf_get_gene_id(const gt_gtf* const gtf, char* const name){
  if(!gt_gtf_contains_gene_id(gtf, name)){
    gt_string* const gene_id = gt_string_set_new(name);
    gt_shash_insert(gtf->gene_ids, name, gene_id, gt_string*);
  }
  return gt_shash_get(gtf->gene_ids, name, gt_string);
}
GT_INLINE bool gt_gtf_contains_gene_id(const gt_gtf* const gtf, char* const name){
	return gt_shash_is_contained(gtf->gene_ids, name);
}

GT_INLINE gt_string* gt_gtf_get_transcript_id(const gt_gtf* const gtf, char* const name){
  if(!gt_gtf_contains_transcript_id(gtf, name)){
    gt_string* const gene_id = gt_string_set_new(name);
    gt_shash_insert(gtf->transcript_ids, name, gene_id, gt_string*);
  }
  return gt_shash_get(gtf->transcript_ids, name, gt_string);
}
GT_INLINE bool gt_gtf_contains_transcript_id(const gt_gtf* const gtf, char* const name){
	return gt_shash_is_contained(gtf->transcript_ids, name);
}

GT_INLINE gt_string* gt_gtf_get_gene_type(const gt_gtf* const gtf, char* const name){
  if(!gt_gtf_contains_gene_type(gtf, name)){
    gt_string* const gene_type = gt_string_set_new(name);
    gt_shash_insert(gtf->gene_types, name, gene_type, gt_string*);
  }
  return gt_shash_get(gtf->gene_types, name, gt_string);
}
GT_INLINE bool gt_gtf_contains_gene_type(const gt_gtf* const gtf, char* const name){
	return gt_shash_is_contained(gtf->gene_types, name);
}

/**
 * Comparator that compares two gtf_entries by starting position
 */
GT_INLINE int gt_gtf_sort_by_start_cmp_(const gt_gtf_entry** a, const gt_gtf_entry** b){
  uint64_t p1 = (*a)->start;
  uint64_t p2 = (*b)->start;
  return p1 < p2 ? -1 : (p1>p2 ? 1 : gt_string_cmp( (*a)->type, (*b)->type ));
}
/**
 * Comparator that compares two gtf_entries by ending position
 */
GT_INLINE int gt_gtf_sort_by_end_cmp_(const gt_gtf_entry** a, const gt_gtf_entry** b){
  uint64_t p1 = (*a)->end;
  uint64_t p2 = (*b)->end;
  return p1 < p2 ? -1 : (p1>p2 ? 1 : gt_string_cmp( (*a)->type, (*b)->type ));
}

/**
 * Sort vector of gt_gtf_entries by starting position
 */
GT_INLINE void gt_gtf_sort_by_start(gt_vector* entries) {
  qsort(gt_vector_get_mem(entries, gt_gtf_entry*),
    gt_vector_get_used(entries),
    sizeof(gt_gtf_entry**),
    (int (*)(const void *,const void *))gt_gtf_sort_by_start_cmp_);
}
/**
 * Sort vector of gt_gtf_entries by ending position
 */
GT_INLINE void gt_gtf_sort_by_end( gt_vector* entries) {
  qsort(gt_vector_get_mem(entries, gt_gtf_entry*),
    gt_vector_get_used(entries),
    sizeof(gt_gtf_entry**),
    (int (*)(const void *,const void *))gt_gtf_sort_by_end_cmp_);
}
GT_INLINE gt_gtf_node* gt_gtf_create_node(gt_vector* entries){
  const uint64_t len = gt_vector_get_used(entries);
  if(len == 0){
    return NULL;
  }
  gt_gtf_node* const node = malloc(sizeof(gt_gtf_node));
  const gt_gtf_entry* mid = *gt_vector_get_elm(entries, len/2, gt_gtf_entry*);
  node->midpoint = mid->start + ((mid->end - mid->start)/2);
  node->entries_by_end = gt_vector_new(16, sizeof(gt_gtf_entry*));
  node->entries_by_start = gt_vector_new(16, sizeof(gt_gtf_entry*));
  gt_vector* to_left = gt_vector_new(16, sizeof(gt_gtf_entry*));
  gt_vector* to_right = gt_vector_new(16, sizeof(gt_gtf_entry*));
  GT_VECTOR_ITERATE(entries, element, counter, gt_gtf_entry*){
    if((*element)->end < node->midpoint){
      gt_vector_insert(to_left, (*element), gt_gtf_entry*);
    }else if((*element)->start > node->midpoint){
      gt_vector_insert(to_right, (*element), gt_gtf_entry*);
    }else{
      gt_vector_insert(node->entries_by_end, (*element), gt_gtf_entry*);
      gt_vector_insert(node->entries_by_start, (*element), gt_gtf_entry*);
    }
  }
  // sort the start and end lists
  gt_gtf_sort_by_start(node->entries_by_start);
  gt_gtf_sort_by_end(node->entries_by_end);

  // delete incoming entry list
  gt_vector_delete(entries);
  if(gt_vector_get_used(to_left) > 0){
    // create left node
    node->left = gt_gtf_create_node(to_left);
  }else{
    node->left = NULL;
    gt_vector_delete(to_left);
  }
  if(gt_vector_get_used(to_right) > 0){
    // create right node
    node->right = gt_gtf_create_node(to_right);
  }else{
    node->right = NULL;
    gt_vector_delete(to_right);
  }
  return node;
}

/*
 * Read next tab separated field from line or return NULL if no such field exists
 */
GT_INLINE char* gt_gtf_read_gtf_field_(char** line){
  char* current = *line;
  GT_READ_UNTIL(line, **line=='\t');
  if(GT_IS_EOL(line)) return NULL;
  **line = EOS;
  GT_NEXT_CHAR(line);
  return current;
}

GT_INLINE gt_status gt_gtf_read_attributes_(char** line, gt_shash* attrs){
  gt_shash_clear(attrs, false);
  while(!GT_IS_EOL(line)){
    while(**line == ' ') GT_NEXT_CHAR(line);
    if(**line == EOL || **line == EOS) return GT_STATUS_OK;
    // get the attribute name
    char* name = *line;
    GT_READ_UNTIL(line, **line==' ')
    if(GT_IS_EOL(line)){
      gt_error_msg("Error parsing GTF attributes. Expected space but found end of line");
      return GT_GTF_INVALID_LINE;
    }
    **line = EOS;
    GT_NEXT_CHAR(line);

    // skip to attribute start
    while(**line == ' ') GT_NEXT_CHAR(line);
    // remove starting quote
    if(**line == '"') GT_NEXT_CHAR(line);
    char* attr = *line;
    // skip until the closing ;
    while(**line != ';') GT_NEXT_CHAR(line);
    if(GT_IS_EOL(line)) return GT_GTF_INVALID_LINE;
    // remove trailing quotes and add EOS
    if(*(*line-1) == '"') *(*line-1) = EOS;
    else **line = EOS;
    GT_NEXT_CHAR(line);

    // add attribute
    if(gt_shash_is_contained(attrs, name)){
      gt_shash_remove(attrs, name, false);
    }
    gt_shash_insert(attrs, name, attr, char*);

    if(gt_shash_is_contained(attrs, "gene_id") &&
       gt_shash_is_contained(attrs, "gene_type") &&
       gt_shash_is_contained(attrs, "transcript_id")){
        return GT_STATUS_OK;
    }

  }
  return GT_STATUS_OK;
}
/**
 * Parse a single GTF line
 */
GT_INLINE gt_status gt_gtf_read_line(char* line, gt_gtf* const gtf, uint64_t counter, gt_shash* attrs){
  // skip comments
  if(line[0] == '#'){
    return GT_STATUS_OK;
  }


  char* ref = NULL;
  char* type = NULL;
  uint64_t start = 0;
  uint64_t end = 0;
  gt_strand strand = UNKNOWN;
  char* current = line;

  ref = gt_gtf_read_gtf_field_(&line);
  if(ref == NULL){
    gt_error_msg("Unable to parse name: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }

  // SKIP source
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL){
    gt_error_msg("Unable to parse source: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }

  // type
  type = gt_gtf_read_gtf_field_(&line);
  if(type == NULL){
    gt_error_msg("Unable to parse type: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }


  // start
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL){
    gt_error_msg("Unable to parse start: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }
  start = atol(current);

  // end
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL){
    gt_error_msg("Unable to parse end: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }
  end = atol(current);

  // SKIP score
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL){
    gt_error_msg("Unable to parse score: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }


  // strand
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL) return GT_GTF_INVALID_LINE;
  if(current == NULL){
    gt_error_msg("Unable to parse strand: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }

  if(*current == '+'){
    strand = FORWARD;
  }else if(*current == '-'){
    strand = REVERSE;
  }

  // SIKP last thing where i can not remember what it was
  current = gt_gtf_read_gtf_field_(&line);
  if(current == NULL){
    gt_error_msg("Unable to parse last: '%s'", line);
    return GT_GTF_INVALID_LINE;
  }


  // WARNING >>> the attribute parser stops after
  // the currently used feels are found. If you want
  // to add a field, also update the attribute parser
  if(gt_gtf_read_attributes_(&line, attrs) != GT_STATUS_OK){
    gt_error_msg("Unable to parse attributes: '%s'", line);
    return GT_GTF_INVALID_ATTRIBUTES;
  }

  // get the type or create it
  gt_string* tp = gt_gtf_get_type(gtf, type);
  gt_gtf_entry* e = gt_gtf_entry_new(start, end, strand, tp);
  e->uid = counter;
  if(gt_shash_is_contained(attrs, "gene_id")){
    e->gene_id = gt_gtf_get_gene_id(gtf, gt_shash_get(attrs, "gene_id", char));
  }
  if(gt_shash_is_contained(attrs, "gene_type")){
    e->gene_type = gt_gtf_get_gene_type(gtf, gt_shash_get(attrs, "gene_type", char));
  }
  if(gt_shash_is_contained(attrs, "transcript_id")){
    e->transcript_id = gt_gtf_get_transcript_id(gtf, gt_shash_get(attrs, "transcript_id", char));
  }
  // get the ref or create it
  gt_gtf_ref* gtref = gt_gtf_get_ref(gtf, ref);
  gt_vector_insert(gtref->entries, e, gt_gtf_entry*);
  return GT_STATUS_OK;
}

bool gt_gtf_hits_junction(gt_map* map, gt_gtf_entry* e){
  uint64_t rs = gt_map_get_begin_mapping_position(map);
  uint64_t re = gt_map_get_end_mapping_position(map);
  bool hit = (rs==e->start) || (rs==e->end) || (re == e->end) || (re == e->start);
  return hit;
}

void gt_gtf_print_entry_(gt_gtf_entry* e, gt_map* map){
  if(map != NULL){
    gt_output_map_fprint_map(stdout, map, NULL);
    printf(" ==> ");
  }
  printf("%s : %"PRIu64" - %"PRIu64" (%c)", e->type->buffer, e->start, e->end, (e->strand==FORWARD?'+':'-') );
  if(e->gene_id != NULL){
    printf(" GID:%s", e->gene_id->buffer);
  }
  if(e->transcript_id != NULL){
    printf(" TID:%s", e->transcript_id->buffer);
  }

  if(e->type != NULL){
    printf(" [%s]", e->type->buffer);
  }

  if(e->gene_type != NULL){
    printf(" [%s]", e->gene_type->buffer);
  }

  if(map != NULL && gt_gtf_hits_junction(map, e)){
    printf(" [Hits JS]");
  }
  printf("\n");
}



GT_INLINE gt_gtf_hit* gt_gtf_hit_new(void){
  gt_gtf_hit* hit = malloc(sizeof(gt_gtf_hit));
  hit->exon_overlap = 0.0;
  hit->intron_length = 0.0;
  hit->is_protein_coding = false;
  hit->junction_hits = 0.0;
  hit->map = NULL;
  hit->num_junctions = 0;
  hit->pairs_transcript = false;
  hit->pairs_splits = false;
  hit->pairs_gene = false;
  hit->num_template_blocks = 0;
  hit->transcripts = NULL;
  hit->genes = NULL;
  hit->hits_exon = false;
  return hit;
}

GT_INLINE void gt_gtf_hit_delete(gt_gtf_hit* hit){
  if(hit->transcripts != NULL){
    gt_shash_delete(hit->transcripts, true);
  }
  if(hit->genes != NULL){
    gt_shash_delete(hit->genes, true);
  }
  free(hit);
}


GT_INLINE gt_status gt_gtf_reload_buffer(gt_buffered_input_file* const buffered_fasta_input) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_fasta_input);
  // Dump buffer if BOF it attached to input, and get new out block (always FIRST)
  gt_buffered_input_file_dump_attached_buffers(buffered_fasta_input->attached_buffered_output_file);
  // Read new input block
  const uint64_t read_lines = gt_buffered_input_file_get_block(buffered_fasta_input, GT_NUM_LINES_50K);
  if (gt_expect_false(read_lines==0)) return GT_INPUT_FILE_EOF;
  // Assign block ID
  gt_buffered_input_file_set_id_attached_buffers(buffered_fasta_input->attached_buffered_output_file,buffered_fasta_input->block_id);
  return GT_STATUS_OK;
}



GT_INLINE gt_status gt_gtf_get_line(gt_buffered_input_file* const buffered_input, gt_string* const line) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input);
  GT_STRING_CHECK(line);

  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_input)) {
    if ((error_code=gt_gtf_reload_buffer(buffered_input))!=GT_IMP_OK) return error_code;
  }
  // Prepare the template
  char* const line_start = buffered_input->cursor;
  gt_string_clear(line);
  GT_INPUT_FILE_SKIP_LINE(buffered_input);
  gt_string_set_nstring_static(line, line_start, (buffered_input->cursor - line_start));
  return GT_IMP_OK;
}

GT_INLINE uint64_t gt_gtf_merge_(const gt_gtf* const target, gt_gtf* source, uint64_t counter){
  // get the type or create it
  GT_SHASH_BEGIN_KEY_ITERATE(source->refs, key){
    gt_gtf_ref* source_ref = gt_gtf_get_ref(source, key);
    gt_gtf_ref* target_ref = gt_gtf_get_ref(target, key);
    GT_VECTOR_ITERATE(source_ref->entries, value, c, gt_gtf_entry*){
      gt_gtf_entry* e = *value;
      e->uid = counter++;
      if(e->gene_id != NULL){
        e->gene_id = gt_gtf_get_gene_id(target, gt_string_get_string(e->gene_id));
      }
      if(e->transcript_id != NULL){
        e->transcript_id = gt_gtf_get_transcript_id(target, gt_string_get_string(e->transcript_id));
      }
      if(e->type != NULL)e->type = gt_gtf_get_type(target, gt_string_get_string(e->type));
      if(e->gene_type != NULL)e->gene_type = gt_gtf_get_gene_type(target, gt_string_get_string(e->gene_type));
      gt_vector_insert(target_ref->entries, e, gt_gtf_entry*);
    }
  }GT_SHASH_END_ITERATE;
  return counter;
}


GT_INLINE gt_gtf* gt_gtf_read_from_stream(FILE* input, uint64_t threads){
  gt_input_file* input_file = gt_input_stream_open(input);
  return gt_gtf_read(input_file, threads);
}
GT_INLINE gt_gtf* gt_gtf_read_from_file(char* input, uint64_t threads){
  gt_input_file* input_file = gt_input_file_open(input, false);
  return gt_gtf_read(input_file, threads);
}

GT_INLINE gt_gtf* gt_gtf_read(gt_input_file* input_file, const uint64_t threads){
  GT_NULL_CHECK(input_file);
  GT_ZERO_CHECK(threads);

  uint64_t counter = 0;
  uint64_t i = 0;
  gt_gtf* const gtf = gt_gtf_new();

  gt_gtf** gtfs = gt_calloc(threads-1, gt_gtf*, true);
  for(i=0; i<threads-1; i++){
    gtfs[i] = gt_gtf_new();
  }

  #pragma omp parallel num_threads(threads)
  {
    uint64_t tid = omp_get_thread_num();
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_string* buffered_line = gt_string_new(GTF_MAX_LINE_LENGTH);
    gt_gtf* thread_gtf;
    if(tid == 0){
      thread_gtf = gtf;
    }else{
      thread_gtf = gtfs[tid-1];
    }
    gt_shash* attrs = gt_shash_new();
    while(gt_gtf_get_line(buffered_input, buffered_line)){
      if(gt_gtf_read_line(buffered_line->buffer, thread_gtf, buffered_input->current_line_num, attrs) != GT_STATUS_OK){
        // raise error
        gt_fatal_error_msg("Failed to parse GTF line '%s'", buffered_line->buffer);
      }
      counter++;
    }
    gt_shash_delete(attrs, false);
    gt_buffered_input_file_close(buffered_input);
    gt_string_delete(buffered_line);
  }
  gt_input_file_close(input_file);

  counter = 0;

  // merge all the thread gtfs into a single one
  for(i=0; i<threads-1; i++){
    counter = gt_gtf_merge_(gtf, gtfs[i], counter);
    gt_gtf_delete(gtfs[i]);
  }
  free(gtfs);


  gt_string* const exon_t = gt_string_set_new("exon");
  gt_string* const intron_t = gt_string_set_new("intron");
  // sort the refs
  GT_SHASH_BEGIN_ELEMENT_ITERATE(gtf->refs,shash_element,gt_gtf_ref) {
    // sort by start position
    gt_gtf_sort_by_start(shash_element->entries);
    uint64_t size = gt_vector_get_used(shash_element->entries);
    uint64_t i = 0;
    gt_shash* last_exons = gt_shash_new();

    for(i=0; i<size; i++){
      gt_gtf_entry* entry = *gt_vector_get_elm(shash_element->entries, i, gt_gtf_entry*);
      if(entry->type != NULL && gt_string_equals(exon_t, entry->type)){
        gt_string* transcript_id = entry->transcript_id;
        if(transcript_id != NULL){
          if(!gt_shash_is_contained(last_exons, gt_string_get_string(transcript_id))){
            gt_shash_insert(last_exons, gt_string_get_string(transcript_id), entry, gt_gtf_entry*);
          }else{
            gt_gtf_entry* prev_exon = gt_shash_get_element(last_exons, gt_string_get_string(transcript_id));
            gt_gtf_entry* intron = gt_gtf_entry_new(prev_exon->end+1,
                                                    entry->start-1,
                                                    prev_exon->strand,
                                                    intron_t);
            intron->transcript_id = transcript_id;
            intron->gene_id = prev_exon->gene_id;
            intron->uid = counter++;
            gt_vector_insert(shash_element->entries, intron, gt_gtf_entry*);
            gt_shash_remove(last_exons, gt_string_get_string(transcript_id),false);
            gt_shash_insert(last_exons, gt_string_get_string(transcript_id), entry, gt_gtf_entry*);
//            gt_gtf_print_entry_(prev_exon);
//            gt_gtf_print_entry_(intron);
//            gt_gtf_print_entry_(entry);
          }
        }

      }
    }
    gt_shash_delete(last_exons, false);


    // create a interval tree node for each ref
    shash_element->node = gt_gtf_create_node(shash_element->entries);
  } GT_SHASH_END_ITERATE
  return gtf;
}

/*
 * Binary search for start position
 */
GT_INLINE uint64_t gt_gtf_bin_search(gt_vector* const entries, const uint64_t t, const uint64_t end){
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
    return l;
  }else{
    return m;
  }
}

GT_INLINE void gt_gtf_search_node_(gt_gtf_node* node, const uint64_t start, const uint64_t end, gt_vector* const target){
  if(node == NULL) return;

  // add overlapping intervals from this node
  GT_VECTOR_ITERATE(node->entries_by_start, element, counter, gt_gtf_entry*){
    if((*element)->start > end){
      break;
    }
    gt_gtf_entry* e = *element;
    //if((*element)->start <= start && (*element)->end >= end){
    if((start < e->end && end > e->start)
      || (start >= e->start && end <=e->end)
      || (start < e->end && end >= e->end)
      || (start < e->start && end > e->end)){
      gt_vector_insert(target, (*element), gt_gtf_entry*);
    }
  }


  if(end < node->midpoint || start < node->midpoint){
    // search left tree
    gt_gtf_search_node_(node->left, start, end, target);
  }
  if (start > node->midpoint || end > node->midpoint){
    gt_gtf_search_node_(node->right, start, end, target);
  }

}




GT_INLINE uint64_t gt_gtf_search(const gt_gtf* const gtf,  gt_vector* const target, char* const ref, const uint64_t start, const uint64_t end, const bool clear_target){
  if(clear_target)gt_vector_clear(target);
  // make sure the target ref is contained
  if (! gt_shash_is_contained(gtf->refs, ref)){
    return 0;
  }
  const gt_gtf_ref* const source_ref = gt_gtf_get_ref(gtf, ref);
  gt_gtf_search_node_(source_ref->node, start, end, target);
  return gt_vector_get_used(target);
}

GT_INLINE void gt_gtf_count_(gt_shash* const table, char* const element){
  if(!gt_shash_is_contained(table, element)){
    uint64_t* v = gt_malloc_uint64();
    *v = 1;
    gt_shash_insert(table, element, v, uint64_t);
  }else{
    uint64_t* v = gt_shash_get(table,element,uint64_t);
    ++(*v);
  }
}

GT_INLINE uint64_t gt_gtf_get_count_(gt_shash* const table, char* const element){
  if(!gt_shash_is_contained(table, element)){
    return 0;
  }
  uint64_t* v = gt_shash_get(table,element,uint64_t);
  return *v;
}

GT_INLINE void gt_gtf_create_hit(gt_vector* search_hits, gt_shash* all_genes, gt_gtf_hits* hits, gt_gtf_hit* template_hit){
  template_hit->transcripts = gt_shash_new();
  template_hit->genes = gt_shash_new();
  template_hit->is_protein_coding = false;
  template_hit->hits_exon = false;

  GT_VECTOR_ITERATE(search_hits, v, c, gt_gtf_entry*){
    gt_gtf_entry* e = *v;
    // count gene id
    if(e->gene_id != NULL){
      gt_gtf_count_(template_hit->genes, gt_string_get_string(e->gene_id));
      gt_gtf_count_(all_genes, gt_string_get_string(e->gene_id));
    }
    if(e->transcript_id != NULL){
      gt_gtf_count_(template_hit->transcripts, gt_string_get_string(e->transcript_id));
    }
    if(strcmp(e->type->buffer, "exon") == 0){
      if(e->gene_type != NULL){
        template_hit->is_protein_coding |= (strcmp(e->gene_type->buffer, "protein_coding"));
        if(!template_hit->hits_exon){
          hits->num_protein_coding++;
        }
      }
      template_hit->hits_exon = true;
    }
  }
  template_hit->pairs_gene = (gt_shash_get_num_elements(template_hit->genes) == 1); // single gene
  template_hit->pairs_transcript = (gt_shash_get_num_elements(template_hit->transcripts) == 1); // single gene

  hits->num_paired_genes += (template_hit->pairs_gene ? 1 : 0);
  gt_vector_insert(hits->exon_hits, template_hit, gt_gtf_hit*);
}

GT_INLINE void gt_gtf_search_template_hits(const gt_gtf* const gtf, gt_gtf_hits* const hits, gt_template* const template_src){
  gt_vector* const search_hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  // reset the hits
  gt_gtf_hits_clear(hits);
  gt_shash* all_genes = gt_shash_new();

  // process paired alignment
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template_src,mmap,mmap_attr) {
    gt_gtf_search_map(gtf, search_hits, mmap[0], true);
    gt_gtf_search_map(gtf, search_hits, mmap[1], false);
    gt_gtf_hit* template_hit = gt_gtf_hit_new();
    template_hit->num_template_blocks = gt_template_get_num_blocks(template_src);
    template_hit->mmap = mmap;
    template_hit->map = NULL;
    template_hit->map_attributes = mmap_attr;
    template_hit->num_junctions = (gt_map_get_num_blocks(mmap[0]) + gt_map_get_num_blocks(mmap[1])) - 2;
    gt_gtf_create_hit(search_hits, all_genes, hits, template_hit);
  }
  hits->num_genes = gt_shash_get_num_elements(all_genes);
  gt_shash_delete(all_genes, true);
  gt_vector_delete(search_hits);
}

GT_INLINE void gt_gtf_search_alignment_hits(const gt_gtf* const gtf, gt_gtf_hits* const hits, gt_alignment* const alignment){
  gt_vector* const search_hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  // reset the hits
  gt_gtf_hits_clear(hits);
  gt_shash* all_genes = gt_shash_new();

  // process paired alignment
  GT_ALIGNMENT_ITERATE(alignment, map){
    gt_gtf_search_map(gtf, search_hits, map, true);
    gt_gtf_hit* template_hit = gt_gtf_hit_new();
    template_hit->map = map;
    template_hit->mmap = NULL;
    template_hit->num_junctions = gt_map_get_num_blocks(map);
    template_hit->num_template_blocks = 1;
    gt_gtf_create_hit(search_hits, all_genes, hits, template_hit);
  }
  hits->num_genes = gt_shash_get_num_elements(all_genes);
  gt_shash_delete(all_genes, true);
  gt_vector_delete(search_hits);
}



GT_INLINE void gt_gtf_count_add_(gt_shash* const source, gt_shash* const target){
  GT_SHASH_BEGIN_ITERATE(source, key, value, uint64_t){
    if(!gt_shash_is_contained(target, key)){
      uint64_t* v = gt_malloc_uint64();
      *v = *value;
      gt_shash_insert(target, key, v, uint64_t);
    }else{
      uint64_t* v = gt_shash_get(target,key,uint64_t);
      *v += (*value);
    }
  }GT_SHASH_END_ITERATE;
}

GT_INLINE void gt_gtf_gene_map_(gt_gtf* const gtf, gt_map* const map, gt_shash* const gene_counts, gt_shash* const gene_type_counts){
  uint64_t start = gt_map_get_begin_mapping_position(map);
  uint64_t end   = gt_map_get_end_mapping_position(map);
  gt_vector* const hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  gt_gtf_search(gtf, hits, gt_map_get_seq_name(map), start, end, true);

  GT_VECTOR_ITERATE(hits, e, i, gt_gtf_entry*){
    gt_gtf_entry* hit = *e;
    // count gene if we find an exonic hit
    if(strcmp(hit->type->buffer, "exon") == 0){
      gt_gtf_count_(gene_counts, gt_string_get_string(hit->gene_id));
      gt_gtf_count_(gene_type_counts, gt_string_get_string(hit->gene_type));
    }
  }
  gt_vector_delete(hits);
}


GT_INLINE void gt_gtf_gene_map(gt_gtf* const gtf, gt_map* const map, gt_shash* const genes, gt_shash* const types){
  uint64_t blocks = gt_map_get_num_blocks(map);
  if(blocks == 0) return;
  GT_MAP_ITERATE(map, map_block){
    gt_gtf_gene_map_(gtf, map_block, genes, types);
  }
}

GT_INLINE uint64_t gt_gtf_count_map_(gt_gtf* const gtf, gt_map* const map, gt_shash* const type_counts, gt_shash* const gene_counts, gt_shash* const gene_type_counts, gt_map* const other_pair){
  uint64_t start = gt_map_get_begin_mapping_position(map);
  uint64_t end   = gt_map_get_end_mapping_position(map);
  gt_vector* const hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  gt_gtf_search(gtf, hits, gt_map_get_seq_name(map), start, end, true);

  gt_shash* const local_type_counts = gt_shash_new();
  gt_shash* local_gene_counts = gt_shash_new();
  gt_shash* local_gene_type_counts = gt_shash_new();

  GT_VECTOR_ITERATE(hits, e, i, gt_gtf_entry*){
    gt_gtf_entry* hit = *e;
    gt_gtf_count_(local_type_counts, gt_string_get_string(hit->type));
    // count gene if we find an exonic hit
    if(strcmp(hit->type->buffer, "exon") == 0){
      gt_gtf_count_(local_gene_counts, gt_string_get_string(hit->gene_id));
      gt_gtf_count_(local_gene_type_counts, gt_string_get_string(hit->gene_type));
    }
  }

  uint64_t num_gene_hits = gt_shash_get_num_elements(local_gene_counts);
  uint64_t hit_by_both_ends = 0;
  // count types
  if(gt_vector_get_used(hits) == 0){
    gt_gtf_count_(type_counts, "na");
  }else if(gt_gtf_get_count_(local_type_counts, "exon") > 0){
    // try to resolve multi gene hits for paired end data
    if(other_pair != NULL){
      gt_shash* const other_gene_counts = gt_shash_new();
      gt_shash* const other_gene_type_counts = gt_shash_new();
      gt_gtf_gene_map(gtf, other_pair, other_gene_counts, other_gene_type_counts);

      uint64_t hits = 0;
      GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
        if(gt_shash_is_contained(other_gene_counts, key)) hits++;
      }GT_SHASH_END_ITERATE;
      if(hits == 1){
        // only a single gene hit that is matched by both pairs
        // replace the count maps
        num_gene_hits = 1;
        hit_by_both_ends = 1;
        gt_shash* const merged_gene_counts = gt_shash_new();
        gt_shash* const merged_gene_type_counts = gt_shash_new();
        GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
          if(gt_shash_is_contained(other_gene_counts, key)){
            gt_gtf_count_(merged_gene_counts, key);
          }
        }GT_SHASH_END_ITERATE;
        GT_SHASH_BEGIN_KEY_ITERATE(local_gene_type_counts, key){
          if(gt_shash_is_contained(other_gene_type_counts, key)){
            gt_gtf_count_(merged_gene_type_counts, key);
          }
        }GT_SHASH_END_ITERATE;
        gt_shash_delete(local_gene_counts, true);
        gt_shash_delete(local_gene_type_counts, true);
        local_gene_counts = merged_gene_counts;
        local_gene_type_counts = merged_gene_type_counts;
      }

//      fprintf(stderr, "HIT: %d\n", hits);
//      if(hits > 1){
//        GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
//          if(gt_shash_is_contained(other_gene_counts, key)){
//            fprintf(stderr, "  SAME GENE HIT %s for ", key);
//            gt_output_map_fprint_map(stderr, map, NULL); fprintf(stderr, "::");gt_output_map_fprint_map(stderr, other_pair, NULL);
//            fprintf(stderr, "\n");
//          }
//        }GT_SHASH_END_ITERATE;
//      }

      gt_shash_delete(other_gene_counts, true);
      gt_shash_delete(other_gene_type_counts, true);
    }

    // count exon in single gene block
    if(num_gene_hits == 1){
      gt_gtf_count_(type_counts, "exon");
    }else{
      gt_gtf_count_(type_counts, "exon_mg");
    }
  }else if(gt_gtf_get_count_(local_type_counts, "intron") > 0){
      gt_gtf_count_(type_counts, "intron");
  }else{
      gt_gtf_count_(type_counts, "unknown");
  }

  // add all gene counts
  GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
    gt_gtf_count_(gene_counts, key);
  }GT_SHASH_END_ITERATE;
  GT_SHASH_BEGIN_KEY_ITERATE(local_gene_type_counts, key){
    gt_gtf_count_(gene_type_counts, key);
  }GT_SHASH_END_ITERATE;

  gt_shash_delete(local_gene_counts, true);
  gt_shash_delete(local_type_counts, true);
  gt_shash_delete(local_gene_type_counts, true);
  gt_vector_delete(hits);
  return hit_by_both_ends;
}


GT_INLINE uint64_t gt_gtf_join_(gt_string* buf, char* base, bool multi_gene, uint64_t blocks){
  if(blocks == 0) return 0;
  uint64_t i = 0;
  uint64_t len = strlen(base);
  for(i=0; i<blocks; i++){
    gt_string_right_append_string(buf, base, len);
    if(multi_gene){
      gt_string_right_append_string(buf, "_mg", 3);
    }
    if(i<blocks-1){
      gt_string_append_char(buf, '^');
    }
  }
  return blocks;
}

GT_INLINE uint64_t gt_gtf_count_map(gt_gtf* const gtf, gt_map* const map, gt_shash* const type_counts, gt_shash* const gene_counts, gt_shash* const gene_type_counts, gt_string* pattern, gt_map* const other_pair){
  // print splitmap blocks
  gt_string_clear(pattern);
  uint64_t blocks = gt_map_get_num_blocks(map);
  if(blocks == 0) return 0;

  gt_shash* const local_type_counts = gt_shash_new();
  gt_shash* const local_gene_counts = gt_shash_new();
  gt_shash* const local_gene_type_counts = gt_shash_new();
  GT_MAP_ITERATE(map, map_block){
    gt_gtf_count_map_(gtf, map_block, local_type_counts, local_gene_counts, local_gene_type_counts, other_pair);
  }
  // count type
  // count types
  uint64_t exons = gt_gtf_get_count_(local_type_counts, "exon") + gt_gtf_get_count_(local_type_counts, "exon_mg");
  uint64_t introns = gt_gtf_get_count_(local_type_counts, "intron");
  uint64_t unknown = gt_gtf_get_count_(local_type_counts, "unknown");
  uint64_t not_annotated = gt_gtf_get_count_(local_type_counts, "na");

  uint64_t num_gene_hits = gt_shash_get_num_elements(local_gene_counts);

  // construct type
  uint64_t total = exons + introns + unknown + not_annotated;
  if(exons > 0){
    total -= gt_gtf_join_(pattern, "exon", num_gene_hits > 1, exons);
    if(total > 0) gt_string_append_char(pattern, '^');
  }
  if(introns > 0){
    total -= gt_gtf_join_(pattern, "intron", false, introns);
    if(total > 0) gt_string_append_char(pattern, '^');
  }
  if(unknown > 0){
    total -= gt_gtf_join_(pattern, "unknown", false, unknown);
    if(total > 0) gt_string_append_char(pattern, '^');
  }
  if(not_annotated > 0){
    total -= gt_gtf_join_(pattern, "na", false, not_annotated);
    if(total > 0) gt_string_append_char(pattern, '^');
  }
  gt_gtf_count_(type_counts, gt_string_get_string(pattern));

  // count genes
  // only count if the split is within a single gene
  if(num_gene_hits == 1){
    //fprintf(stderr, "%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n", exons, introns, unknown, not_annotated);
    // add gene count
    GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
      gt_gtf_count_(gene_counts, key);
      break;
    }GT_SHASH_END_ITERATE;
    GT_SHASH_BEGIN_KEY_ITERATE(local_gene_type_counts, key){
      gt_gtf_count_(gene_type_counts, key);
      break;
    }GT_SHASH_END_ITERATE;
  }

  gt_shash_delete(local_gene_counts, true);
  gt_shash_delete(local_type_counts, true);
  gt_shash_delete(local_gene_type_counts, true);
  return num_gene_hits;
}



GT_INLINE void gt_gtf_count_alignment(gt_gtf* const gtf, gt_alignment* const alignment, gt_shash* const type_count, gt_shash* const gene_counts, gt_shash* const gene_type_counts){
  if(gt_alignment_get_num_maps(alignment) != 1) return; // count only unique
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_string* pattern = gt_string_new(16);
    gt_gtf_count_map(gtf, map,type_count, gene_counts, gene_type_counts, pattern, NULL);
    gt_string_delete(pattern);
  }
}

GT_INLINE void gt_gtf_search_map(const gt_gtf* const gtf, gt_vector* const hits, gt_map* const map, const bool clean_target){
  GT_MAP_ITERATE(map, block){
    uint64_t start = gt_map_get_begin_mapping_position(map);
    uint64_t end   = gt_map_get_end_mapping_position(map);
    gt_gtf_search(gtf, hits, gt_map_get_seq_name(map), start, end, clean_target);
  }
}

GT_INLINE void gt_gtf_search_alignment(const gt_gtf* const gtf, gt_vector* const hits, gt_alignment* const alignment){
  GT_ALIGNMENT_ITERATE(alignment, map){
    gt_gtf_search_map(gtf, hits, map, true);
  }
}

GT_INLINE void gt_gtf_search_template(const gt_gtf* const gtf, gt_vector* const hits, gt_template* const template){
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template, alignment){
    gt_gtf_search_alignment(gtf,hits, alignment);
  }GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_gtf_search_alignment(gtf,hits, gt_template_get_block(template, 0));
  gt_gtf_search_alignment(gtf,hits, gt_template_get_block(template, 1));
}


GT_INLINE uint64_t gt_gtf_count_template(gt_gtf* const gtf, gt_template* const template, gt_shash* const type_counts, gt_shash* const gene_counts, gt_shash* const gene_type_counts, gt_shash* const pair_patterns_counts){
  // process single alignments
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
      gt_gtf_count_alignment(gtf, alignment, type_counts, gene_counts, gene_type_counts);
  } return 0;}

  if(gt_template_get_num_mmaps(template) != 1) return 0;
  // process templates
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attr) {

    gt_string* pattern1 = gt_string_new(16);
    gt_string* pattern2 = gt_string_new(16);
    gt_gtf_count_map(gtf, mmap[0],type_counts, gene_counts, gene_type_counts, pattern1, mmap[1]);
    uint64_t hit_by_both = gt_gtf_count_map(gtf, mmap[1],type_counts, gene_counts, gene_type_counts, pattern2, mmap[0]);

    // append pair pattern
    if(gt_string_cmp(pattern1, pattern2) < 0){
      gt_string_append_char(pattern2, '|');
      gt_string_append_eos(pattern2);
      gt_string_left_append_gt_string(pattern1, pattern2);
    }else{
      gt_string_append_char(pattern1, '|');
      gt_string_append_eos(pattern1);
      gt_string_right_append_gt_string(pattern1, pattern2);
    }
    gt_gtf_count_(pair_patterns_counts, pattern1->buffer);
    gt_string_delete(pattern1);
    gt_string_delete(pattern2);
    return hit_by_both;
  }
  return 0;
}
