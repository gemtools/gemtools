#include "gt_gtf.h"

GT_INLINE gt_gtf_entry* gt_gtf_entry_new(const uint64_t start, const uint64_t end, const gt_strand strand, gt_string* const type){
  gt_gtf_entry* entry = malloc(sizeof(gt_gtf_entry));
  entry->uid = 0;
  entry->start = start;
  entry->end = end;
  entry->type = type;
  entry->strand = strand;
  entry->gene_type = NULL;
  entry->gene_id = NULL;
  entry->transcript_id = NULL;
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
  gtf->genes = gt_shash_new();
  return gtf;
}

GT_INLINE void gt_gtf_delete(gt_gtf* const gtf){
  gt_shash_delete(gtf->refs, true);
  gt_shash_delete(gtf->types, true);
  gt_shash_delete(gtf->gene_ids, true);
  gt_shash_delete(gtf->transcript_ids, true);
  gt_shash_delete(gtf->gene_types, true);
  gt_shash_delete(gtf->genes, false);
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

GT_INLINE gt_gtf_entry* gt_gtf_get_gene_by_id(gt_gtf* const gtf, char* const key){
  if(gt_shash_is_contained(gtf->genes, key)){
    return gt_shash_get_element(gtf->genes, key);
  }
  return NULL;
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
  if(strcmp(e->type->buffer, "gene") == 0){
    gt_shash_insert(gtf->genes, e->gene_id->buffer, e, gt_gtf_entry*);
  }
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
  if(e->type != NULL){
    printf("%s : %"PRIu64" - %"PRIu64" (%c)", e->type->buffer, e->start, e->end, (e->strand==FORWARD?'+':'-') );

  }
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

      if(strcmp(e->type->buffer, "gene") == 0 && !gt_shash_is_contained(target->genes, e->gene_id->buffer)){
        gt_shash_insert(target->genes, e->gene_id->buffer, e, gt_gtf_entry*);
      }

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
GT_INLINE void gt_gtf_count_sum_(gt_shash* const table, char* const element, uint64_t value){
  if(!gt_shash_is_contained(table, element)){
    uint64_t* v = gt_malloc_uint64();
    *v = value;
    gt_shash_insert(table, element, v, uint64_t);
  }else{
    uint64_t* v = gt_shash_get(table,element,uint64_t);
    *v += value;
  }
}

GT_INLINE void gt_gtf_count_weight_(gt_shash* const table, char* const element, double weight){
  if(!gt_shash_is_contained(table, element)){
    double* v = malloc(sizeof(double*));
    *v = weight;
    gt_shash_insert(table, element, v, double);
  }else{
    double* v = gt_shash_get(table,element,double);
    *v += weight;
  }
}

GT_INLINE uint64_t gt_gtf_get_count_(gt_shash* const table, char* const element){
  if(!gt_shash_is_contained(table, element)){
    return 0;
  }
  uint64_t* v = gt_shash_get(table,element,uint64_t);
  return *v;
}

GT_INLINE float gt_gtf_get_count_weight(gt_shash* const table, char* const element){
  if(!gt_shash_is_contained(table, element)){
    return 0.0;
  }
  double* v = gt_shash_get(table,element,double);
  return *v;
}

GT_INLINE void gt_gtf_create_hit(gt_vector* search_hits, gt_shash* all_genes, gt_gtf_hits* hits, gt_gtf_hit* template_hit){
  template_hit->transcripts = gt_shash_new();
  template_hit->genes = gt_shash_new();
  template_hit->is_protein_coding = false;
  template_hit->hits_exon = false;
  bool counted_protein = false;
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
      template_hit->hits_exon = true;
    }
    if(e->gene_type != NULL){
      template_hit->is_protein_coding |= (strcmp(e->gene_type->buffer, "protein_coding") == 0);
      if(!counted_protein){
        hits->num_protein_coding++;
        counted_protein = true;
      }
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

/**
 * This counts a single continuous block and takes the. Note that we do not perform any checks on
 * splits/pairs here and simply count for this single continuous map
 *
 * @param gt_gtf* gtf                the gtf reference
 * @param gt_map*                    continuous map block
 * @param gt_shash* type_counts      the type counts, i.e exon/intron etc
 * @param gt_shash* gene_counts      the gene counts with the gene_id's hit by the map.
 * @return uint64_t num_gene_exons   number of unique gene_ids hit by exons
 */
GT_INLINE uint64_t gt_gtf_count_map_(gt_gtf* const gtf, gt_map* const map,
                                     gt_shash* const type_counts,
                                     gt_shash* const gene_counts){
  // get coordinates
  uint64_t start = gt_map_get_begin_mapping_position(map);
  uint64_t end   = gt_map_get_end_mapping_position(map);
  // store the search hits and search
  gt_vector* const hits = gt_vector_new(32, sizeof(gt_gtf_entry*));
  gt_gtf_search(gtf, hits, gt_map_get_seq_name(map), start, end, true);

  // we do a complete local count for this block
  // and then merge the local count with the global count
  // to be able to resolve genes/gene_types that are
  // through wither the pair information or split information,
  // assuming that the counts for the other pair and/or the other split
  // are already contained in the globally presented count maps
  gt_shash* const local_type_counts = gt_shash_new();
  gt_shash* local_gene_counts = gt_shash_new();
  gt_shash* local_exon_gene_counts = gt_shash_new();
  GT_VECTOR_ITERATE(hits, e, i, gt_gtf_entry*){
    gt_gtf_entry* hit = *e;
    // count type
    gt_gtf_count_(local_type_counts, gt_string_get_string(hit->type));
    // count gene id
    if(hit->gene_id != NULL) gt_gtf_count_(local_gene_counts, gt_string_get_string(hit->gene_id));
    // count gene_id from exons
    if(hit->type != NULL && hit->gene_id != NULL && strcmp("exon", hit->type->buffer) == 0){
      gt_gtf_count_(local_exon_gene_counts, gt_string_get_string(hit->gene_id));
    }
  }
  uint64_t num_gene_hit_exons = gt_shash_get_num_elements(local_exon_gene_counts);
  // count types and merge them with the global
  // counts. NOTE that the order matters here, so
  // we:
  //   1. check for NA hits where nothing is found
  //   2. count exon hits
  //   3. count intron hits
  //   4. count unknown if the hit was neither an intron nor exon hit
  // all counting steps are exclusive, thats why the order matters!
  if(gt_vector_get_used(hits) == 0){
    // count 'NA' type if we did not hit anything
    gt_gtf_count_(type_counts, GT_GTF_TYPE_NA);
  }else if(gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_EXON) > 0){
    gt_gtf_count_(type_counts, GT_GTF_TYPE_EXON);
  }else if(gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_INTRON) > 0){
    gt_gtf_count_(type_counts, GT_GTF_TYPE_INTRON);
  }else{
    gt_gtf_count_(type_counts, GT_GTF_TYPE_UNKNOWN);
  }

  // make gene counts based on exon hits if we found at least one
  if(num_gene_hit_exons > 0){
    GT_SHASH_BEGIN_KEY_ITERATE(local_exon_gene_counts, key){
      gt_gtf_count_(gene_counts, key);
    }GT_SHASH_END_ITERATE;
  }else{
    // add all gene counts
    GT_SHASH_BEGIN_KEY_ITERATE(local_gene_counts, key){
      gt_gtf_count_(gene_counts, key);
    }GT_SHASH_END_ITERATE;
  }


  gt_shash_delete(local_gene_counts, true);
  gt_shash_delete(local_type_counts, true);
  gt_shash_delete(local_exon_gene_counts, true);
  gt_vector_delete(hits);
  return num_gene_hit_exons;
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

GT_INLINE double gt_gtf_count_get_sum_(gt_shash* table){
  double v = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(table, value, uint64_t){
    v += *value;
  }GT_SHASH_END_ITERATE;
  return v;
}
/**
 * Count a map. This respects split maps and unifies gene_id's based on the
 * the split. If the both sides of the split match multiple gene_ids but there is
 * a common gene_id on both side, only that id is counted. Otherwise a count is set
 * for all gene_ids.
 * In addition to the counts, if a pattern string is given, it is filled with the type
 * pattern with respect to split maps. For example:
 *
 *   exon                        -> exon
 *   exon and intron (split map) -> exon^intron
 *   exon in multiple genes      -> exon_mg
 *
 * The function returns the number of gene_ids hit by the map.
 *
 * The first map has to be specified, but the second one is options. If it is set,
 * the second map block is also checked and counted.
 *
 *
 * @param gt_gtf* gtf              the gtf reference
 * @param gt_map* map1             the first map
 * @param gt_map* map2             the scond map
 * @param gt_shash* type_counts    the type counts
 * @param gt_shash* gene_counts    the gene counts
 * @param gt_string pattern        the pattern string filled based on the types
 * @return uint64_t num_gene_hits  the number of gene_ids hit by the map
 */
GT_INLINE uint64_t gt_gtf_count_map(gt_gtf* const gtf, gt_map* const map1, gt_map* const map2,
                                    gt_shash* const pattern_counts, gt_shash* const gene_counts,
                                    gt_string* pattern){
  // clear patterns
  gt_string_clear(pattern);
  // get number of blocks and ensure we have at least one
  uint64_t blocks = gt_map_get_num_blocks(map1);
  if(map2 != NULL){
    blocks += gt_map_get_num_blocks(map2);
  }

  if(blocks == 0) return 0;

  // local counts for all blocks
  // and store the number of multi gene exon hits for each block
  // in addition we create the base pattern per block here
  gt_shash* const local_type_counts = gt_shash_new();
  gt_shash* local_gene_counts = gt_shash_new();
  gt_shash* local_gene_counts_1 = gt_shash_new();
  gt_shash* local_gene_counts_2 = gt_shash_new();
  uint64_t* const local_exon_gene_hits = malloc(blocks * sizeof(uint64_t));
  gt_vector* const local_type_patterns = gt_vector_new(2, sizeof(char*));
  uint64_t exons, introns, unknown, not_annotated;
  exons = introns = unknown = not_annotated = 0;
  uint64_t i = 0;
  GT_MAP_ITERATE(map1, map_block){
    local_exon_gene_hits[i++] = gt_gtf_count_map_(gtf, map_block, local_type_counts, local_gene_counts_1);
    uint64_t _exons = exons + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_EXON);
    uint64_t _introns = introns + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_INTRON);
    uint64_t _unknown = unknown + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_UNKNOWN);
    uint64_t _not_annotated = not_annotated + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_NA);
    // add the pattern string based in the count value that changed
    if(_exons > exons) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_EXON, char*);
    if(_introns > introns) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_INTRON, char*);
    if(_unknown > unknown) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_UNKNOWN, char*);
    if(_not_annotated > not_annotated) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_NA, char*);
    exons = _exons;
    introns = _introns;
    unknown = _unknown;
    not_annotated = _not_annotated;
  }
  // unify the gene counts based on the number of blocks.
  // the gene_counts are reduced to either the ones that are found in
  // all blocks or they are kept as they are
  if(gt_shash_get_num_elements(local_gene_counts_1) > 1){
    gt_shash* merged_counts = gt_shash_new();
    uint64_t blocks1 = gt_map_get_num_blocks(map1);
    GT_SHASH_BEGIN_ITERATE(local_gene_counts_1, gene_id, count, uint64_t){
      if(*count == blocks1){
        gt_gtf_count_sum_(merged_counts, gene_id, blocks1);
      }
    }GT_SHASH_END_ITERATE;
    // if we found some unique ids that are covered by both
    // we flip over to the merged counts
    if(gt_shash_get_num_elements(merged_counts) > 0){
      gt_shash_delete(local_gene_counts_1, true);
      local_gene_counts_1 = merged_counts;
      // we fliped so we reset the exon gene hit counts to ones as well
      for(i=0;i<blocks1;i++){local_exon_gene_hits[i] = 1;}
    }
  }

  if(map2 != NULL){
    GT_MAP_ITERATE(map2, map_block){
      local_exon_gene_hits[i++] = gt_gtf_count_map_(gtf, map_block, local_type_counts, local_gene_counts_2);
      uint64_t _exons = exons + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_EXON);
      uint64_t _introns = introns + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_INTRON);
      uint64_t _unknown = unknown + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_UNKNOWN);
      uint64_t _not_annotated = not_annotated + gt_gtf_get_count_(local_type_counts, GT_GTF_TYPE_NA);
      // add the pattern string based in the count value that changed
      if(_exons > exons) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_EXON, char*);
      if(_introns > introns) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_INTRON, char*);
      if(_unknown > unknown) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_UNKNOWN, char*);
      if(_not_annotated > not_annotated) gt_vector_insert(local_type_patterns, GT_GTF_TYPE_NA, char*);
      exons = _exons;
      introns = _introns;
      unknown = _unknown;
      not_annotated = _not_annotated;
    }
    // unify the gene counts based on the number of blocks.
    // the gene_counts are reduced to either the ones that are found in
    // all blocks or they are kept as they are
    if(gt_shash_get_num_elements(local_gene_counts_2) > 1){
      gt_shash* merged_counts = gt_shash_new();
      uint64_t blocks2 = gt_map_get_num_blocks(map2);
      GT_SHASH_BEGIN_ITERATE(local_gene_counts_2, gene_id, count, uint64_t){
        if(*count == blocks2){
          gt_gtf_count_sum_(merged_counts, gene_id, blocks2);
        }
      }GT_SHASH_END_ITERATE;
      // if we found some unique ids that are covered by both
      // we flip over to the merged counts
      if(gt_shash_get_num_elements(merged_counts) > 0){
        gt_shash_delete(local_gene_counts_2, true);
        local_gene_counts_2 = merged_counts;
        uint64_t blocks1 = gt_map_get_num_blocks(map1);
        // we flipped so we reset the exon gene hit counts to ones as well
        for(i=blocks1;i<(blocks1+blocks2);i++){local_exon_gene_hits[i] = 1;}
      }
    }
  }

  /**
   * Merge everything into a single map and reduce the map if we can resove through pairs
   */
  gt_shash* merged_counts = gt_shash_new();
  uint64_t blocks2 = 0;
  if(map2 != NULL){
    blocks2 = gt_map_get_num_blocks(map2);
  }
  uint64_t blocks1 = gt_map_get_num_blocks(map1);
  GT_SHASH_BEGIN_ITERATE(local_gene_counts_1, gene_id, count, uint64_t){
      gt_gtf_count_sum_(merged_counts, gene_id, *count);
  }GT_SHASH_END_ITERATE;
  if(map2 != NULL){
    GT_SHASH_BEGIN_ITERATE(local_gene_counts_2, gene_id, count, uint64_t){
      gt_gtf_count_sum_(merged_counts, gene_id, *count);
    }GT_SHASH_END_ITERATE;
  }
  uint64_t unique_genes_between_pairs = 0;
  GT_SHASH_BEGIN_ELEMENT_ITERATE(merged_counts, count, uint64_t){
    if(*count == blocks){
      unique_genes_between_pairs++;
    }
  }GT_SHASH_END_ITERATE;

  // we found unique genes through the pair, so we can use
  // the merged map to do the final counts
  if(unique_genes_between_pairs == 1){
    // we flipped so we reset the exon gene hit counts to ones as well
    for(i=0;i<(blocks1+blocks2);i++){local_exon_gene_hits[i] = 1;}

    // merge the gene counts weighted to a single map
    double gene_counts = (double)unique_genes_between_pairs;
    GT_SHASH_BEGIN_ITERATE(merged_counts, gene_id, count, uint64_t){
      if(*count == blocks){
        gt_gtf_count_weight_(local_gene_counts, gene_id, ((double)*count)/gene_counts);
      }
    }GT_SHASH_END_ITERATE;

    // cleanup
    gt_shash_delete(local_gene_counts_1, true);
    gt_shash_delete(local_gene_counts_2, true);
    gt_shash_delete(merged_counts, true);
  }else{
    // merge into the global counts using the two pair maps separately
    // merge the gene counts weighted to a single map
    double gene_counts_1 = gt_gtf_count_get_sum_(local_gene_counts_1);
    GT_SHASH_BEGIN_ITERATE(local_gene_counts_1, gene_id, count, uint64_t){
      gt_gtf_count_weight_(local_gene_counts, gene_id, ((double)*count)/gene_counts_1);
    }GT_SHASH_END_ITERATE;
    gt_shash_delete(local_gene_counts_1, true);
    if(map2 != NULL){
      double gene_counts_2 = gt_gtf_count_get_sum_(local_gene_counts_2);
      GT_SHASH_BEGIN_ITERATE(local_gene_counts_2, gene_id, count, uint64_t){
        gt_gtf_count_weight_(local_gene_counts, gene_id, ((double)*count)/gene_counts_2);
      }GT_SHASH_END_ITERATE;
      gt_shash_delete(local_gene_counts_2, true);
    }
  }

  // get the number of hits of this map
  uint64_t num_gene_hits = gt_shash_get_num_elements(local_gene_counts);
  uint64_t map1_blocks = gt_map_get_num_blocks(map1);
  if(pattern_counts != NULL){
    // now iterate the blocks and construct final pattern
    for(i=0; i<blocks; i++){
      char* p = *(gt_vector_get_elm(local_type_patterns, i, char*));
      // for exons check that in case we have a single gene hit, its exons, in case of a multi-gene hit, append _mg if
      // the multi gene hit comes from the current block
      gt_gtf_join_(pattern, p, (strcmp("exon",p) == 0) ? ((num_gene_hits == 1) ? false : (local_exon_gene_hits[i] > 1)) : false, 1);
      // add paired end spacer
      if(map2 != NULL && i == (map1_blocks-1)){
        gt_string_append_char(pattern, '|');
      }else{
        if(i<blocks-1){
          gt_string_append_char(pattern, '^');
        }
      }
    }
    // count global type based on the constructed pattern
    gt_gtf_count_(pattern_counts, gt_string_get_string(pattern));
//    if(strcmp(pattern->buffer, "exon|exon_mg") == 0){
//      gt_output_map_fprint_map(stdout, map1, NULL); printf("::");gt_output_map_fprint_map(stdout, map2, NULL);printf("\n");
//    }
  }

  if(gene_counts != NULL){
    // count the gene ids
    GT_SHASH_BEGIN_ITERATE(local_gene_counts, key, e, double){
      gt_gtf_count_weight_(gene_counts, key, *e);
    }GT_SHASH_END_ITERATE;
  }

  // cleanup
  gt_vector_delete(local_type_patterns);
  gt_shash_delete(local_gene_counts, true);
  gt_shash_delete(local_type_counts, true);
  free(local_exon_gene_hits);
  return num_gene_hits;
}

GT_INLINE uint64_t gt_gtf_count_alignment(gt_gtf* const gtf, gt_alignment* const alignment, gt_shash* const pattern_count, gt_shash* const gene_counts){
  uint64_t hits = 0;
  gt_string* pattern = gt_string_new(16);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    hits += gt_gtf_count_map(gtf, map, NULL, pattern_count, gene_counts, pattern);
    gt_string_clear(pattern);
  }
  gt_string_delete(pattern);
  return hits;
}

GT_INLINE uint64_t gt_gtf_count_template(gt_gtf* const gtf, gt_template* const template, gt_shash* const pattern_count, gt_shash* const gene_counts){
  uint64_t hits = 0;
  gt_string* pattern = gt_string_new(16);
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attr) {
    hits += gt_gtf_count_map(gtf, mmap[0], mmap[1], pattern_count, gene_counts, pattern);
    gt_string_clear(pattern);
  }
  gt_string_delete(pattern);
  return hits;
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
