/*
 * Module to load and manage GTF annotation
 */

#ifndef GT_GTF_H_
#define GT_GTF_H_

#include "gt_commons.h"
#include "gt_string.h"
#include "gt_vector.h"
#include "gt_shash.h"
#include "gt_template.h"
#include "gt_output_map.h"
#include "gt_input_map_parser.h"
#include <omp.h>


#define GTF_DEFAULT_ENTRIES 1000
#define GTF_MAX_LINE_LENGTH 2048

#define GT_GTF_INVALID_LINE 10
#define GT_GTF_INVALID_ATTRIBUTES 20

#define GT_GTF_IS_EOL(text_line) gt_expect_false((*text_line)==EOL || (*text_line)==EOS)
#define GT_GTF_NEXT_CHAR(text_line) ++(text_line)
#define GT_GTF_READ_UNTIL(text_line,test) \
  while (gt_expect_true(!(test) && !GT_GTF_IS_EOL(text_line))) { \
    GT_GTF_NEXT_CHAR(text_line); \
  }


/*
 * Single gtf entry
 */
typedef struct {
  uint64_t uid;
	uint64_t start; // the start position
	uint64_t end; // the end position
	gt_strand strand; // the strand
	gt_string* type; // the type, i.e. exon, gene
	gt_string* gene_id; // the gene id if it exists
	gt_string* transcript_id; // the transcript id if it exists
	gt_string* gene_type; // the gene id if it exists
} gt_gtf_entry;

typedef struct _gt_gtf_node gt_gtf_node;

struct _gt_gtf_node {
  uint64_t midpoint;
  gt_vector* entries_by_start;
  gt_vector* entries_by_end;
  gt_gtf_node* left;
  gt_gtf_node* right;
} ;

/**
 * Single chromosome reference with
 * all it gtf entries
 */
typedef struct {
	gt_vector* entries; // gt_gtf_entry list
	gt_gtf_node* node; // gt_gtf_entry list
} gt_gtf_ref;

/**
 * GTF file with a map to the chromosome references
 * and all available types (exon, gene, ...).
 */
typedef struct {
	gt_shash* refs; // maps from the ref name to the ref char* -> gt_gtf_ref*
	gt_shash* types; // maps from the type name to the gt_string type ref char* -> gt_string*
	gt_shash* gene_ids; // maps from char* to gt_string* for gene_ids char* -> gt_string*
	gt_shash* transcript_ids; // maps from char* to gt_string* for gene_ids char* -> gt_string*
	gt_shash* gene_types; // maps from char* to gt_string* for gene_types char* -> gt_string*
	gt_shash* genes; // maps from char* to gt_gtf_entry for genes
}gt_gtf;

/**
 * gtf hit that are filled by the template search methods
 */
typedef struct {
  gt_vector* exon_hits; // contains an entry for each map or map pair with the hit target id (gene_id)
  uint64_t num_genes; // number of hit genes
  uint64_t num_paired_genes; // number of hit paired genes
  uint64_t num_protein_coding; // number of protein coding exons
}gt_gtf_hits;


typedef struct {
  gt_map* map;
  gt_map** mmap;
  gt_mmap_attributes* map_attributes;
  gt_shash* transcripts;
  gt_shash* genes;
  float exon_overlap;
  float junction_hits;
  uint64_t num_junctions;
  uint64_t intron_length;
  uint64_t num_template_blocks;
  bool is_protein_coding;
  bool pairs_transcript;
  bool pairs_splits;
  bool pairs_gene;
  bool hits_exon;
}gt_gtf_hit;

GT_INLINE gt_gtf_hit* gt_gtf_hit_new(void);
GT_INLINE void gt_gtf_hit_delete(gt_gtf_hit* hit);

/**
 * Create and delete gtf
 */
GT_INLINE gt_gtf* gt_gtf_new(void);
GT_INLINE void gt_gtf_delete(gt_gtf* const gtf);

/**
 * Create and delete hits
 */
GT_INLINE gt_gtf_hits* gt_gtf_hits_new(void);
GT_INLINE void gt_gtf_hits_delete(gt_gtf_hits* const hits);
GT_INLINE void gt_gtf_hits_clear(gt_gtf_hits* const hits);

/**
 * Create and delete new gtf_entries
 */
GT_INLINE gt_gtf_entry* gt_gtf_entry_new(const uint64_t start, const uint64_t end, const gt_strand strand, gt_string* const type);
GT_INLINE void gt_gtf_entry_delete(gt_gtf_entry* const entry);

/**
 * Create and delete new chromosome refs
 */
GT_INLINE gt_gtf_ref* gt_gtf_ref_new(void);
GT_INLINE void gt_gtf_ref_delete(gt_gtf_ref* const ref);

/**
 * Parse GTF files and return a new gt_gtf*. The ref entries
 * will be sorted by star,end,type
 */
GT_INLINE gt_gtf* gt_gtf_read(gt_input_file* input, uint64_t threads);
GT_INLINE gt_gtf* gt_gtf_read_from_stream(FILE* input, uint64_t threads);
GT_INLINE gt_gtf* gt_gtf_read_from_file(char* input, uint64_t threads);

/**
 * Access the chromosome refs
 */
GT_INLINE gt_gtf_ref* gt_gtf_get_ref(const gt_gtf* const gtf, char* const name);
GT_INLINE bool gt_gtf_contains_ref(const gt_gtf* const gtf, char* const name);

/**
 * Access available types
 */
GT_INLINE gt_string* gt_gtf_get_type(const gt_gtf* const gtf, char* const type);
GT_INLINE bool gt_gtf_contains_type(const gt_gtf* const gtf, char* const name);

/**
 * Access available gene_ids
 */
GT_INLINE gt_string* gt_gtf_get_gene_id(const gt_gtf* const gtf, char* const name);
GT_INLINE bool gt_gtf_contains_gene_id(const gt_gtf* const gtf, char* const name);

/**
 * Access available gene_ids
 */
GT_INLINE gt_string* gt_gtf_get_transcript_id(const gt_gtf* const gtf, char* const name);
GT_INLINE bool gt_gtf_contains_transcript_id(const gt_gtf* const gtf, char* const name);

/**
 * Access available gene_types
 */
GT_INLINE gt_string* gt_gtf_get_gene_type(const gt_gtf* const gtf, char* const name);
GT_INLINE bool gt_gtf_contains_gene_type(const gt_gtf* const gtf, char* const name);

/**
 * Get gt_gtf_entry for gene by its gene_id
 */
GT_INLINE gt_gtf_entry* gt_gtf_get_gene_by_id(gt_gtf* const gtf, char* const key);

GT_INLINE void gt_gtf_count_(gt_shash* const table, char* const element);
GT_INLINE void gt_gtf_count_weight_(gt_shash* const table, char* const element, double weight);


/**
 * Search
 */
/**
 * Search for annotation that overlap with the specified region. The matching entries will be added to the target
 * vector. Note that the target vector is cleared at the beginning of the method!
 */
GT_INLINE uint64_t gt_gtf_search(const gt_gtf* const gtf, gt_vector* const target, char* const ref, const uint64_t start, const uint64_t end, const bool clean_target);
/**
 * Search for exons that overlap with the given template mappings.
 */
GT_INLINE void gt_gtf_search_template_hits(const gt_gtf* const gtf, gt_gtf_hits* const hits, gt_template* const template_src);
GT_INLINE void gt_gtf_search_alignment_hits(const gt_gtf* const gtf, gt_gtf_hits* const hits, gt_alignment* const template_src);

GT_INLINE uint64_t gt_gtf_count_alignment(gt_gtf* const gtf, gt_alignment* const alignment, gt_shash* const type_count, gt_shash* const gene_counts);
GT_INLINE uint64_t gt_gtf_count_template(gt_gtf* const gtf, gt_template* const template, gt_shash* const type_counts, gt_shash* const gene_counts);



/**
 * General searches
 */
GT_INLINE void gt_gtf_search_map(const gt_gtf* const gtf, gt_vector* const hits, gt_map* const map, const bool clean_target);
GT_INLINE void gt_gtf_search_alignment(const gt_gtf* const gtf, gt_vector* const hits, gt_alignment* const alignment);
GT_INLINE void gt_gtf_search_template(const gt_gtf* const gtf, gt_vector* const hits, gt_template* const template);


/*MISC helpers*/
void gt_gtf_print_entry_(gt_gtf_entry* e, gt_map* map);

#endif /* GT_GTF_H_ */
