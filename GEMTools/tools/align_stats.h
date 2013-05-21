#include "uthash.h"

#define MAX_READ_LENGTH 16384 // For sanity checking
#define QUAL_FASTQ 33
#define QUAL_SOLEXA 64
#define DEFAULT_QUAL_OFFSET 33 // FASTQ quality convention
#define MAX_QUAL 64
#define DEFAULT_MAX_INSERT 1000

typedef struct {
  char *input_files[2];
  char *output_file;
  char *dist_file;
  bool mmap_input;
  gt_generic_parser_attributes* parser_attr;
  bool paired_read;
  bool variable_read_length;
  bool ignore_id;
  uint64_t read_length[2]; // Untrimmed read length for both reads
  uint64_t max_read_length; // For sanity checking
  uint64_t min_insert; 
  uint64_t max_insert; 
  int num_threads;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina) 
} as_param;

#define PAIR_TYPE_DS 0 // Reads on different strands, same contig
#define PAIR_TYPE_SS 1 // Reads on same strand, same contig
#define PAIR_TYPE_MM 2 // Reads on different contigs

#define BIS_C2T 0
#define BIS_G2A 1
#define BIS_LAMBDA 2

typedef struct {
  int32_t x;
  uint64_t ct[2];
  UT_hash_handle hh;
} dist_element;

typedef struct {
  u_int16_t loc;
  u_int16_t tile;
  int16_t  dist;
} loc_elem;

#define INIT_LB_SIZE 512

typedef struct {
  size_t n_elem;
  size_t size;
  uint32_t x;
  loc_elem *elem;
  UT_hash_handle hh;
} loc_block;

typedef struct {
  char *ctg;
  loc_block *lblock;
  UT_hash_handle hh;
} loc_hash;

typedef struct {
  uint64_t nreads; // single or paired end
  uint64_t yield[2]; // Total number of non N bases
  uint64_t mapped[2]; // Number of mapped reads for each end
  uint64_t unique[2]; // Uniquely mapping reads (single alignment, none in next strata
  uint64_t ambiguous[2]; // Ambiguously mapping reads (single alignment, none in next strata
  uint64_t paired_mapped; // Number of paired reads
  uint64_t paired_unique;
  uint64_t paired_type[3]; // Strand/contig concordance 
  uint64_t bis_stats[3]; // Bisulphite mapping counts
  uint64_t reads_with_splitmaps[3];
  uint64_t reads_only_with_splitmaps[3];
  uint64_t max_read_length[2];
  uint64_t curr_read_store[2];
  uint64_t *read_length_stats[2];
  uint64_t indel_stats[4];  // [rd*2+x]  x=0 INS, x=1 DEL
  uint64_t *indel_length[4]; // For insertions and deletions
  uint64_t max_indel_length;
  uint64_t **base_counts_by_cycle[2];      // [read][cycle][qual*5+base]
  uint64_t mm_stats[2][MAX_QUAL+1];
  uint64_t tv_stats[2][MAX_QUAL+1];
  uint64_t ts_stats[2][MAX_QUAL+1];
  uint64_t pbc_stats[2][MAX_QUAL+1];
  uint64_t pbn_stats[2][MAX_QUAL+1];
  dist_element* insert_size; // Store insert size distribution
  loc_hash* loc_hash; // Track position and insert sizes to estimate duplicates
  bool paired;
} as_stats;
