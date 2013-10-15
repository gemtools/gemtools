#include "uthash.h"

#define MAX_READ_LENGTH 16384 // For sanity checking
#define QUAL_FASTQ 33
#define QUAL_SOLEXA 64
#define DEFAULT_QUAL_OFFSET 33 // FASTQ quality convention
#define MAX_QUAL 64
#define DEFAULT_MAX_INSERT 1000
#define MAX_INDEL_SIZE 1024
#define DUP_LIMIT 255

#define PAIR_TYPE_DS 0 // Reads on different strands, same contig
#define PAIR_TYPE_SS 1 // Reads on same strand, same contig
#define PAIR_TYPE_MM 2 // Reads on different contigs

#define SG_C2T 1
#define SG_G2A 2
#define SG_LAMBDA 3
#define SG_PHIX174 4
#define SG_MITO 5

#define AS_INSERT_TYPE_PAIRED 0        // Paired mapping declared as valid and unique
#define AS_INSERT_TYPE_ALL_UNIQUE 1    // All pairs of uniquely mapping single reads with consistent orientations
#define AS_INSERT_TYPE_RECOVERED 2     // All pairs (uniquely mapping) where at least 1 end was multi-mapped
#define AS_INSERT_TYPE_SPLIT 3         // All pairs (uniquely mapping) where at least 1 end was split-mapped

#define PHAGE_LAMBDA "NC_001416.1"
#define PHIX174 "NC_001422.1"
#define MITO "chrM"

#define CONV_EST_TYPE_PHIX174 0
#define CONV_EST_TYPE_PHAGE_LAMBDA 1
#define CONV_EST_TYPE_MITO 2
#define CONV_EST_TYPE_NONCPG 3
#define CONV_EST_N 4

typedef struct {
  int64_t x;
  uint64_t ct[4];
  UT_hash_handle hh;
} dist_element;

typedef struct {
  u_int16_t loc;
  u_int16_t tile;
  int16_t  dist;
} loc_elem;

#define INIT_LB_SIZE 128

typedef struct {
  size_t n_elem;
  size_t size;
  uint32_t x;
  loc_elem *elem;
  pthread_mutex_t mutex;
  UT_hash_handle hh;
} loc_block;

typedef struct {
  char *ctg;
  loc_block *lblock;
  pthread_rwlock_t rwlock;
  UT_hash_handle hh;
} loc_hash;

#define ID_END_CHAR 127
#define ID_COLON_CHAR 1
#define ID_HASH_CHAR 2
#define ID_SPACE_CHAR 3
#define ID_SLASH_CHAR 4

#define MAX_LANE_ID 8

#define ID_TAG_OK 0
#define ID_TAG_ERROR_BAD_LANE 1
#define ID_TAG_ERROR_BAD_TILE 2
#define ID_TAG_ERROR_BAD_COORD 3

typedef struct {
	gt_string *instrument_name;
	gt_string *flowcell;
	gt_string *run;
	gt_string *index;
	int lane;
	uint32_t tile;
	uint32_t x;
	uint32_t y;
	int pair_id;
	bool filter;
	uint32_t flags;
} id_tag;

typedef struct {
  uint64_t nreads; // single or paired end
  uint64_t yield[2]; // Total number of non N bases
  uint64_t mapped[2]; // Number of mapped reads for each end
  uint64_t unique[2]; // Uniquely mapping reads (single alignment, none in next strata
  uint64_t ambiguous[2]; // Ambiguously mapping reads (single alignment, none in next strata
  uint64_t paired_mapped; // Number of paired reads
  uint64_t paired_unique;
  uint64_t paired_type[3]; // Strand/contig concordance 
  uint64_t sg_stats[6]; // Special genome counts
  uint64_t reads_with_splitmaps[3];
  uint64_t reads_only_with_splitmaps[3];
  uint64_t max_read_length[2];
  uint64_t curr_read_store[2];
  uint64_t *read_length_stats[2];
  uint64_t *indel_stats[4];  // [rd*2+x][cycle]  x=0 INS, x=1 DEL
  uint64_t *indel_length[4]; // [rd*2+x][indel size]  For insertions and deletions
  uint64_t max_indel_length;
  uint64_t **base_counts_by_cycle[2];      // [read][cycle][qual*5+base]
  uint64_t *mm_stats[2*(MAX_QUAL+1)]; // [read*(MAX_QUAL+1)+qual][cycle]
  uint64_t *qual_stats[2*(MAX_QUAL+1)]; // [read*(MAX_QUAL+1)+qual][cycle]
  uint64_t bs_counts[2][30];
  uint64_t *non_cpg_cytosines[2];
  uint64_t tv_stats[2][MAX_QUAL+1];
  uint64_t ts_stats[2][MAX_QUAL+1];
  uint64_t pbc_stats[2][MAX_QUAL+1];
  uint64_t pbn_stats[2][MAX_QUAL+1];
  uint64_t duplicate_reads_used;
  uint64_t unique_fragment_estimate;
  uint64_t duplicate_counts[5][DUP_LIMIT+1];
  double duplicate_rate[2]; // Overall, optical duplicate fractions
  dist_element* insert_size; // Store insert size distribution
  loc_hash** loc_hash; // Track position and insert sizes to estimate duplicates
  bool paired;
  bool bisulphite;
} as_stats;

typedef struct {
  char *input_files[2];
  char *output_file;
  char *dist_file;
  char *phage_lambda;
  char *phix174;
  char *mito;
  bool mmap_input;
  gt_generic_parser_attributes *parser_attr;
  bool variable_read_length;
  bool ignore_id;
  gt_output_file_compression compress;
  uint64_t read_length[2]; // Untrimmed read length for both reads
  uint64_t max_read_length; // For sanity checking
  uint64_t min_insert;
  uint64_t max_insert;
  int num_threads;
  int qual_offset; // quality offset (33 for FASTQ, 64 for Illumina)
  as_stats **stats;
} as_param;

