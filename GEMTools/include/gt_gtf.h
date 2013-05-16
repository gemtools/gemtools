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

/*
 * GTF strand for x/-/.
 */
typedef enum { GTF_STRAND_FORWARD, GTF_STRAND_REVERSE, GTF_STRAND_UNKNOWN } gt_gtf_strand;

/*
 * GTF entry
 */
typedef struct {
	uint64_t start;
	uint64_t end;
	gt_gtf_strand strand;
	gt_string* type;
	gt_string* gene_id;
} gt_gtf_entry;

typedef struct {
	gt_vector* entries; // the entries of this reference
} gt_gtf_ref;

typedef struct {
	gt_shash* refs; // maps from the ref name to the ref
	gt_shash* types; // maps from the type name to the gt_string type ref	
}gt_gtf;

GT_INLINE bool gt_gtf_search_exact(gt_gtf* gtf, char* const  ref, const uint64_t start, const uint64_t end, char* const type, const gt_gtf_strand strand);

GT_INLINE gt_gtf_entry* gt_gtf_entry_new(const uint64_t start, const uint64_t end, const gt_gtf_strand strand, gt_string* type);
GT_INLINE void gt_gtf_entry_delete(gt_gtf_entry* const entry);

GT_INLINE gt_gtf_ref* gt_gtf_ref_new(void);
GT_INLINE void gt_gtf_ref_delete(gt_gtf_ref* const ref);

GT_INLINE gt_gtf* gt_gtf_new();
GT_INLINE void gt_gtf_delete(gt_gtf* const gtf);

GT_INLINE void gt_gtf_read_line(char* line, gt_gtf* const gtf);

GT_INLINE gt_gtf_ref* gt_gtf_get_ref(gt_gtf* gtf, char* name);

GT_INLINE gt_gtf* gt_gtf_read(FILE* input);

GT_INLINE void gt_gtf_search(gt_gtf* gtf, gt_vector* target, char* const ref, const uint64_t start, const uint64_t end, const gt_gtf_strand strand);


GT_INLINE uint64_t gt_gtf_bin_search(gt_vector* entries, uint64_t t);
#endif /* GT_GTF_H_ */
