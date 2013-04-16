/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.h
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_SAM_H_
#define GT_OUTPUT_SAM_H_

#include "gt_essentials.h"
#include "gt_dna_string.h"
#include "gt_template.h"
#include "gt_sam_data_attributes.h"
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"

/*
 * Error Codes
 */
#define GT_SOE_PRINTING_MISM_STRING 10

/*
 * SAM Output attributes
 */
typedef struct {
  /* Maps */
  bool flag_best_as_primary; // Flags the Map/MMap with more phredQuality as primary
  uint64_t max_printable_maps; // Maximum number of maps printed
  bool print_unpaired_maps;
  bool compact_format; // Compact map representation in SAM via XA field
  /* CIGAR/Mismatch string */
  bool print_misms;
  /* Optional fields */
  bool print_casava; // If available print CASAVA tag as extra field
  gt_vector* sam_attributes; // Optional fields stored as sam_attributes
} gt_output_sam_attributes;
#define GT_OUTPUT_SAM_ATTR_DEFAULT() { \
  /* Maps */ \
  .flag_best_as_primary = true, \
  .max_printable_maps = UINT64_MAX, \
  .print_unpaired_maps = false, \
  .compact_format = false, \
  /* Mismatch/CIGAR string */ \
  .print_misms = false, \
  /* Optional fields */ \
  .print_casava = true, \
  .sam_attributes = NULL \
}
/* Setup */
GT_INLINE gt_output_sam_attributes* gt_output_sam_attributes_new();
GT_INLINE void gt_output_sam_attributes_delete(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_reset_defaults(gt_output_sam_attributes* const attributes);
/* Maps */
GT_INLINE void gt_output_sam_attributes_flag_best_as_primary_on(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_flag_best_as_primary_off(gt_output_sam_attributes* const attributes);
GT_INLINE uint64_t gt_output_sam_attributes_get_max_printable_maps(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_set_max_printable_maps(gt_output_sam_attributes* const attributes,const uint64_t max_printable_maps);
GT_INLINE void gt_output_sam_attributes_print_unpaired_maps_on(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_print_unpaired_maps_off(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_compact_format_on(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_compact_format_off(gt_output_sam_attributes* const attributes);
/* CIGAR/Mismatch string */
GT_INLINE void gt_output_sam_attributes_print_mismatches_on(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_print_mismatches_off(gt_output_sam_attributes* const attributes);
/* Optional fields */
GT_INLINE void gt_output_sam_attributes_print_casava_on(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_print_casava_off(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_set_sam_attributes(gt_output_sam_attributes* const attributes,gt_vector* const sam_attributes);

/*
 * SAM Headers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_headers_sh,gt_sam_headers* const sam_headers);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_headers_sa,gt_sequence_archive* const sequence_archive);

/*
 * SAM QNAME (Tag)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_qname,gt_string* const tag);
/*
 * SAM Flag
 *   Beware of the SAM flags, they might cause severe mind injuries...
 *
 * 0x1   (Bit 0)  => Template having multiple segments in sequencing
 *                   (The read was part of a pair during sequencing) [read paired]
 * 0x2   (Bit 1)  => Each segment properly aligned according to the aligner
 *                   (The read is mapped in a pair) [read mapped in proper pair]
 * 0x4   (Bit 2)  => Segment unmapped. The query sequence is unmapped [read unmapped]
 * 0x8   (Bit 3)  => Next segment in the template unmapped. The mate is unmapped [mate unmapped]
 * 0x10  (Bit 4)  => SEQ being reverse complemented. Strand of query (0=forward 1=reverse) [read reverse strand]
 * 0x20  (Bit 5)  => SEQ of the next segment in the template being reversed [mate reverse strand]
 * 0x40  (Bit 6)  => The first segment in the template [first in pair]
 * 0x80  (Bit 7)  => The last segment in the template [second in pair]
 * 0x100 (Bit 8)  => Secondary alignment [not primary alignment]
 * 0x200 (Bit 9)  => Not passing quality controls [read fails platform/vendor quality checks]
 * 0x400 (Bit 10) => PCR or optical duplicate [read is PCR or optical duplicate]
 *
 * - Bit 0x4 is the only reliable place to tell whether the segment is unmapped. If 0x4 is set,
 *     no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100
 *     and the bit 0x20 of the next segment in the template.
 * - If 0x40 and 0x80 are both set, the segment is part of a linear template, but it is neither
 *     the first nor the last segment. If both 0x40 and 0x80 are unset, the index of the segment
 *     in the template is unknown. This may happen for a non-linear template or the index is
 *     lost in data processing.
 * - Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use
 *     are aware of this bit.
 * - If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
*/
GT_INLINE uint16_t gt_output_sam_calculate_flag(
    const bool paired_end,const bool read_paired,
    const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe(
    const bool read_paired,const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_se(
    const bool read_mapped,const gt_strand read_strand,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_se_map(
    gt_map* const map,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe_map(
    gt_map* const map,gt_map* const mate,const bool is_map_first_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
/*
 * SAM CIGAR
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_cigar,gt_map* const map,gt_output_sam_attributes* const attributes);
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 *   @read and @qualities must be already RC in case of mapping on the reverse strand
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,gt_map* const mate,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
/*
 * SAM Optional fields
 *   - SAM Attributes is a gt_vector of gt_sam_attribute (gt_data_attributes.h)
 *   - SAM Attributes are stored inside general attributes under GT_ATTR_ID_SAM_FLAGS
 *       Thus, they can be found in a template, an alignment and/or a map
 *   - SAM Attributes (gt_data_attributes.h) can be either a value(int,double,string)
 *       or a function of {template,alignment,map} returning a value(int,double,string)
 *   - gt_output_sam_print_optional_fields_values() prints all the values contained in @sam_attributes
 *     gt_output_sam_print_optional_fields() prints all attributes.
 *       Those relying on a function, are generating calling that function with @template,@alignment,@map
 *       as arguments (can be NULL, so the attribute function must be ready for that)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_optional_fields_values,gt_vector* sam_attributes,gt_output_sam_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_optional_fields,
    gt_template* const template,gt_alignment* const alignment,gt_map** const map,
    gt_vector* sam_attributes,gt_output_sam_attributes* const attributes);
/*
 * SAM Record (Full SAM record printers)
 *   - @read and @qualities are NOT reverse-complemented (just printed)
 *   - @tag, @read and @qualities are used to print the SAM field
 *       (@template,@alignment can be null and are only used as to generate the opt-fields, if any)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_sam_pe_record,
    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
    gt_string* const qualities_end1,gt_string* const qualities_end2,
    gt_template* const template,gt_alignment* const alignment_end1,gt_alignment* const alignment_end2,
    gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_sam_se_record,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes);

/*
 * SAM High-level MMAp/Map Printers
 *   - Provides flexibility to print maps individually
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_mmap,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map,
    gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes);
/*
 * SAM High-level Template/Alignment Printers
 *   - Optional fields are generated from the first list of SAM Attributes found in the following order:
 *       1.- @attributes->sam_attributes
 *       2.- template->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       3.- alignment->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       4.- map->attributes{GT_ATTR_ID_SAM_FLAGS}
 *     If multiple SAM-Attributes lists are contained within these levels, only the first found one is output
 *     (completely preventing lower levels of attributes to be output)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_alignment,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_template,gt_template* const template,gt_output_sam_attributes* const output_attributes);

#endif /* GT_OUTPUT_SAM_H_ */
