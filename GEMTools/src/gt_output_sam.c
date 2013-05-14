/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_sam.h"

/*
 * Constants
 */
#define GT_OUTPUT_SAM_FORMAT_VERSION "1.4"

/*
 * Output SAM Attributes
 */
/* Setup */
GT_INLINE gt_output_sam_attributes* gt_output_sam_attributes_new() {
  gt_output_sam_attributes* const attributes = gt_alloc(gt_output_sam_attributes);
  /* Optional fields */
  attributes->sam_attributes=NULL;
  attributes->attribute_func_params=NULL;
  /* Reset defaults */
  gt_output_sam_attributes_clear(attributes);
  return attributes;
}
GT_INLINE void gt_output_sam_attributes_delete(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  if (attributes->sam_attributes!=NULL) gt_sam_attributes_delete(attributes->sam_attributes);
  if (attributes->attribute_func_params!=NULL) gt_sam_attribute_func_params_delete(attributes->attribute_func_params);
  gt_free(attributes);
}
GT_INLINE void gt_output_sam_attributes_clear(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  /* Format */
  attributes->format = GT_SAM;
  /* Read/Qualities */
  attributes->always_output_read__qualities = true;
  attributes->qualities_offset = GT_QUALS_OFFSET_33;
  /* Maps */
  attributes->max_printable_maps = UINT64_MAX;
  attributes->compact_format = true;
  /* Mismatch/CIGAR string */
  attributes->print_mismatches = false;
  /* SAM Optional Fields */
  attributes->print_optional_fields = true;
  if (attributes->sam_attributes!=NULL) {
    gt_sam_attributes_clear(attributes->sam_attributes);
  } else {
    attributes->sam_attributes = gt_sam_attributes_new();
  }
  if (attributes->attribute_func_params!=NULL) {
    gt_sam_attribute_func_params_clear(attributes->attribute_func_params);
  } else {
    attributes->attribute_func_params = gt_sam_attribute_func_params_new();
  }
}

/* Format */
GT_INLINE void gt_output_sam_attributes_set_format(gt_output_sam_attributes* const attributes,gt_output_sam_format_t const format) {
  GT_NULL_CHECK(attributes);
  attributes->format = format;
}
/* Read/Qualities */
GT_INLINE void gt_output_sam_attributes_dump_read__qualities_once(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->always_output_read__qualities = false;
}
GT_INLINE void gt_output_sam_attributes_always_dump_read__qualities(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->always_output_read__qualities = true;
}
/* Maps */
GT_INLINE void gt_output_sam_attributes_set_max_printable_maps(gt_output_sam_attributes* const attributes,const uint64_t max_printable_maps) {
  GT_NULL_CHECK(attributes);
  attributes->max_printable_maps = max_printable_maps;
}
GT_INLINE void gt_output_sam_attributes_set_compact_format(gt_output_sam_attributes* const attributes,const bool compact_format) {
  GT_NULL_CHECK(attributes);
  attributes->compact_format = compact_format;
}
/* CIGAR/Mismatch string */
GT_INLINE void gt_output_sam_attributes_set_print_mismatches(gt_output_sam_attributes* const attributes,const bool print_mismatches) {
  GT_NULL_CHECK(attributes);
  attributes->print_mismatches = print_mismatches;
}
/* SAM Optional Fields */
GT_INLINE void gt_output_sam_attributes_set_print_optional_fields(gt_output_sam_attributes* const attributes,const bool print_optional_fields) {
  GT_NULL_CHECK(attributes);
  attributes->print_optional_fields = print_optional_fields;
}
GT_INLINE void gt_output_sam_attributes_set_reference_sequence_archive(gt_output_sam_attributes* const attributes,gt_sequence_archive* const reference_sequence_archive) {
  GT_NULL_CHECK(attributes);
  attributes->attribute_func_params->sequence_archive = reference_sequence_archive;
}
GT_INLINE gt_sam_attributes* gt_output_sam_attributes_get_sam_attributes(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->sam_attributes;
}

/*
 * // TODO replace with sample
 * SAM Headers
 * @HD  VN:1.0  SO:unsorted
 * @HD  VN:1.3
 *
 * @SQ  SN:chr10   LN:135534747  AS:hg19_ncbi37  SP:human
 * @SQ  SN:chr10   LN:135534747
 * @SQ  SN:chr11   LN:135006516
 *
 * @RG  ID:NOID   PG:tmap  SM:NOSM
 * @RG  ID:0      PG:GEM   PL:ILLUMINA  SM:0
 *
 * @PG  ID:tmap  CL:map4 -f /home/user/references/hsapiens_v37.fa -r /home/user/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq -i fastq -s /home/user/miseq_S/map41.H.Sapiens.1M.S.MiSeq.low.l250_1.fastq.sam   VN:3.0.1
 * @PG  ID:dvtgm PN:stampy     VN:1.0.17_(r1481)  CL:-g /home/devel/user/hsapiens_v37 -h /home/devel/user/hsapiens_v37 --bwaoptions=/home/user/hsapiens_v37.fa --substitutionrate=0.08 --maxbasequal 90 -M /home/user/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq
 * @PG  ID:GEM   PN:gem-2-sam  VN:1.414
 *
 * @CO  TM:Fri, 30 Nov 2012 14:14:13 CET        WD:/home/user/benchmark/DEF/miseq_S/CMD      HN:cn38.bullx   UN:user
 * @CO  BWAVersion: 0.6.1-r104
 *
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_headers
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_headers_sh,gt_sam_headers* const sam_headers);
GT_INLINE gt_status gt_output_sam_gprint_headers_sh(gt_generic_printer* const gprinter,gt_sam_headers* const sam_headers) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // Print all @HD line (Header line)
  if (sam_headers==NULL || gt_string_is_null(sam_headers->header)) {
    gt_gprintf(gprinter,"@HD\tVN:"GT_OUTPUT_SAM_FORMAT_VERSION"\n");
  } else {
    gt_gprintf(gprinter,"@HD\t"PRIgts"\n",PRIgts_content(sam_headers->header));
  }
  // Print all @SQ lines (Reference sequence dictionary)
  if (sam_headers==NULL || sam_headers->sequence_archive!=NULL) {
    gt_output_sam_gprint_headers_sa(gprinter,sam_headers->sequence_archive);
  }
  // Print all @RG lines (Read group)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->read_group)) {
    gt_gprintf(gprinter,"@RG\tID:0\tPG:GTools\tSM:0\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->read_group,rg_line,line_num,gt_string*) {
      gt_gprintf(gprinter,"@RG\t"PRIgts"\n",PRIgts_content(*rg_line));
    }
  }
  // Print all @PG lines (Program)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->program)) {
    gt_gprintf(gprinter,"@PG\tID:GToolsLib\tPN:gt_output_sam\tVN:"GT_VERSION"\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->program,prog_line,line_num,gt_string*) {
      gt_gprintf(gprinter,"@PG\t"PRIgts"\n",PRIgts_content(*prog_line));
    }
  }
  // Print all @CO lines (Comments)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->comments)) {
    // Print Current Date
    gt_gprintf(gprinter,"@CO\tTM:\t");
    time_t current_time=time(0);
    struct tm local_time;
    localtime_r(&current_time,&local_time);
    gt_gprintf(gprinter,"%4d/%d/%d %02d:%02d:%02d CET",
        1900+local_time.tm_year,local_time.tm_mon+1,local_time.tm_mday,
        local_time.tm_hour,local_time.tm_min,local_time.tm_sec);
    // Print GT banner
    gt_gprintf(gprinter,"\tGTools v"GT_VERSION" "GT_GIT_URL"\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->comments,comment,line_num,gt_string*) {
      gt_gprintf(gprinter,"@CO\t"PRIgts"\n",PRIgts_content(*comment));
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sequence_archive
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_headers_sa,gt_sequence_archive* const sequence_archive);
GT_INLINE gt_status gt_output_sam_gprint_headers_sa(gt_generic_printer* const gprinter,gt_sequence_archive* const sequence_archive) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Print all @SQ lines (Reference sequence dictionary)
  //   @SQ  SN:chr10 LN:135534747 AS:hg19_ncbi37  SP:human
  gt_sequence_archive_iterator sequence_archive_it;
  gt_sequence_archive_new_iterator(sequence_archive,&sequence_archive_it);
  gt_segmented_sequence* seq;
  while ((seq=gt_sequence_archive_iterator_next(&sequence_archive_it))) {
    // @SQ  SN:chr10 LN:135534747
    gt_gprintf(gprinter,"@SQ\tSN:"PRIgts"\tLN:%"PRIu64"\n",PRIgts_content(seq->seq_name),seq->sequence_total_length);
  }
  return 0;
}
/*
 * SAM QNAME (Tag)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_qname,gt_string* const tag);
GT_INLINE gt_status gt_output_sam_gprint_qname(gt_generic_printer* const gprinter,gt_string* const tag) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // Print the plain tag (no read pair info, nor extra tag nor nothing)
  char* const tag_buffer = gt_string_get_string(tag);
  int i;
  for (i=0;i<gt_string_get_length(tag);++i) {
    if (tag_buffer[i]==SPACE) break;
  }
  gt_gprintf(gprinter,PRIgts,i,gt_string_get_string(tag));
  return 0;
}
/*
 * SAM Flag
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
 * - RULE1:: Bit 0x4 is the only reliable place to tell whether the segment is unmapped. If 0x4 is set,
 *     no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100
 *     and the bit 0x20 of the next segment in the template.
 * - RULE2:: If 0x40 and 0x80 are both set, the segment is part of a linear template, but it is neither
 *     the first nor the last segment. If both 0x40 and 0x80 are unset, the index of the segment
 *     in the template is unknown. This may happen for a non-linear template or the index is
 *     lost in data processing.
 * - RULE3:: Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use
 *     are aware of this bit.
 * - RULE4:: If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
*/
GT_INLINE uint16_t gt_output_sam_calculate_flag_se_map(
    gt_map* const map,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return (gt_expect_true(map!=NULL)) ?
      gt_output_sam_calculate_flag_se(true,map->strand,secondary_alignment,not_passing_QC,PCR_duplicate): // Mapped
      gt_output_sam_calculate_flag_se(false,FORWARD,secondary_alignment,not_passing_QC,PCR_duplicate); // Unmapped
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe_map(
    gt_map* const map,gt_map* const mate,const bool is_map_first_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return gt_output_sam_calculate_flag_pe(
      map!=NULL && mate!=NULL, /* read_paired */
      map!=NULL,               /* read_mapped */
      mate!=NULL,              /* mate_strand */
      (map!=NULL) ? map->strand : FORWARD,   /* read_strand */
      (mate!=NULL) ? mate->strand : FORWARD, /* mate_strand */
      is_map_first_in_pair,    /* first_in_pair */
      !is_map_first_in_pair,   /* last_in_pair */
      secondary_alignment,not_passing_QC,PCR_duplicate);
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe(
    const bool read_paired,const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  /* 0x1 */
  uint16_t sam_flag = GT_SAM_FLAG_MULTIPLE_SEGMENTS;
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= GT_SAM_FLAG_UNMAPPED;
    /* 0x8  */
    if (!mate_mapped) {
      sam_flag |= GT_SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x20 */if (mate_strand==REVERSE) sam_flag |= GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  } else {
    /* 0x10 */  if (read_strand==REVERSE) sam_flag |= GT_SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= GT_SAM_FLAG_SECONDARY_ALIGNMENT;
    /* 0x8 */
    if (!mate_mapped) {
      sam_flag |= GT_SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x2 Each segment properly aligned can only take place if both ends are mapped */
      if (read_paired) sam_flag |= GT_SAM_FLAG_PROPERLY_ALIGNED;
      /* 0x20 */
      if (mate_strand==REVERSE) sam_flag |= GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  }
  /* 0x40 */  if (first_in_pair)  sam_flag |= GT_SAM_FLAG_FIRST_SEGMENT;
  /* 0x80 */  if (last_in_pair)   sam_flag |= GT_SAM_FLAG_LAST_SEGMENT;
  /* 0x200 */ if (not_passing_QC) sam_flag |= GT_SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */ if (PCR_duplicate)  sam_flag |= GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_se(
    const bool read_mapped,const gt_strand read_strand,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  uint16_t sam_flag = 0; // (**RULE4)
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= GT_SAM_FLAG_UNMAPPED;
  } else {
    /* 0x10 */  if (read_strand==REVERSE) sam_flag |= GT_SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= GT_SAM_FLAG_SECONDARY_ALIGNMENT;
  }
  /* 0x200 */
  if (not_passing_QC) sam_flag |= GT_SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */
  if (PCR_duplicate) sam_flag |= GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
GT_INLINE uint16_t gt_output_sam_calculate_flag(
    const bool paired_end,const bool read_paired,
    const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return (paired_end) ? /* 0x1 */
      gt_output_sam_calculate_flag_pe(
          read_paired,read_mapped,mate_mapped,
          read_strand,mate_strand,first_in_pair,last_in_pair,
          secondary_alignment,not_passing_QC,PCR_duplicate): // PE
      gt_output_sam_calculate_flag_se(read_mapped,read_strand,
          secondary_alignment,not_passing_QC,PCR_duplicate); // SE
}
/*
 * XA maps
 */
GT_INLINE void gt_output_sam_gprint_map_xa(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  // chr12,+91022,101M,0
  gt_map* base_map = map;
  GT_MAP_SEGMENT_ITERATE(map,map_segment) {
    // Print the map
    gt_gprintf(gprinter,PRIgts",%c%lu,",
        PRIgts_content(map_segment->seq_name),
        (map_segment->strand==FORWARD)?'+':'-',
        gt_map_get_position(map_segment));
    gt_output_sam_gprint_cigar(gprinter,map_segment,attributes);
    gt_gprintf(gprinter,",%"PRIu64,gt_map_get_levenshtein_distance(map_segment));
  }
}
/*
 * SAM CIGAR
 */
#define GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"%c",misms_pos-centinel,(attributes->print_mismatches)?'=':'M'); \
    centinel = misms_pos; \
  }
#define GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"%c",centinel-misms_pos,(attributes->print_mismatches)?'=':'M'); \
    centinel = misms_pos; \
  }
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar_reverse(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  // Map auxiliary variables
  const uint64_t map_length = gt_map_get_base_length(map);
  int64_t centinel = map_length-1;
  // Mismatches auxiliary variables
  const uint64_t num_misms = gt_map_get_num_misms(map);
  uint64_t misms_n = num_misms;
  gt_misms* misms;
  while (misms_n > 0) {
    misms = gt_map_get_misms(map,misms_n-1);
    const uint64_t misms_pos = gt_misms_get_position(misms);
    switch (misms->misms_type) {
      case MISMS:
        if (attributes->print_mismatches) {
          GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
          gt_gprintf(gprinter,"1X");
          --centinel;
        }
        break;
      case INS: // SAM Deletion
        GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"D",gt_misms_get_size(misms));
        break;
      case DEL: // SAM Insertion
        GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"I",gt_misms_get_size(misms));
        centinel-=gt_misms_get_size(misms);
        break;
      default:
        gt_error(SELECTION_NOT_VALID);
        return GT_SOE_PRINTING_MISM_STRING;
        break;
    }
    --misms_n;
  }
  if (centinel >= 0) gt_gprintf(gprinter,"%"PRIu64"M",centinel+1);
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar_forward(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  const uint64_t map_length = gt_map_get_base_length(map);
  uint64_t centinel = 0;
  GT_MISMS_ITERATE(map,misms) {
    const uint64_t misms_pos = gt_misms_get_position(misms);
    switch (misms->misms_type) {
      case MISMS:
        if (attributes->print_mismatches) {
          GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
          gt_gprintf(gprinter,"1X");
          ++centinel;
        }
        break;
      case INS: // SAM Deletion
        GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"D",gt_misms_get_size(misms));
        break;
      case DEL: // SAM Insertion
        GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"I",gt_misms_get_size(misms));
        centinel+=gt_misms_get_size(misms);
        break;
      default:
        gt_error(SELECTION_NOT_VALID);
        return GT_SOE_PRINTING_MISM_STRING;
        break;
    }
  }
  if (centinel < map_length) gt_gprintf(gprinter,"%"PRIu64"M",map_length-centinel);
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar(
    gt_generic_printer* const gprinter,gt_map* const map_block,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_block);
  gt_status error_code = 0;
  if (gt_map_get_strand(map_block)==REVERSE) {
    // Check following map blocks
    gt_map* const next_map_block = gt_map_get_next_block(map_block);
    if (next_map_block!=NULL && GT_MAP_IS_SAME_SEGMENT(map_block,next_map_block)) { // SplitMap (Otherwise is a quimera)
      error_code = gt_output_sam_gprint_map_block_cigar(gprinter,next_map_block,attributes);
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map_block));
    }
    // Print CIGAR for current map block
    gt_output_sam_gprint_map_block_cigar_reverse(gprinter,map_block,attributes);
  } else {
    // Print CIGAR for current map block
    gt_output_sam_gprint_map_block_cigar_forward(gprinter,map_block,attributes);
    // Check following map blocks
    gt_map* const next_map_block = gt_map_get_next_block(map_block);
    if (next_map_block!=NULL && GT_MAP_IS_SAME_SEGMENT(map_block,next_map_block)) { // SplitMap (Otherwise is a quimera)
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map_block));
      error_code = gt_output_sam_gprint_map_block_cigar(gprinter,next_map_block,attributes);
    }
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_map_cigar(
    gt_generic_printer* const gprinter,gt_map* const map_segment,gt_output_sam_attributes* const attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_segment);
  gt_status error_code = 0;
  // Check strandness
  if (gt_map_get_strand(map_segment)==FORWARD) {
    if (hard_left_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_left_trim_read);
    error_code=gt_output_sam_gprint_map_block_cigar(gprinter,map_segment,attributes);
    if (hard_right_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_right_trim_read);
  } else {
    if (hard_right_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_right_trim_read);
    error_code=gt_output_sam_gprint_map_block_cigar(gprinter,map_segment,attributes);
    if (hard_left_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_left_trim_read);
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map_segment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_cigar,gt_map* const map_segment,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_cigar(gt_generic_printer* const gprinter,gt_map* const map_segment,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_segment);
  return gt_output_sam_gprint_map_cigar(gprinter,map_segment,attributes,0,0);
}
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 *   Don't handle quimeras (just print one record out of the first map segment)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map,position,phred_score, \
    hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_core_fields_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_gprintf(gprinter,"\t%"PRId16,gt_output_sam_calculate_flag_se_map(map,secondary_alignment,not_passing_QC,PCR_duplicate));
  // Is mapped?
  if (gt_expect_true(map!=NULL)) {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRId8"\t",PRIgts_content(map->seq_name),position,phred_score);
    // (6) Print CIGAR
    gt_output_sam_gprint_map_cigar(gprinter,map,attributes,hard_left_trim_read,hard_right_trim_read);
  } else {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    // (6) Print CIGAR
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  //  (7) Print RNEXT
  //  (8) Print PNEXT
  //  (9) Print TLEN
  // (10) Print SEQ
  // (11) Print QUAL
  if (read!=NULL && qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t"PRIgts,
        PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read),
        PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  } else if (read!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t*",PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read));
  } else if (qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t*\t"PRIgts,PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities, \
    map,position,phred_score,mate,mate_position,template_length, \
    hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    gt_map* const mate,const uint64_t mate_position,const uint64_t template_length,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_core_fields_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    gt_map* const mate,const uint64_t mate_position,const uint64_t template_length,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_gprintf(gprinter,"\t%"PRId16,gt_output_sam_calculate_flag_pe_map(
      map,mate,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate));
  // (3) Print RNAME
  // (4) Print POS
  // (5) Print MAPQ
  // (6) Print CIGAR
  if (map!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRId8"\t",PRIgts_content(map->seq_name),position,phred_score);
    gt_output_sam_gprint_map_cigar(gprinter,map,attributes,hard_left_trim_read,hard_right_trim_read); // CIGAR
  } else {
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  // (7) Print RNEXT
  // (8) Print PNEXT
  // (9) Print TLEN
  if (mate!=NULL) {
    if (gt_string_equals(map->seq_name,mate->seq_name)) {
      gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRIu64,PRIgts_content(mate->seq_name),mate_position,template_length);
    } else {
      gt_gprintf(gprinter,"\t=\t%"PRIu64"\t%"PRIu64,mate_position,template_length);
    }
  } else {
    gt_gprintf(gprinter,"\t*\t0\t0");
  }
  // (10) Print SEQ
  // (11) Print QUAL
  if (read!=NULL && qualities!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t"PRIgts,
        PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read),
        PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  } else if (read!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t*",PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read));
  } else if (qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t"PRIgts,PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_segment,hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map_segment,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_core_fields_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map_segment,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  return gt_output_sam_gprint_core_fields_se(gprinter,tag,read,qualities,
      map_segment,gt_map_get_position(map_segment),gt_map_get_phred_score(map_segment),
      hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_segment,mate_segment,mmap_attributes, \
    hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map_segment,gt_map* const mate_segment,gt_mmap_attributes* const mmap_attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_core_fields_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map_segment,gt_map* const mate_segment,gt_mmap_attributes* const mmap_attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  const int64_t observed_template_size = gt_map_get_observed_template_size(map_segment,mate_segment);
  return gt_output_sam_gprint_core_fields_pe(gprinter,tag,read,qualities,
      map_segment,gt_map_get_position(map_segment),gt_map_get_phred_score(map_segment),
      mate_segment,gt_map_get_position(map_segment),observed_template_size,
      hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_placeholder,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_placeholder,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map_placeholder* const map_placeholder,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map_placeholder* const map_placeholder,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  if (map_placeholder->type==GT_MAP_PLACEHOLDER) {
    return gt_output_sam_gprint_map_core_fields_se(gprinter,tag,read,qualities,map_placeholder->map,
        map_placeholder->hard_trim_left,map_placeholder->hard_trim_right,
        map_placeholder->secondary_alignment,map_placeholder->not_passing_QC,map_placeholder->PCR_duplicate,
        output_attributes);
  } else {
    return gt_output_sam_gprint_map_core_fields_pe(gprinter,tag,read,qualities,
        map_placeholder->map,map_placeholder->paired_end.mate,map_placeholder->paired_end.mmap_attributes,
        map_placeholder->paired_end.paired_end_position==0,
        map_placeholder->hard_trim_left,map_placeholder->hard_trim_right,
        map_placeholder->secondary_alignment,map_placeholder->not_passing_QC,map_placeholder->PCR_duplicate,output_attributes);
  }
}
/*
 * SAM Optional fields
 *   - SAM Attributes is a shash of gt_sam_attribute (gt_sam_data_attributes.h)
 *   - SAM Attributes (gt_sam_data_attributes.h) can be either a value(int,double,string)
 *       or a function -> f(gt_sam_attribute_func_params* params) returning a value(int,double,string)
 *   - gt_output_sam_print_optional_fields_values() prints all the values contained in @sam_attributes
 *     gt_output_sam_print_optional_fields() prints all attributes.
 *       Those relying on a function, are generating calling that function with @gt_sam_attribute_func_params
 *       as argument (some fields can be NULL, so the attribute function must be ready to deal with that)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields_values,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields_values(gt_generic_printer* const gprinter,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  if (sam_attributes==NULL) sam_attributes = output_attributes->sam_attributes;
  if (sam_attributes!=NULL) {
    GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
    GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,sam_attribute) {
      if (sam_attribute->attribute_type == SAM_ATTR_INT_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_value);
      } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->f_value);
      } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(sam_attribute->s_value));
      }
    } GT_SAM_ATTRIBUTES_END_ITERATE;
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields(gt_generic_printer* const gprinter,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  if (sam_attributes==NULL) sam_attributes = output_attributes->sam_attributes;
  if (sam_attributes!=NULL) {
    GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
    GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,sam_attribute) {
      // Values
      if (sam_attribute->attribute_type == SAM_ATTR_INT_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_value);
      } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->f_value);
      } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_VALUE) {
        gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(sam_attribute->s_value));
      } else
      // Functions
      if (sam_attribute->attribute_type == SAM_ATTR_INT_FUNC) {
        if (sam_attribute->i_func(output_attributes->attribute_func_params)==0) { // Generate i-value
          gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,
              output_attributes->attribute_func_params->return_i);
        }
      } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_FUNC) {
        if (sam_attribute->f_func(output_attributes->attribute_func_params)==0) { // Generate f-value
          gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,
              output_attributes->attribute_func_params->return_f);
        }
      } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_FUNC) {
        if (sam_attribute->s_func(output_attributes->attribute_func_params)==0) { // Generate s-value
          gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,
              PRIgts_content(output_attributes->attribute_func_params->return_s));
        }
      }
    } GT_SAM_ATTRIBUTES_END_ITERATE;
  }
  return 0;
}
///*
// * SAM Record (Full SAM record printers)
// *   - @read and @qualities are NOT reverse-complemented (just printed)
// *   - Prints properly quimeric reads
// */
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_se,
//    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const attributes);
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_sam_pe,
//    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
//    gt_string* const qualities_end1,gt_string* const qualities_end2,
//    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const output_attributes);
///*
// * SAM High-level MMap/Map Printers
// */
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_mmap,
//    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const output_attributes);
//GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map,
//    gt_alignment* const alignment,gt_map* const map,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const output_attributes);


///*
// * SAM Record (Full SAM record printers)
// *   - @read and @qualities are NOT reverse-complemented (just printed)
// *   - Prints properly quimeric reads
// */
//#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
//#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
//GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_se,
//    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const attributes);
//GT_INLINE gt_status gt_output_sam_gprint_map_se(gt_generic_printer* const gprinter,
//    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const attributes) {
//  GT_GENERIC_PRINTER_CHECK(gprinter);
//  GT_STRING_CHECK(tag);
//  gt_status error_code = 0;
//  if (map==NULL) { // Unmapped
//    // Print Core Fields
//    error_code|=gt_output_sam_gprint_core_fields_se(
//        gprinter,tag,read,qualities,false,NULL,0,FORWARD,0,0,0,
//        secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
//    // Print Optional Fields
//    if (attributes->print_optional_fields) {
//      error_code|=gt_output_sam_gprint_optional_fields(gprinter,NULL,attributes);
//    }
//    gt_gprintf(gprinter,"\n");
//  } else {
//    // Traverse all connected segments of the map
//    uint64_t num_segments_printed = 0; // There can only be 1 primary alignment (in case of quimeras)
//    GT_MAP_SEGMENT_ITERATE(map,map_segment_iterator) {
//      gt_map* const base_map = map_segment_iterator.begin_map_segment;
//      const int64_t accumulated_offset = gt_map_segments_iterator_get_accumulated_offset(&map_segment_iterator);
//      const int64_t remaining_bases = gt_map_segments_iterator_get_remaining_bases(&map_segment_iterator,read);
//      if (base_map->strand==FORWARD) {
//        error_code|=gt_output_sam_gprint_core_fields_se(gprinter,tag,read,qualities,true,base_map->seq_name,
//            gt_map_get_position(base_map),base_map->strand,map->phred_score,
//            accumulated_offset,remaining_bases,
//            secondary_alignment||num_segments_printed>0,not_passing_QC,PCR_duplicate,attributes);
//      } else {
//        error_code|=gt_output_sam_gprint_core_fields_se(gprinter,tag,read,qualities,true,base_map->seq_name,
//            gt_map_get_position(map_segment_iterator.end_map_segment),base_map->strand,map->phred_score,
//            remaining_bases,accumulated_offset,
//            secondary_alignment||num_segments_printed>0,not_passing_QC,PCR_duplicate,attributes);
//      }
//      // Print Optional Fields
//      if (attributes->print_optional_fields) {
//        error_code|=gt_output_sam_gprint_optional_fields(gprinter,gt_attributes_get_sam_attributes(base_map->attributes),attributes);
//      }
//      gt_gprintf(gprinter,"\n"); ++num_segments_printed;
//    }
//  }
//  return error_code;
//}
//GT_INLINE gt_status gt_output_sam_gprint_map_pe_record(gt_generic_printer* const gprinter,
//    gt_string* const tag,gt_string* const read,gt_string* const qualities,
//    gt_map* const map,gt_map* const mate,gt_mmap_attributes* const mmap_attributes,
//    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const attributes) {
//  GT_GENERIC_PRINTER_CHECK(gprinter);
//  GT_STRING_CHECK(tag);
//  gt_status error_code = 0;
//  if (map==NULL) { // Map unmapped
//    if (mate==NULL) { // Mate unmapped
//
//    } else {
//      // Get the first segment of the mate
//      gt_map_segment_iterator mate_segment_iterator;
//      gt_map_new_segments_iterator(mate,&mate_segment_iterator);
//      gt_map_next_segment(&mate_segment_iterator);
//      return gt_output_sam_gprint_map_core_fields_pe(
//          gprinter,tag,read,qualities,false,NULL,0,FORWARD,0,0,
//          true,mate->seq_name,gt_map_segments_iterator_get_position(&mate_segment_iterator),FORWARD,0,0,
//          is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
//    }
//  } else {
//    gt_map* const fist_block_mate = mate_segment_iterator.begin_map_segment;
//    gt_map* const last_block_mate = mate_segment_iterator.end_map_segment;
//    // Traverse all connected segments of the map
//    //   TODO: secondary alignment in case of quimeras
//    //   TODO: mate/next segment in case of quimeras
//    GT_MAP_SEGMENT_ITERATE(map,map_segment_iterator) {
//      gt_map* const fist_block_map = map_segment_iterator.begin_map_segment;
//      gt_map* const last_block_map = map_segment_iterator.end_map_segment;
//      const int64_t accumulated_offset = gt_map_segments_iterator_get_accumulated_offset(&map_segment_iterator);
//      const int64_t remaining_bases = gt_map_segments_iterator_get_remaining_bases(&map_segment_iterator,read);
//      if (base_map->strand==FORWARD) {
//        error_code=gt_output_sam_gprint_map_core_fields_pe(gprinter,tag,read,qualities,true,base_map->seq_name,
//            gt_map_get_position(base_map),base_map->strand,map->phred_score,
//
//
//            accumulated_offset,(remaining_bases>0)?remaining_bases:0,
//            secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
//      } else {
//        error_code=gt_output_sam_gprint_map_core_fields_pe(gprinter,tag,read,qualities,true,base_map->seq_name,
//            gt_map_get_position(map_segment_iterator->end_map_segment),base_map->strand,map->phred_score,
//            (remaining_bases>0)?remaining_bases:0,map_segment_iterator->accumulated_offset,
//            secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
//      }
//      if (error_code) return error_code;
//    }
//  }
//  return 0;
//}
//#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
//#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read_end1,read_end2,qualities_end1,qualities_end2, \
//    map_end1,map_end2,mmap_attributes,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes
//GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_pe,
//    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
//    gt_string* const qualities_end1,gt_string* const qualities_end2,
//    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const output_attributes);
//GT_INLINE gt_status gt_output_sam_gprint_map_pe(gt_generic_printer* const gprinter,
//    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
//    gt_string* const qualities_end1,gt_string* const qualities_end2,
//    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
//    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
//    gt_output_sam_attributes* const output_attributes) {
//  // TODO
//}
/*
 * SAM Placeholders Printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map_placeholder,primary_position,primary_position
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_alignment_map_placeholder_xa,
    gt_vector* const map_placeholder,const uint64_t primary_position,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_alignment_map_placeholder_se_xa(gt_generic_printer* const gprinter,
    gt_vector* const map_placeholder,const uint64_t primary_position,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder);
  GT_NULL_CHECK(attributes);
  uint64_t num_maps_printed = 0;
  if (attributes->max_printable_maps == 0) return 0;
  if (gt_vector_get_used(map_placeholder) > 1) {
    gt_gprintf(gprinter,"\tXA:Z:");
    GT_VECTOR_ITERATE(map_placeholder,map_ph,map_placeholder_position,gt_map_placeholder) {
      // Filter PH
      if (num_maps_printed > attributes->max_printable_maps) break;
      if (map_ph->type!=GT_MAP_PLACEHOLDER || map_placeholder_position==primary_position) continue;
      // Print the map (XA)
      gt_output_sam_gprint_map_xa(gprinter,map_ph->map,attributes);
      ++num_maps_printed;
    }
  }
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_alignment_map_placeholder_se(gt_generic_printer* const gprinter,
    gt_vector* const map_placeholder,const uint64_t primary_position,
    const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  // Produce RC of the mapping's read/qualities
  gt_string *read_f=alignment->read, *qualities_f=alignment->qualities;
  gt_string *read_rc=NULL, *qualities_r=NULL;
  // Get primary map
  gt_map* const map = gt_vector_get_elm(mmap_placeholder,primary_position,gt_map_placeholder)->map;
  if (map!=NULL && attributes->compact_format) {
    // Fetch sam attributes
    gt_vector* const sam_attributes = gt_attribute_sam_fetch_attributes(map->attributes);
    gt_vector* const current_sam_attributes = (sam_attributes!=NULL) ? sam_attributes : attributes->sam_attributes;
    // Print primary alignment
    if (gt_map_get_strand(primary_map)==FORWARD) {
      error_code |= gt_output_sam_gprint_map_placeholder(gprinter,);
    } else {
      if (gt_expect_false(read_rc==NULL)) { // Check RC
        read_rc = gt_dna_string_reverse_complement_dup(read_f);
        qualities_r = gt_string_reverse_dup(qualities_f);
      }
      error_code |= gt_output_sam_gprint_map_placeholder();
    }
    // Print XA:Z field
    gt_output_sam_gprint_alignment_map_placeholder_xa(gprinter,maps_sort,alignment,primary_map,output_attributes);
    // Print Optional Fields
    gt_output_sam_gprint_optional_fields(gprinter,NULL,alignment,(gt_map** const)&primary_map,current_sam_attributes,output_attributes);
    // Print Xtras
    if (output_attributes->print_casava) gt_output_sam_gprint_op_casava(gprinter,alignment,output_attributes);
    gt_gprintf(gprinter,"\n");
  } else {
    uint64_t num_maps_printed = 0;
    GT_VECTOR_ITERATE(maps_sort,map_placeholder,map_placeholder_it,gt_map_placeholder) {
      if (map_placeholder->type != GT_MAP_PLACEHOLDER) continue;
      if (num_maps_printed > output_attributes->max_printable_maps) break;
      gt_map* const map = gt_alignment_get_map(alignment,map_placeholder->mmap_position);
      // Choose sam attributes
      gt_vector* const map_sam_attributes = gt_attribute_sam_fetch_attributes(map->attributes);
      gt_vector* const current_sam_attributes = (map_sam_attributes!=NULL) ? map_sam_attributes : sam_attributes;
      // Print MAP
      if (gt_map_get_strand(map)==FORWARD) {
        error_code |= gt_output_sam_gprint_sam_se_record(gprinter,
            alignment->tag,read_f,qualities_f,alignment,map,
            map!=primary_map,not_passing_QC,PCR_duplicate,current_sam_attributes,output_attributes);
      } else {
        if (gt_expect_false(read_rc==NULL)) { // Check RC
          read_rc = gt_dna_string_reverse_complement_dup(read_f);
          qualities_r = gt_string_reverse_dup(qualities_f);
        }
        error_code |= gt_output_sam_gprint_sam_se_record(gprinter,
            alignment->tag,read_rc,qualities_r,alignment,map,
            map!=primary_map,not_passing_QC,PCR_duplicate,current_sam_attributes,output_attributes);
      }
      ++num_maps_printed;
    }
  }
  // Free
  if (read_rc!=NULL) {
    gt_string_delete(read_rc);
    gt_string_delete(qualities_r);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_template_map_placeholder(gt_generic_printer* const gprinter,
    gt_vector* const maps_sort,const bool compact_format,gt_template* const template,gt_map** const primary_mmap,
    const bool not_passing_QC,const bool PCR_duplicate,gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes) {
  // TODO
  return 0;
}
/*
 * SAM High-level MMap/Map Printers
 */
#define GT_OUTPUT_SAM_SETUP_READ__QUALITIES(alignment,map,read_f,qualities_f,read_rc,qualities_r,rc_complement) \
  /* Produce RC of the mapping's read/qualities */ \
  rc_complement = map!=NULL && gt_map_get_strand(map)==REVERSE; \
  if (alignment!=NULL) { \
    if (alignment->read!=NULL) { \
      read_f = alignment->read; \
      if (rc_complement) read_rc = gt_dna_string_reverse_complement_dup(read_f); \
    } \
    if (alignment->qualities!=NULL) { \
      qualities_f = alignment->qualities; \
      if (rc_complement) qualities_r = gt_string_reverse_dup(qualities_f); \
    } \
  }
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,map_end1,map_end2,mmap_attributes,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_mmap,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_mmap(gt_generic_printer* const gprinter,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  // Get alignments
  gt_alignment* const alignment_end1 = (gt_template_get_num_blocks(template)>0) ? gt_template_get_block(template,0): NULL;
  gt_alignment* const alignment_end2 = (gt_template_get_num_blocks(template)>0) ? gt_template_get_block(template,0): NULL;
  // Produce RC of the mapping's read/qualities
  bool rc_complement_end1, rc_complement_end2;
  gt_string *read_f_end1=NULL, *qualities_f_end1=NULL, *read_rc_end1=NULL, *qualities_r_end1=NULL;
  gt_string *read_f_end2=NULL, *qualities_f_end2=NULL, *read_rc_end2=NULL, *qualities_r_end2=NULL;
  GT_OUTPUT_SAM_SETUP_READ__QUALITIES(alignment_end1,map_end1,read_f_end1,qualities_f_end1,read_rc_end1,qualities_r_end1,rc_complement_end1);
  GT_OUTPUT_SAM_SETUP_READ__QUALITIES(alignment_end2,map_end2,read_f_end2,qualities_f_end2,read_rc_end2,qualities_r_end2,rc_complement_end2);
  // Print SAM record
  gt_output_sam_gprint_sam_pe_record(gprinter,template->tag,
      ((rc_complement_end1) ? read_rc_end1 : read_f_end1),
      ((rc_complement_end2) ? read_rc_end2 : read_f_end2),
      ((rc_complement_end1) ? qualities_r_end1 : qualities_f_end1),
      ((rc_complement_end2) ? qualities_r_end2 : qualities_f_end2),
      template,map_end1,map_end2,
      secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  if (read_rc_end1) gt_string_delete(read_rc_end1);
  if (qualities_r_end1) gt_string_delete(qualities_r_end1);
  if (read_rc_end2) gt_string_delete(read_rc_end2);
  if (qualities_r_end2) gt_string_delete(qualities_r_end2);
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,map,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map,
    gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_map(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  // Produce RC of the mapping's read/qualities
  bool rc_complement;
  gt_string *read_f=NULL, *qualities_f=NULL, *read_rc=NULL, *qualities_r=NULL;
  GT_OUTPUT_SAM_SETUP_READ__QUALITIES(alignment,map,read_f,qualities_f,read_rc,qualities_r,rc_complement);
  // Print SAM record
  gt_output_sam_gprint_sam_se_record(gprinter,alignment->tag,
      ((rc_complement) ? read_rc : read_f),((rc_complement) ? qualities_r : qualities_f),
      alignment,map,secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  if (read_rc) gt_string_delete(read_rc);
  if (qualities_r) gt_string_delete(qualities_r);
  return 0;
}
/*
 * SAM High-level Template/Alignment Printers
 *   - Optional fields are generated from the first SAM-Attributes object found in the following order:
 *       1.- map->attributes{GT_ATTR_ID_SAM_FLAGS} / mmap_attributes->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       2.- @output_attributes->sam_attributes
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_alignment,gt_alignment* const alignment,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_attributes);
  gt_status error_code;
  // Get passingQC and PCRDuplicate flags
  bool passing_QC = true, PCR_duplicate = false, *aux;
  aux = gt_attributes_get(alignment->attributes,GT_ATTR_ID_SAM_PASSING_QC);
  if (aux!=NULL) passing_QC = *aux;
  aux = gt_attributes_get(alignment->attributes,GT_ATTR_ID_SAM_PCR_DUPLICATE);
  if (aux!=NULL) PCR_duplicate = *aux;
  // Check unmapped
  const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  if (num_maps==0 || output_attributes->max_printable_maps==0) {
    return gt_output_sam_gprint_sam_se_record(gprinter,
        alignment->tag,alignment->read,alignment->qualities,alignment,NULL,
        false,!passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  }
  // Create a sorted vector with scores to output
  gt_vector* const maps_sort = gt_vector_new(num_maps,sizeof(gt_map_placeholder));
  gt_map_placeholder_create_from_alignment(alignment,maps_sort,PHREAD_SCORE);
  // Print maps !!
  gt_map* primary_map = gt_alignment_get_map_primary(alignment);
  if (primary_map==NULL) {
    primary_map = gt_alignment_get_map(alignment,gt_vector_get_elm(maps_sort,0,gt_map_placeholder)->mmap_position);
  }
  error_code = gt_output_sam_gprint_alignment_map_placeholder(gprinter,
      maps_sort,output_attributes->compact_format,alignment,primary_map,
      !passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  gt_vector_delete(maps_sort);
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_template,gt_template* const template,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_attributes);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_output_sam_gprint_alignment(gprinter,alignment,output_attributes);
  } GT_TEMPLATE_END_REDUCTION;
  gt_status error_code;
  // Get Attributes
  gt_vector* const template_sam_attributes = gt_attribute_sam_fetch_attributes(template->attributes);
  gt_vector* const sam_attributes =
      (template_sam_attributes!=NULL) ? template_sam_attributes : output_attributes->sam_attributes;
  // Get passingQC and PCRDuplicate flags
  bool passing_QC = true, PCR_duplicate = false, *aux;
  aux = gt_attribute_get(template->attributes,GT_ATTR_ID_SAM_PASSING_QC);
  if (aux!=NULL) passing_QC = *aux;
  aux = gt_attribute_get(template->attributes,GT_ATTR_ID_SAM_PCR_DUPLICATE);
  if (aux!=NULL) PCR_duplicate = *aux;
  // Check unmapped
  const uint64_t num_maps = gt_template_get_num_mmaps(template);
  if (num_maps==0 || output_attributes->max_printable_maps==0) {
    gt_alignment* const alignment_end1 = gt_template_get_block(template,0);
    gt_alignment* const alignment_end2 = gt_template_get_block(template,1);
    return gt_output_sam_gprint_sam_pe_record(gprinter,
        template->tag,alignment_end1->read,alignment_end2->read,alignment_end1->qualities,alignment_end2->qualities,
        template,NULL,NULL,
        false,!passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  }
  // Create a sorted vector with scores to output
  gt_vector* const maps_sort = gt_vector_new(num_maps,sizeof(gt_map_placeholder));
  gt_map_placeholder_create_from_template(template,maps_sort,PHREAD_SCORE,output_attributes->print_unpaired_maps,NULL);
  // Find primary mmap
  gt_map** primary_mmap = gt_template_get_mmap_primary(template);
  if (primary_mmap==NULL) {
    // Find first paired mmap (in the vector of placeholders)
    uint64_t first_paired_map_placeholder = 0;
    GT_VECTOR_ITERATE(maps_sort,map_placeholder,map_placeholder_it,gt_map_placeholder) {
      if (map_placeholder->type == GT_MMAP_PLACEHOLDER_PAIRED) {
        first_paired_map_placeholder = map_placeholder_it; break;
      }
    }
    primary_mmap = gt_template_get_mmap(template,
        gt_vector_get_elm(maps_sort,first_paired_map_placeholder,gt_map_placeholder)->mmap_position,NULL);
  }
  error_code = gt_output_sam_gprint_template_map_placeholder(gprinter,
      maps_sort,output_attributes->compact_format,template,primary_mmap,
      !passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  gt_vector_delete(maps_sort);
  return error_code;
  return 0;
}


