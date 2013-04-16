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
  gt_output_sam_attributes_reset_defaults(attributes);
  return attributes;
}
GT_INLINE void gt_output_sam_attributes_delete(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  gt_free(attributes);
}
GT_INLINE void gt_output_sam_attributes_reset_defaults(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  /* Maps */
  attributes->max_printable_maps = UINT64_MAX;
  attributes->print_unpaired_maps = false;
  attributes->compact_format = false;
  /* Mismatch/CIGAR string */
  attributes->print_misms = false;
  /* Optional fields */
  attributes->print_casava = true;
  attributes->sam_attributes = NULL;
}
/* Maps */
GT_INLINE void gt_output_sam_attributes_flag_best_as_primary_on(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->flag_best_as_primary = true;
}
GT_INLINE void gt_output_sam_attributes_flag_best_as_primary_off(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->flag_best_as_primary = false;
}
GT_INLINE uint64_t gt_output_sam_attributes_get_max_printable_maps(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->max_printable_maps;
}
GT_INLINE void gt_output_sam_attributes_set_max_printable_maps(gt_output_sam_attributes* const attributes,const uint64_t max_printable_maps) {
  GT_NULL_CHECK(attributes);
  attributes->max_printable_maps = max_printable_maps;
}
GT_INLINE void gt_output_sam_attributes_print_unpaired_maps_on(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_unpaired_maps = true;
}
GT_INLINE void gt_output_sam_attributes_print_unpaired_maps_off(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_unpaired_maps = false;
}
GT_INLINE void gt_output_sam_attributes_compact_format_on(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->compact_format = true;
}
GT_INLINE void gt_output_sam_attributes_compact_format_off(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->compact_format = false;
}
/* CIGAR/Mismatch string */
GT_INLINE void gt_output_sam_attributes_print_mismatches_on(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_misms = true;
}
GT_INLINE void gt_output_sam_attributes_print_mismatches_off(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_misms = false;
}
/* Optional fields */
GT_INLINE void gt_output_sam_attributes_print_casava_on(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_casava = true;
}
GT_INLINE void gt_output_sam_attributes_print_casava_off(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->print_casava = false;
}
GT_INLINE void gt_output_sam_attributes_set_sam_attributes(gt_output_sam_attributes* const attributes,gt_vector* const sam_attributes) {
  GT_NULL_CHECK(attributes);
  attributes->sam_attributes = sam_attributes;
}

/*
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
 * @PG  ID:tmap  CL:map4 -f /home/devel/mdabad/scratch/references/HSA/TMAP/hsapiens_v37.fa -r /home/devel/mdabad/scratch/READS/DEF/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq -i fastq -s /home/devel/mdabad/scratch/benchmark/DEF/miseq_S/OUT/TMAP.map41.H.Sapiens.1M.S.MiSeq.low.l250_1.fastq.sam   VN:3.0.1
 * @PG  ID:dvtgm PN:stampy     VN:1.0.17_(r1481)  CL:-g /home/devel/mdabad/scratch/references/HSA/STAMPY/hsapiens_v37 -h /home/devel/mdabad/scratch/references/HSA/STAMPY/hsapiens_v37 --bwaoptions=/home/devel/mdabad/scratch/references/HSA/BWA/hsapiens_v37.fa --substitutionrate=0.08 --maxbasequal 90 -M /home/devel/mdabad/scratch/READS/DEF/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq
 * @PG  ID:GEM   PN:gem-2-sam  VN:1.414
 *
 * @CO  TM:Fri, 30 Nov 2012 14:14:13 CET        WD:/scratch/devel/mdabad/benchmark/DEF/miseq_S/CMD      HN:cn38.bullx   UN:mdabad
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
      gt_output_sam_calculate_flag_se(true,map_end1->strand,secondary_alignment,not_passing_QC,PCR_duplicate): // Mapped
      gt_output_sam_calculate_flag_se(false,NULL,secondary_alignment,not_passing_QC,PCR_duplicate); // Unmapped
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe_map(
    gt_map* const map,gt_map* const mate,const bool is_map_first_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return gt_output_sam_calculate_flag_pe(
      map!=NULL && mate!=NULL, /* read_paired */
      map!=NULL,               /* read_mapped */
      mate!=NULL,              /* mate_strand */
      (map!=NULL) ? map->strand : NULL,   /* read_strand */
      (mate!=NULL) ? mate->strand : NULL, /* mate_strand */
      is_map_first_in_pair,    /* first_in_pair */
      !is_map_first_in_pair,   /* last_in_pair */
      secondary_alignment,not_passing_QC,PCR_duplicate);
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe(
    const bool read_paired,const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
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
  GT_GENERIC_PRINTER_CHECK(gprinter);
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
  uint16_t sam_flag = 0;
  return (paired_end) ? /* 0x1 */
      gt_output_sam_calculate_flag_pe(
          read_paired,read_mapped,mate_mapped,
          read_strand,mate_strand,first_in_pair,last_in_pair,
          secondary_alignment,not_passing_QC,PCR_duplicate): // PE
      gt_output_sam_calculate_flag_se(read_mapped,read_strand,
          secondary_alignment,not_passing_QC,PCR_duplicate); // SE
}
/*
 * SAM CIGAR
 */
#define GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"M",misms_pos-centinel); \
    centinel = misms_pos; \
  }
#define GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"M",centinel-misms_pos); \
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
        if (attributes->print_misms) {
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
  GT_MAP_MISMS_ITERATOR(map,misms,misms_n) {
    const uint64_t misms_pos = gt_misms_get_position(misms);
    switch (misms->misms_type) {
      case MISMS:
        if (attributes->print_misms) {
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
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_cigar,gt_map* const map,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_cigar(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  if (gt_map_get_strand(map)==REVERSE) {
    if (gt_map_has_next_block(map)) {
      gt_output_sam_gprint_cigar(gprinter,gt_map_get_next_block(map),attributes);
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map));
    }
    gt_output_sam_gprint_map_block_cigar_reverse(gprinter,map,attributes);
  } else {
    gt_output_sam_gprint_map_block_cigar_forward(gprinter,map,attributes);
    if (gt_map_has_next_block(map)) {
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map));
      gt_output_sam_gprint_cigar(gprinter,gt_map_get_next_block(map),attributes);
    }
  }
  return 0;
}
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 *   @read and @qualities must be already RC in case of mapping on the reverse strand
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_output_sam_gprint_flag(gprinter,"\t%"PRId16,gt_output_sam_calculate_flag_se_map(map,secondary_alignment,not_passing_QC,PCR_duplicate));
  // Is mapped?
  if (gt_expect_true(map!=NULL)) {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRId8"\t",
        PRIgts_content(gt_map_get_string_seq_name(map)),gt_map_get_global_position(map),gt_map_get_phred_score(map));
    // (6) Print CIGAR
    gt_output_sam_gprint_cigar(gprinter,map,attributes);
  } else {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    // (6) Print CIGAR
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  // (7) Print RNEXT
  // (8) Print PNEXT
  // (9) Print TLEN
  // (10) Print SEQ
  // (11) Print QUAL
  if (read!=NULL && qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t"PRIgts,PRIgts_content(read),PRIgts_content(qualities));
  } else if (read!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t*",PRIgts_content(read));
  } else if (qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t0\t0\t*\t"PRIgts,PRIgts_content(qualities));
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map,mate,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,gt_map* const mate,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map,gt_map* const mate,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_output_sam_gprint_flag(gprinter,"\t%"PRId16,
      gt_output_sam_calculate_flag_pe_map(map,mate,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate));
  // (3) Print RNAME
  // (4) Print POS
  // (5) Print MAPQ
  // (6) Print CIGAR
  if (map!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRId8"\t");
    gt_output_sam_gprint_cigar(gprinter,map,attributes); // CIGAR
  } else {
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  // Calculate Template Length
  const bool same_sequence = gt_string_equals(gt_map_get_string_seq_name(map),gt_map_get_string_seq_name(mate));
  int64_t template_length = 0;
  if (same_sequence && map!=NULL && mate!=NULL) {
    uint64_t error;
    gt_map* map[2] = {map, mate};
    template_length = gt_template_get_insert_size(map,&error);
  }
  // (7) Print RNEXT
  // (8) Print PNEXT
  // (9) Print TLEN
  if (mate!=NULL) {
    if (same_sequence) {
      gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRIu64,PRIgts_content(gt_map_get_string_seq_name(map)),gt_map_get_global_position(map),template_length);
    } else {
      gt_gprintf(gprinter,"\t=\t%"PRIu64"\t%"PRIu64,gt_map_get_global_position(map),template_length);
    }
  } else {
    gt_gprintf(gprinter,"\t*\t0\t0");
  }
  // (10) Print SEQ
  // (11) Print QUAL
  if (read!=NULL && qualities!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t"PRIgts,PRIgts_content(read),PRIgts_content(qualities));
  } else if (read!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t*",PRIgts_content(read));
  } else if (qualities!=NULL) {
    gt_gprintf(gprinter,"\t*\t"PRIgts,PRIgts_content(qualities));
  }
  return 0;
}

/*
 * Print Optional fields
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
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attributes,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields_values,
    gt_vector* sam_attributes,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields_values(gt_generic_printer* const gprinter,
    gt_vector* sam_attributes,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(sam_attributes);
  GT_ATTRIBUTES_SAM_ITERATE(attributes,sam_attribute) {
    if (sam_attribute->attribute_type == SAM_ATTR_INT_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_value);
    } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->d_value);
    } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(sam_attribute->s_value));
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,alignment,map,sam_attributes,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields,
    gt_template* const template,gt_alignment* const alignment,gt_map** const map,
    gt_vector* sam_attributes,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields(gt_generic_printer* const gprinter,
    gt_template* const template,gt_alignment* const alignment,gt_map** const map,
    gt_vector* sam_attributes,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(sam_attributes);
  GT_ATTRIBUTES_SAM_ITERATE(attributes,sam_attribute) {
    // Values
    if (sam_attribute->attribute_type == SAM_ATTR_INT_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_value);
    } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->d_value);
    } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_VALUE) {
      gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(sam_attribute->s_value));
    } else
    // Functions
    if (sam_attribute->attribute_type == SAM_ATTR_INT_FUNC) {
      gt_gprintf(gprinter,"\t%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_func(template,alignment,map));
    } else if (sam_attribute->attribute_type == SAM_ATTR_FLOAT_FUNC) {
      gt_gprintf(gprinter,"\t%c%c:%c:%3.2f",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->d_func(template,alignment,map));
    } else if (sam_attribute->attribute_type == SAM_ATTR_STRING_FUNC) {
      gt_string* const aux = sam_attribute->s_func(template,alignment,map);
      gt_gprintf(gprinter,"\t%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(aux));
      gt_string_delete(aux);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_op_casava(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  gt_string* const casava_string = gt_attribute_get(alignment->attributes,GT_ATTR_ID_TAG_CASAVA);
  if (casava_string!=NULL) gt_gprintf(gprinter,"\tcs:Z"PRIgts,PRIgts_content(casava_string));
  return 0;
}
/*
 * SAM Record (Full SAM record printers)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read_end1,read_end2,qualities_end1,qualities_end2,\
  template,alignment_end1,alignment_end2,map_end1,map_end2,secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_sam_pe_record,
    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
    gt_string* const qualities_end1,gt_string* const qualities_end2,
    gt_template* const template,gt_alignment* const alignment_end1,gt_alignment* const alignment_end2,
    gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_sam_pe_record(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
    gt_string* const qualities_end1,gt_string* const qualities_end2,
    gt_template* const template,gt_alignment* const alignment_end1,gt_alignment* const alignment_end2,
    gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  /*
   * Print SAM Record (end1)
   */
  if (map_end1!=NULL) {
    gt_map* mmap[2] = {map_end1,map_end2};
    // Print CORE fields
    gt_output_sam_gprint_map_pe(gprinter,tag,read_end1,qualities_end1,map_end1,map_end2,
        true,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes);
    // Print Optional Fields
    gt_output_sam_gprint_optional_fields(gprinter,template,alignment_end1,mmap,sam_attributes,output_attributes);
    // Print Xtras
    if (output_attributes->print_casava && alignment_end1!=NULL) gt_output_sam_gprint_op_casava(gprinter,alignment_end1,output_attributes);
    gt_gprintf(gprinter,"\n");
  }
  /*
   * Print SAM Record (end2)
   */
  if (map_end2!=NULL) {
    gt_map* mmap[2] = {map_end2,map_end1};
    // Print CORE fields
    gt_output_sam_gprint_map_pe(gprinter,tag,read_end2,qualities_end2,map_end2,map_end1,
        true,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes);
    // Print Optional Fields
    gt_output_sam_gprint_optional_fields(gprinter,template,alignment_end2,mmap,sam_attributes,output_attributes);
    // Print Xtras
    if (output_attributes->print_casava && alignment_end2!=NULL) gt_output_sam_gprint_op_casava(gprinter,alignment_end2,output_attributes);
    gt_gprintf(gprinter,"\n");
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,template,alignment,map,secondary_alignment,\
  not_passing_QC,PCR_duplicate,sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_sam_se_record,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_sam_se_record(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  /*
   * Print SAM Record
   */
  gt_output_sam_gprint_map_se(gprinter,tag,read,qualities,map,
      secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes); // Print CORE fields
  // Print Optional Fields
  if (sam_attributes) gt_output_sam_gprint_optional_fields(gprinter,template,alignment_end1,&map,sam_attributes,output_attributes);
  // Print Xtras
  if (output_attributes->print_casava && alignment_end1!=NULL) gt_output_sam_gprint_op_casava(gprinter,alignment,output_attributes);
  gt_gprintf(gprinter,"\n");
  return 0;
}

/*
 * SAM High-level MMAp/Map Printers
 *   - Provides flexibility to print maps individually
 *   - @tag, @read and @qualities are used to print the SAM field
 *       (@template,@alignment can be null and are only used as to generate the opt-fields, if any)
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
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,map_end1,map_end2,secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_mmap,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_mmap(gt_generic_printer* const gprinter,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* const sam_attributes,gt_output_sam_attributes* const output_attributes) {
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
      alignment_end1,alignment_end2,map_end1,map_end2,
      secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  if (read_rc_end1) gt_string_delete(read_rc_end1);
  if (qualities_r_end1) gt_string_delete(qualities_r_end1);
  if (read_rc_end2) gt_string_delete(read_rc_end2);
  if (qualities_r_end2) gt_string_delete(qualities_r_end2);
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,alignment,map,secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map,
    gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_map(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_vector* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  // Produce RC of the mapping's read/qualities
  bool rc_complement;
  gt_string *read_f=NULL, *qualities_f=NULL, *read_rc=NULL, *qualities_r=NULL;
  GT_OUTPUT_SAM_GENERATE_RC_READ__QUALITIES(alignment,map,read_f,qualities_f,read_rc,qualities_r,rc_complement);
  // Print SAM record
  gt_output_sam_gprint_sam_se_record(gprinter,template->tag,
      ((rc_complement) ? read_rc : read_f),((rc_complement) ? qualities_r : qualities_f),
      NULL,alignment,map,secondary_alignment,not_passing_QC,PCR_duplicate,sam_attributes,output_attributes);
  // Free
  if (read_rc) gt_string_delete(read_rc);
  if (qualities_r) gt_string_delete(qualities_r);
  return 0;
}

/*
 * SAM High-level Template/Alignment Printers
 *   - Depending on @output_attributes->flag_best_as_primary
 *       - True. The Map/MMap with best phredScore is output as primary
 *       - False. Maps/MMaps with attribute
 *   - Optional fields are generated from the first list of SAM Attributes found in the following order:
 *       1.- @attributes->sam_attributes
 *       2.- template->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       3.- alignment->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       4.- map->attributes{GT_ATTR_ID_SAM_FLAGS}
 *     If multiple SAM-Attributes lists are contained within these levels, only the first found one is output
 *     (completely preventing lower levels of attributes to be output)
 */
int gt_output_sam_cmp_maps_phred_scores(const void *sa, const void *sb) {
  const uint8_t score_a = ((gt_map*)sa)->phred_score;
  const uint8_t score_b = ((gt_map*)sb)->phred_score;
  // No scored
  if (score_a==GT_MAP_NO_PHRED_SCORE || score_b==GT_MAP_NO_PHRED_SCORE) {
    if (score_a!=GT_MAP_NO_PHRED_SCORE) return -1;
    if (score_b!=GT_MAP_NO_PHRED_SCORE) return 1;
    return 0;
  }
  // Cmp regular scores
  if (score_a>score_b) {
    return -1;
  } else if (score_b>score_a) {
    return 1;
  } else {
    return 0;
  }
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_alignment,gt_alignment* const alignment,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_attributes);
  const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  // Find attributes
  gt_vector* sam_attributes;
  if (output_attributes->sam_attributes!=NULL) {
    sam_attributes=output_attributes->sam_attributes;
  } else {
    sam_attributes=gt_attribute_sam_fetch_attributes(alignment->attributes);
  }
  // Check unmapped
  if (num_maps==0) {
    return gt_output_sam_gprint_sam_se_record(gprinter,
        alignment->tag,alignment->read,alignment->qualities,NULL,alignment,NULL,
        true,false,false,sam_attributes,output_attributes);
  }
  // TODO Compact


  // Produce RC of the mapping's read/qualities
  gt_string *read_f=alignment->read, *qualities_f=alignment->qualities;
  gt_string *read_rc=NULL, *qualities_r=NULL;
  // Get passingQC and PCRDuplicate flags
  bool passing_QC = true, PCR_duplicate = false, aux;
  aux = gt_attribute_get(alignment->attributes,GT_ATTR_ID_SAM_PASSING_QC);
  if (aux!=NULL) passing_QC = *aux;
  aux = gt_attribute_get(alignment->attributes,GT_ATTR_ID_SAM_PCR_DUPLICATE);
  if (aux!=NULL) PCR_duplicate = *aux;
  // Create a vector with scores as to sort and output
  gt_vector* const maps_order = gt_vector_new(num_maps,sizeof(gt_map*));
  gt_vector_set_used(maps_order,num_maps);
  GT_VECTOR_ITERATE(maps_order,map,pos,gt_map*) {
    (*map) = gt_alignment_get_map(alignment,pos);
  }
  qsort(gt_vector_get_mem(maps_order,void),num_maps,sizeof(gt_map*),gt_output_sam_cmp_maps_phred_scores); // Sort !
  // Traverse all maps
  gt_status error_code;
  const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  uint64_t map_it = 0;
  for (map_it=0;map_it<num_maps;++map_it) {
    if (num_map > output_attributes->max_printable_maps) break;
    gt_map* const map = *gt_vector_get_elm(maps_order,map_it,gt_map*);
    // Set primary/secondary
    bool primary_alignment;
    if (output_attributes->flag_best_as_primary) {
      primary_alignment = (num_map==0);
    } else {
      aux = gt_attribute_get(alignment->attributes,GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT);
      if (aux!=NULL) primary_alignment = *aux;
    }
    // Print SAM record
    if (gt_map_get_strand(map)==FORWARD) {
      error_code=gt_output_sam_gprint_sam_se_record(gprinter,
          alignment->tag,read_f,qualities_f,NULL,alignment,map,
          !primary_alignment,passing_QC,PCR_duplicate,sam_attributes,output_attributes);
    } else {
      if (gt_expect_false(read_rc==NULL)) { // Check RC
        read_rc = gt_dna_string_reverse_complement_dup(read_f);
        qualities_r = gt_string_reverse_dup(qualities_f);
      }
      error_code=gt_output_sam_gprint_sam_se_record(gprinter,
          alignment->tag,read_rc,qualities_r,NULL,alignment,map,
          !primary_alignment,passing_QC,PCR_duplicate,sam_attributes,output_attributes);
    }
    if (error_code) return error_code;
  }
  // Free
  gt_vector_delete(maps_order);
  if (read_rc!=NULL) {
    gt_string_delete(read_rc);
    gt_string_delete(qualities_r);
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_template,gt_template* const template,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_attributes);
  // TODO
  return 0;
}

