/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.c
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Output printers to FASTA/FASTQ/MULTIFASTA
 */

#include "gt_output_fasta.h"

/*
 * FASTA building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fasta,gt_string* const tag,gt_string* const read,gt_shash* attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fastq,gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_shash* attributes);


#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_fasta,gt_string* const tag,gt_string* const read,gt_shash* attributes);
GT_INLINE gt_status gt_output_fasta_gprint_fasta(
    gt_generic_printer* const gprinter,gt_string* const tag,gt_string* const read,gt_shash* attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  GT_STRING_CHECK(read);
  gt_gprintf(gprinter,">"PRIgts"\n",PRIgts_content(tag));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(read));
  // TODO attributes
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_fastq,gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_shash* attributes);
GT_INLINE gt_status gt_output_fasta_gprint_fastq(
    gt_generic_printer* const gprinter,gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_shash* attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  GT_STRING_CHECK(read);
  gt_gprintf(gprinter,"@"PRIgts"\n",PRIgts_content(tag));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(read));
  if (!gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"+\n"PRIgts"\n",PRIgts_content(qualities));
  } else { // Print dummy qualities
    register const uint64_t read_length = gt_string_get_length(read);
    register uint64_t i;
    for (i=0;i<read_length;++i) gt_gprintf(gprinter,"X");
    gt_gprintf(gprinter,"\n");
  }
  // TODO attributes
  return 0;
}


/*
 * FASTA High-level Printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS fasta_format,dna_read
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_dna_read,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
GT_INLINE gt_status gt_output_fasta_gprint_dna_read(
    gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_DNA_READ_CHECK(dna_read);
  switch (fasta_format) {
    case F_FASTA:
      return gt_output_fasta_gprint_fasta(gprinter,dna_read->tag,dna_read->read,dna_read->attributes);
      break;
    case F_FASTQ:
      return gt_output_fasta_gprint_fastq(gprinter,dna_read->tag,dna_read->read,dna_read->qualities,dna_read->attributes);
      break;
    default:
      return -1;
  }
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS fasta_format,alignment
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_alignment,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
GT_INLINE gt_status gt_output_fasta_gprint_alignment(
    gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_alignment* const alignment) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  switch (fasta_format) {
    case F_FASTA:
      return gt_output_fasta_gprint_fasta(gprinter,alignment->tag,alignment->read,alignment->attributes);
      break;
    case F_FASTQ:
      return gt_output_fasta_gprint_fastq(gprinter,alignment->tag,alignment->read,alignment->qualities,alignment->attributes);
      break;
    default:
      return -1;
  }
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS fasta_format,template
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_template,gt_file_fasta_format fasta_format,gt_template* const template);
GT_INLINE gt_status gt_output_fasta_gprint_template(
    gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_template* const template) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  register gt_status error_code;
  GT_TEMPLATE_ALIGNMENT_ITERATE(template,alignment) {
    error_code=gt_output_fasta_gprint_alignment(gprinter,fasta_format,alignment);
    if (error_code) return error_code;
  }
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sequence_archive,column_width
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_sequence_archive,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
GT_INLINE gt_status gt_output_fasta_gprint_sequence_archive(
    gt_generic_printer* const gprinter,gt_sequence_archive* const sequence_archive,const uint64_t column_width) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // Dump the content of the reference file (using Iterators)
  register gt_segmented_sequence* seq;
  gt_sequence_archive_iterator seq_arch_it;
  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
    gt_gprintf(gprinter,">"PRIgts"\n",seq->seq_name);
    gt_segmented_sequence_iterator sequence_iterator;
    gt_segmented_sequence_new_iterator(seq,0,GT_ST_FORWARD,&sequence_iterator);
    register uint64_t chars_written = 0;
    if (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
      while (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
        gt_gprintf(gprinter,"%c",gt_segmented_sequence_iterator_next(&sequence_iterator));
        if ((++chars_written)%column_width==0) gt_gprintf(gprinter,"\n");
      }
      gt_gprintf(gprinter,"\n");
    }
  }
  return 0;
}
