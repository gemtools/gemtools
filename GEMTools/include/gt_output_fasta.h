/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_fasta.h
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Output printers to FASTA/FASTQ/MULTIFASTA
 */

#ifndef GT_OUTPUT_FASTA_H_
#define GT_OUTPUT_FASTA_H_

#include "gt_commons.h"
#include "gt_dna_read.h"
#include "gt_sequence_archive.h"
#include "gt_input_file.h"

#include "gt_alignment.h"
#include "gt_template.h"

#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"


/*
 * FASTA building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fasta,gt_string* const tag,gt_string* const read,gt_shash* attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fastq,gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_shash* attributes);

/*
 * FASTA High-level Printers
 */

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_dna_read,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_gprint_dna_read(gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_fprint_dna_read(FILE* file,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_sprint_dna_read(gt_string* const string,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_bprint_dna_read(gt_output_buffer* const output_buffer,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_ofprint_dna_read(gt_output_file* const output_file,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);
inline gt_status gt_output_fasta_bofprint_dna_read(gt_buffered_output_file* const buffered_output_file,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_alignment,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_gprint_alignment(gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_fprint_alignment(FILE* file,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_sprint_alignment(gt_string* const string,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_bprint_alignment(gt_output_buffer* const output_buffer,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_ofprint_alignment(gt_output_file* const output_file,gt_file_fasta_format fasta_format,gt_alignment* const alignment);
inline gt_status gt_output_fasta_bofprint_alignment(gt_buffered_output_file* const buffered_output_file,gt_file_fasta_format fasta_format,gt_alignment* const alignment);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_template,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_gprint_template(gt_generic_printer* const gprinter,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_fprint_template(FILE* file,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_sprint_template(gt_string* const string,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_bprint_template(gt_output_buffer* const output_buffer,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_ofprint_template(gt_output_file* const output_file,gt_file_fasta_format fasta_format,gt_template* const template);
inline gt_status gt_output_fasta_bofprint_template(gt_buffered_output_file* const buffered_output_file,gt_file_fasta_format fasta_format,gt_template* const template);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_sequence_archive,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_gprint_sequence_archive(gt_generic_printer* const gprinter,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_fprint_sequence_archive(FILE* file,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_sprint_sequence_archive(gt_string* const string,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_bprint_sequence_archive(gt_output_buffer* const output_buffer,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_ofprint_sequence_archive(gt_output_file* const output_file,gt_sequence_archive* const sequence_archive,const uint64_t column_width);
inline gt_status gt_output_fasta_bofprint_sequence_archive(gt_buffered_output_file* const buffered_output_file,gt_sequence_archive* const sequence_archive,const uint64_t column_width);

#endif /* GT_OUTPUT_FASTA_H_ */
