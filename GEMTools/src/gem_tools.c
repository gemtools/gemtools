/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

/*
 * gt.filter menu options
 */
gt_option gt_filter_options[] = {
  /* I/O */
  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (MultiFASTA/FASTA)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (GEM2-Index)" , "" },
  { 200, "annotation", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (GTF Annotation)" , "" },
  { 201, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , false, "" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },
  { 202, "output-format", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "'FASTA'|'MAP'|'SAM' (default='InputFormat')" , "" },
  { 203, "discarded-output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "" , "" },
  { 204, "no-output", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },
  { 205, "check-duplicates", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "Check for duplicated mappings" },
  /* Filter Read/Qualities */
  { 300, "hard-trim", GT_OPT_REQUIRED, GT_OPT_FLOAT, 3 , true, "<left>,<right>" , "" },
  { 301, "quality-trim", GT_OPT_REQUIRED, GT_OPT_FLOAT, 3 , false, "<quality-threshold>,<min-read-length>" , "" },
  { 302, "restore-trim", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "(Previously annotated in the read)" , "" },
  { 303, "uniform-read", GT_OPT_OPTIONAL, GT_OPT_STRING, 3 , true, "['strict']" , "" },
  { 304, "qualities-to-offset-33", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 305, "qualities-to-offset-64", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 306, "remove-qualities", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 307, "add-qualities", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  /* Filter Template/Alignments */
  { 400, "mapped", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , true, "" , "" },
  { 401, "unmapped", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , true, "" , "" },
  { 402, "unique-level", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "<number>|<float>" , "" },
  { 403, "min-length", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<number>" , "" },
  { 404, "max-length", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<number>" , "" },
  { 405, "min-maps", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<number>" , "" },
  { 406, "max-maps", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<number>" , "" },
  { 407, "allow-alignment", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<MapPattern>[,...] (Eg 'Chr1','Chr2:*:1-100')" , "" },
  { 408, "forbid-alignment", GT_OPT_REQUIRED, GT_OPT_INT, 4 , true, "<MapPattern>[,...] (Eg 'Chr1','Chr2')" , "" },
  /* Filter SE-Maps */
  { 500, "first-map", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 'k', "keep-first-map", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 'u', "keep-unique", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 'd', "max-decoded-matches", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number> (stratum-wise)" , "" },
  { 'D', "min-decoded-strata", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number> (stratum-wise)" , "" },
  { 501, "max-output-matches", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number> (to be output, NOT-stratum-wise)" , "" },
  { 502, "max-input-matches", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number> (to be read, stratum-wise)" , "" },
  { 503, "max-strata-after-map", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "" , "" },
  { 504, "make-counters", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 505, "min-strata", GT_OPT_REQUIRED, GT_OPT_FLOAT, 5 , true, "<number>|<float>" , "" },
  { 506, "max-strata", GT_OPT_REQUIRED, GT_OPT_FLOAT, 5 , true, "<number>|<float>" , "" },
  { 507, "min-levenshtein-error", GT_OPT_REQUIRED, GT_OPT_FLOAT, 5 , true, "<number>|<float>" , "" },
  { 508, "max-levenshtein-error", GT_OPT_REQUIRED, GT_OPT_FLOAT, 5 , true, "<number>|<float>" , "" },
  { 509, "map-id", GT_OPT_REQUIRED, GT_OPT_STRING, 5 , true, "<SequenceId>[,...] (Eg 'Chr1','Chr2')" , "" },
  { 510, "strandedness", GT_OPT_REQUIRED, GT_OPT_STRING, 5 , true, "'R'|'F' (default='F,R')" , "" },
  { 511, "filter-quality", GT_OPT_REQUIRED, GT_OPT_STRING, 5 , true, "<min-quality>,<max-quality>" , "" },
  { 512, "reduce-to-unique-strata", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number>" , "" },
  { 513, "reduce-by-quality", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<number> (Quality difference)" , "" },
  { 514, "reduce-by-gene-id", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 515, "reduce-to-max-maps", GT_OPT_REQUIRED, GT_OPT_INT, 5 , true, "<maps>" , "Reduce reads with > <maps> mappings to unmapped" },
  { 516, "reduce-to-pairs", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "Set all unpaired templates to unmapped" },
  { 517, "reduce-to-protein-coding", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 518, "reduce-by-junctions", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "Reduce mappings by checking junction sites. Annotated junctions are preferred." },
  /* Filter RNA-Maps */
  { 600, "no-split-maps", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 6 , true, "" , "" },
  { 601, "only-split-maps", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 6 , true, "" , "" },
  { 's', "no-penalty-for-splitmaps", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 6 , true, "" , "" },
  { 603, "min-intron-length", GT_OPT_REQUIRED, GT_OPT_INT, 6 , true, "<number>" , "" },
  { 604, "min-block-length", GT_OPT_REQUIRED, GT_OPT_INT, 6 , true, "<number>" , "" },
  /* Filter PE-Maps */
  { 700, "pair-strandedness", GT_OPT_REQUIRED, GT_OPT_STRING, 7 , true, "<STRAND>[,...] ('FR'|'RF'|'FF'|'RR')" , "" },
  { 701, "min-inss", GT_OPT_REQUIRED, GT_OPT_STRING, 7 , true, "<number>" , "" },
  { 702, "max-inss", GT_OPT_REQUIRED, GT_OPT_STRING, 7 , true, "<number>" , "" },
  /* Realign/Check */
  { 800, "mismatch-recovery", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , true, "" , "" },
  { 801, "hamming-realign", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , true, "" , "" },
  { 802, "levenshtein-realign", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , true, "" , "" },
  { 'c', "check", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , true, "" , "" },
  { 'C', "check-only", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , false, "(check only, no output)" , "" },
  { 803, "check-format", GT_OPT_REQUIRED, GT_OPT_STRING, 8 , true, "" , "" },
  /* Split/Grouping */
  { 900, "split-reads", GT_OPT_REQUIRED, GT_OPT_NONE, 9 , true, "<number>[,'lines'|'files'] (default=files)" , "" },
  { 901, "sample-read", GT_OPT_REQUIRED, GT_OPT_STRING, 9 , true, "<chunk_size>,<step_size>,<left_trim>,<right_trim>[,<min_remainder>]" , "" },
  { 902, "group-read-chunks", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 9 , true, "" , "" },
  /* Display/Information */
  { 1000, "error-plot", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 10 , false, "" , "" },
  { 1001, "insert-size-plot", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 10 , false, "" , "" },
  { 1002, "sequence-list", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 10 , true, "" , "" },
  { 1003, "display-pretty", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 10 , true, "" , "" },
  /* Misc */
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 11 , true, "" , "" },
#endif
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 11 , true, "" , "" },
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 11 , true, "" , "" },
  { 'H', "help-full", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 11 , false, "" , "" },
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 11 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_filter_options_short = "i:o:r:I:pd:D:Ckst:hHv";
char* gt_filter_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Filter Read/Qualities",
  /*  4 */ "Filter Alignments",
  /*  5 */ "Filter SE-maps",
  /*  6 */ "Filter RNA-Maps",
  /*  7 */ "Filter PE-maps",
  /*  8 */ "Realign/Check",
  /*  9 */ "Split/Grouping",
  /* 10 */ "Display/Information",
  /* 11 */ "Misc"
};

/*
 * gt.stats menu options
 */
gt_option gt_stats_options[] = {
  /* I/O */
  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 200, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , false, "" , "" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , false, "<file> (MultiFASTA/FASTA)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , false, "<file> (GEM2-Index)" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },
  { 'n', "num-reads", GT_OPT_REQUIRED, GT_OPT_INT, 2 , true, "<number>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'f', "output-format", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "'report'|'json'|'both' (default='report')" , "" },
  /* Analysis */
  { 300, "first-map", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'a', "all-tests", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'M', "maps-profile", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'T', "mismatch-transitions", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'Q', "mismatch-quality", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'R', "rna-profile", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  { 'P', "population-profile", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  // { 'D', "indel-profile", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", ""},
  /* MAP Specific */
  { 400, "use-only-decoded-maps", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "(instead of counters)", ""},
  /* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 5, true, "", ""},
#endif
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
  { 'H', "help-full", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , false, "" , "" },
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_stats_options_short = "i:r:I:pn:o:f:aMTQRPDvt:hH";
char* gt_stats_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Analysis",
  /*  4 */ "Misc"
};

/*
 * gt.mapset menu options
 */
gt_option gt_mapset_options[] = {
  /* Operations */
  { 'C', "operation", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true,
      "<operation>\n"
      "     [Set Operators]\n"
      "        union\n"
      "        intersection\n"
      "        difference\n"
      "     [Compare/Display Files]\n"
      "        compare\n"
      "        join\n"
      "        display-compact\n"
      "     [Map Specific]\n"
      "        merge-map\n" , "" },
  /* I/O */
  { 300, "i1", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 301, "i2", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 302, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , false, "" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  /* Compare Function */
  { 's', "files-with-same-reads", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , true, "" , "" },
  { 400, "eq-th", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "<integer>|<float> (Difference tolerated between positions)" , "" },
  { 401, "strict", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , true, "(Strict comparison of mappings)" , "" },
  /* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 5, true, "", ""},
#endif
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_mapset_options_short = "C:po:svt:h";
char* gt_mapset_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "Operations",
  /*  3 */ "I/O",
  /*  4 */ "Compare Function",
  /*  5 */ "Misc"
};

/*
 * gt.scorereads menu options
 */
gt_option gt_scorereads_options[] = {
  /* Operations */
   /* I/O */
  { 300, "i1", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 301, "i2", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 'i', "insert-dist", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 'z', "gzip", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 'j', "bzip2", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 'Z', "no-compress", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 302, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , false, "" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file> (MultiFASTA/FASTA)" , "" },
  { 's', "header", GT_OPT_REQUIRED, GT_OPT_STRING, 3, true, "<file (SAM Header)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file> (GEM2-Index)" , "" },

  /* Score Function */
  { 'q', "quality-format", GT_OPT_REQUIRED, GT_OPT_STRING, 4 , true, "'offset-33'|'offset-64'" , "" },
    { 401, "min-insert", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },
  { 402, "max-insert", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },
  { 403, "indel-penalty", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },
  { 'm', "mismatch-limit", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },
  { 'M', "mapping-quality-cutoff", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },
  { 'S', "split-penalty", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "" , "" },

  /* Optional Fields */
  { 500, "NH", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 503, "XS", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 504, "md", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },

  /* Headers */
  { 600, "read-group-id", GT_OPT_REQUIRED, GT_OPT_STRING, 6 , true, "<read group id>" , "" },

  /* Format */
  { 704, "output-format", GT_OPT_REQUIRED, GT_OPT_STRING, 7 , true, "'MAP'|'SAM' (default='MAP')" , "" },
  { 'c', "compact", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 7 , false, "" , "" },

  /* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8, true, "", ""},
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 8, true, "", ""},
#endif
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8, true, "", ""},
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 8 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_scorereads_options_short = "i:q:pzjcZo:s:S:I:r:M:m:vt:hJ";
char* gt_scorereads_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "Operations",
  /*  3 */ "I/O",
  /*  4 */ "Score Function",
  /*  5 */ "Optional Fields",
  /*  6 */ "Headers",
  /*  7 */ "Format",
  /*  8 */ "Misc",
};

/*
 * gt.map2sam menu options
 */
gt_option gt_map2sam_options[] = {
  /* I/O */
	  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (MultiFASTA/FASTA)" , "" },
  { 's', "header", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file (SAM Header)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (GEM2-Index)" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },
   { 200, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , false, "" , "" },
  /* Headers */
   { 300, "read-group-id", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<read group id>" , "" },
     /* Alignments */
   { 'q', "quality-format", GT_OPT_REQUIRED, GT_OPT_STRING, 4 , true, "'offset-33'|'offset-64'" , "" },
  /* Optional Fields */
  { 500, "NH", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 501, "NM", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 502, "XT", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 503, "XS", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 504, "md", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  { 'Q', "calc-mapq", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
//  { 500, "", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , true, "" , "" },
  /* Format */
  { 'c', "compact", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 6 , false, "" , "" },
  /* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 7, true, "", ""},
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 7, true, "", ""},
#endif
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 7, true, "", ""},
  { 'H', "help-full", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 7 , false, "" , "" },
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 7 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_map2sam_options_short = "i:o:r:I:pq:ct:s:QhHv";
char* gt_map2sam_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Headers",
  /*  4 */ "Alignments",
  /*  5 */ "Optional Fields",
  /*  6 */ "Format",
  /*  7 */ "Misc",
};
/*
 * gt.gtfcount menu options
 */
gt_option gt_gtfcount_options[] = {
  /* I/O */
  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'g', "counts", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "Output file for the gene counts" },
  { 'a', "annotation", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "GTF annotation" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },
  { 'f', "output-format", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "'report'|'json'|'both' (default='report')" , "" },
  /*Counts*/
  { 'w', "weighted", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "Count multi-gene hits (and multi-maps if non unique counts are on) weighted"},
  { 'm', "multi-maps", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "Count multi-maps"},
  //{ 's', "multi-genes", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "Count hits to multiple genes"},
  { 's', "count-single-ends", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "Treat paired reads as single ends for counting (not for matching)"},
  { 'e', "exon-overlap", GT_OPT_REQUIRED, GT_OPT_FLOAT, 3, true, "<overlap>" , "Fraction (0<=overlap<=1) of overlap of the fragment with exon to be counted (default 0, disabled)" },
  /* Misc */
  { 500, "shell", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", "Interactive shell to query the annotation"},
  { 'c', "coverage", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", "Compute coverage profiles (stored in JSON output)"},
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 4, true, "", ""},
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  { 'H', "help-full", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , false, "" , "" },
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_gtfcount_options_short = "i:o:a:p:t:mwhHv";
char* gt_gtfcount_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Counts",
  /*  4 */ "Misc",
};

/*
 * gt.gtfcount menu options
 */
gt_option gt_region_options[] = {
  /* I/O */
  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'a', "annotation", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "GTF annotation" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2 , true, "" , "" },

  /* Misc */
  { 'g', "gene-id", GT_OPT_REQUIRED, GT_OPT_NONE, 5 , true, "" , "" },
  { 'r', "ref-id", GT_OPT_REQUIRED, GT_OPT_NONE, 5 , true, "" , "" },
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 5, true, "", ""},
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
  { 'H', "help-full", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5 , false, "" , "" },
  { 'J', "help-json", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , false, "" , "" },
  {  0, 0, 0, 0, 0, false, "", ""}
};
char* gt_region_options_short = "i:o:a:p:t:hH";
char* gt_region_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Optional Fields",
  /*  4 */ "Format",
  /*  5 */ "Misc",
};



GT_INLINE uint64_t gt_options_get_num_options(const gt_option* const options) {
  uint64_t num_options = 0, i = 0;
  while (options[i++].option_id != 0) ++num_options;
  return num_options;
}
GT_INLINE struct option* gt_options_adaptor_getopt(const gt_option* const options) {
  const uint64_t num_options = gt_options_get_num_options(options);
  struct option* menu_options = gt_malloc(sizeof(struct option)*(1+num_options));
  // Adapt all the records
  uint64_t i = 0;
  for (i=0;i<=num_options;++i) {
    menu_options[i].name = options[i].long_option;
    menu_options[i].has_arg = options[i].option_type;
    menu_options[i].flag = 0;
    menu_options[i].val = options[i].option_id;
  }
  return menu_options;
}
GT_INLINE gt_string* gt_options_adaptor_getopt_short(const gt_option* const options) {
  const uint64_t num_options = gt_options_get_num_options(options);
  gt_string* const options_short = gt_string_new(2*num_options);
  // Adapt all the short options
  uint64_t i = 0;
  for (i=0;i<num_options;++i) {
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      gt_string_append_char(options_short,short_option);
      if (options[i].option_type==GT_OPT_REQUIRED || options[i].option_type==GT_OPT_OPTIONAL) {
        gt_string_append_char(options_short,COLON);
      }
    }
  }
  gt_string_append_eos(options_short);
  return options_short;
}
GT_INLINE void gt_options_fprint_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive) {
  const uint64_t num_options = gt_options_get_num_options(options);
  int64_t i, last_group = -1;
  for (i=0;i<num_options;++i) {
    if (!print_inactive && !options[i].active) continue;
    // Print group (if not printed yet)
    if (last_group!=options[i].group_id) {
      fprintf(stream,"    [%s]\n",groups[options[i].group_id]);
      last_group=options[i].group_id;
    }
    // Print Long Option
    fprintf(stream,"      --%s",options[i].long_option);
    // Print Short Option (if it has)
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      fprintf(stream,"|-%c",short_option);
    }
    // Print extra command line syntax info
    fprintf(stream," %s\n",options[i].command_info);
    // Print description (@print_description)
    if (print_description && !gt_streq(options[i].description,"")) {
      fprintf(stream,"%s",options[i].description);
    }
  }
}
GT_INLINE void gt_options_fprint_json_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive) {
  const uint64_t num_options = gt_options_get_num_options(options);
  int64_t i;
  bool at_least_one_printed = false;
  fprintf(stream,"{ \n"); // Begin JSON record
  fprintf(stream,"\"numOptions\": %"PRIu64",\n",num_options);
  fprintf(stream,"\"options\": [ \n");
  for (i=0;i<num_options;++i) {
    if (!print_inactive && !options[i].active) continue;
    if(at_least_one_printed){
      fprintf(stream,",\n");
    }
    at_least_one_printed = true;
    fprintf(stream,"\t{ \n");
    // Print ID/Short Option
    fprintf(stream,"\t  \"ID\": %d,\n",options[i].option_id);
    // Print Long Option
    fprintf(stream,"\t  \"longOption\": \"%s\",\n",options[i].long_option);
    // Print Short Option
    const char short_option = options[i].option_id;
    if (options[i].option_id<128 && gt_is_alphanumeric(short_option)) {
      fprintf(stream,"\t  \"shortOption\": \"%c\",\n",short_option);
    } else {
      fprintf(stream,"\t  \"shortOption\": null,\n");
    }
    // Group
    fprintf(stream,"\t  \"group\": \"%s\",\n",groups[options[i].group_id]);
    // Option Type
    switch (options[i].option_type) {
      case GT_OPT_NO_ARGUMENT: fprintf(stream,"\t  \"optionType\": \"noArgument\",\n"); break;
      case GT_OPT_REQUIRED: fprintf(stream,"\t  \"optionType\": \"required\",\n"); break;
      case GT_OPT_OPTIONAL: fprintf(stream,"\t  \"optionType\": \"optional\",\n"); break;
    }
    // Argument Type
    switch (options[i].argument_type) {
      case GT_OPT_NONE: fprintf(stream,"\t  \"argumentType\": null,\n"); break;
      case GT_OPT_INT: fprintf(stream,"\t  \"argumentType\": \"int\",\n"); break;
      case GT_OPT_FLOAT: fprintf(stream,"\t  \"argumentType\": \"float\",\n"); break;
      case GT_OPT_CHAR: fprintf(stream,"\t  \"argumentType\": \"char\",\n"); break;
      case GT_OPT_STRING: fprintf(stream,"\t  \"argumentType\": \"string\",\n"); break;
      case GT_OPT_BOOL: fprintf(stream,"\t  \"argumentType\": \"bool\",\n"); break;
    }
    // Print extra command line syntax info
    fprintf(stream,"\t  \"commandInfo\": \"%s\"",options[i].command_info);
    // Print description (@print_description)
    if (print_description && !gt_streq(options[i].description,"")) {
      fprintf(stream,",\n\t  \"description\": \"%s\"",options[i].description);
    }
    fprintf(stream,"\n\t}");
  }
  fprintf(stream,"\n    ]\n");
  fprintf(stream,"}\n");
}


