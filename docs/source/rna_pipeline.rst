.. _rna_pipeline:

RNA Pipeline
============

General remarks
---------------
To map things during this tutorial, we are going to use the GEM (GEnome Multi-tool)
programs for short-read processing (see http://gemlibrary.sourceforge.net). Several
programs are provided (continuous mapper, split mapper, etc.) all based on a common
very optimized alignment library. You find more information about the mapper in
http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.2221.html.
In general, alignment is performed by querying an FM-index. This requires the index to
be pre-generated (command gem-indexer) and then loaded in memory (takes 1-2 times
the length of the sequence to be indexed depending on the sampling rate). After that
you map reads with either the gem-mapper or the gem-split-mapper, or both. After
mapping you can re-reprocess the output —which is produced in a format proprietary to
GEM, as SAM is not general enough and much bulkier— to score/select/combine/verify
matches and create pipelines (command gem-map-2-map). Finally, you would convert to
SAM (command gem-2-sam).
Typical common options for all programs are -I (index file), -i (input file) and -o
(output file). Input/output file can usually be omitted, as all the mappers work as filters
(accepting input from stdin and writing the results on stdout — commands can be
piped). In most cases programs are multi-threaded to improve wall-clock time, you just
need to specify the option -T.
In this tutorial we examine in some detail the GEM commends. However, if you are using
the latest GEM distributions you do not actually need to learn all the details: a toplevel
pipeline named gemtools is provided which performs for you all the steps needed to do
RNA-seq mapping. You will just need to learn how to configure it.


Output Format
-------------
The typical output format of GEM looks like this:

    HW I - S T 6 6 1 : 1 3 0 : C 0 3 7 K A C X X : 8 : 1 1 0 1 : 1 5 4 6 : 2 1 6 2 t a b
    N G T T N A C A C C A A C T T G C G C A A A A A C A A C A G A C A G C C C T A T G C T G T C A G T G AA T T A G C A G G T C A T C A G A C T A G T G Cs p a c e
    C T C G A A C T C T G G G A A T T C G A G C C A C A G C T C T G C C A G T A C C C C A A G A C T C A GC A C T A G T C T G A T G A C C T G C T A A T Tt a b
    # 0 ; ? # 2 8 > ? @ ? > ? ? @ @ ? ? ? = < @ @ ? ? ? < < < > > ? 7 > ? ? ? ? ? > < ? ? ? > ? ? ? ? ? ?? ? ? ? ? ? ? @ @ > ? ? = ? ? > ; = ; = > > ? >s p a c e
    @ C @ F D F F F H G H H D F D H I I I G H I J J > G C F G I I I J J I J I I I I I I J I G I G I G G BH < F F G D @ F G I H J > D G H H G G G E @ ? Ct a b
    0 : 0 : 0 : 0 : 0 : 0 : 0 : 1 : 1 t a b
    c h r 3 : - : 1 8 5 1 3 6 3 7 6 : T 3 C 3 4 C 4 C 2 4 C 5 : : c h r 3 : + : 1 8 5 1 3 6 3 2 6 : 2 4 T3 0 G 1 9 : : : 3 2 5 7 6 ,
    c h r 1 5 : - : 6 6 7 9 5 4 8 7 : T 3 C 5 4 > 1 9 3 * 1 6 : : c h r 1 5 : + : 6 6 7 9 5 4 3 7 : 6 8 AA A 1 T 1 G : : : 1 3 8 8 8

There are 5 tab-separated fields:

Read name. When paired alignments found, the trailing /1 and /2 are removed
Read sequence. When more than one end present and paired alignments found, the
sequences are presented one after the other, separated by spaces
Read qualities. When more than one end present and paired alignments found, the
sequences are presented one after the other, separated by spaces

Match summary. Counts of the alignments found for the read. A list of colonseparated numbers. The first number describes how many matches have been
found having error 0, the second number how many alignments have been found
having error 1, and so on. Occasionally a + can be found in place of a colon.
Alignments. Either - (means ``no alignment found'') or a list of comma-separated
alignments.

Each alignment has the following structure:
Matching blocks. A list of one or more matching blocks, separated by a doublecolon sign `::'. Each matching block describes one or more alignment operations
performed on the same (sequence,strand) couple.
Annotations. An optional sequence of annotations, starting by a triple-colon
sign `:::', and separated by colons. They encode additional global information for the alignment, typically its quality.
Each matching block has the following fields, separated by colon signs:
Matching sequence name. The name of the sequence in the reference the block is
matching to.
Matching sequence strand. The strand (`+' for forward, and `-' for reverse) of the
sequence in the reference the block is matching to.
Position. The lowest position (in 1-based strand-independent format) in the reference
sequence the block is matching to.
``GIGAR'' string. The ``GEM-CIGAR'', a description of the alignment in terms of basic
operations. Such operations are written consecutively without any separator, and
can be:
Match. A number. For instance, 75 would mean that 75 bases have been
matched exactly.
Substitution. A letter. Indicates that in the reference that letter is present
instead of the corresponding one present in the read.
Skip. A sign `>' followed by a number n, and by one of the symbols `+', `-', `*'
and `/'. Such four possible kinds of skips are:
`+'. Skip of n positions in the genome, according to the natural direction
of each strand (forward or left-to-right if the match is on the + strand,
and backward or right-to-left if the match is on the - strand). Note
that the skip may well be negative (for instance in the case of two
partially overlapping sub-blocks).
`-'. Skip of n positions in the read, walking the read from left to right.
`*'. Splice of n positions (that is, as a skip of n positions in the genome,
but with the additional constraint that a suitable splice-site consensus is obeyed by the sequences flanking the skipped block).
One should note that, differently from the ordinary notion of ``insertion/deletion'', this notation is unambiguous, as no confusion is possible
about where the indel is present (either read or sequence).
Trim. The notation `(' followed by a number n followed by a symbol `)' represents a trim, that is a `-'-skip in the read with the additional constraint that
the skip happens either at the beginning or at the end of the read.
2An example of aligned block represented as a GIGAR string is the following:
T3C54>193*16
means ``replace a T at the beginning of the read to get the reference, match 3 bases
exactly, replace the next one with a C, match 54 bases, skip 193 bases in the
genome (it is a splice, and the consensus is obeyed), and finally match 16 bases
exactly.
Annotations. An optional sequence of annotations, separated by colons. They encode
additional information for the alignment block.
An example of alignment might thus be:
chr15:-:66795487:T3C54>193*16::chr15:+:66795437:68AAA1T1G:::13888
that is: the first read fragment (the first read end in this case) matches to the negative
strand of chromosome 15 at position 66795487, with the alignment seen in the example
above; the second read fragment (the second read end in this case) matches to the
positive strand of the same chromosome, at position 66795437 (68 bases match exactly,
then several mismatches are present). The overall quality for the paired end alignment
is 13888.
One should note that:
• the GEM mapper always returns correct counters for all the matches up to some
given stratum as per the mapping parameters specified by the user, although
not all the matches might have been printed (depends on the command-line
configuration)
• the matches are sorted by increasing error (first come all the matches with error=0,
if any, the all the matches with error=1, and so on)
• any alignment can be represented in this format, without limitations (in particular,
the GEM format is much more powerful than the SAM format). Hence, one should
delay the transformation to SAM as much as possible, as it implies a loss of
information
• when a few matches are present (as typical in the case of longer reads) this format
is several times more compact that SAM
• tasks like SNP-calling can be performed without consulting the reference, thanks
to the definition of the GIGAR string.
Installation
Nothing special is required, just unpack your GEM binaries into a directory of your
choice. Then set the PATH so as to include that directory:
export PATH=/my/GEM/directory/:$PATH
Several kinds of binaries exist, some (the core_i3 series) providing better performance
than others (the core_2 series). However, the most optimized ones might crash on your
machine if your CPU does not provide support for some hardware instructions.
Generating the index
You have to run gem-indexer, which in turn will automatically call other programs to
do the job. The input format is FASTA, the output a single file containing a GEM archive
(extension .gem). In our case

g e m - i n d e x e r - i d m e l - a l l - c h r o m o s o m e - r 5 . 4 8 . f a s t a - o d m e l - a l l - c h r o m o s o m e - r 5 . 4 8
with no special options will do the job. The only thing one might wish to tune is whether
both strands or just one are indexed.
Continuous mapping
Once the index has been generated, running the mapper is pretty much straightforward.
The mapper can perform both paired-end and single-end alignment, depending on your
input. A typical command line might be (for single-end mapping):
g e m - m a p p e r - I d m e l - a l l - c h r o m o s o m e - r 5 . 4 8 . g e m - i T E S T . f a s t q - q o f f s e t - 3 3 - m 4 - s 1 - e 8
- - f a s t - m a p p i n g = 0 - - m i n - d e c o d e d - s t r a t a 2 - T 8
Most relevant parameters are:
-q. quality type, necessary when FASTQ input
-m. number or fraction of mismatches
-e. number or fraction of errors
-s. number of match strata to be reported
-T. number of threads.
In the case of paired-end mapping things go pretty much the same way, but for a few
parameters more to specify how the mapper should handle pairing-related issues:
g e m - m a p p e r - I d m e l - a l l - c h r o m o s o m e - r 5 . 4 8 . g e m - i T E S T . f a s t q - q o f f s e t - 3 3 - m 4 - s 1 - e 8
- - f a s t - m a p p i n g = 0 - - m i n - d e c o d e d - s t r a t a 2 - T 8 - p - - m a x - i n s e r t - s i z e = 6 0 0 . . .
The fine print
There are many more subtleties about how you can tune the behavior of the GEM
mapper, in particular with respect to paired-end mapping. As a general concept,
the GEM mapper can map a pair either by first finding matches for one of the ends
and subsequently extending them to the other end, or by first mapping the two ends
independently and then recombining the matches for the two ends. In our case (RNA
mapping) we will map the two ends of each pair separately, which is the simplest strategy
(more control, less tuning — slower in the case of DNA mapping, the only possibility in
the case of RNA).
Hence we will use the gem-mapper in several invocations:
1. as a single-end mapper (two times) to map end 1 and end 2 of each read
2. as a pairer, giving it the final results of the mapping pipeline for both ends
(the mapper decides what to do depending on the format of the input).
General remarks about how GEMmaps reads
GEM is quite different with respect to other mappers, in that it always reports all
matches that exist within the alignment parameters specified by the user. No arbitrary
choice is imposed on you — in the basic output of the mapper there is no such a thing
as a ``probabilistic score'' or a ``best match'', as those concepts depend on priors which
might not be the correct ones for your situation. You can rescore the matches by your
score of choice later on, though.

Also, you can tune virtually all the alignment parameters as you like best. But remember:
privilege implies responsibility! You will get better results, but you have to understand
how the basics of the mapper work. We believe this is the correct approach to get the
most out of your analysis (the field is full of biases generated by the fact that mappers
work as ``black boxes'').
RNAmapping
Mapping RNA data is definitely more complicated, as it requires
1. either a good knowledge of the annotation (splice junctions) — you can then
perform continuous mapping to the annotated transcriptome
2. or the ability of performing spliced mapping — you can then reconstruct/find denovo splice junctions out of splitting reads which have a sufficient support, and
go back to step (1).
To map our example dataset I have used a pipeline similar to the one that mapped the
Geuvadis data. Actually it does several complicated things one after the other:

A few words of caution. How this is done is not that relevant, what is done is. Inflexible
setups where all the decisions have already been taken for you by somebody else are
dangerous, so one should understand and master the data flow — the good data analyst
should always be in control. This one is not necessarily the best analysis pipeline for all
cases: some steps require decisions (pairing: what to discard? scoring: what is best?)
and the pipeline might need some rewiring depending on the parameters of the run (for
instance, the read length).

Preparation
First of all, you need to have an annotation —as good as possible— of the organism you
intend to map to. Typically such annotation will be made available to you as a .gff/.gtf
file.

    2 L F l y B a s e e x o n 1 1 9 2 2 5 6 1 1 9 2 6 2 5 . + . t r a n s c r i p t _ i d F B g n 0 0 3 1 3 2 2
    2 L F l y B a s e e x o n 1 1 9 2 2 6 8 1 1 9 2 6 2 5 . + . t r a n s c r i p t _ i d F B g n 0 0 3 1 3 2 2
    2 L F l y B a s e e x o n 1 1 9 6 2 1 6 1 1 9 6 8 0 7 . + . t r a n s c r i p t _ i d F B g n 0 0 3 1 3 2 2
    2 L F l y B a s e e x o n 1 1 9 8 4 6 5 1 1 9 9 2 0 4 . - . t r a n s c r i p t _ i d F B g n 0 0 3 1 3 2 3
    2 L F l y B a s e e x o n 1 1 9 9 4 1 2 1 1 9 9 9 2 4 . + . t r a n s c r i p t _ i d F B g n 0 2 6 3 0 8 0
    2 L F l y B a s e e x o n 1 1 9 9 4 1 2 1 1 9 9 7 0 9 . + . t r a n s c r i p t _ i d F B g n 0 2 6 3 0 8 0
    2 L F l y B a s e e x o n 1 2 0 0 2 0 5 1 2 0 0 6 6 3 . + . t r a n s c r i p t _ i d F B g n 0 2 6 3 0 8 1
    2 L F l y B a s e e x o n 1 2 0 0 2 0 5 1 2 0 0 5 8 0 . + . t r a n s c r i p t _ i d F B g n 0 2 6 3 0 8 1

The mapper actually does not understand this format, but a much simpler one defining
only splice junctions:

    2 L + 1 0 0 0 6 5 8 5 2 L + 1 0 0 0 6 9 0 8 3
    2 L + 1 0 0 0 7 1 7 8 2 L + 1 0 0 0 7 2 5 2 1 7
    2 L + 1 0 0 0 7 1 7 8 2 L + 1 0 0 0 7 2 5 5 2 2
    2 L + 1 0 0 0 7 3 4 1 2 L + 1 0 0 0 7 4 1 4 1
    2 L + 1 0 0 0 7 5 0 3 2 L + 1 0 0 0 7 5 6 7 4
    2 L + 1 0 0 0 9 5 3 6 2 L + 1 0 0 1 0 3 6 0 2 0
    2 L + 1 0 0 0 9 5 3 6 2 L + 1 0 0 1 0 7 3 7 1 0
    2 L + 1 0 0 1 0 5 5 7 2 L + 1 0 0 1 0 7 3 7 8
    2 L + 1 0 0 1 1 5 5 2 2 L + 1 0 0 1 1 6 0 7 7

In fact, all the needed conversions are performed for you internally by the GEM RNA
mapping pipeline provided with the GEM distribution.
The input can be supplied in several ways (in case of paired-end data, as two separate
files containing one end each, or as a single file containing interleaved ends, i.e. a file
where for each read the two ends are presented one after the other).
If you say gemtools, the pipeline will show you a (very long) list of options that you can
use to modify its behavior. You can also use configuration files to drive it.