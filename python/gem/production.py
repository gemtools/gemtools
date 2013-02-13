#!/usr/bin/env python
"""Production pipelines"""
#!/usr/bin/env python
import os
import logging
import sys

import gem
import gem.commands

from gem.pipeline import MappingPipeline, PipelineError
from gem.utils import Command, CommandException, Timer


class Merge(Command):
    title = "Merge .map files"
    description = """Merge two .map files. The first file has to
    be the master file that contains all the reads, the second file can
    contain a subset of the reads with the same ID tags and the same order.
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", nargs="+", help='List of files to merge', required=True)
        parser.add_argument('-t', '--threads', dest="threads", default=1, help='Number of threads')
        parser.add_argument('-o', '--output', dest="output", help='Output file name, prints to stdout if nothing is specified')
        parser.add_argument('--same', dest="same", action="store_true", default=False, help="File contain the same reads")

    def run(self, args):
        if len(args.input) < 2:
            logging.error("You have to specify at least 2 files")

        files = []
        for f in args.input:
            files.append(gem.files.open(f))

        gem.merger(files[0], files[1:]).merge(args.output, threads=int(args.threads), same_content=args.same)


class Junctions(Command):
    title = "Extract junctions from GTF"
    description = """Specify an input GTF to extract the junctions
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Input GTF file', required=True)
        parser.add_argument('-o', '--output', dest="output", help='Output file name, prints to stdout if nothing is specified')

    def run(self, args):
        infile = args.input
        junctions = set([x for x in gem.junctions.from_gtf(infile)])
        logging.info("%d Junctions loaded from file" % (len(junctions)))
        gem.junctions.write_junctions(junctions, args.output)


class Index(Command):
    title = "Index genomes"
    description = """This command can be used to index genomes
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Path to a single uncompressed fasta file with the genome', required=True)
        parser.add_argument('-o', '--output', dest="output", help='Output file name (has to end in .gem), defaults to input file name + .gem extension')
        parser.add_argument('-t', '--threads', dest="threads", help='Number of threads', default=2)

    def run(self, args):
        input = args.input
        output = os.path.basename(input)
        if args.output is not None:
            output = args.output
        else:
            output = output[:output.rfind(".")] + ".gem"
        if not output.endswith(".gem"):
            raise CommandException("Output file name has to end in .gem")
        if not os.path.exists(input):
            raise CommandException("Input file not found : %s" % input)
        if input.endswith(".gz"):
            raise CommandException("Compressed input is currently not supported!")

        gem.index(input, output, threads=args.threads)


class TranscriptIndex(Command):
    title = "Create and index transcriptomes"
    description = """This command creates a transcriptome and its index from a gem
    index and a GTF annotation.

    Currently the output name is striclty taken from the name of the annotation
    file given. The command creates a set of files:

        <gtf>.junctions       -- the junction sites of the GTF
        <gtf>.junctions.fa    -- the transcriptome sequences
        <gtf>.junctions.keys  -- the translation table from transcriptome to genome coordinates
        <gtf>.gem             -- the GEM transcriptome index
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--index', dest="index", help='Path to the GEM genome index', required=True)
        parser.add_argument('-a', '--annotation', dest="annotation", help='Path to the GTF annotation', required=True)
        parser.add_argument('-m', '--max-length', dest="maxlength", help='Maximum read length, defaults to 150', default=150)
        parser.add_argument('-t', '--threads', dest="threads", help='Number of threads', default=2)

    def run(self, args):
        if not args.index.endswith(".gem"):
            raise CommandException("No valid GEM index specified, the file has to end in .gem")
        if not os.path.exists(args.index):
            raise CommandException("GEM index not found")
        if not os.path.exists(args.annotation):
            raise CommandException("Annotation not found")

        junctions_out = os.path.basename(args.annotation) + ".junctions"
        index_out = os.path.basename(args.annotation) + ".gem"

        print "Loading Junctions"
        junctions = set(gem.junctions.from_gtf(args.annotation))
        print "%d Junctions loaded from annotation " % (len(junctions))
        gem.junctions.write_junctions(junctions, junctions_out)
        print "Junctions writen to %s " % (junctions_out)

        print "Computing transcriptome..."
        (transcriptome, keys) = gem.compute_transcriptome(args.maxlength, args.index, junctions_out)

        print "Indexing transcriptome"
        gem.index(transcriptome, index_out, threads=args.threads)
        print "Done"


class RnaPipeline(Command):
    description = """The RNASeq pipeline alignes reads against a reference genome as well as
    agains a specified transcriptome. The transcriptome can be generated from an annotation.
    In addition, the pipeline performes a denovo-junction detection to find unknown junctions.

    Input file detection: If you do not specify --single to disable read pairing, we look automatically for
    the second pair file if you only specify one file. For that to work, the second file has to end with
    wither .2 or _2, with the file extension .fastq or .txt (+ .gz for compressed files). For example,

    gemtools rna-pipeline -f myinput_1.fastq.gz ...

    will sarch for a file myinput_2.fastq.gz and use it as the second pair file.
    """
    title = "GEMTools RNASeq Pipeline"

    def register(self, parser):

        general_group = parser.add_argument_group('General')
        ## general pipeline paramters
        general_group.add_argument('-f', '--files', dest="input", nargs="+", metavar="input",
            help='''Single fastq input file or both files for a paired-end run separated by space.
            Note that if you specify only one file, we will look for the pair counter-part
            automatically and start a paired-end run. Add the --single-end parameter to disable
            pairing and file search. The file search for the second pair detects pairs
            ending in [_|.][1|2].[fastq|txt][.gz].''',
            required=True)
        general_group.add_argument('-q', '--quality', dest="quality", metavar="quality",
            default=33, help='Quality offset. 33, 64 or "ignore" to disable qualities. Default 33')
        general_group.add_argument('-i', '--index', dest="index", metavar="index", help='Path to the .gem genome index')
        general_group.add_argument('-a', '--annotation', dest="annotation", metavar="gtf", help='''Path to the GTF annotation. If specified the transcriptome generated from teh annotation is
            used in addition to denovo junctions.''')
        general_group.add_argument('-r', '--transcript-index', dest="transcript_index", help='''GTF Transcriptome index. If not specified and an annotation is given,
            it is assumed to be <gtf>.gem. ''')
        general_group.add_argument('-k', '--transcript-keys', dest="transcript_keys", help='''Transcriptome .keys file. If not specified and an annotation is given,
            it is assumed to be <gtf>.junctions.keys''')

        general_group.add_argument('-n', '--name', dest="name", metavar="name", help="""Name used for the results. If not specified, the name is inferred from
            the input files""")
        general_group.add_argument('-o', '--output-dir', dest="output_dir", metavar="dir", help='Optional output folder. If not specified the current working directory is used.')
        general_group.add_argument('-t', '--threads', dest="threads", default=2, metavar="threads", type=int, help="Number of threads to use")
        general_group.add_argument('--max-read-length', dest="max_read_length", default=150, help='''The maximum read length. This is used to create the denovo
            transcriptom and acts as an upper bound. The default is 150.''')
        general_group.add_argument('-s', '--scoring', dest="scoring_scheme", metavar="scheme", help='''The default scoring scheme to use''')
        general_group.add_argument('--filter', dest="filter", metavar="filter", help='''The final result filter, defaults to 1,2,25''')
        general_group.add_argument('--no-bam', dest="bam_create", action="store_false", default=True, help="Do not create bam file")
        general_group.add_argument('--no-bam-sort', dest="bam_sort", action="store_false", default=True, help="Do not sort bam file")
        general_group.add_argument('--no-bam-index', dest="bam_index", action="store_false", default=True, help="Do not index the bam file")
        general_group.add_argument('--map-quality', dest="bam_mapq", default=0, help="Filter resulting bam for minimum map quality, Default 0 (no filtering)")
        general_group.add_argument('-g', '--no-gzip', dest="compress", action="store_false", default=True, help="Do not compress result mapping")
        general_group.add_argument('--keep-temp', dest="remove_temp", action="store_false", default=True, help="Keep temporary files")
        general_group.add_argument('--single-end', dest="single_end", action="store_true", default=False, help="Single end reads")

        # genome mapping parameter
        mapping_group = parser.add_argument_group('General Mapping Parameter')
        mapping_group.add_argument('-m', '--mismatches', dest="genome_mismatches", metavar="mm",
            default=None, help='Set the allowed mismatch rate as 0 < mm < 1')
        mapping_group.add_argument('--quality-threshold', dest="genome_quality_threshold", metavar="qth",
            default=None, help='The quality threshold, Default 26')
        mapping_group.add_argument('--max-decoded-matches', dest="genome_max_decoded_matches", metavar="mdm",
            default=None, help='Maximum decoded matches. Default 20')
        mapping_group.add_argument('--min-decoded-strata', dest="genome_min_decoded_strata", metavar="mds",
            default=None, help='Minimum decoded strata. Default to 1')
        mapping_group.add_argument('--min-matched-bases', dest="genome_min_matched_bases", metavar="mmb",
            default=None, help='Minimum ratio of bases that must be matched. Default 0.8')
        mapping_group.add_argument('--max-big-indel-length', dest="genome_max_big_indel_length", metavar="mbi",
            default=None, help='Maximum length of a single indel. Default 15')
        mapping_group.add_argument('--max-edit-distance', dest="genome_max_edit_distance", metavar="med",
            default=None, help='Maximum edit distance (ratio) allowed for an alignment. Default 0.2')
        mapping_group.add_argument('--mismatch-alphabet', dest="genome_mismatch_alphabet", metavar="alphabet",
            default=None, help='The mismatch alphabet. Default ACGT')
        mapping_group.add_argument('--strata-after-best', dest="genome_strata_after_best", metavar="strata",
            default=None, help='The number of strata examined after the best one. Default 1')

        # transcript mapping parameter
        transcript_mapping_group = parser.add_argument_group('Transcript Mapping Parameter')
        transcript_mapping_group.add_argument('-tm', '--transcript-mismatches', dest="transcript_mismatches", metavar="mm",
            default=None, help='Set the allowed mismatch rate as 0 < mm < 1')
        transcript_mapping_group.add_argument('--transcript-quality-threshold', dest="transcript_quality_threshold", metavar="qth",
            default=None, help='The quality threshold, Default 26')
        transcript_mapping_group.add_argument('--transcript-max-decoded-matches', dest="transcript_max_decoded_matches", metavar="mdm",
            default=None, help='Maximum decoded matches. Default 20')
        transcript_mapping_group.add_argument('--transcript-min-decoded-strata', dest="transcript_min_decoded_strata", metavar="mds",
            default=None, help='Minimum decoded strata. Default to 1')
        transcript_mapping_group.add_argument('--transcript-min-matched-bases', dest="transcript_min_matched_bases", metavar="mmb",
            default=None, help='Minimum ratio of bases that must be matched. Default 0.8')
        transcript_mapping_group.add_argument('--transcript-max-big-indel-length', dest="transcript_max_big_indel_length", metavar="mbi",
            default=None, help='Maximum length of a single indel. Default 15')
        transcript_mapping_group.add_argument('--transcript-max-edit-distance', dest="transcript_max_edit_distance", metavar="med",
            default=None, help='Maximum edit distance (ratio) allowed for an alignment. Default 0.2')
        transcript_mapping_group.add_argument('--transcript-mismatch-alphabet', dest="transcript_mismatch_alphabet", metavar="alphabet",
            default=None, help='The mismatch alphabet. Default ACGT')
        transcript_mapping_group.add_argument('--transcript-strata-after-best', dest="transcript_strata_after_best", metavar="strata",
            default=None, help='The number of strata examined after the best one. Default 1')

        # junction detection parameter
        junctions_group = parser.add_argument_group('Junction detection Parameter')
        junctions_group.add_argument('-jm', '--junction-mismatches', dest="junction_mismatches", metavar="jmm",
            default=None, help='Set the allowed mismatch rate for junction detection as 0 < mm < 1')
        junctions_group.add_argument('--junction-max-matches', dest="junctions_max_junction_matches", metavar="mm",
            default=None, help='Number of matches allowed for a junction. Default 5')  # todo : fix description
        junctions_group.add_argument('--min-split-length', dest="junctions_min_split_length", metavar="msl",
            default=None, help='Minimum length of a split. Default 4')
        junctions_group.add_argument('--max-split-length', dest="junctions_max_split_length", metavar="msl",
            default=None, help='Maximum length of a split. Default 500000')
        junctions_group.add_argument('--refinement-step', dest="junctions_refinement_step_size", metavar="r",
            default=None, help='Refinement step. Default 2')  # todo: fix description
        junctions_group.add_argument('--min-split-size', dest="junctions_min_split_size", metavar="mss",
            default=None, help='Minimum split size. Default 15')  # todo: fix description and check with min-split-length
        junctions_group.add_argument('--matched-threshold', dest="junctions_matches_threshold", metavar="mt",
            default=None, help='Matches threshold. Default 75')  # todo: fix description
        junctions_group.add_argument('--junction-coverage', dest="junctions_coverage", metavar="jc",
            default=None, help='Maximum allowed junction converage. Defautl 2')  # todo: fix description

        # pair alignment parameter
        # self.pairing_max_edit_distance = 0.30
        # self.pairing_min_matched_bases = 0.80
        # self.pairing_max_extendable_matches = 0
        #self.pairing_max_matches_per_extension = 1  # todo: do we need to expose these ?
        pairing_group = parser.add_argument_group('Pairing Parameter')
        pairing_group.add_argument('--pairing-quality-threshold', dest="pairing_quality_threshold", metavar="pq",
            default=None, help='Quality threshold for pairing. Defaults to general quality threshold')
        pairing_group.add_argument('--pairing-max-decoded-matches', dest="pairing_max_decoded_matches", metavar="pdm",
            default=None, help='Maximum decoded matches. Default 20')
        pairing_group.add_argument('--pairing-min-decoded-strata', dest="pairing_min_decoded_strata", metavar="pds",
            default=None, help='Minimum decoded strata. Default to 1')
        pairing_group.add_argument('--pairing-min-insert-size', dest="pairing_min_insert_size", metavar="is",
            default=None, help='Minimum insert size allowed for pairing. Default 0')
        pairing_group.add_argument('--pairing-max-insert-size', dest="pairing_max_insert_size", metavar="is",
            default=None, help='Maximum insert size allowed for pairing. Default to maximum split length in junctions settings')

    def run(self, args):
        ## parsing command line arguments
        pipeline = MappingPipeline()

        ## update parameter
        pipeline.update(vars(args))

        ## initialize pipeline and check values
        try:
            pipeline.initialize()
        except PipelineError, e:
            sys.stderr.write("\n" + e.message + "\n")
            exit(1)

        # define pipeline steps
        map_initial = pipeline.map(name="initial")
        map_gtf = pipeline.transcripts_annotation(name="annotation_mapping")
        map_denovo = pipeline.transcripts_denovo(name="denovo_mapping")
        merged = pipeline.merge(name="merge", dependencies=[map_initial, map_gtf, map_denovo], final=pipeline.single_end)
        last = merged

        if not pipeline.single_end:
            paired = pipeline.pair(name="pair", dependencies=[merged], final=True)
            last = paired

        if pipeline.bam_create:
            pipeline.bam(name="bam", dependencies=[last], final=True)

        try:
            pipeline.run()
        except PipelineError, e:
            sys.stderr.write("\n" + e.message + "\n")
            exit(1)
