#!/usr/bin/env python
"""Production pipelines"""
#!/usr/bin/env python
import os
import logging
import json
import sys
from sys import exit
import subprocess

import gem
import gem.commands
import gem.reports
import gem.gemtools as gt

from gem.pipeline import MappingPipeline, PipelineError
from gem.utils import Command, CommandException


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
        parser.add_argument('-s', '--same', dest="same", action="store_true", default=False, help="File contain the same reads")
        parser.add_argument('-p', '--paired', dest="paired", action="store_true", default=False,
                            help="Content is paired")
        parser.add_argument('-c', '--compress', dest='compress',
                            action="store_true", default=False,
                            help="Write gzip compressed output if a output "
                            "file is specified")

    def run(self, args):
        if len(args.input) < 2:
            logging.error("You have to specify at least 2 files")
        files = args.input
        output = args.output
        if output is None:
            output = sys.stdout
        gem.merge(files[0], files[1:], output,
                  threads=int(args.threads), same_content=args.same,
                  paired=args.paired, compress=args.compress)


class Convert(Command):
    title = "Convert .map to .bam"
    description = """Take a .map file or reads from stdin and converts
    to .bam.
    """

    def register(self, parser):
        parser.add_argument("-i", "--input", dest="input",
                            help="The .map file. Defaults to stdin")
        parser.add_argument("-I", "--index", dest="index",
                            help="The GEM index", required=True)
        parser.add_argument("-o", "--output", dest="output",
                            help="Output .bam file", required=True)
        parser.add_argument("-q", "--quality", dest="quality",
                            help="Quality offset (33,64,ignore)")
        parser.add_argument("-t", "--threads", dest="threads", type=int,
                            default=1,
                            help="Number of threads")
        parser.add_argument("-m", "--memory", dest="sort_memory",
                            default="768M",
                            help="Memory to use per sorting threads. "
                            "Default 768M")
        parser.add_argument("--no-xs", action="store_true", default=False,
                            dest="no_xs",
                            help="Disable computation of XS field")
        parser.add_argument("-p", "--paired", action="store_true", default=False,
                            dest="paired",
                            help="Paired-end reads")
        parser.add_argument("--no-sort", action="store_true", default=False,
                            dest="no_sort",
                            help="Disable sorting")
        parser.add_argument("--no-index", action="store_true", default=False,
                            dest="no_index",
                            help="Disable indexing the bam file")


    def run(self, args):
        quality = gem._prepare_quality_parameter(args.quality)
        raw = False
        if args.input is not None:
            map_file = gem.files.open(args.input, quality=quality)
        else:
            map_file = gem.files.open(sys.stdin, quality=quality)
            raw = True
        cons = gem.extended_splice_consensus
        if args.no_xs:
            cons = None
        sam = gem.gem2sam(map_file,
                            index=args.index,
                            threads=args.threads,
                            quality=args.quality,
                            consensus=cons,
                            raw=raw
                            )
        gem.sam2bam(sam, output=args.output,
                    sorted=not args.no_sort,
                    threads=args.threads,
                    sort_memory=str(args.sort_memory))
        if not args.no_index:
            gem.bamIndex(args.output)


class Stats(Command):
    title = "Create .map stats"
    description = """Calculate stats on a map file"""

    def register(self, parser):
        self.add_options('gt.stats', parser)

    def run(self, args):
        cmd = self.get_command(args)
        p = subprocess.Popen(cmd)
        sys.exit(p.wait())


class GtfCount(Command):
    title = "Create gene counts and gtf statistics"
    description = """This tools can be used to create GTF statistics and
    simple gene read counts. The assumtion here is that the given annotation
    containes a gene model (different transcript_ids belonging to the same
    gene_id, i.e. gencode or ensemble). In addition, the GTF entries should
    contain a gene_type attribute to count different types, for example rRNA.
    """

    def register(self, parser):
        self.add_options('gt.gtfcount', parser)

    def run(self, args):
        cmd = self.get_command(args)
        p = subprocess.Popen(cmd)
        sys.exit(p.wait())


class Filter(Command):
    title = "Filter .map files"
    description = """Filter .map files"""

    def register(self, parser):
        self.add_options('gt.filter', parser)

    def run(self, args):
        """Run gt filter"""
        cmd = self.get_command(args)
        p = subprocess.Popen(cmd)
        sys.exit(p.wait())


class StatsReport(Command):
    title = "Create a stats report from a .json stats file"
    description = """Takes a .json stats file and create a HTML report including
plots of the main statistics.
"""

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Input map file', required=True)
        parser.add_argument('-o', '--output', dest="output", help='The output name. Defaults to the input name + _stats')
        parser.add_argument('-p', '--paired', dest="paired", action="store_true", default=False, help="Paired end reads")
        parser.add_argument('-e', '--extract', dest="extract", action="store_true", default=False, help="Keep the directory next to the .zip archive")

    def run(self, args):
        output = args.output
        if output is None:
            output = os.path.abspath(args.input) + "_stats"
        gem.reports.create_report(args.input, output, paired=args.paired, extract=args.extract)


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

        logging.gemtools.gt("Creating index")
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
        parser.add_argument('-o', '--output', dest="name", help='Optional output prefix. If this is not set, the annotation name will be used', default=None)

    def run(self, args):
        if not args.index.endswith(".gem"):
            raise CommandException("No valid GEM index specified, the file has to end in .gem")
        if not os.path.exists(args.index):
            raise CommandException("GEM index not found")
        if not os.path.exists(args.annotation):
            raise CommandException("Annotation not found")

        name = args.name
        if name is None:
            name = os.path.basename(args.annotation)

        junctions_out = name + ".junctions"
        index_out = name + ".junctions.gem"

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
        pipeline = MappingPipeline()
        pipeline.register_parameter(parser)

    def run(self, args):
        ## parsing command line arguments
        try:
            ## initialize pipeline and check values
            pipeline = MappingPipeline(args=args)
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)

        # check if we want to do a preparation step
        input_dep = []
        if not pipeline.direct_input and (pipeline.input is not None and ((len(pipeline.input) > 1 or len(filter(lambda x: x.endswith(".gz"), pipeline.input)) > 0))):
            input_dep.append(pipeline.prepare_input(name="prepare"))

        # basic pipeline steps
        map_initial = pipeline.map(name="initial", description="Map to index", dependencies=input_dep)
        map_gtf = pipeline.transcripts_annotation(name="annotation-mapping", dependencies=input_dep, description="Map to transcript-index")
        map_denovo = pipeline.transcripts_denovo(name="denovo-mapping", dependencies=input_dep, description="Map to denovo transcript-index")

        # for single end, just merge, otherwise merge and pair
        merged = -1
        if not pipeline.single_end:
            merged = pipeline.merge_and_pair(name="merge_and_pair", dependencies=[map_initial, map_gtf, map_denovo], final=True)
        else:
            merged = pipeline.merge(name="merge", dependencies=[map_initial, map_gtf, map_denovo], final=True)

        # add stats
        if pipeline.stats_create:
            pipeline.create_stats(name="stats", dependencies=[merged], final=True)

        # add the bam step
        if pipeline.bam_create:
            bam = pipeline.bam(name="bam", dependencies=[merged], final=True)
            if pipeline.bam_index:
                pipeline.index_bam(name="index-bam", dependencies=[bam], final=True)

        # add filter step
        if pipeline.filtered_create:
            filtered = pipeline.filtered_map(name="filtered", dependencies=[merged], final=True)
            if pipeline.counts_create and pipeline.annotation is not None:
                pipeline.create_gtfcounts(name="filtered.counts", suffix=".filtered", dependencies=[filtered], final=True)
            if pipeline.stats_create:
                pipeline.create_stats(name="filtered.stats", dependencies=[filtered], suffix=".filtered", final=True)
            if pipeline.bam_create:
                filtered_bam = pipeline.bam(name="filtered.bam", suffix=".filtered", dependencies=[filtered], final=True)
                if pipeline.bam_index:
                    pipeline.index_bam(name="filtered.index-bam", suffix=".filtered", dependencies=[filtered_bam], final=True)

        if pipeline.counts_create and pipeline.annotation is not None:
            pipeline.create_gtfcounts(name="counts", dependencies=[merged], final=True)



        # show parameter and step configuration
        pipeline.log_parameter()

        # run the pipeline
        try:
            pipeline.run()
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)


class JunctionExtraction(Command):
    description = """Run the split mapper to extract junctions"""
    title = "Junction Extraction"

    def register(self, parser):
        pipeline = MappingPipeline()

        general_group = parser.add_argument_group('General')
        ## general pipeline paramters
        general_group.add_argument('-f', '--files', dest="input", nargs="+", metavar="input",
            help='''Single fastq input file or both files for a paired-end run separated by space.
            Note that if you specify only one file, we will look for the file containing the other pairs
            automatically and start a paired-end run. Add the --single-end parameter to disable
            pairing and file search. The file search for the second pair detects pairs
            ending in [_|.|-][0|1|2].[fq|fastq|txt][.gz].''')
        general_group.add_argument('--single-end', dest="single_end", action="store_true", default=None, help="Single end reads")
        general_group.add_argument('-q', '--quality', dest="quality", metavar="quality",
            default=pipeline.quality, help='Quality offset. 33, 64 or "ignore" to disable qualities.')
        general_group.add_argument('-i', '--index', dest="index", metavar="index", help='Path to the .gem genome index')
        general_group.add_argument('-a', '--annotation', dest="gtf", metavar="gtf", help='Path to the .gtf annotation that will be used to keep junctions even if the site does not have enough coverage.')
        general_group.add_argument('--dry', dest="dry", action="store_true", default=None, help="Print and write configuration but do not start the pipeline")
        general_group.add_argument('-t', '--threads', dest="threads", metavar="threads", type=int, help="Number of threads to use. Default %d" % pipeline.threads)

        pipeline.register_junctions(parser)

    def run(self, args):
        ## parsing command line arguments
        try:
            ## initialize pipeline and check values
            pipeline = MappingPipeline(args=args)
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)

        # add annotation to the game
        if args.gtf is not None:
            pipeline.update({"annotation": args.gtf})

        pipeline.extract_junctions("extract", description="Extract Junctions", final=True)

        # show parameter and step configuration
        pipeline.log_parameter()

        # run the pipeline
        try:
            pipeline.run()
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)
