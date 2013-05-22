#!/usr/bin/env python
"""Production pipelines"""
#!/usr/bin/env python
import os
import logging
import json
import sys
from sys import exit

import gem
import gem.commands
import gem.stats
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
        parser.add_argument('--same', dest="same", action="store_true", default=False, help="File contain the same reads")

    def run(self, args):
        if len(args.input) < 2:
            logging.error("You have to specify at least 2 files")

        files = []
        for f in args.input:
            files.append(gem.files.open(f))

        gem.merger(files[0], files[1:]).merge(args.output, threads=int(args.threads), same_content=args.same)


class Stats(Command):
    title = "Create .map stats"
    description = """Calculate stats on a map file"""

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Input map file', required=True)
        parser.add_argument('-t', '--threads', dest="threads", type=int, default=1, help='Number of threads')
        parser.add_argument('-o', '--output', dest="output", help='Output file name, prints to stdout if nothing is specified')
        parser.add_argument('-p', '--paired', dest="paired", action="store_true", default=False, help="Paired end reads")
        parser.add_argument('-b', '--best', dest="best", action="store_true", default=False, help="Calculates stats only for the best maps")
        parser.add_argument('--json', dest="json", default=None, help='Output stats as json to specified file')

    def run(self, args):
        infile = gem.files.open(args.input)
        stats = gt.Stats(args.best, args.paired)
        stats.read(infile, args.threads)

        if args.json is not None:
            with open(args.json, 'w') as f:
                json.dump(stats.__dict__, f)

        of = sys.stdout
        if args.output is not None:
            of = open(args.output, 'w')
        stats.write(of)
        of.close()


class Filter(Command):
    title = "Filter .map files"
    description = """Filter .map files"""

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input",
                            help='Input map file. If no specified, stdin'
                            ' is used to read input',
                            )
        parser.add_argument('-o', '--output', dest="output",
                            help='Output file prefix')
        parser.add_argument('-t', '--threads', dest="threads",
                            type=int, default=1,
                            help='Number of threads')

        filter_group = parser.add_argument_group("Filter")
        filter_group.add_argument('--max-alignments', dest="max_matches",
                                  type=int, default=0,
                                  help="Reduce the maximum number of "
                                  "alignments reported")
        filter_group.add_argument('--min-strata', dest="min_event_distance",
                                  type=int, default=0,
                                  help="Minimum number of strata (default 0)")
        filter_group.add_argument('--max-strata',
                                  dest="max_event_distance", type=int,
                                  default=None,
                                  help="Maximum number of "
                                  "strata (default all)")
        filter_group.add_argument('--min-levenshtein-error',
                                  dest="min_levenshtein_distance", type=int,
                                  default=0,
                                  help="Minimum levenshtein "
                                  "distance (default 0)")
        filter_group.add_argument('--max-levenshtein-error',
                                  dest="max_levenshtein_distance", type=int,
                                  default=None,
                                  help="Maximum levenshtein "
                                  "distance (default all)")
        filter_group.add_argument('--max-insert-size', dest="max_inss",
                                  type=int, default=None,
                                  help="Maximum insert "
                                  "size for paired reads")
        filter_group.add_argument('--min-insert-size', dest="min_inss",
                                  type=int, default=None,
                                  help="Minimum insert size for paired reads")
        filter_group.add_argument('--filter-strand', dest="filter_strand",
                                  action="store_true", default=False,
                                  help="Filter strand and "
                                  "allow only F-R or R-F")
        filter_group.add_argument('--keep-unique', dest="keep_unique",
                                  action="store_true", default=False,
                                  help="Always keep unique"
                                  " mappings as they are")
        filter_group.add_argument('--min-score', dest="min_score",
                                  default=None,
                                  help="Filter by score")
        filter_group.add_argument('--group', dest="include_groups",
                                  nargs="*",
                                  choices=["I", "II", "III", "IV"],
                                  help="Include only the best mappings for "
                                  "groups I,II, and III and all "
                                  "mappings for IV")
        filter_group.add_argument("--annotation", dest="annotation",
                                 default=None,
                                 help="Apply annotation filtering with the "
                                 "given annotation"
                                 )

        rescore_group = parser.add_argument_group("SAM/BAM and rescoring")
        rescore_group.add_argument("--rescore", dest="rescore",
                                   action="store_true", default=False,
                                   help="Rescore the alignments. You have to "
                                   "specify the index and the quality offset "
                                   "to apply rescoring.")
        rescore_group.add_argument("--create-bam", dest="create_bam",
                                   action="store_true",
                                   default=False,
                                   help="Create BAM file. You need to specify "
                                   "the index and the quality offset as well "
                                   "as the output name. This can not be "
                                   "used with stdout output")
        rescore_group.add_argument("--no-sort", dest="no_sort",
                                   action="store_true",
                                   default=False,
                                   help="Do not sort the resulting bam file")
        rescore_group.add_argument("--no-index", dest="no_index",
                                   action="store_true",
                                   default=False,
                                   help="Do not index the resulting bam file")
        rescore_group.add_argument('--index', dest="index",
                                   help='Index to support rescoring')
        rescore_group.add_argument('-q', '--quality', dest="quality",
                                   choices=["33", "64", "ignore"],
                                   help='Quality Offset 33|64|ignore')

    def run(self, args):
        ## check rescore params
        if args.rescore or args.create_bam:
            if args.index is None:
                args.error("You have to specify an index to do rescoring""")
                return False
            if args.quality is None:
                args.error("You have to specify the quality offset")
                return False
            if args.create_bam and args.output is None:
                args.error("You have to specify an output name prefix to "
                           "create a BAM file")
                return False

        # prepare input
        infile = None
        if args.input:
            infile = gem.files.open(args.input)
        else:
            infile = gem.files.open(sys.stdin)

        name = args.output
        threads = int(args.threads)

        outfile = None
        if name is None:
            outfile = gt.OutputFile(sys.stdout)
        else:
            outfile = gt.OutputFile("%s.map" % name)

        params = {
            "max_matches": args.max_matches,
            "min_event_distance": args.min_event_distance,
            "min_levenshtein_distance": args.min_levenshtein_distance,
            "filter_strand": args.filter_strand,
            "keep_unique": args.keep_unique,
        }

        if args.max_event_distance is not None:
            params["max_event_distance"] = args.max_event_distance
        if args.max_levenshtein_distance is not None:
            params["max_levenshtein_distance"] = args.max_levenshtein_distance
        if args.min_inss is not None:
            params["min_inss"] = args.min_inss
        if args.max_inss is not None:
            params["max_inss"] = args.max_inss
        if args.min_score is not None:
            params["min_score"] = args.min_score
        if args.include_groups is not None:
            params["filter_groups"] = True
            if "I" in args.include_groups:
                params["group_1"] = True
            if "II" in args.include_groups:
                params["group_2"] = True
            if "III" in args.include_groups:
                params["group_3"] = True
            if "IV" in args.include_groups:
                params["group_4"] = True
        else:
            params["filter_groups"] = False

        if args.annotation is not None:
            params["annotation"] = args.annotation

        scored = infile
        if args.rescore:
            scored = gem.score(infile, args.index, quality=args.quality,
                               threads=threads)


        gt.filter_map(scored, outfile, params, threads=threads,
                      background_process=False)

        if name is not None:
            outfile.close()


        if args.create_bam:
            map_file = gem.files.open("%s.map" % name, quality=args.quality)
            sam = gem.gem2sam(map_file,
                              index=args.index,
                              threads=threads,
                              quality=args.quality,
                              consensus=gem.extended_splice_consensus,
                              )
            gem.sam2bam(sam, output=("%s.bam" % name),
                        sorted=not args.no_sort,
                        threads=threads, sort_memory="2G")
            if not args.no_index:
                gem.bamIndex("%s.bam" % name)


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
        gem.stats.create_report(args.input, output, paired=args.paired, extract=args.extract)


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
        parser.add_argument('--no-hash', dest="create_hash", default=True, action="store_false", help='Don not create the .hash file')
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

        if args.create_hash:
            logging.gemtools.gt("Creating hash")
            hash_name = output[:-4] + ".hash"
            gem.hash(input, hash_name)


class Hash(Command):
    title = "Hash genomes"
    description = """Create a .hash file out of a fasta genome
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Path to a single uncompressed fasta file with the genome', required=True)
        parser.add_argument('-o', '--output', dest="output", help='Output file name (has to end in .hash), defaults to input file name + .hash extension')

    def run(self, args):
        input = args.input
        output = os.path.basename(input)
        if args.output is not None:
            output = args.output
        else:
            output = output[:output.rfind(".")] + ".hash"

        if not output.endswith(".hash"):
            raise CommandException("Output file name has to end in .hash")

        if not os.path.exists(input):
            raise CommandException("Input file not found : %s" % input)

        if input.endswith(".gz"):
            raise CommandException("Compressed input is currently not supported!")

        gem.hash(input, output)


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


class SamConverter(Command):
    description = """Convert a .map or .map.gz file to SAM/BAM and index/sort it
    """
    title = "Convert to SAM/BAM"

    def register(self, parser):
        pipeline = MappingPipeline()
        pipeline.register_general(parser)
        pipeline.register_bam(parser)
        pipeline.register_execution(parser)

    def run(self, args):
        ## parsing command line arguments
        try:
            ## initialize pipeline and check values
            pipeline = MappingPipeline(args=args)
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)

        # add the bam step
        bam = pipeline.bam(name="bam", dependencies=[], final=True)
        if pipeline.bam_index:
            pipeline.index_bam(name="index-bam", dependencies=[bam], final=True)

        # show parameter and step configuration
        pipeline.log_parameter()

        # run the pipeline
        try:
            pipeline.run()
        except PipelineError, e:
            sys.stderr.write("\nERROR: " + e.message + "\n")
            exit(1)

