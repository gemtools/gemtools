#!/usr/bin/env python
"""Production pipelines"""
import os
import logging
import sys
from sys import exit
import subprocess

from gem.commands import cli
import gem
import gem.commands
import gem.gemtools as gt

from gem.utils import Command, CommandException


def _check(msg, statement=None):
    if statement is None or statement:
        raise CommandException(msg)


@cli("merge",
     inputs=['--input'],
     outputs=['--output'])
class Merge(Command):
    title = "Merge .map files"
    description = """Merge two .map files. The first file has to
    be the master file that contains all the reads, the second file can
    contain a subset of the reads with the same ID tags and the same order.
    """

    def register(self, parser):
        parser.description = Merge.description
        ## required parameters
        parser.add_argument('-i', '--input',
                            nargs="+",
                            help='List of files to merge',
                            required=True)
        parser.add_argument('-o', '--output',
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified')
        parser.add_argument('-t', '--threads',
                            dest="threads",
                            type=int,
                            default=1,
                            help='Number of threads')
        parser.add_argument('-s', '--same',
                            action="store_true",
                            default=False,
                            help="File contain the same reads")
        parser.add_argument('-p', '--paired',
                            action="store_true",
                            default=False,
                            help="Content is paired")
        parser.add_argument('-c', '--compress', dest='compress',
                            action="store_true", default=False,
                            help="Write gzip compressed output if a output "
                            "file is specified")

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        if len(args.input) < 2:
            logging.error("You have to specify at least 2 files")
            exit(1)
        files = args.input
        output = args.output
        if output is None:
            output = sys.stdout
        gem.merge(files[0], files[1:], output,
                  threads=int(args.threads), same_content=args.same,
                  paired=args.paired, compress=args.compress)


@cli("score",
     inputs=['--input'],
     outputs=['--output'])
class Score(Command):
    title = "Score .map files"
    description = """Apply the gemtools default scoring scheme to the
    given mappings
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-I', '--index', dest="index",
                            metavar="<index>",
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality', dest="quality",
                            metavar="<quality>",
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-i', '--input',
                            default=sys.stdin,
                            help='Input .map file. Reads from stdin if '
                            'not specified')
        parser.add_argument('-o', '--output',
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified')
        parser.add_argument('-f', '--filter',
                            help='Apply an additional filtering step where. '
                            'The expectex value is s string with: '
                            '<max_strata>,<max_distance>,<max_alignments>')
        parser.add_argument('-t', '--threads',
                            dest="threads",
                            type=int,
                            default=1,
                            help='Number of threads')
        parser.add_argument('-c', '--compress', dest='compress',
                            action="store_true", default=False,
                            help="Write gzip compressed output if a output "
                            "file is specified")

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        gem.score(args.input,
                  args.index,
                  args.output,
                  quality=args.quality,
                  filter=args.filter,
                  threads=args.threads,
                  compress=args.compress)


@cli("prepare",
     inputs=['--input'],
     outputs=['--output'])
class PrepareInput(Command):
    title = "Prepare input files for a pipeline run"
    description = """Decompresses/Interleaves/Converts input files
    into a single fasta/fastq file that can be passed to the gem-mapper.
    """

    def register(self, parser):
        parser.description = PrepareInput.description
        ## required parameters
        parser.add_argument('-i', '--input', dest="input",
                            nargs="+",
                            default=[],
                            metavar="<input>",
                            help='List of files to prepare', required=True)
        parser.add_argument('-o', '--output', dest="output",
                            metavar="<output>",
                            default=sys.stdout,
                            help='Output file name, prints to stdout "\
                            "if nothing is specified')
        parser.add_argument('-t', '--threads', dest="threads", default=1,
                            metavar="<threads>",
                            type=int,
                            help='Number of threads')
        parser.add_argument('-d', '--discover', dest="discover", default=False,
                            action="store_true",
                            help="Search for a second paird-end reads file")

    def validate(self):
        input = self.options['input'].raw()
        _check("No input specified", input is None or len(input) == 0)
        _check("You can not perpare more than 2 files", len(input) > 2)
        if len(input) == 1 and self.options['discover'].raw():
            # search for second file
            (n, p) = gem.utils.find_pair(input[0])
            if p is not None:
                self.options['input'].append(p)
        return True

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        threads = args.threads
        input = args.input

        output = args.output if not isinstance(args.output, file) \
            else args.output
        outfile = gt.OutputFile(output, clean_id=True, append_extra=False)

        if len(input) == 1:
            infile = gem.files.open(input[0])
        else:
            infile = gem.filter.interleave(
                [gem.files.open(f) for f in input],
                threads=max(1, threads))
        infile.write_stream(outfile, write_map=False)
        infile.close()
        outfile.close()


@cli("mapper",
     inputs=['--input'],
     outputs=['--output'])
class GemMapper(Command):
    title = "Run the GEM-mapper"
    description = """Runs the GEM-mapper"""

    def register(self, parser):
        parser.description = GemMapper.description
        ## required parameters
        parser.add_argument('-I', '--index', dest="index",
                            metavar="<index>",
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality', dest="quality",
                            metavar="<quality>",
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-i', '--input', dest="input",
                            metavar="<input>",
                            default=sys.stdin,
                            help='Input File. Reads from stdin if nothing '
                            'is specified.')
        parser.add_argument('-o', '--output', dest="output",
                            metavar="<output>",
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified')
        parser.add_argument('-t', '--threads', dest="threads",
                            metavar="<threads>",
                            default=1,
                            type=int,
                            help='Number of threads')
        parser.add_argument('-s', '--strata-after-best',
                            dest="strata_after_best", type=int,
                            metavar="strata",
                            default=0,
                            help='The number of strata examined after '
                            'the best one. Default 0')
        parser.add_argument('-m', '--mismatches',
                            dest="mismatches", metavar="mm",
                            type=float,
                            default=0.04,
                            help='Set the allowed mismatch ratio as '
                            '0 < mm < 1. Default 0.04')
        parser.add_argument('--quality-threshold', dest="quality_threshold",
                            type=int, metavar="qth",
                            default=26,
                            help='Good quality threshold. Bases with a '
                            'quality score >= threshold are considered '
                            'good. Default 26')
        parser.add_argument('--max-decoded-matches',
                            dest="max_decoded_matches",
                            type=int,
                            default=20,
                            metavar="mdm",
                            help='Maximum decoded matches. Default 20')
        parser.add_argument('--min-decoded-strata',
                            dest="min_decoded_strata",
                            type=int,
                            default=1,
                            metavar="mds",
                            help='Minimum decoded strata. '
                            'Default to 1')
        parser.add_argument('--min-big-indel-length',
                            dest="max_big_indel_length",
                            type=int,
                            default=15,
                            metavar="mil",
                            help='Maximum big indel length. '
                            'Default to 15')
        parser.add_argument('--min-matched-bases',
                            dest="min_matched_bases",
                            type=float,
                            default=0.80,
                            metavar="mmb",
                            help='Minimum ratio of bases that '
                            'must be matched. Default 0.80')
        parser.add_argument('--max-edit-distance',
                            dest="max_edit_distance",
                            metavar="med",
                            default=0.20,
                            type=float,
                            help='Maximum edit distance (ratio) '
                            'allowed for an alignment. Default 0.20')
        parser.add_argument('--keys', '-k',
                            dest="keys",
                            metavar="keys",
                            help='Key file to translate result back to '
                            'genomic coordinates')
        parser.add_argument('--mismatch-alphabet',
                            dest="mismatch_alphabet", metavar="alphabet",
                            default="ACGT",
                            help='The mismatch alphabet. Default "ACGT"')
        parser.add_argument('-C', '--compress', dest="compress",
                            default=False,
                            action="store_true",
                            help="Compress the output file")

    def validate(self):
        _check("No index specified", not self.options['index'].get())
        return True

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        outfile = args.output if not isinstance(args.output, file) \
            else args.output
        infile = gem.files.open(args.input) \
            if not isinstance(args.input, file) else args.input

        gem.mapper(
            infile,
            args.index,
            outfile,
            quality=args.quality,
            mismatches=args.mismatches,
            quality_threshold=args.quality_threshold,
            max_decoded_matches=args.max_decoded_matches,
            min_decoded_strata=args.min_decoded_strata,
            min_matched_bases=args.min_matched_bases,
            max_big_indel_length=args.max_big_indel_length,
            max_edit_distance=args.max_edit_distance,
            mismatch_alphabet=args.mismatch_alphabet,
            delta=args.strata_after_best,
            key_file=args.keys,
            #trim=args.trim,
            threads=args.threads,
            compress=args.compress
        )


@cli("pairalign",
     inputs=['--input'],
     outputs=['--output'])
class GemPairalign(Command):
    title = "Run the GEM-Mapper to find pairs"
    description = """Runs the GEM-mapper in pair align mode"""

    def register(self, parser):
        parser.description = GemPairalign.description
        ## required parameters
        parser.add_argument('-I', '--index', dest="index",
                            metavar="<index>",
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality', dest="quality",
                            metavar="<quality>",
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-i', '--input', dest="input",
                            metavar="<input>",
                            default=sys.stdin,
                            help='Input File. Reads from stdin if nothing '
                            'is specified.')
        parser.add_argument('-o', '--output', dest="output",
                            metavar="<output>",
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified')
        parser.add_argument('-t', '--threads', dest="threads",
                            metavar="<threads>",
                            default=1,
                            type=int,
                            help='Number of threads')
        parser.add_argument('-s', '--strata-after-best',
                            dest="strata_after_best", type=int,
                            metavar="strata",
                            default=0,
                            help='The number of strata examined after '
                            'the best one. Default 0')
        parser.add_argument('--min-insert-size',
                            type=int,
                            default=0,
                            help='The minimum insert size for a pair')
        parser.add_argument('--max-insert-size',
                            type=int,
                            default=1000,
                            help='The maximum insert size for a pair')
        parser.add_argument('--quality-threshold', dest="quality_threshold",
                            type=int, metavar="qth",
                            default=26,
                            help='Good quality threshold. Bases with a '
                            'quality score >= threshold are considered '
                            'good. Default 26')
        parser.add_argument('--max-decoded-matches',
                            dest="max_decoded_matches",
                            type=int,
                            default=20,
                            metavar="mdm",
                            help='Maximum decoded matches. Default 20')
        parser.add_argument('--min-decoded-strata',
                            dest="min_decoded_strata",
                            type=int,
                            default=1,
                            metavar="mds",
                            help='Minimum decoded strata. '
                            'Default to 1')
        parser.add_argument('--min-matched-bases',
                            dest="min_matched_bases",
                            type=float,
                            default=0.80,
                            metavar="mmb",
                            help='Minimum ratio of bases that '
                            'must be matched. Default 0.80')
        parser.add_argument('--max-edit-distance',
                            dest="max_edit_distance",
                            metavar="med",
                            default=0.30,
                            type=float,
                            help='Maximum edit distance (ratio) '
                            'allowed for an alignment. Default 0.20')
        parser.add_argument('--max-extendable-matches',
                            default=0,
                            type=int,
                            help='Maximal extendable matches. 0 disables '
                            'extension mode.')
        parser.add_argument('--max-matches-per-extension',
                            default=0,
                            type=int,
                            help='Maximum matches per extension')
        parser.add_argument('-C', '--compress', dest="compress",
                            default=False,
                            action="store_true",
                            help="Compress the output file")

    def validate(self):
        _check("No index specified", not self.options['index'].get())
        return True

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        outfile = args.output if not isinstance(args.output, file) \
            else args.output
        infile = gem.files.open(args.input) \
            if not isinstance(args.input, file) else args.input

        gem.pairalign(
            infile,
            args.index,
            outfile,
            quality=args.quality,
            quality_threshold=args.quality_threshold,
            max_decoded_matches=args.max_decoded_matches,
            min_decoded_strata=args.min_decoded_strata,
            min_insert_size=args.min_insert_size,
            max_insert_size=args.max_insert_size,
            min_matched_bases=args.min_matched_bases,
            max_edit_distance=args.max_edit_distance,
            max_extendable_matches=args.max_extendable_matches,
            max_matches_per_extension=args.max_matches_per_extension,
            threads=args.threads,
            compress=args.compress
        )


@cli("convert",
     inputs=['--input'],
     outputs=['--output'])
class Convert(Command):
    title = "Convert .map to .bam"
    description = """Take a .map file or reads from stdin and converts
    to .bam.
    """

    def register(self, parser):
        parser.description = Convert.description
        parser.add_argument("-i", "--input", dest="input",
                            default=sys.stdin,
                            help="The .map file. Defaults to stdin")
        parser.add_argument("-o", "--output", dest="output",
                            help="Output .bam file", required=True)
        parser.add_argument("-I", "--index", dest="index",
                            help="The GEM index")
        parser.add_argument("-q", "--quality", dest="quality",
                            help="Quality offset (33,64,ignore)")
        parser.add_argument("-t", "--threads", dest="threads", type=int,
                            default=1,
                            help="Number of threads")
        parser.add_argument("-m", "--memory",
                            default="768M",
                            help="Memory to use per sorting threads. "
                            "Default 768M")
        parser.add_argument("--no-xs", action="store_true", default=False,
                            dest="no_xs",
                            help="Disable computation of XS field")
        parser.add_argument("-p", "--paired", action="store_true",
                            default=False,
                            dest="paired",
                            help="Paired-end reads")
        parser.add_argument("--no-sort", action="store_true", default=False,
                            dest="no_sort",
                            help="Disable sorting")
        parser.add_argument("--no-index", action="store_true", default=False,
                            dest="no_index",
                            help="Disable indexing the bam file")

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        quality = gem._prepare_quality_parameter(args.quality)
        map_file = gem.files.open(args.input, quality=quality)
        raw = isinstance(args.input, file)
        cons = gem.extended_splice_consensus
        if args.no_xs:
            cons = None
        sam = gem.gem2sam(map_file,
                          index=args.index,
                          threads=args.threads,
                          quality=args.quality,
                          consensus=cons,
                          raw=raw)
        gem.sam2bam(sam, output=args.output,
                    sorted=not args.no_sort,
                    threads=args.threads,
                    sort_memory=str(args.memory))
        if not args.no_index:
            gem.bamIndex(args.output)


@cli("stats",
     inputs=['--input'],
     outputs=['--output'])
class Stats(Command):
    title = "Create .map stats"
    description = """Calculate stats on a map file"""

    def register(self, parser):
        parser.description = Stats.description
        self.add_options('gt.stats', parser, stream_in="input",
                         stream_out="output")
        parser.add_argument("--json", help="Write json stats to a file")

    def run(self, args):
        stderr = None
        if args["json"]:
            args["output_format"] = "both"
            stderr = open(args["json"], 'wb')
        cmd = self.get_command(args)
        p = subprocess.Popen(cmd, stderr=stderr)
        p.wait()
        if stderr is not None:
            stderr.close()
        sys.exit(p.wait())


@cli("gtf-stats",
     inputs=['--input'],
     outputs=['--output'])
class GtfCount(Command):
    title = "Create gene counts and gtf statistics"
    description = """This tools can be used to create GTF statistics and
    simple gene read counts. The assumtion here is that the given annotation
    containes a gene model (different transcript_ids belonging to the same
    gene_id, i.e. gencode or ensemble). In addition, the GTF entries should
    contain a gene_type attribute to count different types, for example rRNA.
    """

    def register(self, parser):
        self.add_options('gt.gtfcount', parser, stream_in="input",
                         stream_out="output")
        parser.add_argument("--json", help="Write json stats to a file")

    def run(self, args):
        stderr = None
        if args["json"]:
            args["output_format"] = "both"
            stderr = open(args["json"], 'wb')
        cmd = self.get_command(args)
        p = subprocess.Popen(cmd, stderr=stderr)
        p.wait()
        if stderr is not None:
            stderr.close()


@cli("filter",
     inputs=['--input'],
     outputs=['--output'])
class Filter(Command):
    title = "Filter .map files"
    description = """Filter .map files"""

    def register(self, parser):
        parser.description = Filter.description
        self.add_options('gt.filter', parser, stream_in="input",
                         stream_out="output")
        parser.add_argument("--compress", action="store_true",
                            default=False,
                            help="Compress output. Only works if "
                            "an output file is specified.")

    def run(self, args):
        """Run gt filter"""
        outstream = sys.stdout
        threads = args.get('threads', 1)
        if not threads:
            threads = 1
        compressor = None
        if not isinstance(args['output'], file) and args['compress']:
            outstream, out, compressor = gem.filter.create_output_stream(
                args['output'],
                compress=True,
                threads=threads
            )
            args['output'] = None

        cmd = self.get_command(args)
        p = subprocess.Popen(cmd, stdout=outstream)
        r = p.wait()
        if compressor:
            compressor.stdin.close()
            compressor.wait()
        sys.exit(r)


@cli("gtf-junctions",
     title="Extract junctions from GTF",
     description="Specify an input GTF to extract the junctions")
def gtf_junctions(args):
    """\
    Extracts junctions sites from a given GTF file. If no file is specified,
    stdin is used as input.

    Usage:
        gemtools gtf-junctions [-h] [-i <input>] [-o <output>]

    Inputs:
        -i, --input <input>    The input gtf file
                               [default: stdin]
    Outputs:
        -o, --output <output>  The output junctions file
                               [default: stdout]
    Options:
        -h, --help             Show this help message
    """
    infile = args["input"]
    junctions = set([x for x in gem.junctions.from_gtf(infile)])
    logging.info("%d Junctions loaded from file" % (len(junctions)))
    gem.junctions.write_junctions(junctions, args["output"])


@cli("index",
     inputs=['--input'],
     outputs=['--output'])
class Index(Command):
    title = "Index genomes"
    description = """This command can be used to index genomes
    """

    def register(self, parser):
        parser.description = Index.description
        ## required parameters
        parser.add_argument('-i', '--input',
                            help='Path to a single uncompressed '
                            'fasta file with the genome',
                            required=True)
        parser.add_argument('-o', '--output',
                            help='Output file name (has to end in .gem), '
                            'defaults to input file name + .gem extension')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=2,
                            help='Number of threads')

    def validate(self):
        input = self.options['input'].get()
        output = os.path.basename(input)
        if self.options['output'].raw() is not None:
            output = self.options['output'].get()
        else:
            output = output[:output.rfind(".")] + ".gem"
        self.options['output'] = output

        if not output.endswith(".gem"):
            raise CommandException("Output file name has to end in .gem")
        if not self.options['input'].is_dependency() and \
                not os.path.exists(input):
            raise CommandException("Input file not found : %s" % input)
        if input.endswith(".gz"):
            raise CommandException("Compressed input is currently "
                                   "not supported!")

    def run(self, args):
        gem.index(args['input'], args['output'], threads=args['threads'])
        logfile = args['output'][:-4] + ".log"
        if os.path.exists(logfile):
            os.remove(logfile)


@cli("t-index",
     add_outputs=[
         "index_out",
     ])
class TranscriptIndex(Command):
    title = "Create and index transcriptomes"
    description = """This command creates a transcriptome and its index
    from a gem index and a GTF annotation or a junctions file.

    The output name is created from the name of the annotation
    file given if its not explicitly set. The command creates a set of files:


        <name>.junctions.fa    -- the transcriptome sequences
        <name>.junctions.keys  -- the translation table from transcriptome
                                  to genome coordinates
        <name>.junctions.gem   -- the GEM transcriptome index

    If a GTF annotation is specified, the junctions found are reported as

        <name>.junctions       -- the junction sites of the GTF
    """

    def register(self, parser):
        parser.description = TranscriptIndex.description
        ## required parameters
        parser.add_argument('-i', '--index',
                            help='Path to the GEM genome index',
                            required=True)
        parser.add_argument('-a', '--annotation',
                            dest="annotation",
                            help='Path to the GTF annotation')
        parser.add_argument('-j', '--junctions',
                            dest="junctions",
                            help='Path to a .junctions file')
        parser.add_argument('-n', '--name',
                            help='Optional output prefix. '
                            'If this is not set, the annotation/junctions '
                            'file name will be used',
                            default=None)
        parser.add_argument('-m', '--max-length',
                            type=int,
                            help='Maximum read length, defaults to 150',
                            default=150)
        parser.add_argument('-t', '--threads',
                            type=int,
                            help='Number of threads',
                            default=2)

    def validate(self):
        trans = ComputeTranscriptome()
        trans.validate(oargs)
        args = gem.utils.dict_2_tuple(oargs)
        name = args.name
        index_out = name + ".junctions.gem"
        oargs['index_out'] = index_out
        return True

    def run(self, args):
        trans = ComputeTranscriptome()
        trans.run(args)
        print "Indexing transcriptome"
        gem.index(args['fasta_out'], args['index_out'],
                  threads=args['threads'])
        logfile = args['index_out'][:-4] + ".log"
        if os.path.exists(logfile):
            os.remove(logfile)
        print "Done"


@cli("compute-transcriptome",
     add_outputs=[
         "fasta-out",
         "keys-out",
         "junctions-out"
     ])
class ComputeTranscriptome(Command):
    title = "Create a transcriptome from a genome index"
    description = """This command creates a transcriptome from a gem index
    and a GTF annotation or a junctions file.

    The output name is created from the name of the annotation/junctions
    file given if its not explicitly set. The command creates a set of files:

        <name>.junctions.fa    -- the transcriptome sequences
        <name>.junctions.keys  -- the translation table from transcriptome
                                  to genome coordinates

    If a GTF annotation is specified, the junctions found are reported as

        <name>.junctions       -- the junction sites of the GTF
    """

    def register(self, parser):
        parser.description = ComputeTranscriptome.description
        ## required parameters
        parser.add_argument('-i', '--index',
                            help='Path to the GEM genome index',
                            required=True)
        parser.add_argument('-a', '--annotation',
                            dest="annotation",
                            help='Path to the GTF annotation')
        parser.add_argument('-j', '--junctions',
                            dest="junctions",
                            help='Path to a .junctions file')
        parser.add_argument('-n', '--name',
                            help='Optional output prefix. '
                            'If this is not set, the annotation/junctions '
                            'file name will be used',
                            default=None)
        parser.add_argument('-m', '--max-length',
                            type=int,
                            help='Maximum read length, defaults to 150',
                            default=150)

    def validate(self):
        args = self.options.to_dict()
        if not args['index'].endswith(".gem"):
            raise CommandException("No valid GEM index specified, "
                                   "the file has to end in .gem")
        self.options['index'].check_files()
        self.options['annotation'].check_files()
        self.options['junctions'].check_files()
        if not args['junctions'] and not args['annotation']:
            raise CommandException("You have to specify either the "
                                   "junctions file or a reference annotation")
        name = args['name']
        if name is None:
            fname = args['annotation']
            if not fname:
                fname = args['junctions']
            name = os.path.basename(fname)
            i = name.index('.')
            name = name[:i] if i > 0 else name

        junctions_out = name + ".junctions"
        keys_out = name + ".junctions.keys"
        fasta_out = name + ".junctions.fa"

        self.options['name'] = name
        self.options.add_output('junctions_out',
                                junctions_out if args['annotation']
                                else args['junctions'])
        self.options.add_output('keys_out', keys_out)
        self.options.add_output('fasta_out', fasta_out)
        return True

    def run(self, args):
        print "Loading Junctions"
        if args['junctions']:
            junctions = set(gem.junctions.from_junctions(args['junctions']))
        else:
            junctions = set(gem.junctions.from_gtf(args['annotation']))
        print "%d Junctions loaded" % (len(junctions))
        if args['annotation']:
            gem.junctions.write_junctions(junctions, args['junctions_out'])
            print "Junctions writen to %s " % (args['junctions_out'])

        print "Computing transcriptome..."
        (transcriptome, keys) = gem.compute_transcriptome(
            args['max_length'],
            args['index'],
            args['junctions_out'],
            output_name=args['name'] + ".junctions")
        print "Done"


@cli("rna-pipeline",
     pipeline='pipeline',
     add_outputs=[
     ])
class RnaPipeline(Command):
    description = """The RNASeq pipeline alignes reads against a reference
    genome as well as agains a specified transcriptome. The transcriptome can
    be generated from an annotation.
    In addition, the pipeline performes a denovo-junction detection
    to find unknown junctions.

    Input file detection: If you do not specify --single to disable read
    pairing, we look automatically for the second pair file if you only
    specify one file. For that to work, the second file has to end with
    either .2 or _2, with the file extension .fastq or .txt (+ .gz for
    compressed files). For example,

    gemtools rna-pipeline -f myinput_1.fastq.gz ...

    will search for a file myinput_2.fastq.gz and use it as the second
    pair file.
    """
    title = "GEMTools RNASeq Pipeline"

    def register(self, parser):
        parser.description = RnaPipeline.description
        input_group = parser.add_argument_group('Input and names')
        ## general pipeline paramters
        input_group.add_argument(
            '-f', '--files', nargs="+",
            default=[],
            help='''Single fastq input file or both files for a paired-end run
            separated by space. Note that if you specify only one file, we
            will look for the file containing the other pairs automatically
            and start a paired-end run. Add the --single-end parameter to
            disable pairing and file search. The file search for the second
            pair detects pairs ending in [_|.|-][0|1|2].[fq|fastq|txt][.gz].
            ''')
        input_group.add_argument('-q', '--quality',
                                 default="33",
                                 help='Quality offset. 33, 64 or '
                                 '"ignore" to disable qualities.',
                                 required=True)
        input_group.add_argument('-i', '--index',
                                 help='Path to the .gem genome index',
                                 required=True)
        input_group.add_argument('-a', '--annotation',
                                 help='Path to the .gtf annotation file '
                                 'used to create the transcriptome and '
                                 'denovo transcriptome index. This '
                                 'is optionan and if not specified, only '
                                 'denovo-transcript mapping will be performed')
        input_group.add_argument('-r', '--transcript-index',
                                 help='Path to the .gem transcript index '
                                 'that was generated from the given '
                                 'annotation. ')
        input_group.add_argument('-k', '--transcript-keys',
                                 help='Path to the .keys file that '
                                 'was generated from the given transcript '
                                 'index. ')
        input_group.add_argument('-n', '--name',
                                 help="Specify a prefix for all result files. "
                                 "If no name is specified, the name of the "
                                 "primary input file will be used (without "
                                 "any file extensions")
        input_group.add_argument('-t', '--threads',
                                 default=1,
                                 type=int,
                                 help="Number of threads")
        input_group.add_argument('--single-end',
                                 action="store_true",
                                 default=False,
                                 help="Single end reads")
        input_group.add_argument('--direct-input',
                                 default=False,
                                 action="store_true",
                                 help="Skip preparation step and pipe the "
                                 "input directly into the first mapping step")
        input_group.add_argument('--compress-all',
                                 default=False,
                                 action="store_true",
                                 help="Compress also intermediate files")
        ######################################################################
        # Job control
        #####################################################################
        ctrl_group = parser.add_argument_group('Job and Pipeline controls')
        ctrl_group.add_argument("--dry", action="store_true", default=False,
                                help="Show the pipeline configuraiton but "
                                "do not execute the pipeline")

    def _find_second_pair(self, args):
        ## try to guess the second file
        (n, p) = gem.utils.find_pair(args['files'][0])
        if p is not None:
            args['files'].append(p)
        # we guessed a name
        if args['name'] is None:
            args['name'] = n

    def _guess_name(self, args):
        name = os.path.basename(args['files'][0])
        if name.endswith(".gz"):
            name = name[:-3]
        idx = name.rfind(".")
        if idx > 0:
            name = name[:idx]
        args['name'] = name

    def _file_name(self, suffix=None, args=None, compress=None):
        if compress is None:
            compress = args["compress_all"]
        return "%s%s%s" % (args['name'], "" if suffix is None else suffix,
                           ".gz" if compress else "")

    def _find_transcript_index(self, args):
        from os.path import basename, dirname, join, exists
        from jip.utils import rreplace
        # try to detect the transcript index based on the annotation
        # file name. We remove the .gtf[.gz] extension and check for
        # a <name>.gem and <name>.junctions.gem
        name = basename(args['annotation'])
        base = dirname(args['annotation'])
        i = name.rindex(".")
        if i > 0:
            name = name[:i]

        if not args['transcript_index']:
            t_index = join(base, "%s.gem" % (name))
            if not exists(t_index):
                t_index = join(base, "%s.junctions.gem" % (name))
            args['transcript_index'] = t_index

        if not args['transcript_keys']:
            # we look for the key file just next to the transcript_index
            # replacing the .gem extension with .keys. If that fails,
            # we try .junctions.keys
            t_keys = rreplace(args['transcript_index'], ".gem", ".keys", 1)
            if not exists(t_keys):
                t_keys = rreplace(args['transcript_index'], ".gem",
                                  ".junctions.keys", 1)
            args['transcript_keys'] = t_keys

    def validate(self):
        from os.path import exists
        args = self.options.to_dict()
        ###########################################################
        # General input data checks
        ###########################################################
        # check for at least one input file
        _check("No input files specified!", len(args['files']) == 0)
        # check input files for single end alignments
        _check("Single end runs take only one input file!",
               args['single_end'] and len(args['files']) > 1)
        if not args['single_end']:
            # check paired end
            if len(args['files']) == 1:
                self._find_second_pair(args)
            _check("Paired end runs take no more than two "
                   "input files!", len(args['files']) > 2)
        ## check for a name or guess it
        if not args['name']:
            self._guess_name(args)
        ## check that all input files exists
        self.options['files'].check_files()

        ## check the .gem index
        _check("Genome GEM index not found: %s" % args['index'],
               not exists(args['index']))
        _check("Genome GEM does not end in .gem: %s" % args['index'],
               not args['index'].endswith(".gem"))

        ## check quality value
        _check("Quality offset value is not valid! Supported are 33|64|ignore",
               not args['quality'] in ("33", "64", "ignore"))

        ###########################################################
        # Transcriptome mapping with annotation
        ###########################################################
        if args['annotation']:
            _check("GTF annotation not found: %s" % args['annotation'],
                   not exists(args['annotation']))
            # check the transcript index and keys
            self._find_transcript_index(args)
            _check("No Transcript index found: %s" %
                   args['transcript_index'],
                   not exists(args['transcript_index']))
            _check("No Transcript keys found: %s" %
                   args['transcript_keys'],
                   not exists(args['transcript_keys']))

        ###########################################################
        # Output file configuration
        ###########################################################
        def fn(name, c=None):
            return self._file_name(name, args, compress=c)

        opts = self.options
        opts.add_output('prepare_out', fn("_prepare.fq", c=False))
        opts.add_output('initial_map_out', fn("_initial.map"))
        opts.add_output('transcript_map_out', fn("_transcripts.map"))
        opts.add_output('denovo_junctions_out', fn("_denovo.junctions",
                                                   c=False))
        opts.add_output('denovo_map_out', fn("_denovo.map"))
        opts.add_output('final_out', fn(".map", c=True))
        opts.add_output('bam_out', fn(".bam", c=False))
        opts.add_output('filtered_out', fn(".filtered.map", c=True))
        opts.add_output('filtered_bam_out', fn(".filtered.bam", c=False))
        opts.add_output('stats_out', fn(".stats.txt", c=False))
        opts.add_output('stats_json_out', fn(".stats.json", c=False))
        opts.add_output('filtered_stats_out',
                        fn(".filtered.stats.txt", c=False))
        opts.add_output('filtered_stats_json_out',
                        fn(".filtered.stats.json", c=False))
        opts.add_output('gtf_stats_out', fn(".gtf.stats.txt", c=False))
        opts.add_output('gtf_stats_json_out', fn(".gtf.stats.json", c=False))
        opts.add_output('gtf_counts_out', fn(".gtf.counts.txt", c=False))
        opts.add_output('filtered_gtf_stats_out',
                        fn(".filtered.gtf.stats.txt", c=False))
        opts.add_output('filtered_gtf_stats_json_out',
                        fn(".filtered.gtf.stats.json", c=False))
        opts.add_output('filtered_gtf_counts_out',
                        fn(".filtered.gtf.counts.txt", c=False))
        return True

    def pipeline(self):
        from jip import Pipeline
        args = self.options.to_dict()
        threads = int(args['threads'])
        quality = args['quality']
        index = args['index']

        p = Pipeline()
        job = p.job(threads=threads)
        mapping_inputs = args['files']
        prepare_step = None
        if not args['direct_input']:
            ## run the prepare step
            prepare_step = job('Prepare').run('gemtools_prepare',
                                              input=args['files'],
                                              output=args['prepare_out'],
                                              threads=threads)
            mapping_inputs = prepare_step

        # we collect all mapping steps here for mergin
        all_mappings = []

        ###################################################################
        # Initial mapping step
        ###################################################################
        mapper_cfg = dict(input=mapping_inputs,
                          threads=threads,
                          index=index,
                          output=args['initial_map_out'],
                          compress=args['compress_all'],
                          quality=quality)
        initial_mapping = job('Initial.Mapping').run('gemtools_mapper',
                                                     **mapper_cfg)
        all_mappings.append(initial_mapping)

        ###################################################################
        # Transcriptome mapping
        ###################################################################
        if args['annotation']:
            t_map_cfg = dict(input=mapping_inputs,
                             threads=threads,
                             index=args['transcript_index'],
                             keys=args['transcript_keys'],
                             quality=quality)
            t_map_filter_cfg = dict(only_split_maps=True,
                                    threads=threads,
                                    compress=args['compress_all'],
                                    output=args['transcript_map_out'])
            t_mapping = job('GTF.Mapping').run('gemtools_mapper', **t_map_cfg)
            # filter for splitmaps
            transcript_mapping = t_mapping | \
                job('GTF.Filter').run('gemtools_filter',
                                      **t_map_filter_cfg)
            all_mappings.append(transcript_mapping)

        ###################################################################
        # Denovo transcriptome mapping
        ###################################################################
        junctions = job('Denovo.Junctions').run(
            'gemtools_denovo_junctions',
            index=args['index'],
            threads=threads,
            input=mapping_inputs,
            quality=quality,
            output=args["denovo_junctions_out"]
        )
        denovo_transcriptome = job('Denovo.Transcriptome').run(
            'gemtools_compute_transcriptome',
            index=args['index'],
            junctions=junctions.output,
            max_length=150)  # todo add auto calc

        junctions_index = job('Denovo.Index').run(
            'gemtools_index',
            threads=threads,
            input=denovo_transcriptome.fasta_out)
        mapper_cfg = dict(input=mapping_inputs,
                          index=junctions_index,
                          threads=threads,
                          keys=denovo_transcriptome.keys_out,
                          output=args['denovo_map_out'],
                          compress=args['compress_all'],
                          quality=quality)
        d_map_filter_cfg = dict(only_split_maps=True,
                                threads=threads,
                                compress=args['compress_all'],
                                output=args['transcript_map_out'])
        denovo_mapping = job('Denovo.Mapping').run(
            'gemtools_mapper', **mapper_cfg) | job('Denovo.Filter').run(
                'gemtools_filter', **d_map_filter_cfg)

        all_mappings.append(denovo_mapping)

        merge = job('Merge').run(
            'gemtools_merge',
            same=True,
            threads=threads,
            input=all_mappings)
        pair = job('Pair').run(
            'gemtools_pairalign',
            quality=quality,
            threads=threads,
            index=args['index'])
        score = job('Score').run(
            'gemtools_score',
            index=args['index'],
            threads=threads,
            quality=quality,
            compress=True,
            output=args['final_out'])
        merge | pair | score

        # create sorted bam
        bam_cfg = dict(
            threads=threads,
            index=args['index'],
            quality=quality,
            memory="768M",
            paired=True,
            no_sort=False,
            no_xs=False,
            no_index=False
        )
        job('Covert.BAM').run(
            'gemtools_convert',
            input=score,
            output=args['bam_out'],
            **bam_cfg)
        # create filtered output
        filtered = job('Filter').run(
            'gemtools_filter',
            input=score,
            threads=threads,
            output=args['filtered_out'],
            annotation=args['annotation'],
            no_penalty_for_splitmaps=True,
            paired_end=True,  # paired
            keep_unique=True,
            min_intron_length=20,
            min_block_length=5,
            reduce_by_gene_id=False,
            reduce_by_junctions=False,
            reduce_to_unique_strata=0,
            reduce_to_max_maps=5,
            max_strata=0,  # max error events
            compress=True)
        job('Filter.Convert.BAM').run(
            'gemtools_convert',
            input=filtered,
            output=args['filtered_bam_out'],
            **bam_cfg)

        # create stats
        job('Stats').run(
            'gemtools_stats',
            input=score,
            threads=threads,
            output=args['stats_out'],
            json=args['stats_json_out'],
            paired_end=True,
            all_tests=True)

        job('Filter.Stats').run(
            'gemtools_stats',
            input=filtered,
            threads=threads,
            output=args['filtered_stats_out'],
            json=args['filtered_stats_json_out'],
            paired_end=True,
            all_tests=True)

        # create counts
        job('GTF.Stats').run(
            'gemtools_gtf_stats',
            input=score,
            threads=threads,
            annotation=args['annotation'],
            output=args['gtf_stats_out'],
            json=args['gtf_stats_json_out'],
            paired_end=True,
            exon_overlap=1,
            counts=args['gtf_counts_out'],
            multi_maps=True,
            weighted=True
        )
        job('Filter.GTF.Stats').run(
            'gemtools_gtf_stats',
            input=filtered,
            threads=threads,
            annotation=args['annotation'],
            output=args['filtered_gtf_stats_out'],
            json=args['gtf_stats_json_out'],
            paired_end=True,
            exon_overlap=1,
            counts=args['filtered_gtf_counts_out'],
            multi_maps=True,
            weighted=True
        )
        clean = job('Cleanup').run('cleanup', files=all_mappings)
        clean << junctions
        clean << denovo_transcriptome.keys_out
        clean << denovo_transcriptome.fasta_out
        clean << denovo_transcriptome.junctions_out
        clean << junctions_index
        if prepare_step:
            clean << prepare_step
        return p

    def run(self, args):
        import jip
        jip.run(self.tool_instance, force=False, dry=self.options['dry'].raw())


@cli("denovo-junctions",
     inputs=['--input'],
     outputs=['--output'])
class JunctionExtraction(Command):
    description = """Run the split mapper to extract junctions"""
    title = "Junction Extraction"

    def register(self, parser):
        parser.add_argument('-I', '--index', dest="index",
                            metavar="<index>",
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality', dest="quality",
                            metavar="<quality>",
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-a', '--annotation', dest="gtf",
                            metavar="<gtf>",
                            help='Path to the .gtf annotation that will '
                            'be used to keep junctions even if the '
                            'site does not have enough coverage.')
        parser.add_argument('-i', '--input', dest="input",
                            metavar="<input>",
                            help='Input File. Reads from stdin if nothing '
                            'is specified.')
        parser.add_argument('-o', '--output', dest="output",
                            metavar="<output>",
                            help='Output file name, prints to stdout '
                            'if nothing is specified')
        parser.add_argument('-t', '--threads', dest="threads",
                            metavar="<threads>",
                            default=1,
                            type=int,
                            help='Number of threads')
        parser.add_argument('-s', '--strata-after-best',
                            dest="strata_after_best", type=int,
                            metavar="strata",
                            default=0,
                            help='The number of strata examined after '
                            'the best one. Default 0')
        parser.add_argument('-m', '--mismatches',
                            dest="mismatches", metavar="mm",
                            type=float,
                            default=0.04,
                            help='Set the allowed mismatch ratio as '
                            '0 < mm < 1. Default 0.04')
        parser.add_argument("--filter",
                            default="ordered,non-zero-distance",
                            help="Split-mapper filtering options")
        parser.add_argument("--consensus",
                            help="Consensus Sequence. Default "
                            "'(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'")
        parser.add_argument("--min-coverage",
                            dest="coverage",
                            type=int,
                            default=2,
                            help="The minimum coverage for a junction to "
                            "be taken into account. Default ")
        parser.add_argument('--min-intron-length',
                            type=int, metavar="<mil>",
                            default=4,
                            help='Minimum intron length. Default 4')
        parser.add_argument('--max-intron-length',
                            type=int, metavar="<mil>",
                            default=500000,
                            help='Maximum intron length. Default 500000')
        parser.add_argument('--refinement-step', dest="refinement_step_size",
                            metavar="<r>",
                            type=int,
                            default=2,
                            help='Refine the minimum split size when '
                            'constraints on number of candidates are '
                            'not met. Default 2')
        parser.add_argument('--min-split-size', dest="min_split_size",
                            default=15,
                            type=int,
                            metavar="<mss>",
                            help='Minimum split length. Default 15')
        parser.add_argument('--matches-threshold',
                            metavar="<mt>",
                            type=int,
                            default=75,
                            help='Maximum number canidates considered '
                            'when splitting the read. Default 75')
        parser.add_argument('--max-matches',
                            metavar="<mm>",
                            type=int,
                            default=5,
                            help='Maximum number of multi-maps '
                            'allowed for a junction. Default 5')

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        from gem.junctions import filter_by_distance
        infile = gem.files.open(args.input) \
            if args.input is not None else sys.stdin
        denovo_junctions = gem.extract_junctions(
            infile,
            args.index,
            filter=args.filter,
            splice_consensus=args.consensus,
            mismatches=args.mismatches,
            threads=args.threads,
            strata_after_first=args.strata_after_best,
            coverage=args.coverage,
            min_split=args.min_intron_length,
            max_split=args.max_intron_length,
            refinement_step_size=args.refinement_step_size,
            min_split_size=args.min_split_size,
            matches_threshold=args.matches_threshold,
            max_junction_matches=args.max_matches,
            annotation=args.gtf,
        )

        logging.gemtools.gt("Found de-novo Junctions %d with coverage >= %s" %
                            (len(denovo_junctions), str(args.coverage)))
        filtered = set(filter_by_distance(denovo_junctions,
                                          args.min_intron_length,
                                          args.max_intron_length))

        logging.gemtools.gt("de-novo junction passing distance "
                            "filter (min: %s max: %s): %d (%s removed)" %
                            (str(args.min_intron_length),
                             str(args.max_intron_length),
                             len(filtered),
                             (len(denovo_junctions) - len(filtered))))

        gem.junctions.write_junctions(filtered, args.output, args.index)
