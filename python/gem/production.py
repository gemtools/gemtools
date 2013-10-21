#!/usr/bin/env python
"""Production pipelines"""
import os
import sys
import subprocess
import re

from gem.commands import cli
import gem
import gem.commands
import gem.gemtools as gt
import logging
from gem.utils import Command, CommandException

from jip import Pipeline

log = logging.getLogger("gem.commands")


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

    def validate(self):
        self.options['input'].check_files()

    def run(self, args):
        if len(args['input']) < 2:
            args['input'].append(args['input'][0])
            args['same'] = True
        files = args['input']
        output = args['output']
        if output is None:
            output = sys.stdout
        p = gem.merge(files[0], files[1:], output,
                      threads=int(args['threads']), same_content=args['same'],
                      paired=args['paired'], compress=args['compress'])
        if p.process:
            if p.process.wait() != 0:
                raise Exception("Merger execution failed!")


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
            if p is not None and os.path.exists(p):
                self.options['input'].append(p)
        return True

    def run(self, args):
        args = gem.utils.dict_2_tuple(args)
        threads = args.threads
        input = args.input

        output = args.output if not isinstance(args.output, file) \
            else args.output
        outfile = gt.OutputFile(output, clean_id=True, append_extra=False)
        input = filter(lambda f: os.path.exists(f), input)
        if len(input) == 1:
            infile = gem.files.open(input[0])
        else:
            infile = gem.filter.interleave(
                [gem.files.open(f) for f in input],
                threads=max(1, threads))
        try:
            infile.write_stream(outfile, write_map=False)
        except ValueError:
            # catch the value error delegate that comes from flush on
            # closed stdout
            pass
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
        parser.add_argument('-I', '--index',
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality',
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-i', '--input', dest="input",
                            default=sys.stdin,
                            help='Input File. Reads from stdin if nothing '
                            'is specified. (default: stdin)')
        parser.add_argument('-o', '--output',
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified. (default: stdout)')
        parser.add_argument('-t', '--threads', dest="threads",
                            default=1,
                            type=int,
                            help='Number of threads')
        parser.add_argument('-s', '--strata-after-best',
                            type=int,
                            default=1,
                            help=""" A stratum is a set of matches all having
                            the same string distance from the query. The GEM
                            mapper is able to find not only the matches
                            belonging to the best stratum (i.e., the best
                            matches having minimum string distance from the
                            query) but also additional sets of matches (the
                            next-to-best matches, the next-to-next-to-best
                            matches, and so on) having alignment score worse
                            than that of the best matches. This parameter
                            determines how many strata should be explored after
                            the best one (i.e., --strata-after-best 1 will list
                                          all the best and all the second best
                                          matches).""")
        parser.add_argument('-m', '--mismatches',
                            type=float,
                            default=0.04,
                            help="""The maximum number of nucleotide
                            substitutions allowed while mapping each read.  It
                            is always guaranteed that, however other options
                            are chosen, all the matches up to the specified
                            number of substitutions will be found by the
                            program. In case qualities are being taken into
                            account (multi-FASTQ input), this parameter assumes
                            the meaning of the maximum number of mismatches
                            which can be found in high-quality bases.""")
        parser.add_argument('--quality-threshold',
                            type=int,
                            default=26,
                            help="""In case of a mapping which takes qualities
                            into account (as deduced from a FASTQ input in
                            Phred or Solexa format), consider as good-quality
                            the bases having a quality score above this
                            threshold, and as bad-quality those having a
                            quality score below this threshold. It should be
                            noted that the interpretation of the specified
                            quality threshold in terms of an error probability
                            differs for Phred and Solexa scales when q<7, or
                            error>=0.17; aside from that, the program takes
                            care to automatically convert this value to the
                            correct ASCII encoding (which again differs
                                                    depending on whether Phred
                                                    or Solexa conventions are
                                                    being used in the
                                                    input).""")
        parser.add_argument('-d', '--max-decoded-matches',
                            type=int,
                            default=1000,
                            help="""The GEM mapper always provides a complete
                            count of all the existing matches up to the
                            selected number of mismatches; however, not all
                            matches are printed, since only a few will be
                            needed for the typical application. This options
                            allows to fine-tune this behaviour. You should
                            specify all only if due to some reason you already
                            know that the maximum number of matches has a
                            reasonable bound (which is not the case for typical
                                              mammalian genomes).""")
        parser.add_argument('-D', '--min-decoded-strata',
                            type=int,
                            default=0,
                            help="""In some occasions (when maximum sensitivity
                            is desirable) it might be useful to be sure that
                            all the matches belonging to a number of strata are
                            always output, irrespectively of their number. By
                            default, the first stratum is always printed in
                            full.  If max_decoded_matches is greater than the
                            number of matches belonging to the strata that
                            should be printed mandatorily, additional strata
                            are possibly printed.""")
        parser.add_argument('--max-big-indel-length',
                            type=int,
                            default=15,
                            help="""The GEM mapper implements a special
                            algorithm that, in addition to ordinary matches, is
                            sometimes able to find a single long indel (in
                            particular, a long insertion in the read).  This
                            option specifies the maximum allowed size for
                            such long indel.""")
        parser.add_argument('--min-matched-bases',
                            type=float,
                            default=0.80,
                            help="""This parameter limits the number of
                            deletions that can occur in the read (if there are
                            too many deletions, the quality of the alignment
                            will be questionable). The default says that at
                            least the 80%% of the bases must be mapped, i.e.
                            there cannot be more than 20%% of the bases
                            deleted.""")
        parser.add_argument('-e', '--max-edit-distance',
                            default=0.20,
                            type=float,
                            help="""The maximum number of edit operations
                            allowed while verifying candidate matches by
                            dynamic programming. Saying --e 0 disables the
                            possibilities of finding alignments with indels
                            (however, "big" indels might still be found, see
                            option --max-big-indel-length).""")
        parser.add_argument('--fast-mapping',
                            const=True,
                            nargs="?",
                            help="""Activates fast mapping modes, whereby the
                            aligner does not align "hard" reads (that is, reads
                            which would require too large a computational
                            budget, usually a few).  Other reads are aligned as
                            in the normal modes.  The parameter number defines
                            the computational budget (and hence --fast-mapping
                            0 will be the cheapest fast mode, --fast mapping 1
                            the next-to-cheapest, and so on). If you just put
                            --fast-mapping, the aligner will use an adaptive
                            mode to figure out a limit.""")
        parser.add_argument('--keys', '-k',
                            dest="keys",
                            metavar="keys",
                            help='Key file to translate result back to '
                            'genomic coordinates')
        parser.add_argument('--mismatch-alphabet',
                            default="ACGT",
                            help="""Specifies the set of characters which are
                            valid replacements in case of mismatch. Note that
                            if you would like to consider Ns in the reference
                            as wildcards, you should specify ACGNT here;
                            otherwise, the mapper will never return positions
                            in the reference containing Ns.""")
        parser.add_argument('-C', '--compress', dest="compress",
                            default=False,
                            action="store_true",
                            help="Compress the output file")

    def validate(self):
        _check("No index specified", not self.options['index'].get())
        if float(self.options['min_matched_bases'].get()) > 1.0:
            self.options['min_matched_bases'].value = int(
                float(self.options['min_matched_bases'].get())
            )
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
            fast_mapping=args.fast_mapping[0] if args.fast_mapping else False,
            #trim=args.trim,
            threads=args.threads,
            compress=args.compress
        )


@cli("split-mapper",
     inputs=['--input'],
     outputs=['--output'])
class GemSplitMapper(Command):
    title = "Run the GEM-split-mapper"
    description = """Runs the GEM-split-mapper"""

    def register(self, parser):
        parser.description = GemSplitMapper.description
        ## required parameters
        parser.add_argument('-I', '--index',
                            help='The GEM index',
                            required=True)
        parser.add_argument('-q', '--quality',
                            help='Quality offset (33/64/ignore)',
                            required=True)
        parser.add_argument('-i', '--input', dest="input",
                            default=sys.stdin,
                            help='Input File. Reads from stdin if nothing '
                            'is specified. (default: stdin)')
        parser.add_argument('-o', '--output',
                            default=sys.stdout,
                            help='Output file name, prints to stdout '
                            'if nothing is specified. (default: stdout)')
        parser.add_argument('-t', '--threads', dest="threads",
                            default=1,
                            type=int,
                            help='Number of threads')
        parser.add_argument('-s', '--strata-after-best',
                            type=int,
                            default=1,
                            help=""" A stratum is a set of matches all having
                            the same string distance from the query. The GEM
                            mapper is able to find not only the matches
                            belonging to the best stratum (i.e., the best
                            matches having minimum string distance from the
                            query) but also additional sets of matches (the
                            next-to-best matches, the next-to-next-to-best
                            matches, and so on) having alignment score worse
                            than that of the best matches. This parameter
                            determines how many strata should be explored after
                            the best one (i.e., --strata-after-best 1 will list
                                          all the best and all the second best
                                          matches).""")
        parser.add_argument('-m', '--mismatches',
                            type=float,
                            default=0.04,
                            help="""The maximum number of nucleotide
                            substitutions allowed while mapping each read.  It
                            is always guaranteed that, however other options
                            are chosen, all the matches up to the specified
                            number of substitutions will be found by the
                            program. In case qualities are being taken into
                            account (multi-FASTQ input), this parameter assumes
                            the meaning of the maximum number of mismatches
                            which can be found in high-quality bases.""")
        parser.add_argument('--quality-threshold',
                            type=int,
                            default=26,
                            help="""In case of a mapping which takes qualities
                            into account (as deduced from a FASTQ input in
                            Phred or Solexa format), consider as good-quality
                            the bases having a quality score above this
                            threshold, and as bad-quality those having a
                            quality score below this threshold. It should be
                            noted that the interpretation of the specified
                            quality threshold in terms of an error probability
                            differs for Phred and Solexa scales when q<7, or
                            error>=0.17; aside from that, the program takes
                            care to automatically convert this value to the
                            correct ASCII encoding (which again differs
                                                    depending on whether Phred
                                                    or Solexa conventions are
                                                    being used in the
                                                    input).""")
        parser.add_argument('--mismatch-alphabet',
                            default="ACGT",
                            help="""Specifies the set of characters which are
                            valid replacements in case of mismatch. Note that
                            if you would like to consider Ns in the reference
                            as wildcards, you should specify ACGNT here;
                            otherwise, the mapper will never return positions
                            in the reference containing Ns.""")
        parser.add_argument('--min-split-size',
                            default=15,
                            type=int,
                            help="""Minimum split length""")
        parser.add_argument('--matches-threshold',
                            default=100,
                            type=int,
                            help="""Maximum number of matches allowed for a
                            split map.""")
        parser.add_argument('--refinement-step-size',
                            type=int,
                            default=2,
                            help="Refine the minimum split size when "
                            "constraints on number of candidates "
                            "are not met.")
        parser.add_argument('-f', '--filter',
                            type=str,
                            help="""Filter for allowed splits. Possible
                            filters: 'same-chromosome' 'same-strand'
                            'minimum-distance='<distance>
                            'maximum-distance='<distance> 'non-zero-distance'
                            'ordered'. You can specify a comma separated list
                            of filters.
                            """)
        parser.add_argument(
            '-c', '--consensus',
            default=(",".join(["%s+%s" % (c[0], c[1])
                               for c in gem.extended_splice_consensus])),
            help="Consensus used to detect junction sites.")
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
        gem.splitmapper(
            infile,
            args.index,
            outfile,
            mismatches=args.mismatches,
            splice_consensus=args.consensus,
            refinement_step_size=args.refinement_step_size,
            min_split_size=args.min_split_size,
            matches_threshold=args.matches_threshold,
            strata_after_first=args.strata_after_best,
            mismatch_alphabet=args.mismatch_alphabet,
            quality=args.quality,
            filter=args.filter,
            trim=None,
            filter_splitmaps=True,
            post_validate=True,
            threads=args.threads,
            compress=args.compress,
            extra=None
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
        parser.add_argument('--filter-max-matches',
                            default=0,
                            type=int,
                            help='Reduce the number of maximum matches')
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
            max_edit_distance=args.max_edit_distance,
            max_extendable_matches=args.max_extendable_matches,
            max_matches_per_extension=args.max_matches_per_extension,
            filter_max_matches=args.filter_max_matches,
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
     outputs=['--output', '--json'])
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
        log.info("Starting:\n\t%s" % (" ".join(cmd)))
        p = subprocess.Popen(cmd, stderr=stderr)
        p.wait()
        if stderr is not None:
            stderr.close()
        sys.exit(p.wait())


@cli("gtf-stats",
     inputs=['--input'],
     outputs=['--output', '--json'])
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
        log.info("Starting:\n\t%s" % (" ".join(cmd)))
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
        log.info("Starting:\n\t%s" % (" ".join(cmd)))
        p = subprocess.Popen(cmd, stdout=outstream)
        r = p.wait()
        if compressor:
            compressor.stdin.close()
            compressor.wait()
        sys.exit(r)


@cli("map2sam",
     inputs=['--input'],
     outputs=['--output'])
class Map2Sam(Command):
    title = "Convert .map to .sam using the gt.map2sam"
    description = """Take a .map file or reads from stdin and converts
    to .sam.
    """

    def register(self, parser):
        parser.description = Map2Sam.description
        self.add_options('gt.map2sam', parser, stream_in="input",
                         stream_out="output")
        parser.add_argument("-m", "--memory",
                            default="768M",
                            help="Memory to use per sorting threads. "
                            "Default 768M")
        parser.add_argument("--no-bam", action="store_true", default=False,
                            help="Disable sorting")
        parser.add_argument("--no-sort", action="store_true", default=False,
                            help="Disable sorting")
        parser.add_argument("--no-index", action="store_true", default=False,
                            help="Disable indexing the bam file")

    def run(self, args):
        original_out = args['output']
        if not args['no_bam']:
            # reset to stdout to pipe to samtools
            args['output'] = sys.stdout

        cmd = self.get_command(args)
        log.info("Starting:\n\t%s" % (" ".join(cmd)))
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE if not args['no_bam'] else None
        )
        if not args['no_bam']:
            gem.sam2bam(gem._prepare_output(p),
                        output=original_out,
                        threads=args['threads'],
                        sorted=not args['no_sort'],
                        sort_memory=args['memory'])
            if not args['no_index'] and isinstance(original_out, basestring):
                # index the bam file if we actually wrote a file
                gem.bamIndex(original_out)
        sys.exit(p.wait())


@cli("scorereads",
     inputs=['--i1', '--i2'],
     outputs=['--output'])
class Scorereads(Command):
    title = "Calculate mapq scores and pair reads"
    description = """Calculate scores and pairs reads

    """

    def register(self, parser):
        parser.description = Scorereads.description
        self.add_options('gt.scorereads', parser, stream_in="i1",
                         stream_out="output")

    def run(self, args):
        cmd = self.get_command(args)
        log.info("Starting:\n\t%s" % (" ".join(cmd)))
        p = subprocess.Popen(cmd)
        sys.exit(p.wait())


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
     pipeline='pipeline')
class TranscriptIndex(Command):
    title = "Create and index transcriptomes"
    description = """This command creates a transcriptome and its index
    from a gem index and a GTF annotation or a junctions file.

    The output name is created from the name of the annotation
    file given if its not explicitly set. The command creates a set of files:


        <name>.junctions       -- the junction sites of the GTF
        <name>.junctions.fa    -- the transcriptome sequences
        <name>.junctions.keys  -- the translation table from transcriptome
                                  to genome coordinates
        <name>.junctions.gem   -- the GEM transcriptome index
    """

    def register(self, parser):
        parser.description = TranscriptIndex.description
        ## required parameters
        parser.add_argument('-i', '--index',
                            help='Path to the GEM genome index',
                            required=True)
        parser.add_argument('-a', '--annotation',
                            dest="annotation",
                            required=True,
                            help='Path to the GTF annotation')
        parser.add_argument('-n', '--name',
                            help='Optional output prefix. '
                            'If this is not set, the annotation'
                            'file name will be used',
                            default=None)
        parser.add_argument('-m', '--max-length',
                            help='Maximum read length, defaults to 150',
                            type=int,
                            default=150)
        parser.add_argument('-t', '--threads',
                            type=int,
                            help='Number of threads',
                            default=2)

    def validate(self):
        args = self.options.to_dict()
        self.options['index'].check_files()
        if not args['index'].endswith(".gem"):
            raise CommandException("No valid GEM index specified, "
                                   "the file has to end in .gem")
        self.options['annotation'].check_files()
        name = self.options['name'].raw()
        if name is None:
            fname = self.options['annotation'].get()
            name = os.path.basename(fname)
            if name.endswith(".gz"):
                name = name[:-3]
            if name.endswith(".gtf"):
                name = name[:-4]
            self.options['name'] = name
        self.options.add_output('junctions_out', name + ".junctions")

    def pipeline(self):
        p = Pipeline()
        gtf_junctions = p.job('GTF-Junctions').run(
            'gemtools_gtf_junctions',
            input=self.options['annotation'].get(),
            output=self.options['junctions_out'].get()
        )
        trans = p.job('Compute-Transcriptome', 1).run(
            'gemtools_compute_transcriptome',
            index=self.options['index'].raw(),
            name=self.options['name'].raw(),
            max_length=self.options['max_length'].raw(),
            junctions=gtf_junctions
        )
        p.job('Index-Transcriptome', self.options['threads'].raw()).run(
            'gemtools_index',
            input=trans.fasta_out,
            threads=self.options['threads'].raw()
        )
        return p

    def run(self, args):
        import jip
        jobs = jip.create_jobs(self.tool_instance)
        # group the jobs
        for group in jip.group(jobs):
            job = group[0]
            name = ", ".join(str(j) for j in group)
            if job.state == jip.STATE_DONE and not args['force']:
                print "Skipping jobs: {name:30} Done".format(name=name)
            else:
                t = gem.utils.Timer()
                print "Running jobs: {name:30} Running".format(name=name)
                success = jip.run_job(job)
                if success:
                    print "Finished jobs: {name:29} {time}" \
                        .format(name=name, time=str(t.stop()))
                else:
                    print "Execution failed for:", name
                    sys.exit(1)


@cli("compute-transcriptome",
     inputs=['--junctions'],
     outputs=['annotation_junctions']
     )
class ComputeTranscriptome(Command):
    title = "Create a transcriptome from a genome index"
    description = """This command creates a transcriptome from a gem index
    and a junctions file. If an annotation is given, the junctions
    in the annotation are substracted from the set of junctions given
    in the junctions file. The latter is used to create denovo-transcriptomes.


    The output name is created from the name of the annotation/junctions
    file given if its not explicitly set. The command creates a set of files:

        <name>.junctions.fa    -- the transcriptome sequences
        <name>.junctions.keys  -- the translation table from transcriptome
                                  to genome coordinates
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
        parser.add_argument('-J', '--annotation-junctions',
                            help='Path to the GTF annotation junctions')
        parser.add_argument('-j', '--junctions',
                            dest="junctions",
                            required=True,
                            help='Path to a .junctions file')
        parser.add_argument('-n', '--name',
                            help='Optional output prefix. '
                            'If this is not set, the junctions '
                            'file name will be used',
                            default=None)
        parser.add_argument('-m', '--max-length',
                            help='Maximum read length, defaults to 150',
                            default=150)

    def validate(self):
        args = self.options.to_dict()
        if not args['index'].endswith(".gem"):
            raise CommandException("No valid GEM index specified, "
                                   "the file has to end in .gem")
        self.options['index'].check_files()
        self.options['annotation'].check_files()
        if not args['junctions']:
            raise CommandException("You have to specify the junctions file")
        self.options['junctions'].check_files()

        name = args['name']
        if name is None:
            fname = args['junctions']
            name = os.path.basename(fname)
            i = name.index('.')
            name = name[:i] if i > 0 else name

        keys_out = name + ".junctions.keys"
        fasta_out = name + ".junctions.fa"

        self.options['name'] = name
        self.options.add_output('keys_out', keys_out)
        self.options.add_output('fasta_out', fasta_out)
        if args['annotation']:
            junctions_out = name + ".gtf.junctions"
            self.options['annotation_junctions'].value = junctions_out
        return True

    def run(self, args):
        substract = None
        if args['annotation']:
            gtf_junctions = set(gem.junctions.from_gtf(args['annotation']))
            gem.junctions.write_junctions(gtf_junctions,
                                          args['annotation_junctions'])
            substract = args['annotation_junctions']
        elif args['annotation_junctions']:
            substract = args['annotation_junctions']

        max_len = 150
        try:
            max_len = int(args['max_length'])
        except:
            # see if this is a json file
            import json
            with open(args['max_length']) as f:
                data = json.load(f)
                max_len = data['general']['read_lenght_max']

        gem.compute_transcriptome(
            max_len,
            args['index'],
            args['junctions'],
            substract=substract,
            output_name=args['name'] + ".junctions"
        )


@cli("rna-pipeline",
     pipeline='pipeline',
     inputs=['--first'])
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
            '-f', '--first',
            required=True,
            help='''Primary input file'''
        )
        input_group.add_argument(
            '-s', '--second',
            help='''Secodary input file'''
        )
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
        input_group.add_argument('-l', '--read-length',
                                 default=None,
                                 type=int,
                                 help="Specify the maximum reads length."
                                 "This is used to compute the "
                                 "denovo-transcriptome and will be calclulated"
                                 " dynamically if not specified")
        input_group.add_argument('--single-end',
                                 action="store_true",
                                 default=False,
                                 help="Single end reads")
        input_group.add_argument('--direct-input',
                                 default=False,
                                 action="store_true",
                                 help="Skip preparation step and pipe the "
                                 "input directly into the first mapping step")
        input_group.add_argument('--no-pair-search',
                                 default=False,
                                 action="store_true",
                                 help="Skip searching for the second pair file"
                                 "")
        input_group.add_argument('--compress-all',
                                 default=False,
                                 action="store_true",
                                 help="Compress also intermediate files")

        ######################################################################
        # Global mapping parameter
        ######################################################################
        initial_group = parser.add_argument_group(
            'Global Mapping Paramter',
            description="The global mapping parameters are applied to all "
            "mapping steps. You can refine the parameters for the transcript "
            "and de-novo mapping steps independently"
        )
        initial_group.add_argument('-m', '--mismatches',
                                   type=float,
                                   default=0.06,
                                   help="Allowed mismatched for initial "
                                   "alignment to the genome index")
        initial_group.add_argument('--quality-threshold',
                                   type=int,
                                   default=26,
                                   help="Read quality threshold. Bases with a "
                                   "quality less than the threshold are "
                                   "treated as bad bases and more mismatches "
                                   "are allowed.")
        initial_group.add_argument("--max-decoded-matches",
                                   type=int,
                                   default=25,
                                   help="Maximum decoded matches. NOTE that "
                                   "this is taken into account only for "
                                   "additional strata. min-decoded-strata "
                                   "are always fully printed.")
        initial_group.add_argument("--min-decoded-strata",
                                   type=int,
                                   default=1,
                                   help="Minimum number of strata that are "
                                   "always fully decoded.")
        initial_group.add_argument("--min-matched-bases",
                                   type=float,
                                   default=0.80,
                                   help="This parameter limits the number of "
                                   "deletions that can occur in the read (if "
                                   "there are too many deletions, the quality "
                                   "of the alignment will be questionable). "
                                   "The default says that at least the 80%% "
                                   "of the bases must be mapped, i.e. there "
                                   "cannot be more than 20%% of the bases "
                                   "deleted.")
        initial_group.add_argument("--max-big-indel-length",
                                   type=int,
                                   default=15,
                                   help="The GEM mapper implements a special "
                                   "algorithm that, in addition to ordinary "
                                   "matches, is sometimes able to find a "
                                   "single long indel (in particular, a long "
                                   "insertion in the read).  This option "
                                   "specifies the maximum allowed size "
                                   "for such long indel.")
        initial_group.add_argument('-e', "--max-edit-distance",
                                   type=float,
                                   default=0.20,
                                   help="The maximum number of edit "
                                   "operations allowed while verifying "
                                   "candidate matches by dynamic "
                                   "programming. Saying --e 0 disables the "
                                   "possibilities of finding alignments with "
                                   "indels (however, 'big' indels might still "
                                   "be found, see option "
                                   "--max-big-indel-length).")
        initial_group.add_argument("--mismatch-alphabet",
                                   type=str,
                                   default='ACTG',
                                   help="Specifies the set of characters "
                                   "which are valid replacements in case of "
                                   "mismatch. Note that if you would like "
                                   "to consider Ns in the reference as "
                                   "wildcards, you should specify ACGNT "
                                   "here; otherwise, the mapper will never "
                                   "return positions in the reference "
                                   "containing Ns.")
        initial_group.add_argument('-S', "--strata-after-best",
                                   type=int,
                                   default=1,
                                   help="A stratum is a set of matches all "
                                   "having the same string distance from the "
                                   "query. The GEM mapper is able to find not "
                                   "only the matches belonging to the best "
                                   "stratum (i.e., the best matches having "
                                   "minimum string distance from the query) "
                                   "but also additional sets of matches "
                                   "(the next-to-best matches, the "
                                   "next-to-next-to-best matches, and so on) "
                                   "having alignment score worse than that of "
                                   "the best matches. This parameter "
                                   "determines how many strata should be "
                                   "explored after the best one")

        initial_group.add_argument('--output-max-matches',
                                   type=int,
                                   default=25,
                                   help='Maximum number of printed matches. '
                                   'Default 25')
        initial_group.add_argument('--output-min-strata',
                                   type=int,
                                   default=1,
                                   help='Minimum number of printed strata. '
                                   'Default 1')
        initial_group.add_argument('--output-max-strata',
                                   type=int,
                                   default=2,
                                   help='Maximum number of printed strata. '
                                   'Default 2')
        ######################################################################
        # Pairing parameter
        ######################################################################
        pair_group = parser.add_argument_group(
            'Pairing parameter'
        )
        pair_group.add_argument('--pairing-quality-threshold',
                                type=int,
                                help="Override read quality threshold for "
                                "pairing. Default to global setting")
        pair_group.add_argument("--pairing-max-decoded-matches",
                                type=int,
                                default=25,
                                help="The number of decoded matches for "
                                "pairing. Default 25")
        pair_group.add_argument("--pairing-min-decoded-strata",
                                type=int,
                                default=1,
                                help="Minimum number of strata that are "
                                "always fully decoded.")
        pair_group.add_argument("--pairing-min-insert-size",
                                type=int,
                                default=0,
                                help="Minimum insert size.")
        pair_group.add_argument("--pairing-max-insert-size",
                                type=int,
                                default=500000,
                                help="Minimum insert size.")
        pair_group.add_argument("--pairing-min-matched-bases",
                                type=float,
                                default=0.80,
                                help="Minimum number of matched bases.")
        pair_group.add_argument("--pairing-max-edit-distance",
                                type=float,
                                default=0.30,
                                help="Maximum edit distance.")
        ######################################################################
        # Transcript mapping parameters
        ######################################################################
        trans_group = parser.add_argument_group(
            'Transcript Mapping Parameters'
        )
        trans_group.add_argument('--transcript-mismatches',
                                 type=float,
                                 help="Allowed mismatched for transcript "
                                 "alignment. Defaults to global settings.")
        trans_group.add_argument('--transcript-quality-threshold',
                                 type=int,
                                 help="Override read quality threshold for "
                                 "transcript mapping")
        trans_group.add_argument("--transcript-max-decoded-matches",
                                 type=int,
                                 default=150,
                                 help="The number of decoded matches for "
                                 "the transcript mapping is increased to a "
                                 "default of 150 to catch all possible "
                                 "mappings.")
        trans_group.add_argument("--transcript-min-decoded-strata",
                                 type=int,
                                 help="Minimum number of strata that are "
                                 "always fully decoded.")
        trans_group.add_argument("--transcript-min-matched-bases",
                                 type=float,
                                 help="Override min matched bases for "
                                 "transcript mapping. Defaults to global "
                                 "setting")
        trans_group.add_argument("--transcript-max-big-indel-length",
                                 type=int,
                                 help="Override big indel length for "
                                 "transcript mapping. Defaults to global "
                                 "setting")
        trans_group.add_argument("--transcript-max-edit-distance",
                                 type=float,
                                 help="Override max edit distance for "
                                 "transcript mapping. Defaults to global "
                                 "setting")
        trans_group.add_argument("--transcript-strata-after-best",
                                 type=int,
                                 help="Override strata after best for "
                                 "transcript mapping. Defaults to global "
                                 "setting")
        ######################################################################
        # Denovo mapping parameters
        ######################################################################
        denovo_group = parser.add_argument_group(
            'De-Novo Mapping Parameters'
        )
        denovo_group.add_argument('--denovo-mismatches',
                                  type=float,
                                  help="Allowed mismatched for denovo "
                                  "alignment. Defaults to global settings.")
        denovo_group.add_argument('--denovo-quality-threshold',
                                  type=int,
                                  help="Override read quality threshold for "
                                  "denovo mapping")
        denovo_group.add_argument("--denovo-max-decoded-matches",
                                  type=int,
                                  default=150,
                                  help="The number of decoded matches for "
                                  "the denovo mapping is increased to a "
                                  "default of 150 to catch all possible "
                                  "mappings.")
        denovo_group.add_argument("--denovo-min-decoded-strata",
                                  type=int,
                                  help="Minimum number of strata that are "
                                  "always fully decoded.")
        denovo_group.add_argument("--denovo-min-matched-bases",
                                  type=float,
                                  help="Override min matched bases for "
                                  "denovo mapping. Defaults to global "
                                  "setting")
        denovo_group.add_argument("--denovo-max-big-indel-length",
                                  type=int,
                                  help="Override big indel length for "
                                  "denovo mapping. Defaults to global "
                                  "setting")
        denovo_group.add_argument("--denovo-max-edit-distance",
                                  type=float,
                                  help="Override max edit distance for "
                                  "denovo mapping. Defaults to global "
                                  "setting")
        denovo_group.add_argument("--denovo-strata-after-best",
                                  type=int,
                                  help="Override strata after best for "
                                  "denovo mapping. Defaults to global "
                                  "setting")
        ######################################################################
        # Split map junction detection paramters
        ######################################################################
        junctions_group = parser.add_argument_group(
            'Split mapper and junction detection parameters'
        )
        junctions_group.add_argument('--junction-mismatches',
                                     type=float,
                                     default=0.04,
                                     help="Allowed mismatched for junction "
                                     "detection alignments.")
        junctions_group.add_argument('--junction-max-matches',
                                     type=int,
                                     default=5,
                                     help="Maximum number of multi-maps "
                                     "allowed for a junction.")
        junctions_group.add_argument('--junction-strata-after-best',
                                     type=int,
                                     default=0,
                                     help="Maximum number of strata to "
                                     "examine after best.")
        junctions_group.add_argument('--min-denovo-intron-length',
                                     type=int,
                                     default=4,
                                     help="Minimum intron length.")
        junctions_group.add_argument('--max-denovo-intron-length',
                                     type=int,
                                     default=500000,
                                     help="Maximum intron length.")
        junctions_group.add_argument('--refinement-step',
                                     type=int,
                                     default=2,
                                     help="Refine the minimum split size when "
                                     "constraints on number of candidates "
                                     "are not met.")
        junctions_group.add_argument('--min-split-size',
                                     type=int,
                                     default=15,
                                     help="Minimum split length.")
        junctions_group.add_argument('--matches-threshold',
                                     type=int,
                                     default=75,
                                     help="Maximum number candidates "
                                     "considered when splitting the read.")
        junctions_group.add_argument('--junction-coverage',
                                     type=int,
                                     default=2,
                                     help="Minimum allowed junction coverage.")
        junctions_group.add_argument(
            '--junction-consensus',
            default=(",".join(["%s+%s" % (c[0], c[1])
                               for c in gem.extended_splice_consensus])),
            help="Consensus used to detect junction sites")
        ######################################################################
        # filter parameters
        ######################################################################
        filter_group = parser.add_argument_group(
            'Filter parameters'
        )
        filter_group.add_argument('--filter-intron-length',
                                  type=int,
                                  default=20,
                                  help="Filter alignment for intron length")
        filter_group.add_argument('--filter-block-length',
                                  type=int,
                                  default=5,
                                  help="Filter alignment for minimum exon "
                                  "overlap")
        filter_group.add_argument('--filter-level',
                                  type=int,
                                  default=0,
                                  help="Reduce multi-maps by uniqueness level")
        filter_group.add_argument('--filter-max-maps',
                                  type=int,
                                  default=5,
                                  help="Exclude alignment with more than "
                                  "max-maps multi-maps.")
        filter_group.add_argument('--filter-max-errors',
                                  type=int,
                                  default=0,
                                  help="Exclude alignment with more than "
                                  "max-error error events.")
        filter_group.add_argument('--filter-all',
                                  default=False,
                                  action="store_true",
                                  help="Do not keep unique mappings but apply "
                                  "the filter to all alignments.")
        filter_group.add_argument('--no-annotation-filter',
                                  action="store_true",
                                  default=False,
                                  help="Do not filter by annotation. "
                                  "The annotation filter checks that pairs "
                                  "and splits fall into the same gene "
                                  "(assuming the gene_id is set in the "
                                  "annotation)")
        ######################################################################
        # filter parameters
        ######################################################################
        counts_group = parser.add_argument_group(
            'Gene counts'
        )
        counts_group.add_argument('--count-no-multi-maps',
                                  action="store_true",
                                  help="Do not count multi-maps")
        counts_group.add_argument('--count-no-weights',
                                  action="store_true",
                                  help="Do not weight multi-maps and "
                                  "multi-gene hits.")
        counts_group.add_argument('--count-exon-threshold',
                                  type=float,
                                  default=1.0,
                                  help="Minimum overlap with exon to be "
                                  "considered.")

        ######################################################################
        # Split map junction detection parameters
        ######################################################################
        bam_group = parser.add_argument_group(
            'SAM/BAM conversion'
        )
        bam_group.add_argument("--no-sort",
                               action="store_true",
                               default=False,
                               help="Do not sort the BAM file")
        bam_group.add_argument("--sort-mem",
                               default="768M",
                               help="Memory passed to samtools for sorting")
        bam_group.add_argument("--no-index",
                               action="store_true",
                               default=False,
                               help="Do not create a BAM index")
        bam_group.add_argument("--no-xs",
                               action="store_true",
                               default=False,
                               help="Do not calculate the XS field for splits")
        ######################################################################
        # Job control
        ######################################################################
        ctrl_group = parser.add_argument_group('Job and Pipeline controls')
        ctrl_group.add_argument("--dry", action="store_true", default=False,
                                help="Show the pipeline configuraiton but "
                                "do not execute the pipeline")
        ctrl_group.add_argument("--force", action="store_true", default=False,
                                help="Force execution")
        ctrl_group.add_argument("--keep", action="store_true", default=False,
                                help="Keep temporary files")
        ctrl_group.add_argument("--load",
                                help="Load a configuration form a file. "
                                "You can use --load/--save to easily persist "
                                "more complex configurations and then load "
                                "and modify them from the command line.")
        ctrl_group.add_argument("--save",
                                help="Save a configuration to a file. "
                                "You can use --load/--save to easily persist "
                                "more complex configurations and then load "
                                "and modify them from the command line.")
        ctrl_group.add_argument("--skip",
                                nargs="*",
                                help="Exclude steps from the pipeline. This "
                                "parameter takes a space separated list of "
                                "job names that will be excluded. Try "
                                "the pipeline with --dry to see a list of "
                                "active jobs")

    def _find_second_pair(self, args):
        ## try to guess the second file
        (n, p) = gem.utils.find_pair(args['first'])
        if p is not None and os.path.exists(p):
            args['second'] = p
            self.options['second'] = p
        # we guessed a name
        if args['name'] is None:
            args['name'] = n

    def _guess_name(self, args):
        if not "first" in args or not args['first']:
            args['name'] = 'unknown'
            return
        name = os.path.basename(args['first'])
        if name.endswith(".gz"):
            name = name[:-3]
        idx = name.rfind(".")
        if idx > 0:
            name = name[:idx]
        args['name'] = name
        if not args['no_pair_search']:
            args['name'] = re.sub(r'[_\.-]\d$', '', name)

    def _file_name(self, suffix=None, args=None, compress=None):
        if compress is None:
            compress = args["compress_all"]
        return "%s%s%s" % (args['name'], "" if suffix is None else suffix,
                           ".gz" if compress else "")

    def _find_transcript_index(self, args, options):
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
            options['transcript_index'] = t_index
            args['transcript_index'] = t_index

        if not args['transcript_keys']:
            # we look for the key file just next to the transcript_index
            # replacing the .gem extension with .keys. If that fails,
            # we try .junctions.keys
            t_keys = rreplace(args['transcript_index'], ".gem", ".keys", 1)
            if not exists(t_keys):
                t_keys = rreplace(args['transcript_index'], ".gem",
                                  ".junctions.keys", 1)
            options['transcript_keys'] = t_keys
            args['transcript_keys'] = t_keys

    def validate(self):
        from os.path import exists
        ###########################################################
        # Check --load and load config from file
        ###########################################################
        if self.options['load']:
            import json
            try:
                with open(self.options['load'].get()) as f:
                    data = json.load(f)
                    for k, v in data.iteritems():
                        opt = self.options[k]
                        if opt is None:
                            raise CommandException(
                                "Configuration key not found in options: %s" %
                                k)
                        if not opt.user_specified:
                            opt.set(v)
            except Exception as err:
                raise CommandException("Error while loading config : %s" % err)

        args = self.options.to_dict()
        ## check for a name or guess it
        if not args['name']:
            self._guess_name(args)

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
        opts.add_output('read_length_stats_out',
                        fn("_rl.stats.json", c=False))

        ###########################################################
        # set transcript and de-novo defaults
        ###########################################################
        def _set_default(prefix, glob):
            name = '%s_%s' % (prefix, glob)
            if opts[name].raw() is None:
                opts[name].value = opts[glob].get()

        _set_default('transcript', 'mismatches')
        _set_default('transcript', 'quality_threshold')
        _set_default('transcript', 'min_decoded_strata')
        _set_default('transcript', 'min_matched_bases')
        _set_default('transcript', 'max_big_indel_length')
        _set_default('transcript', 'max_edit_distance')
        _set_default('transcript', 'strata_after_best')

        _set_default('denovo', 'mismatches')
        _set_default('denovo', 'quality_threshold')
        _set_default('denovo', 'min_decoded_strata')
        _set_default('denovo', 'min_matched_bases')
        _set_default('denovo', 'max_big_indel_length')
        _set_default('denovo', 'max_edit_distance')
        _set_default('denovo', 'strata_after_best')

        _set_default('pairing', 'quality_threshold')

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
            self._find_transcript_index(args, self.options)
            _check("No Transcript index found: %s" %
                   args['transcript_index'],
                   not exists(args['transcript_index']))
            _check("No Transcript keys found: %s" %
                   args['transcript_keys'],
                   not exists(args['transcript_keys']))

        if self.options['save']:
            import json
            try:
                ex_opts = ['help', 'dry', 'load', 'save']
                values = {}
                for o in filter(lambda o: o.name not in ex_opts,
                                self.options):
                    if o.raw() != o.default:
                        values[o.name] = o.raw()
                with open(self.options['save'].get(), 'w') as f:
                    json.dump(values, f, indent=4)
                print "Saved configuration in %s" % self.options['save']
            except Exception as err:
                raise CommandException("Error while saving config : %s" % err)

        ###########################################################
        # General input data checks
        ###########################################################
        # check for at least one input file
        _check("No input file specified!", not args['first'])
        # check input files for single end alignments
        _check("Single end runs take only one input file!",
               args['single_end'] and args['second'])
        if not args['single_end'] and not args['no_pair_search']:
            # check paired end
            if not args['second']:
                self._find_second_pair(args)
        ## check that all input files exists
        self.options['first'].check_files()
        self.options['second'].check_files()
        return True

    def pipeline(self):
        args = self.options.to_dict()
        threads = int(args['threads'])
        quality = args['quality']
        index = args['index']

        p = Pipeline()
        name = args['name']
        if not name:
            self._guess_name(args)
        name = args['name']
        p.name(name)
        job = p.job(threads=threads)
        prepare_step = None

        # small helper to create
        # a prepare stream quickly
        def create_prepare(name):
            infiles = [args['first']]
            if args['second']:
                infiles.append(args['second'])
            return job(name).run(
                'gemtools_prepare',
                input=infiles,
                threads=threads
            )

        # set the input of a mapping job
        def set_mapping_input(map_job):
            if prepare_step is not None:
                map_job.input = prepare_step
            else:
                create_prepare("Prepare") | map_job
            return map_job

        ###################################################################
        # Prepare
        ###################################################################
        if not args['direct_input']:
            infiles = [args['first']]
            if args['second']:
                infiles.append(args['second'])
            ## run the prepare step
            prepare_step = job('Prepare', temp=True).run(
                'gemtools_prepare',
                input=infiles,
                output=args['prepare_out'],
                threads=threads
            )

        # we collect all mapping steps here for mergin
        all_mappings = []

        ###################################################################
        # Initial mapping step
        ###################################################################
        initial_mapping = job('Initial.Mapping', temp=True).run(
            'gemtools_mapper',
            threads=threads,
            index=index,
            output=args['initial_map_out'],
            compress=args['compress_all'],
            quality=quality,
            mismatches=args['mismatches'],
            quality_threshold=args['quality_threshold'],
            max_decoded_matches=args['max_decoded_matches'],
            min_decoded_strata=args['min_decoded_strata'],
            min_matched_bases=args['min_matched_bases'],
            max_big_indel_length=args['max_big_indel_length'],
            max_edit_distance=args['max_edit_distance'],
            mismatch_alphabet=args['mismatch_alphabet'],
            strata_after_best=args['strata_after_best']
        )
        set_mapping_input(initial_mapping)
        all_mappings.append(initial_mapping)

        ###################################################################
        # Transcriptome mapping
        ###################################################################
        if args['annotation']:
            t_mapping = set_mapping_input(job('GTF.Mapping').run(
                'gemtools_mapper',
                threads=threads,
                index=args['transcript_index'],
                keys=args['transcript_keys'],
                quality=quality,
                mismatch_alphabet=args['mismatch_alphabet'],
                mismatches=args['transcript_mismatches'],
                quality_threshold=args['transcript_quality_threshold'],
                max_decoded_matches=args['transcript_max_decoded_matches'],
                min_decoded_strata=args['transcript_min_decoded_strata'],
                min_matched_bases=args['transcript_min_matched_bases'],
                max_big_indel_length=args['transcript_max_big_indel_length'],
                max_edit_distance=args['transcript_max_edit_distance'],
                strata_after_best=args['transcript_strata_after_best']
            ))
            t_filter = job('GTF.Filter', temp=True).run(
                'gemtools_filter',
                only_split_maps=True,
                threads=threads,
                compress=args['compress_all'],
                output=args['transcript_map_out']
            )

            transcript_mapping = t_mapping | t_filter
            all_mappings.append(transcript_mapping)

        ###################################################################
        # Denovo transcript mapping
        ###################################################################
        junctions = set_mapping_input(job('Denovo.Junctions', temp=True).run(
            'gemtools_denovo_junctions',
            index=args['index'],
            threads=threads,
            quality=quality,
            output=args["denovo_junctions_out"],
            mismatches=args['junction_mismatches'],
            max_matches=args['junction_max_matches'],
            strata_after_best=args['junction_strata_after_best'],
            min_intron_length=args['min_denovo_intron_length'],
            max_intron_length=args['max_denovo_intron_length'],
            consensus=args['junction_consensus'],
            coverage=args['junction_coverage'],
            refinement_step_size=args['refinement_step'],
            min_split_size=args['min_split_size'],
            matches_threshold=args['matches_threshold']
        ))
        if not args['read_length']:
            # calculate simple stats on the initial mapping
            # to get the max_read length
            rl_stats = job('ReadLength.Stats', temp=True).run(
                'gemtools_stats',
                input=initial_mapping,
                output_format='json',
                output=args['read_length_stats_out'],
                threads=threads
            )
            args['read_length'] = rl_stats.output
        denovo_transcriptome = job('Denovo.Transcriptome', 1, temp=True).run(
            'gemtools_compute_transcriptome',
            index=args['index'],
            annotation=args['annotation'] if args['annotation'] else None,
            junctions=junctions.output,
            max_length=args['read_length'])
        # this is a workaround for a jip bug! Somehow
        # the tool is not valudated and therefore the output
        # is not added ?
        #denovo_transcriptome._tool.validate()

        junctions_index = job('Denovo.Index', temp=True).run(
            'gemtools_index',
            threads=threads,
            input=denovo_transcriptome.fasta_out
        )
        d_map = set_mapping_input(job('Denovo.Mapping').run(
            'gemtools_mapper',
            index=junctions_index,
            threads=threads,
            keys=denovo_transcriptome.keys_out,
            compress=args['compress_all'],
            quality=quality,
            mismatch_alphabet=args['mismatch_alphabet'],
            mismatches=args['denovo_mismatches'],
            quality_threshold=args['denovo_quality_threshold'],
            max_decoded_matches=args['denovo_max_decoded_matches'],
            min_decoded_strata=args['denovo_min_decoded_strata'],
            min_matched_bases=args['denovo_min_matched_bases'],
            max_big_indel_length=args['denovo_max_big_indel_length'],
            max_edit_distance=args['denovo_max_edit_distance'],
            strata_after_best=args['denovo_strata_after_best']
        ))
        d_filter = job('Denovo.Filter', temp=True).run(
            'gemtools_filter',
            only_split_maps=True,
            threads=threads,
            compress=args['compress_all'],
            output=args['denovo_map_out']
        )
        denovo_mapping = d_map | d_filter
        all_mappings.append(denovo_mapping)

        ###################################################################
        # Merge and pair if not single end
        ###################################################################
        merge = job('Merge').run(
            'gemtools_merge',
            same=True,
            threads=threads,
            input=all_mappings)

        score = job('Score').run(
            'gemtools_score',
            index=args['index'],
            filter=",".join(str(i) for i in (args['output_min_strata'],
                                             (args['output_max_strata'] -
                                              args['output_min_strata']),
                                             args['output_max_matches'])),
            threads=threads,
            quality=quality,
            compress=True,
            output=args['final_out'])

        if not args['single_end']:
            pair = job('Pair').run(
                'gemtools_pairalign',
                quality=quality,
                threads=threads,
                index=args['index'],
                max_decoded_matches=args['pairing_max_decoded_matches'],
                min_decoded_strata=args['pairing_min_decoded_strata'],
                min_insert_size=args['pairing_min_insert_size'],
                max_insert_size=args['pairing_max_insert_size'],
                max_edit_distance=args['pairing_max_edit_distance'],
                max_extendable_matches=0,
                max_matches_per_extension=0,
                filter_max_matches=0,
            )
            merge | pair | score
        else:
            merge | score

        # create sorted bam
        bam_cfg = dict(
            threads=threads,
            index=args['index'],
            quality=quality,
            memory=args['sort_mem'],
            paired=not args['single_end'],
            no_sort=args['no_sort'],
            no_xs=args['no_xs'],
            no_index=args['no_index']
        )
        job('BAM').run(
            'gemtools_convert',
            input=score,
            output=args['bam_out'],
            **bam_cfg)
        # create filtered output
        filtered = job('Filtered').run(
            'gemtools_filter',
            input=score,
            threads=threads,
            output=args['filtered_out'],
            annotation=args['annotation'],
            no_penalty_for_splitmaps=True,
            paired_end=not args['single_end'],  # paired
            keep_unique=not args['filter_all'],
            min_intron_length=args['filter_intron_length'],
            min_block_length=args['filter_block_length'],
            reduce_by_gene_id=not args['no_annotation_filter'],
            reduce_by_junctions=not args['no_annotation_filter'],
            reduce_to_unique_strata=args['filter_level'],
            reduce_to_max_maps=args['filter_max_maps'],
            max_strata=args['filter_max_errors'],  # max error events
            compress=True)
        job('Filtered.BAM').run(
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
            paired_end=not args['single_end'],
            all_tests=True)

        job('Filtered.Stats').run(
            'gemtools_stats',
            input=filtered,
            threads=threads,
            output=args['filtered_stats_out'],
            json=args['filtered_stats_json_out'],
            paired_end=not args['single_end'],
            all_tests=True)

        # create counts
        job('GTF.Stats').run(
            'gemtools_gtf_stats',
            input=score,
            threads=threads,
            annotation=args['annotation'],
            output=args['gtf_stats_out'],
            json=args['gtf_stats_json_out'],
            paired_end=not args['single_end'],
            exon_overlap=args['count_exon_threshold'],
            counts=args['gtf_counts_out'],
            multi_maps=not args['count_no_multi_maps'],
            weighted=not args['count_no_weights']
        )
        job('Filtered.GTF.Stats').run(
            'gemtools_gtf_stats',
            input=filtered,
            threads=threads,
            annotation=args['annotation'],
            output=args['filtered_gtf_stats_out'],
            json=args['filtered_gtf_stats_json_out'],
            paired_end=not args['single_end'],
            exon_overlap=args['count_exon_threshold'],
            counts=args['filtered_gtf_counts_out'],
            multi_maps=not args['count_no_multi_maps'],
            weighted=not args['count_no_weights']
        )
        if args['keep']:
            p.excludes.append('cleanup')
        if args['skip']:
            p.excludes.extend(args['skip'])
        return p


@cli("vc-pipeline",
     pipeline='pipeline',
     inputs=['--first'])
class VcPipeline(Command):
    description = """Variant calling pipeline.

    Input file detection: If you do not specify --single to disable read
    pairing, we look automatically for the second pair file if you only
    specify one file. For that to work, the second file has to end with
    either .2 or _2, with the file extension .fastq or .txt (+ .gz for
    compressed files). For example,

    gemtools vc-pipeline -f myinput_1.fastq.gz ...

    will search for a file myinput_2.fastq.gz and use it as the second
    pair file.
    """
    title = "GEMTools VCSeq Pipeline"

    def register(self, parser):
        parser.description = RnaPipeline.description
        input_group = parser.add_argument_group('Input and names')
        ## general pipeline paramters
        input_group.add_argument(
            '-f', '--first',
            required=True,
            help='''Primary input file'''
        )
        input_group.add_argument(
            '-s', '--second',
            help='''Secodary input file'''
        )
        input_group.add_argument('-q', '--quality',
                                 default="33",
                                 help='Quality offset. 33, 64 or '
                                 '"ignore" to disable qualities.',
                                 required=True)
        input_group.add_argument('-i', '--index',
                                 help='Path to the .gem genome index',
                                 required=True)
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
        input_group.add_argument('--no-pair-search',
                                 default=False,
                                 action="store_true",
                                 help="Skip searching for the second pair file"
                                 "")
        input_group.add_argument('--compress-all',
                                 default=False,
                                 action="store_true",
                                 help="Compress also intermediate files")

        ######################################################################
        # Global mapping parameter
        ######################################################################
        initial_group = parser.add_argument_group(
            'Global Mapping Paramter',
            description="The global mapping parameters are applied to all "
            "mapping steps. You can refine the parameters for the transcript "
            "and de-novo mapping steps independently"
        )
        initial_group.add_argument('-m', '--mismatches',
                                   type=float,
                                   default=0.06,
                                   help="Allowed mismatched for initial "
                                   "alignment to the genome index")
        initial_group.add_argument('--quality-threshold',
                                   type=int,
                                   default=26,
                                   help="Read quality threshold. Bases with a "
                                   "quality less than the threshold are "
                                   "treated as bad bases and more mismatches "
                                   "are allowed.")
        initial_group.add_argument("--max-decoded-matches",
                                   type=int,
                                   default=25,
                                   help="Maximum decoded matches. NOTE that "
                                   "this is taken into account only for "
                                   "additional strata. min-decoded-strata "
                                   "are always fully printed.")
        initial_group.add_argument("--min-decoded-strata",
                                   type=int,
                                   default=1,
                                   help="Minimum number of strata that are "
                                   "always fully decoded.")
        initial_group.add_argument("--min-matched-bases",
                                   type=float,
                                   default=0.80,
                                   help="This parameter limits the number of "
                                   "deletions that can occur in the read (if "
                                   "there are too many deletions, the quality "
                                   "of the alignment will be questionable). "
                                   "The default says that at least the 80%% "
                                   "of the bases must be mapped, i.e. there "
                                   "cannot be more than 20%% of the bases "
                                   "deleted.")
        initial_group.add_argument("--max-big-indel-length",
                                   type=int,
                                   default=15,
                                   help="The GEM mapper implements a special "
                                   "algorithm that, in addition to ordinary "
                                   "matches, is sometimes able to find a "
                                   "single long indel (in particular, a long "
                                   "insertion in the read).  This option "
                                   "specifies the maximum allowed size "
                                   "for such long indel.")
        initial_group.add_argument('-e', "--max-edit-distance",
                                   type=float,
                                   default=0.20,
                                   help="The maximum number of edit "
                                   "operations allowed while verifying "
                                   "candidate matches by dynamic "
                                   "programming. Saying --e 0 disables the "
                                   "possibilities of finding alignments with "
                                   "indels (however, 'big' indels might still "
                                   "be found, see option "
                                   "--max-big-indel-length).")
        initial_group.add_argument("--mismatch-alphabet",
                                   type=str,
                                   default='ACTG',
                                   help="Specifies the set of characters "
                                   "which are valid replacements in case of "
                                   "mismatch. Note that if you would like "
                                   "to consider Ns in the reference as "
                                   "wildcards, you should specify ACGNT "
                                   "here; otherwise, the mapper will never "
                                   "return positions in the reference "
                                   "containing Ns.")
        initial_group.add_argument('-S', "--strata-after-best",
                                   type=int,
                                   default=1,
                                   help="A stratum is a set of matches all "
                                   "having the same string distance from the "
                                   "query. The GEM mapper is able to find not "
                                   "only the matches belonging to the best "
                                   "stratum (i.e., the best matches having "
                                   "minimum string distance from the query) "
                                   "but also additional sets of matches "
                                   "(the next-to-best matches, the "
                                   "next-to-next-to-best matches, and so on) "
                                   "having alignment score worse than that of "
                                   "the best matches. This parameter "
                                   "determines how many strata should be "
                                   "explored after the best one")
        ######################################################################
        # Split map junction detection paramters
        ######################################################################
        split_group = parser.add_argument_group(
            'Split mapper parameter'
        )
        split_group.add_argument('--split-errors',
                                 type=int,
                                 default=0,
                                 help="""If set to > 0, not only unmapped reads
                                 but also reads mapped with split-error error
                                 are passed on to the split mapper.""")
        split_group.add_argument('--junction-mismatches',
                                 type=float,
                                 default=0.04,
                                 help="Allowed mismatched for junction "
                                 "detection alignments. Default 0.04")
        split_group.add_argument('--junction-max-matches',
                                 type=int,
                                 default=5,
                                 help="Maximum number of multi-maps "
                                 "allowed for a junction. Default 5")
        split_group.add_argument('--junction-strata-after-best',
                                 type=int,
                                 default=0,
                                 help="Maximum number of strata to "
                                 "examine after best. Default 0")
        split_group.add_argument('--min-denovo-intron-length',
                                 type=int,
                                 default=4,
                                 help="Minimum intron length. Default 4")
        split_group.add_argument('--max-denovo-intron-length',
                                 type=int,
                                 default=500000,
                                 help="Maximum intron length. "
                                 "Default 500000")
        split_group.add_argument('--refinement-step',
                                 type=int,
                                 default=2,
                                 help="Refine the minimum split size when "
                                 "constraints on number of candidates "
                                 "are not met. Default 2")
        split_group.add_argument('--min-split-size',
                                 type=int,
                                 default=15,
                                 help="Minimum split length. Default 15")
        split_group.add_argument('--matches-threshold',
                                 type=int,
                                 default=75,
                                 help="Maximum number candidates "
                                 "considered when splitting the read. "
                                 "Default 75")
        split_group.add_argument('--junction-coverage',
                                 type=int,
                                 default=2,
                                 help="Minimum allowed junction coverage. "
                                 "Default 2")

        ######################################################################
        # SAM/BAM parameter
        ######################################################################
        bam_group = parser.add_argument_group(
            'SAM/BAM conversion'
        )
        bam_group.add_argument("--no-sort",
                               action="store_true",
                               default=False,
                               help="Do not sort the BAM file")
        bam_group.add_argument("--sort-mem",
                               default="768M",
                               help="Memory passed to samtools for sorting")
        bam_group.add_argument("--no-index",
                               action="store_true",
                               default=False,
                               help="Do not create a BAM index")
        ######################################################################
        # Job control
        ######################################################################
        ctrl_group = parser.add_argument_group('Job and Pipeline controls')
        ctrl_group.add_argument("--dry", action="store_true", default=False,
                                help="Show the pipeline configuraiton but "
                                "do not execute the pipeline")
        ctrl_group.add_argument("--force", action="store_true", default=False,
                                help="Force execution")
        ctrl_group.add_argument("--keep", action="store_true", default=False,
                                help="Keep temporary files")
        ctrl_group.add_argument("--load",
                                help="Load a configuration form a file. "
                                "You can use --load/--save to easily persist "
                                "more complex configurations and then load "
                                "and modify them from the command line.")
        ctrl_group.add_argument("--save",
                                help="Save a configuration to a file. "
                                "You can use --load/--save to easily persist "
                                "more complex configurations and then load "
                                "and modify them from the command line.")
        ctrl_group.add_argument("--skip",
                                nargs="*",
                                help="Exclude steps from the pipeline. This "
                                "parameter takes a space separated list of "
                                "job names that will be excluded. Try "
                                "the pipeline with --dry to see a list of "
                                "active jobs")

    def _find_second_pair(self, args):
        ## try to guess the second file
        (n, p) = gem.utils.find_pair(args['first'])
        if p is not None and os.path.exists(p):
            args['second'] = p
            self.options['second'] = p
        # we guessed a name
        if args['name'] is None:
            args['name'] = n

    def _guess_name(self, args):
        if not "first" in args or not args['first']:
            args['name'] = 'unknown'
            return
        name = os.path.basename(args['first'])
        if name.endswith(".gz"):
            name = name[:-3]
        idx = name.rfind(".")
        if idx > 0:
            name = name[:idx]
        args['name'] = name
        if not args['no_pair_search']:
            args['name'] = re.sub(r'[_\.-]\d$', '', name)

    def _file_name(self, suffix=None, args=None, compress=None):
        if compress is None:
            compress = args["compress_all"]
        return "%s%s%s" % (args['name'], "" if suffix is None else suffix,
                           ".gz" if compress else "")

    def validate(self):
        from os.path import exists
        ###########################################################
        # Check --load and load config from file
        ###########################################################
        if self.options['load']:
            import json
            try:
                with open(self.options['load'].get()) as f:
                    data = json.load(f)
                    for k, v in data.iteritems():
                        opt = self.options[k]
                        if opt is None:
                            raise CommandException(
                                "Configuration key not found in options: %s" %
                                k)
                        if not opt.user_specified:
                            opt.set(v)
            except Exception as err:
                raise CommandException("Error while loading config : %s" % err)

        args = self.options.to_dict()
        ## check for a name or guess it
        if not args['name']:
            self._guess_name(args)

        ###########################################################
        # Output file configuration
        ###########################################################
        def fn(name, c=None):
            return self._file_name(name, args, compress=c)

        opts = self.options
        opts.add_output('prepare_out', fn("_prepare.fq", c=False))
        opts.add_output('initial_map_out', fn("_initial.map"))
        opts.add_output('split_out', fn(".map", c=False))
        opts.add_output('final_out', fn(".map", c=True))
        opts.add_output('bam_out', fn(".bam", c=False))
        opts.add_output('stats_out', fn(".stats.txt", c=False))
        opts.add_output('stats_json_out', fn(".stats.json", c=False))

        ###########################################################
        # set transcript and de-novo defaults
        ###########################################################
        #def _set_default(prefix, glob):
            #name = '%s_%s' % (prefix, glob)
            #if opts[name].raw() is None:
                #opts[name].value = opts[glob].get()

        #_set_default('pairing', 'quality_threshold')

        ## check the .gem index
        _check("Genome GEM index not found: %s" % args['index'],
               not exists(args['index']))
        _check("Genome GEM does not end in .gem: %s" % args['index'],
               not args['index'].endswith(".gem"))

        ## check quality value
        _check("Quality offset value is not valid! Supported are 33|64|ignore",
               not args['quality'] in ("33", "64", "ignore"))

        if self.options['save']:
            import json
            try:
                ex_opts = ['help', 'dry', 'load', 'save', 'name']
                values = {}
                for o in filter(lambda o: o.name not in ex_opts,
                                self.options):
                    if o.raw() != o.default:
                        values[o.name] = o.raw()
                with open(self.options['save'].get(), 'w') as f:
                    json.dump(values, f, indent=4)
                print "Saved configuration in %s" % self.options['save']
                return False

            except Exception as err:
                raise CommandException("Error while saving config : %s" % err)

        ###########################################################
        # General input data checks
        ###########################################################
        # check for at least one input file
        _check("No input file specified!", not args['first'])
        # check input files for single end alignments
        _check("Single end runs take only one input file!",
               args['single_end'] and args['second'])
        if not args['single_end'] and not args['no_pair_search']:
            # check paired end
            if not args['second']:
                self._find_second_pair(args)
        ## check that all input files exists
        self.options['first'].check_files()
        self.options['second'].check_files()
        return True

    def pipeline(self):
        args = self.options.to_dict()
        threads = int(args['threads'])
        quality = args['quality']
        index = args['index']

        p = Pipeline()
        name = args['name']
        if not name:
            self._guess_name(args)
        name = args['name']
        p.name(name)
        job = p.job(threads=threads)
        prepare_step = None

        # small helper to create
        # a prepare stream quickly
        def create_prepare(name):
            infiles = [args['first']]
            if args['second']:
                infiles.append(args['second'])
            return job(name).run(
                'gemtools_prepare',
                input=infiles,
                threads=threads
            )

        # set the input of a mapping job
        def set_mapping_input(map_job):
            if prepare_step is not None:
                map_job.input = prepare_step
            else:
                create_prepare("Prepare") | map_job
            return map_job

        ###################################################################
        # Prepare
        ###################################################################
        if not args['direct_input']:
            infiles = [args['first']]
            if args['second']:
                infiles.append(args['second'])
            ## run the prepare step
            prepare_step = job('Prepare', temp=True).run(
                'gemtools_prepare',
                input=infiles,
                output=args['prepare_out'],
                threads=threads
            )

        # we collect all mapping steps here for mergin
        all_mappings = []

        ###################################################################
        # Initial mapping step
        ###################################################################
        initial_mapping = job('Initial.Mapping', temp=True).run(
            'gemtools_mapper',
            threads=threads,
            index=index,
            output=args['initial_map_out'],
            compress=args['compress_all'],
            quality=quality,
            mismatches=args['mismatches'],
            quality_threshold=args['quality_threshold'],
            max_decoded_matches=args['max_decoded_matches'],
            min_decoded_strata=args['min_decoded_strata'],
            min_matched_bases=args['min_matched_bases'],
            max_big_indel_length=args['max_big_indel_length'],
            max_edit_distance=args['max_edit_distance'],
            mismatch_alphabet=args['mismatch_alphabet'],
            strata_after_best=args['strata_after_best']
        )
        set_mapping_input(initial_mapping)
        all_mappings.append(initial_mapping)

        ###################################################################
        # Run the split mapper un unmapped or max errors
        ###################################################################
        #filter_maps = job('Filter.Initial', temp=True).run(
            #'gemtools_filter',
            #unmapped=True,
        #)
        ###################################################################
        # Merge and pair if not single end
        ###################################################################
        final_mapping = None
        if len(all_mappings) == 1:
            final_mapping = all_mappings[0]
        else:
            merge = job('Merge').run(
                'gemtools_merge',
                same=True,
                threads=threads,
                output=args['final_out'] if args['single_end'] else None,
                input=all_mappings,
                compress=args['single_end'])
            final_mapping = merge

        pair = job('Score').run(
            'gemtools_scorereads',
            threads=threads,
            output=args['final_out'],
            paired_end=not args['single_end'],
            gzip=True
        )
        final_mapping >> pair
        final_mapping = pair

        job('BAM').run(
            'gemtools_map2sam',
            input=final_mapping,
            output=args['bam_out'],
            gem_index=args['index'],
            paired_end=not args['single_end'],
            threads=threads,
            memory=args['sort_mem'],
            no_sort=args['no_sort'],
            no_index=args['no_index']
        )

        # create stats
        job('Stats').run(
            'gemtools_stats',
            input=final_mapping,
            threads=threads,
            output=args['stats_out'],
            json=args['stats_json_out'],
            paired_end=not args['single_end'],
            all_tests=True)

        if args['keep']:
            p.excludes.append('cleanup')
        if args['skip']:
            p.excludes.extend(args['skip'])
        return p


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
                            default=sys.stdin,
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
        if isinstance(args.input, file):
            infile = args.input
        else:
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

        filtered = set(filter_by_distance(denovo_junctions,
                                          args.min_intron_length,
                                          args.max_intron_length))
        gem.junctions.write_junctions(filtered, args.output, args.index)
