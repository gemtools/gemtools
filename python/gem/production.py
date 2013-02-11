#!/usr/bin/env python
"""Production pipelines"""
#!/usr/bin/env python
import sys
import os
import time
import argparse
import logging
import gem
import gem.commands

from gem.filter import interleave, unmapped
from gem.filter import filter as gf
from gem.pipeline import MappingPipeline
from gem.utils import Command, CommandException
import gem.gemtools as gt


class Merge(Command):
    title = "Merge .map files"
    description = """Merge two .map files. The first file has to
    be the master file that contains all the reads, the second file can
    contain a subset of the reads with the same ID tags and the same order.
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Master file with all the reads', required=True)
        parser.add_argument('-s', '--second', dest="second", help='Second file with a subset of the reads in the same order', required=True)
        parser.add_argument('-t', '--threads', dest="threads", default=1, help='Number of threads')
        parser.add_argument('-o', '--output', dest="output", help='Output file name, prints to stdout if nothing is specified')
        parser.add_argument('--same', dest="same", action="store_true", default=False, help="File contain the same reads")

    def run(self, args):
        i1 = args.input
        i2 = args.second
        gem.merger(gem.files.open(i1), [gem.files.open(i2)]).merge(args.output, threads=int(args.threads), same_content=args.same)


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
        input_group = parser.add_argument_group('Input')

        ## required parameters
        input_group.add_argument('-f', '--files', dest="file", nargs="+", help='''Single fastq input file or both files for a paired-end run separated by space.
            Note that if you specify only one file, we will look for the pair counter-part
            automatically and start a paired-end run. Add the --single-end parameter to disable
            pairing and file search. The file search for the second pair detects pairs
            ending in [_|.][1|2].[fastq|txt][.gz].''', required=True)
        input_group.add_argument('-q', '--quality', dest="quality", default=33, help='Quality offset. 33, 64 or "ignore". Default 33')
        # the index
        input_group.add_argument('-i', '--index', dest="index", help='Path to the .gem genome index', required=True)
        # annotation parameters
        input_group.add_argument('-a', '--annotation', dest="annotation", help='''Path to the GTF annotation. If specified the transcriptome generated from teh annotation is
            used in addition to denovo junctions.''')
        input_group.add_argument('-r', '--transcript-index', dest="trans_index", help='''GTF Transcriptome index. If not specified and an annotation is given,
            it is assumed to be <gtf>.gem. ''')
        input_group.add_argument('-k', '--keys', dest="transcript_keys", help='''Transcriptome .keys file. If not specified and an annotation is given,
            it is assumed to be <gtf>.junctions.keys''')

        output_group = parser.add_argument_group('Output')
        output_group.add_argument('-n', '--name', dest="name", help="Name used for the results")
        output_group.add_argument('-o', '--output-dir', dest="output", help='Optional output folder. If not specified the current working directory is used.')
        output_group.add_argument('--no-sam', dest="nosam", action="store_true", default=False, help="Do not create sam/bam file")
        output_group.add_argument('-g', '--no-gzip', dest="gzip", action="store_false", default=True, help="Do not compress result mapping")
        output_group.add_argument('--keep-temp', dest="rmtemp", action="store_false", default=True, help="Keep temporary files")

        mapping_group = parser.add_argument_group('Mapping')
        mapping_group.add_argument('-m', '--max-mismatches', dest="max_mismatches", help='Number of mismatches allowed in genome mapping step. Default 0.06', default=0.06)

        transcript_group = parser.add_argument_group('Transcript Mapping')
        transcript_group.add_argument('--junction-coverage', dest="junctioncoverage", help='A denovo junction must be covered by > coverage reads to be taken into account, 0 to disable. Default 2', default=2)
        transcript_group.add_argument('--junction-length', dest="junctionlength", help='Maximum length of a denovo junction. Default 500000', default=500000)
        transcript_group.add_argument('-tm', '--transcript-max-mismatches', dest="transcript_max_mismatches", help='Number of mismatches allowed in transcriptome mapping step. Defaults to allowed mismatchs in genome mapping', default=None)
        transcript_group.add_argument('--max-read-length', dest="maxlength", default=150, help='''The maximum read length. This is used to create the denovo
            transcriptom and acts as an upper bound. The default is 150.''', required=True)

        general_group = parser.add_argument_group('General')
        general_group.add_argument('-s', '--strata-after-best', dest="delta", help='Number of strata that are examined after the best one. Default 1', default=1)
        general_group.add_argument('-t', '--threads', dest="threads", default=2, type=int, help="Number of threads to use")

    def run(self, args):
        ## parsing command line arguments
        input_file = os.path.abspath(args.file[0])
        input_file2 = None
        name = None
        if len(args.file) > 1:
            input_file2 = os.path.abspath(args.file[1])
            name = os.path.basename(input_file)[:input_file.rfind(".")]
        else:
            (name, input_file2) = gem.utils.find_pair(input_file)

        if args.name:
            name = args.name

        if args.extendname:
            name = "%s_%d_%d" % (name, args.delta, args.junctioncoverage)
        logging.info("Using dataset name: %s" % (name))

        pipeline = MappingPipeline(
            name=name,
            index=args.index,
            output_dir=args.output,
            annotation=args.annotation,
            threads=int(args.threads),
            quality=int(args.quality),
            junctioncoverage=int(args.junctioncoverage),
            maxlength=int(args.maxlength),
            transcript_index=args.trans_index,
            transcript_keys=args.transcript_keys,
            delta=int(args.delta),
            remove_temp=args.rmtemp
        )

        timer = gem.utils.Timer()

        # main_input = gem.files.open(input_file)
        # main_input2 = gem.files.open(input_file2)

        pipeline.add_step("Initial-Mapping", "mapping")
        pipeline.add_step("Initial-Transcript-Mapping", "transcript_mapping")
        pipeline.add_step("Merge", "merge")
        pipeline.add_step("Pairalign", "pairaling")


        # # initial mapping
        # initial_mapping = pipeline.mapping_step(interleave([main_input, main_input2]), "initial")
        # # create denovo transcriptome
        # pipeline.create_denovo_transcriptome(initial_mapping)

        # ## run initial transcript mapping
        # pipeline.transcript_mapping_step(initial_mapping.clone(), "initial")




        # # merge
        # #pipeline.merge("step_1")

        # ## trim 20 mappings
        # #pipeline.mapping_step(gf(initial_split_mapping, unmapped), "trim_20", trim=(0, 20))
        # #pipeline.transcript_mapping_step(gf(initial_split_mapping, unmapped), "trim_20", trim=(0, 20))

        # ## trium 5 mappings
        # #pipeline.mapping_step(gf(initial_split_mapping, unmapped), "trim_5", trim=(5, 20))
        # #pipeline.transcript_mapping_step(gf(initial_split_mapping, unmapped), "trim_5", trim=(5, 20))

        # ## merge everything
        # merged = pipeline.merge("merged", same_content=True)

        # paired_mapping = pipeline.pair_align(merged, compress=args.gzip)

        # if not args.nosam:
        #     pipeline.create_bam(paired_mapping, sort=True)

        # pipeline.cleanup()

        timer.stop("Completed job in %s")
