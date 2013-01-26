#!/usr/bin/env python
import sys
import os
import time
import argparse
import logging
import gem

from gem.filter import interleave, unmapped
from gem.filter import filter as gf
from gem.pipeline import MappingPipeline


class RnaPipeline(gem.command.Command):
    def register(parser):
        ## required parameters
        parser.add_argument('-f', '--file', dest="file", help='The fastq input file', required=True)
        parser.add_argument('-i', '--index', dest="index", nargs=2, help='First genome index sencond transcript index', required=True)
        parser.add_argument('-k', '--keys', dest="transcript_keys", help='transcriptome keys file', required=True)
        parser.add_argument('-a', '--annotation-file', dest="annotation", help='The annotation file', required=True)
        parser.add_argument('-m', '--max-read-length', dest="maxlength", help='The maximum read length', required=True)

        ## optional parameters
        parser.add_argument('-o', '--output-dir', dest="output", help='The output folder. If not specified the current working directory is used.')
        parser.add_argument('--unmapped-threshold', dest="unmappedthreshold",
            help='Number of mismatches for a read to be treated as unmapped. If 0 < i < 1, this is interpreted as percentage with respect to the read length.', default=2)
        parser.add_argument('--junction-coverage', dest="junctioncoverage",
            help='A denovo junction must be covered by > coverage reads to be taken into account, 0 to disable', default=4)
        parser.add_argument('-s', '--strata-after-best', dest="delta",
            help='Number of strata that are examined after the best one', default=1)
        parser.add_argument('-g', '--no-gzip', dest="gzip", action="store_false", default=True, help="Do not compress result mapping")
        parser.add_argument('--keep-temp', dest="rmtemp", action="store_false", default=True, help="Keep temporary files")
        parser.add_argument('-t', '--threads', dest="threads", default=8, type=int, help="Number of threads to use")
        parser.add_argument('--no-sam', dest="nosam", action="store_true", default=False, help="Do not create sam/bam file")
        parser.add_argument('--extend-name', dest="extendname", action="store_true", default=False, help="Extend the name by prefixing it with the parameter combination")

    def run(args):
        ## parsing command line arguments
        args = parser.parse_args()
        if args.loglevel is not None:
            gem.loglevel(args.loglevel)

        input_file = os.path.abspath(args.file)
        input_file2 = input_file.replace("0.f", "1.f")
        name = os.path.splitext(os.path.basename(input_file))[0].replace(".0", "")
        if args.extendname:
            name = "%s_%d_%.2f_%d" % (name, args.delta, args.unmappedthreshold, args.junctioncoverage)

        pipeline = MappingPipeline(
            name=name,
            index=args.index[0],
            output_dir=args.output,
            annotation=args.annotation,
            threads=int(args.threads),
            junctioncoverage=int(args.junctioncoverage),
            maxlength=int(args.maxlength),
            transcript_index=args.index[1],
            transcript_keys=args.transcript_keys,
            delta=int(args.delta)
        )

        start_time = time.time()

        main_input = gem.files.open(input_file)
        main_input2 = gem.files.open(input_file2)

        # initial mapping
        initial_mapping = pipeline.mapping_step(interleave([main_input, main_input2], add_id=False), "initial")
        # create denovo transcriptome
        pipeline.create_denovo_transcriptome(initial_mapping)

        ## run initial transcript mapping
        initial_split_mapping = pipeline.transcript_mapping_step(initial_mapping.clone(), "initial_transcripts")

        # merge
        #pipeline.merge("step_1")

        ## trim 20 mappings
        pipeline.mapping_step(gf(initial_split_mapping, unmapped), "trim_20", trim=(0, 20))
        pipeline.transcript_mapping_step(gf(initial_split_mapping, unmapped), "trim_20", trim=(0, 20))

        ## trium 5 mappings
        pipeline.mapping_step(gf(initial_split_mapping, unmapped), "trim_5", trim=(5, 20))
        pipeline.transcript_mapping_step(gf(initial_split_mapping, unmapped), "trim_5", trim=(5, 20))

        ## merge everything
        merged = pipeline.merge("merged")

        paired_mapping = pipeline.pair_align(merged, compress=True)

        pipeline.create_bam(paired_mapping, sort=True)

        pipeline.cleanup()

        end_time = (time.time() - start_time) / 60

        print "Completed job in %0.2f mins" % end_time


if __name__ == "__main__":
    main()

