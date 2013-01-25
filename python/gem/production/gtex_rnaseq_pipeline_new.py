#!/usr/bin/env python
import sys
import os
import time
import argparse
import logging
import gem

from gem.filter import interleave, unmapped
from gem.filter import filter as gf


class MappingPipeline(object):

    def __init__(self, name=None, index=None,
        output_dir=None,
        annotation=None,
        threads=2,
        junctioncoverage=4,
        maxlength=100,
        transcript_index=None,
        transcript_keys=None,
        delta=1
        ):
        self.name = name
        self.index = index
        self.output_dir = output_dir
        self.annotation = annotation
        self.threads = threads
        self.junctioncoverage = junctioncoverage
        self.maxlength = maxlength
        self.transcript_index = transcript_index
        self.transcript_keys = transcript_keys
        self.denovo_index = None
        self.denovo_keys = None
        self.delta = delta
        self.mappings = []
        self.files = []

        if self.output_dir is None:
            self.output_dir = os.getcwd()

        ## check if output is specified and create folder in case it does not exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if self.transcript_keys is None and self.transcript_index is not None:
            self.transcript_keys = self.transcript_index[:-4] + ".keys"

    def cleanup(self):
        """Delete all remaining temporary and intermediate files
        """
        for f in self.files:
            if os.path.exists(f):
                logging.info("Removing intermediate file : %s" % (f))
                os.remove(f)

    def create_file_name(self, suffix, file_suffix="map", final=False):
        """Create a result file name"""
        file = ""
        if suffix is not None and len(suffix) > 0:
            file = "%s/%s_%s.%s" % (self.output_dir, self.name, suffix, file_suffix)
        else:
            file = "%s/%s.%s" % (self.output_dir, self.name, file_suffix)
        if not final:
            self.files.append(file)
        return file

    def _guess_input_name(self, input):
        try:
            return input.filename
        except Exception:
            return "STREAM"

    def merge(self, suffix):
        """Merge current set of mappings and delete last ones"""
        out = self.create_file_name(suffix)
        logging.info("Merging %d mappings into %s" % (len(self.mappings), out))
        merged = gem.merger(
            self.mappings[0].clone(),
            [m.clone() for m in self.mappings[1:]]
        ).merge(out)

        for m in self.mappings:
            logging.info("Removing temporary mapping %s ", m.filename)
            os.remove(m.filename)
        self.mappgins = []
        self.mappings.append(merged)
        return merged

    def mapping_step(self, input, suffix, trim=None):
        """Single mapping step where the input is passed on
        as is to the mapper. The output name is created based on
        the dataset name in the parameters the suffix.

        The function returns an opened ReadIterator on the resulting
        mapping

        input -- mapper input
        suffix -- output name suffix
        trim -- optional trimming options, i.e. '0,20'
        """
        mapping_out = self.create_file_name(suffix)
        input_name = self._guess_input_name(input)
        logging.debug("Mapping from %s to %s" % (input_name, mapping_out))
        mapping = gem.mapper(input,
                self.index,
                mapping_out,
                mismatches=0.06,
                delta=self.delta,
                trim=trim,
                threads=self.threads)
        self.mappings.append(mapping)
        return mapping

    def transcript_mapping_step(self, input, suffix, trim=None):
        """Single transcript mapping step where the input is passed on
        as is to the mapper. The output name is created based on
        the dataset name in the parameters the suffix.

        The function returns an opened ReadIterator on the resulting
        mapping

        input -- mapper input
        suffix -- output name suffix
        trim -- optional trimming options, i.e. '0,20'
        """
        mapping_out = self.create_file_name(suffix + "_transcript")
        input_name = self._guess_input_name(input)
        logging.debug("Mapping transcripts from %s to %s" % (input_name, mapping_out))

        mapping = gem.transcript_mapper(
                                input,
                                [self.transcript_index, self.denovo_index],
                                [self.transcript_keys, self.denovo_keys],
                                mapping_out,
                                mismatches=0.06,
                                trim=trim,
                                delta=self.delta,
                                threads=self.threads
                                )
        self.mappings.append(mapping)
        return mapping

    def gtf_junctions(self):
        """check if there is a .junctions file for the given annotation, if not,
        create it. Returns a tuple of the set of junctions and the output file.
        """
        logging.info("Loading junctions from %s" % (self.annotation))
        junctions = set(gem.junctions.from_gtf(self.annotation))
        logging.info("%d Junctions from GTF" % (len(junctions)))
        out = self.create_file_name("gtf", file_suffix="junctions")
        gem.junctions.write_junctions(junctions, out, self.index)
        return (junctions, out)

    def create_denovo_transcriptome(self, input):
        """Squeeze input throw the split mapper to extract denovo junctions
        and create a transcriptome out of the junctions"""
        (junctions, junctions_gtf_out) = self.gtf_junctions()
        ## get de-novo junctions
        logging.info("Getting de-novo junctions")
        denovo_out = self.create_file_name("denovo")
        junctions_out = self.create_file_name("all", file_suffix="junctions")

        ## add .fa and .keys to list of files to delete
        self.files.append(junctions_out + ".fa")
        self.files.append(junctions_out + ".keys")

        index_denovo_out = self.create_file_name("denovo_transcripts", file_suffix="gem")
        (denovo_mapping, junctions) = gem.extract_junctions(
            input,
            self.index,
            denovo_out,
            mismatches=0.04,
            threads=self.threads,
            strata_after_first=0,
            coverage=self.junctioncoverage,
            merge_with=junctions)
        logging.info("Total Junctions %d" % (len(junctions)))

        gem.junctions.write_junctions(gem.junctions.filter_by_distance(junctions, 500000), junctions_out, self.index)

        logging.info("Computing denovo transcriptome")
        (denovo_transcriptome, denovo_keys) = gem.compute_transcriptome(self.maxlength, self.index, junctions_out, junctions_gtf_out)
        logging.info("Indexing denovo transcriptome")

        idx = gem.index(denovo_transcriptome, index_denovo_out, threads=self.threads)
        # add indexer log file to list of files
        self.files.append(idx[:-4] + ".log")
        self.denovo_keys = denovo_keys
        self.denovo_index = idx
        self.mappings.append(denovo_mapping)
        return denovo_mapping

    def pair_align(self, input, compress=False):
            logging.info("Running pair aligner")
            paired_out = self.create_file_name("", final=True)
            paired_mapping = gem.pairalign(input, self.index, None, max_insert_size=100000, threads=self.threads)
            scored = gem.score(paired_mapping, self.index, paired_out, threads=self.threads)
            if compress:
                logging.info("Compressing final mapping")
                gem.utils.gzip(paired_out, threads=self.threads)

            return scored

    def create_bam(self, input, sort=True):
        logging.info("Converting to sam/bam")
        sam_out = self.create_file_name("", file_suffix="bam", final=True)
        sam = gem.gem2sam(input, self.index, threads=max(1, int(self.threads / 2)))
        gem.sam2bam(sam, sam_out, sorted=sort)


def main():
    parser = argparse.ArgumentParser()

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
    parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
    parser.add_argument('--extend-name', dest="extendname", action="store_true", default=False, help="Extend the name by prefixing it with the parameter combination")

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

