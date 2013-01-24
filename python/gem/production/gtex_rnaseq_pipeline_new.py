#!/usr/bin/env python
import sys
import os
import time
import argparse

import gem


class MappingParameters(object):
    def __init__(self, index, name, output_dir, delta, threads):
        self.index = index
        self.name = name
        self.delta = delta
        self.output_dir = output_dir
        self.threads = threads
        self.transcript_index = None
        self.transcript_index_keys = None
        self.transcript_index_denovo = None
        self.transcript_index_denovo_keys = None


def map_step(input, suffix, params, trim=None, do_transcript_mapping=True):
    mapping = None
    transcript_mapping = None
    mapping_out = "%s/%s_%s.map" % (params.output_dir, suffix, params.name)
    transcript_out = "%s/%s_%s_transcript.map" % (params.output_dir, suffix, params.name)

    ## create initial mapping
    if not os.path.exists(mapping_out):
        mapping = gem.mapper(input,
                params.index,
                mapping_out,
                mismatches=0.06,
                delta=params.delta,
                trim=trim,
                threads=params.threads)
    else:
        mapping = gem.files.open(mapping_out)

    if do_transcript_mapping:
        if not os.path.exists(transcript_out):
            print "Running transcriptome mapping"
            transcript_mapping = gem.transcript_mapper(
                                input.clone(),
                                [params.transcript_index, params.transcript_index_denovo],
                                [params.transcript_index_keys, params.transcript_index_denovo_keys],
                                transcript_out,
                                mismatches=0.06,
                                trim=trim,
                                delta=params.delta,
                                threads=params.threads)
        else:
            transcript_mapping = gem.files.open(transcript_out)
        return (mapping, transcript_mapping)
    else:
        return mapping


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
    parser.add_argument('-e', '--exclusive', dest="exclusive", action="store_true", default=False, help="Exclusive merge")
    parser.add_argument('-g', '--no-gzip', dest="gzip", action="store_false", default=True, help="Do not compress result mapping")
    parser.add_argument('--keep-temp', dest="rmtemp", action="store_false", default=True, help="Keep temporary files")
    parser.add_argument('-t', '--threads', dest="threads", default=8, type=int, help="Number of threads to use")
    parser.add_argument('--no-sam', dest="nosam", action="store_true", default=False, help="Do not create sam/bam file")
    parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
    parser.add_argument('--extend-name', dest="extendname", action="store_true", default=False, help="Extend the name by prefixing it with the parameter combination")

    ## parsing command line arguments
    args = parser.parse_args()
    annotation = os.path.abspath(args.annotation)
    index = os.path.abspath(args.index[0])
    transcript_index = os.path.abspath(args.index[1])

    THREADS = args.threads
    rmFiles = args.rmtemp
    gzipOut = args.gzip
    exclusive = args.exclusive
    maxlength = args.maxlength
    nosam = args.nosam
    delta = int(args.delta)
    unmappedthreshold = float(args.unmappedthreshold)
    junctioncoverage = int(args.junctioncoverage)
    transcript_keys = args.transcript_keys
    if args.loglevel is not None:
        gem.loglevel(args.loglevel)

    input_file = os.path.abspath(args.file)
    name = os.path.splitext(os.path.basename(input_file))[0].replace(".0", "")

    input_file2 = input_file.replace("0.f", "1.f")

    if args.extendname:
        name = "%s_%d_%.2f_%d" % (name, delta, unmappedthreshold, junctioncoverage)

    output_dir = os.getcwd()
    ## check if output is specified and create folder in case it does not exists
    if args.output is not None:
        output_dir = os.path.abspath(args.output)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    print "Loading %s with data set name %s" % (input_file, name)
    print ""
    print "Index file", index
    print "Annotation file", annotation
    print "Output folder", output_dir
    print ""
    print "****** Parameters ******"
    print "Keep Temp", (not rmFiles)
    print "GZIP Map", gzipOut
    print "Threads", THREADS
    print "Exclusive Merge", exclusive
    print "Prefix name", name
    print "Unmapped Threshold", unmappedthreshold
    print "Junction Coverage", junctioncoverage
    print "Strata after best", delta
    print "************************"
    print ""

    initial_out = "%s/%s_initial.map" % (output_dir, name)
    denovo_out = "%s/%s_denovo.map" % (output_dir, name)
    index_denovo_out = "%s/%s_denovo.gem" % (output_dir, name)
    junctions_out = "%s/%s.junctions" % (output_dir, name)
    junctions_gtf_out = "%s/%s.gtf.junctions" % (output_dir, name)
    initial_split_out = "%s/%s_initial_split.map" % (output_dir, name)
    trim_20_out = "%s/%s_trim_20.map" % (output_dir, name)
    trim_20_split_out = "%s/%s_trim_20_split.map" % (output_dir, name)
    trim_5_out = "%s/%s_trim_5.map" % (output_dir, name)
    trim_5_split_out = "%s/%s_trim_5_split.map" % (output_dir, name)
    final_out = "%s/%s_final.map" % (output_dir, name)
    scored_out = "%s/%s.map" % (output_dir, name)
    paired_out = "%s/%s_paired.map" % (output_dir, name)
    sam_out = "%s/%s.bam" % (output_dir, name)

    start_time = time.time()

    ## create initial mapping
    if not os.path.exists(initial_out):
        main_input = gem.files.open(input_file)
        main_input2 = gem.files.open(input_file2)
        print "Running initial mapping"
        initial_mapping = gem.mapper(gem.filter.interleave([main_input, main_input2], add_id=False), index, initial_out,
                mismatches=0.06,
                delta=delta,
                threads=THREADS)
    else:
        initial_mapping = gem.files.open(initial_out)

    ## get junctions from the annotation

    if not os.path.exists(denovo_out):
        print "Loading GTF junctions from %s" % annotation
        junctions = set(gem.junctions.from_gtf(annotation))
        print "%d Junctions from GTF" % (len(junctions))
        gem.junctions.write_junctions(junctions, junctions_gtf_out, index)

        ## get de-novo junctions
        print "Getting de-novo junctions"
        (denovo_mapping, junctions) = gem.extract_junctions(
            gem.filter.unmapped(initial_mapping, unmappedthreshold),
            index,
            denovo_out,
            mismatches=0.04,
            threads=THREADS,
            strata_after_first=0,
            coverage=junctioncoverage,
            merge_with=junctions)
        print "Total Junctions %d" % (len(junctions))
        ## merge de-novo and known junctions
        print "Writing junctions file"
        gem.junctions.write_junctions(gem.junctions.filter_by_distance(junctions, 500000), junctions_out, index)
    else:
        denovo_mapping = gem.files.open(denovo_out)

    print "Computing denovo transcriptome"
    (denovo_transcriptome, denovo_keys) = gem.compute_transcriptome(maxlength, index, junctions_out, junctions_gtf_out)
    print "Indexing denovo transcriptome"
    gem.index(denovo_transcriptome, index_denovo_out, threads=THREADS)

    ## map to transcriptome
    ## create initial split-map
    if not os.path.exists(initial_split_out):
        print "Running transcriptome mapping"
        initial_split_mapping = gem.transcript_mapper(
                            gem.filter.filter(denovo_mapping, gem.filter.unmapped),
                            [transcript_index, index_denovo_out],
                            [transcript_keys, denovo_keys],
                            initial_split_out,
                            mismatches=0.06,
                            delta=delta,
                            threads=THREADS)
    else:
        initial_split_mapping = gem.files.open(initial_split_out)

    params = MappingParameters(index, name, output_dir, delta, THREADS)
    params.transcript_index = transcript_index
    params.transcript_index_keys = transcript_keys
    params.transcript_index_denovo = index_denovo_out
    params.transcript_index_denovo_keys = denovo_keys

    (trim_20_mapping, trim_20_split_mapping) = map_step(initial_split_mapping, "trim_20", params, trim=(0, 20))
    (trim_5_mapping, trim_5_split_mapping) = map_step(initial_split_mapping, "trim_5", params, trim=(5, 20))

    ## merge
    print "Merging results"
    merged = gem.merger(
        initial_mapping.clone(),
        [denovo_mapping.clone(),
         initial_split_mapping.clone(),
         trim_20_mapping.clone(),
         trim_20_split_mapping.clone(),
         trim_5_mapping.clone(),
         trim_5_split_mapping.clone()],
         exclusive=exclusive
    ).merge(final_out)

    ## remove files
    if rmFiles:
        print "Removing files"
        os.remove(initial_mapping.filename)
        os.remove(initial_split_mapping.filename)
        os.remove(denovo_mapping.filename)
        os.remove(junctions_out)
        os.remove(trim_20_mapping.filename)
        os.remove(trim_20_split_mapping.filename)
        os.remove(trim_5_mapping.filename)
        os.remove(trim_5_split_mapping.filename)

    ## pair align
    print "Running pair aligner"
    #validation = gem.validate(merged, index, threads=8)
    paired_mapping = gem.pairalign(merged, index, paired_out, max_insert_size=100000, threads=THREADS)
    if rmFiles:
        os.remove(final_out)

    ## score the alignments
    print "Validating and scoring alignment"
    scored = gem.score(paired_mapping, index, scored_out, threads=THREADS)
    if rmFiles:
        os.remove(paired_out)

    ## create a sorted bamfile
    if not nosam:
        print "Converting to sam"
        sam = gem.gem2sam(scored, index, threads=max(1, int(THREADS / 2)))
        bam = gem.sam2bam(sam, sam_out, sorted=True)

    ## create a gzipped map file
    if gzipOut:
        print "Compressing map"
        gem.utils.gzip(scored_out, threads=THREADS)

    end_time = (time.time() - start_time) / 60
    print "Done!"
    print "Completed job in %0.2f mins" % end_time


if __name__ == "__main__":
    main()

