#!/usr/bin/env python
import sys
import os
import time
import argparse

import gem

def main():
    parser = argparse.ArgumentParser()

    ## required parameters
    parser.add_argument('-f', '--file', dest="file", help='The fastq input file', required=True)
    parser.add_argument('-i', '--index-file', dest="index", help='The index file', required=True)
    parser.add_argument('-a', '--annotation-file', dest="annotation", help='The annotation file', required=True)

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
    parser.add_argument('-k', '--keep-temp', dest="rmtemp", action="store_false", default=True, help="Keep temporary files")
    parser.add_argument('-t', '--threads', dest="threads", default=8, type=int, help="Number of threads to use")
    parser.add_argument('--no-sam', dest="nosam", action="store_true", default=False, help="Do not create sam/bam file")
    parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
    parser.add_argument('--extend-name', dest="extendname", action="store_true", default=False, help="Extend the name by prefixing it with the parameter combination")

    ## parsing command line arguments
    args = parser.parse_args()
    annotation = os.path.abspath(args.annotation)
    index = os.path.abspath(args.index)
    THREADS = args.threads
    rmFiles = args.rmtemp
    gzipOut = args.gzip
    exclusive = args.exclusive
    nosam = args.nosam
    delta = int(args.delta)
    unmappedthreshold = float(args.unmappedthreshold)
    junctioncoverage = int(args.junctioncoverage)

    if args.loglevel is not None:
        gem.loglevel(args.loglevel)

    input_file = os.path.abspath(args.file)
    name = os.path.splitext(os.path.basename(input_file))[0].replace(".0","")

    print name

    input_file2 = input_file.replace("0.f","1.f")

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
    junctions_out = "%s/%s.junctions" % (output_dir, name)
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

        ## get de-novo junctions
        print "Getting de-novo junctions"
        (denovo_mapping, junctions) = gem.extract_junctions(
            gem.filter.unmapped(initial_mapping, unmappedthreshold),
            index,
            denovo_out,
            mismatches=0.04,
            threads=THREADS,
            coverage=junctioncoverage,
            merge_with=junctions)
        print "Total Junctions %d" % (len(junctions))
        ## merge de-novo and known junctions
        print "Writing junctions file"
        gem.junctions.write_junctions(gem.junctions.filter_by_distance(junctions, 500000), junctions_out, index)
    else:
        denovo_mapping = gem.files.open(denovo_out)

    ## create initial split-map
    if not os.path.exists(initial_split_out):
        print "Running initial split-map"
        initial_split_mapping = gem.splitmapper(
            gem.filter.unmapped(denovo_mapping),
            index,
            initial_split_out,
            junctions_file=junctions_out,
            threads=THREADS,
            mismatches=0.06)
    else:
        initial_split_mapping = gem.files.open(initial_split_out)

    ## create trim-20 mapping
    if not os.path.exists(trim_20_out):
        print "Running trim 20 mapping"
        trim_20_mapping = gem.mapper(
            gem.filter.unmapped(initial_split_mapping),
            index,
            trim_20_out,
            mismatches=0.07,
            delta=delta,
            trim=(0, 20),
            threads=THREADS)
    else:
        trim_20_mapping = gem.files.open(trim_20_out)

    ## create trim-20 split-map
    if not os.path.exists(trim_20_split_out):
        print "Running trim 20 split-map"
        trim_20_split_mapping = gem.splitmapper(
            gem.filter.unmapped(trim_20_mapping),
            index,
            trim_20_split_out,
            junctions_file=junctions_out,
            trim=(0, 20),
            threads=THREADS,
            mismatches=0.06)
    else:
        trim_20_split_mapping = gem.files.open(trim_20_split_out)

    if not os.path.exists(trim_5_out):
        print "Running trim 5 mapping"
        trim_5_mapping = gem.mapper(
            gem.filter.unmapped(trim_20_split_mapping),
            index,
            trim_5_out,
            mismatches=0.07,
            delta=delta,
            trim=(5, 20),
            threads=THREADS)
    else:
        trim_5_mapping = gem.files.open(trim_5_out)

    ## create trim-20 split-map
    if not os.path.exists(trim_5_split_out):
        print "Running trim 5 split-map"
        trim_5_split_mapping = gem.splitmapper(
            gem.filter.unmapped(trim_5_mapping),
            index,
            trim_5_split_out,
            junctions_file=junctions_out,
            trim=(5, 20),
            threads=THREADS,
            mismatches=0.06)
    else:
        trim_5_split_mapping = gem.files.open(trim_5_split_out)

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
        os.remove(initial_out)
        os.remove(initial_split_out)
        os.remove(denovo_out)
        os.remove(junctions_out)
        os.remove(trim_20_out)
        os.remove(trim_20_split_out)
        os.remove(trim_5_out)
        os.remove(trim_5_split_out)

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
