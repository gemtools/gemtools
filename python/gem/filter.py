#!/usr/bin/env python
"""Default filters for to filter reads and mappings"""
from __future__ import print_function
import gem
import gem.gemtools as gt
import gem.files as gf

import sys
import subprocess


def create_output_stream(output, compress=False, threads=1):
    """Create a writable output stream based on the given output.
    If the given output is None, this returns a subprocess.PIPE. If
    the output is a filename, a file descriptor for the file is opened
    and returned. Compression can be activated only for files but is
    ignored for stdout, stderr streams.

    The threads are used for the compression.

    The method returns a tuple of the target output stream and the modified
    output. The latter can be used to pass the output to gem._prepare_output().

    Parameter
    ---------
    output   - the target output, None, a file name or an already opend
               stream are allowed
    compress - if set to True and the output is not stdout/err, the output
               stream will be gzip compressed.
    threads  - if compression is enabled and installed compressors support
               multiple threads, they are passed on to the compressor
    """
    output_stream = subprocess.PIPE
    if output is not None:
        # its a file name
        if isinstance(output, basestring):
            output_stream = open(output, 'wb')
        elif isinstance(output, file):
            output_stream = file
            output = None  # disable for streams
        else:
            raise ValueError("Unable to open an output stream in %s" % output)

    # if the output goes to stderr or stdin, we reset
    # outout and disable compression
    if output_stream in [sys.stdout, sys.stderr]:
        compress = False
        output = None

    input_stream = output_stream

    p = None
    if compress and output is not None:
        p = subprocess.Popen(gem._compressor(threads=threads),
                             stdout=output_stream,
                             stdin=subprocess.PIPE)
        input_stream = p.stdin

    return input_stream, output, p


def only_split_maps(mappings, output=None, threads=1, paired=False,
                    compress=False):
    """Filter the mappings and remove all non-split maps"""
    args = ["--only-split-maps"]
    return run_filter(mappings, args, output=output, threads=threads,
                      paired=paired, compress=compress)


def rnaseq_filter(mappings, output=None, threads=1, paired=False,
                  min_intron=0,
                  min_block=0,
                  level=-1,
                  max_strata=0,
                  max_multi_maps=0,
                  gene_pairing=False,
                  junction_filter=False,
                  keep_unique=True,
                  annotation=None,
                  compress=False):
    """Filter rna seq mappings"""
    if (gene_pairing or junction_filter) and annotation is None:
        raise ValueError("You have to specify an annotation to perform "
                         "junction or gene-pairing filtering!")

    args = ["--min-intron-length", str(min_intron),
            "--min-block-length", str(min_block)]
    if keep_unique:
        args.append("-u")
    if gene_pairing:
        args.append("--reduce-by-gene-id")
    if junction_filter:
        args.append("--reduce-by-junctions")
    if level >= 0:
        args.extend(["--reduce-to-unique-strata", str(level)])
    if max_multi_maps >= 1:
        args.extend(["--reduce-to-max-maps", str(max_multi_maps)])
    if max_strata > 0:
        args.extend(["--max-strata", str(max_strata)])


    return run_filter(mappings, args, output=output, threads=threads,
                      paired=paired, compress=compress, rna_seq=True)


def run_filter(input, args, output=None,
               threads=1,
               paired=False,
               rna_seq=False,
               compress=False):
    """General purpose gt.filter execution"""
    quality = None
    # set qualitu if its an InputFile
    in_stream = gf.get_stream(input)
    if isinstance(input, gt.InputFile):
        quality = input.quality

    output_stream, output, compressor = create_output_stream(output,
                                                             compress=compress,
                                                             threads=threads)
    pa = [gem.executables['gt.filter'], '-t', str(threads)]
    if paired:
        pa.append('-p')
    if rna_seq:
        pa.append('--no-penalty-for-splitmaps')
    # add general argument
    pa.extend(args)
    # run the process
    process = subprocess.Popen(pa, stdin=in_stream, stdout=output_stream)
    return gem._prepare_output([process, compressor],
                               output=output, quality=quality)


def unmapped(reads, exclude=-1):
    """Yield only unmapped reads and reads
    that have only mappings with mismatches >= exclude
    if exclude is a nubmer between 0 and 1 we use it as a
    percentage of mismatches
    """
    return gt.unmapped(reads.__iter__(), exclude)


def length(input, min=-1, max=65536):
    """Filter reads by sequence length"""
    for read in input:
        l = read.length
        if l >= min and l <= max:
            yield read


class filter(object):
    """General clonable filter. Pass a read iterator and a filter function
    with arguments. This filter is clonable and can be used if you want
    to iterate ofer a base iterator more than once with the same filter.

    read_iterator   -- the read itereator
    filter_function -- the filter function
    args            -- filter function arguments
    kwargs          -- filter functions kw arguments
    """
    def __init__(self, read_iterator, filter_function, *args, **kwargs):
        """Clones the read iterator and fields output from the filter"""
        self.itereator = read_iterator.__iter__().clone()
        self.filter_function = filter_function
        self.args = args
        self.kwargs = kwargs
        self.fun = filter_function(self.itereator, *args, **kwargs)

    def clone(self):
        self.itereator = self.itereator.clone()
        self.fun = self.filter_function(self.itereator, *self.args, **self.kwargs)
        return self

    def __iter__(self):
        return self

    def next(self):
        return self.fun.next()


def interleave(reads, exclude=-1, threads=1):
    return gt.interleave(reads, threads=threads)


def cat(reads):
    return gt.cat(reads)


def trim(reads, left_trim=0, right_trim=0, min_length=10, append_label=False):
    """Trim reads from left and/pr right side.
    min_length is set to 10 by default and is the minimum
    length of the resulting read. If read length < min_length,
    the read is not trimmed and returned as is
    """
    return gt.trim(reads.__iter__(), left_trim, right_trim, min_length, append_label)


def _iterate_iterators(input):
    """Generator function that iterates over a list of iterables
    and return their values
    """
    for i in input:
        for j in i:
            yield j
