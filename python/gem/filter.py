#!/usr/bin/env python
"""Default filters for to filter reads and mappings"""
from __future__ import print_function
import sys
import re
import gem
import math

class BasicStats(object):
    """Basic statistic computation for the reads that
    pass through this filter.
    Statistics are gathered by calling to basic_stats, and they
    are printed when a print_stats call is done at the end.
    """

    def __init__(self, output=None):
        """Initialize the statistics fields"""
        self.mapped     = 0
        self.unmapped   = 0
        self.outputFile = output

    def basic_stats(self, reads):
        """Yield all the reads and count how many
        are mapped and how many are unmapped.
        In short, compute the desired statistics.
        """
        for read in reads:
            if read.mappings is None or len(read.mappings) == 0:
                # This check gets rid of reads coming from not mapped input line fasta files
                self.unmapped += 1
            else:
                # Here we have an alignment try, which can work or not
                if read.min_mismatches() == -1:
                    self.unmapped += 1
                else:
                    self.mapped += 1

            yield read

    def print_stats(self):
        """Compute the final statistics and write them to a file"""
        # Compute basic statistics
        total        = self.mapped   + self.unmapped
        pct_mapped   = self.mapped   * 100.0 / total
        pct_unmapped = self.unmapped * 100.0 / total

        # Open output file
        if self.outputFile == None or self.outputFile == sys.stdout:
            output = None
        else:
            output = open(self.outputFile, "w")

        # Write statistics
        print("ALL:      {0:>9d}".format(total), file=output)
        print("Mapped:   {0:>9d} ({1:5.2f}%)".format(self.mapped, pct_mapped), file=output)
        print("Unmapped: {0:>9d} ({1:5.2f}%)".format(self.unmapped, pct_unmapped), file=output)

        # Close file, if needed
        if self.outputFile != None and self.outputFile != sys.stdout:
            output.close()


def unmapped(reads, exclude=-1):
    """Yield only unmapped reads and reads
    that have only mappings with mismatches >= exclude
    if exclude is a nubmer between 0 and 1 we use it as a percentage of mismatches
    """
    for read in reads:
        if read.mappings is None or len(read.mappings) == 0:
            yield read
        else:
            mis = read.min_mismatches()
            excl = exclude
            if excl > 0 and excl < 1:
                excl = math.floor(read.length() * exclude)
            if mis < 0:
                yield read
            elif excl > 0 and mis >= excl:
                yield read


def length(input, min=-1, max=65536):
    """Filter reads by sequence length"""
    for read in input:
        l = len(read.sequence)
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
        self.itereator = read_iterator.clone()
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


class interleave(object):
    """Interleaving iterator that takes a sequence
    of Reads and interleaves them. By default the read
    ids are checked for /1 /2 endings and the endings
    are appended if not found.
    """

    def __init__(self, streams, add_id=True):
        self.streams = streams
        self.idx = 0
        self.counter = 0
        self.len = len(streams)
        self.add_id = add_id

    def __iter__(self):
        return self

    def next(self):
        while self.counter < self.len:
            if self.idx >= self.len:
                self.idx = 0

            ret = self.streams[self.idx].next()
            self.idx += 1
            if ret is not None:
                if self.add_id:
                    id = ret.id.split(" ")[0]
                    if not re.search(".*/\d+$", id):
                        ret.id = "%s/%d" % (id, self.idx)
                self.counter = 0
                return ret
            else:
                self.counter += 1
        raise StopIteration()


class cat(object):
    """Concatenate the content of a list of
    read iterators
    """

    def __init__(self, streams):
        self.streams = streams
        self.idx = 0
        self.len = len(streams)

    def __iter__(self):
        return self

    def next(self):
        while self.idx < self.len:
            try:
                ret = self.streams[self.idx].next()
                if ret is not None:
                    return ret
                else:
                    self.idx += 1
            except StopIteration:
                self.idx += 1
        raise StopIteration()


def trim(reads, left_trim=0, right_trim=0, min_length=10, append_label=False):
    """Trim reads from left and/pr right side.
    min_length is set to 10 by default and is the minimum
    length of the resulting read. If read length < min_length,
    the read is not trimmed and returned as is
    """
    seq = None
    quals = None
    for read in reads:
        if right_trim > 0:
            seq = read.sequence[left_trim:-right_trim]
            if read.qualities is not None:
                quals = read.qualities[left_trim:-right_trim]
        elif left_trim > 0:
            seq = read.sequence[left_trim:]
            if  read.qualities is not None:
                quals = read.qualities[left_trim:]
        if len(seq) > min_length:
            if append_label:
                left_qual = ""
                right_qual = ""
                left_seq = read.sequence[:left_trim]
                right_seq = read.sequence[-right_trim:]
                if read.qualities:
                    left_qual = read.qualities[:left_trim]
                    right_qual = read.qualities[-right_trim:]
                read.id = "%s B T %s %s %s %s" % (read.id, left_seq, right_seq, left_qual, right_qual)

            read.sequence = seq
            read.qualities = quals
        yield read


def keep_short_indel(read):
    """
    Takes a read and removes all splitmaps that are not short indels

    @param read: the read
    @type read: Read
    """
    (sums, maps) = read.get_maps()
    accs = []
    for s in sums:
        accs.append(0)
    map_idx = 0
    matches = []
    for i, sum in enumerate(sums):
        if sum == gem._max_mappings:
            sum = 0
        for j in range(sum):
            match = re.search(">([0-9]+)\*", maps[map_idx])
            if match:
                hit = int(match.groups()[0])
                if 0 < hit and hit < 4:
                    accs[i] += 1
                    matches.append(maps[map_idx])
            map_idx += 1
    if len(matches) > 0:
        read.mappings = ",".join(matches)
        read.summary = ":".join(str(x) for x in accs)
    else:
        read.mappings = "-"
        read.summary = "0"


def keep_short_indels(geminput):
    """Filter gem output and remove all split-mappings
    (the mapping, not the read!) where the
    splice site is >= 4 bases. Only short indel
    split-maps are kept.
    """
    for read in geminput:
        keep_short_indel(read)
        yield read


def _iterate_iterators(input):
    """Generator function that iterates over a list of iterables
    and return their values
    """
    for i in input:
        for j in i:
            yield j
