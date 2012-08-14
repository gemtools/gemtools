#!/usr/bin/env python
"""Default filters for to filter reads and mappings"""
import re


def sam_2_map(reads):
    """Convert sam reads to simple maps"""
    for read in reads:
        read.sam_2_gem()
        yield read


def unmapped(reads, exclude=-1):
    """Yield only unmapped reads and reads
    that have only mappings with mismatches >= exclude
    """
    for read in reads:
        if read.mappings is None or len(read.mappings) == 0:
            yield read
        else:
            mis = read.min_mismatches()
            if mis < 0:
                yield read
            elif exclude > 0 and mis >= exclude:
                yield read


def length(input, min=-1, max=65536):
    """Filter reads by sequence length"""
    for read in input:
        l = len(read.sequence)
        if l >= min and l <= max:
            yield read


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


def trim(reads, left_trim=0, right_trim=0, min_length=10):
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
