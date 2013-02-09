#!/usr/bin/env python
"""Default filters for to filter reads and mappings"""
from __future__ import print_function
import gem.gemtools as gt


def unmapped(reads, exclude=-1):
    """Yield only unmapped reads and reads
    that have only mappings with mismatches >= exclude
    if exclude is a nubmer between 0 and 1 we use it as a percentage of mismatches
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


def interleave(reads, exclude=-1):
    return gt.interleave(reads)


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
