#!/usr/bin/env python
import subprocess

from . import utils as utils
from . import multiplexer

"""Default filter im[plementations
"""


def multiplex(files, labels=False, filtered=False, sort=True):
    """Return a multiplexing generator. The generator takes
    a lit of files and multiplexes them. The input has to be in the same
    order, no checks are done. You can add conanical labels to the
    read ids and include reads that are filtered by the illumina filter.
    By default, the input files are sorted by name to identify read 1 and 2
    files.
    Note that the multiplexer splits the id by space and only adds the first
    part. The classical /1 and /2 are appended to the ids if not present.
    """
    return multiplexer.multiplex_generator(files, labels, filtered, sort)


def zcat(input):
    """Uses zcat to read teh input data. Input can e s single
    file or a list of files, but the string name. No other generator
    and no file handle is supported here."""

    inputs = []
    if isinstance(input, basestring):
        inputs.append(input)
    else:
        inputs = input

    input_fd = None
    for file in inputs:
        if file.endswith(".gz"):
            zcat = subprocess.Popen(["zcat", "-f", file],
                                    stdout=subprocess.PIPE)
            input_fd = zcat.stdout
        else:
            input_fd = open(file, 'r')

        for line in input_fd:
            yield line
        input_fd.close()


def gemoutput(input, exclude=None):
    """Generator that reads gem output and splits the
    string by tab."""

    if input is None:
        raise ValueError("No input specified")
    inputs = input
    if not isinstance(input, (list, tuple)):
        inputs = []
        inputs.append(input)

    for i in inputs:
        input_fd = None
        if isinstance(i, basestring):
            if i.endswith(".gz"):
                input_fd = zcat(i)
            else:
                input_fd = open(i, 'r')
        for line in input_fd:
            line = line.rstrip()
            e = line.split("\t")
            idx = 3
            if len(e) < 5:
                idx = 2
            if e[idx] in ["+", "-", "*"]:
                    e[idx] = "0"
            yield e
        try:
            input_fd.close()
        except:
            pass


def length(input, min=-1, max=65536):
    """Filter read length"""
    for gem in input:
        l = len(gem[1])
        if l >= min and l <= max:
            yield gem


def unmapped(input, exclude=None):
    """Generator that reads splitted gem output and returns the unmapped gem lines.
    If exclude is specified, also mappings with >= exclude missmatches are
    returned.
    """
    for gem in input:
        mismatches = -1
        if len(gem) >= 5:
            mismatches = utils.get_mismatches(gem[3])
        else:
            mismatches = utils.get_mismatches(gem[2])
        if exclude is None:
            if mismatches < 0:
                yield gem
        else:
            if mismatches < 0 or mismatches >= exclude:
                yield gem


def trim(input, right_trim=None, left_trim=None):
    """Reads gem output and trims the sequence and the qualities and returns
    gem output

    right_trim -- trim output read by read_trim characters from the right side
    left_trim -- trim output read by left_trim characters from the left side
    """

    if right_trim is None:
        right_trim = 0
    if left_trim is None:
        left_trim = 0

    for gem in input:
        if right_trim > 0:
            gem[1] = gem[1][left_trim:-right_trim]
            if len(gem) >= 5:
                gem[2] = gem[2][left_trim:-right_trim]  # trim qualities
        elif left_trim > 0:
            gem[1] = gem[1][left_trim:]
            if len(gem) >= 5:
                gem[2] = gem[2][left_trim:]
        yield gem


def to_input(input):
    """Reads gem output and converts it to gem gem input"""
    for gem in input:
        id = gem[0]
        if len(gem) >= 5:
            if not id.startswith("@"):
                id = "".join(["@", id])
            # empty string for final newline
            yield "\n".join([id, gem[1], "+", gem[2], ""])

        else:
            if not id.startswith(">"):
                id = "".join(">", id)
            yield "\n".join([id, gem[1], ""])


def to_gem(input):
    """Reads gem output and converts it back to a gem line"""
    for gem in input:
        yield "\t".join(gem)


def prepare_input(input):
    """Takes possible mapper input and converts it
    to an iterable.
    This is able to deal with string input, generators
    and file descriptors
    """
    if isinstance(input, (list, tuple)):
        ## deal with list or tupel input
        all = []
        for i in input:
            all.append(_prepare_single_input(i))
        return _iterate_iterators(all)
    else:
        ## handle single input
        return _prepare_single_input(input)


def _prepare_single_input(input):
    """Take a single input and convert it to something iterable"""

    if input is None:
        raise ValueError("No input specified")

    result = input
    if isinstance(input, basestring):
        ## handle files
        if input.endswith(".gz"):
            # wrap in zcat
            result = zcat(input)
            if input.endswith(".map.gz"):
                result = to_input(gemoutput(result))
        else:
            result = _iterate_file(input)
            if input.endswith(".map"):
                result = to_input(gemoutput(result))
        return result
    else:
        ## assume something iterable
        return result


def _iterate_iterators(input):
    """Generator function that iterates over a list of iterables
    and return their values
    """
    for i in input:
        for j in i:
            yield j


def _iterate_file(file):
    """Generator opens the given files and yields the lines"""
    fd = open(file, 'r')
    for line in fd:
        yield line
    fd.close()
