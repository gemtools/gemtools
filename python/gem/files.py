#!/usr/bin/env python
"""
gem.files handles opening files and streams
"""
import subprocess
import __builtin__
import gem.gemtools as gt
from utils import which
import shutil
import os

__author__ = 'Thasso Griebel'
__zcat_path = None

__open_iterators = []

delete_on_exit = []

def open(input, quality=None):
    """
    Open the given file and return on iterator
    over Reads.

    The file parameter has to be a string.

    @param input: string of the file name or an open stream
    @type input: string or stream
    @param type: the type of the input file or none for auto-detection
    @type type: string
    @param process: optional process associated with the input
    @type process: Process
    @param remove_after_iteration: if set to True, the input file is removed after a completed iteration of the content
    @type remove_after_iteration: boolean
    @param quality: the gem quality parameter
    @type quality: string
    """
    is_string = isinstance(input, basestring)
    stream = None
    if is_string:
        type = _guess_type(input)
        if type == "bam":
            stream = open_bam(input)
    else:
        stream = input

    it = None
    if stream is not None:
        it = gt.InputFile(stream, quality=quality)
        __open_iterators.append(it)
    else:
        it = gt.InputFile(input, quality=quality)
    return it


def open_bam(input):
    if isinstance(input, basestring):
       return subprocess.Popen(["samtools", "view", "-h", input],
            stdout=subprocess.PIPE,
            stderr=__builtin__.open("/dev/null", 'w'),
            close_fds=True).stdout
    return subprocess.Popen(["samtools", "view", "-h", "-"],
            stdout=subprocess.PIPE,
            stdin=input,
            stderr=__builtin__.open("/dev/null", 'w'),
            close_fds=True).stdout


def open_bam_stream(input):
    return subprocess.Popen(["samtools", "view", "-h", "-"],
        stdout=subprocess.PIPE,
        stdin=input,
        #stderr=__builtin__.open("/dev/null", 'w'),
        close_fds=True).stdout


def _cleanup():
    for r in __open_iterators:
        r.close()
    if delete_on_exit is not None:
        for f in delete_on_exit:
            if os.path.exists(f):
                if os.path.isfile(f):
                    os.remove(f)
                else:
                    shutil.rmtree(f)


def get_stream(input):
    """Checks the input and returns a raw stream. The input can be either
    a path to a file, a gemtools.IntputFile, or an already opend stream
    """
    if isinstance(input, basestring):
        return open_file(input)
    if isinstance(input, gt.InputFile):
        return input.raw_stream()
    if hasattr(input, "reads"):
        return input
    raise ValueError("Unable to open a stream on the given input")


def open_file(file):
    """
    Open the given file and return a stream
    of the file content. The method checks if the file
    ends with .gz and open a gzip stream in that case.

    The file parameter has to be a string.

    @param file: string of the file name
    @type file: string
    """
    if file.endswith(".gz"):
        return open_gzip(file)
    else:
        return __builtin__.open(file, 'r')


def _guess_type(name):
    """
    guess the type based on the given file name
    and returns one of teh following:

    fasta, fastq, map

    or None if no type could be detected

    @param name: the name of the file
    @type name: string
    @return: one of fasta, fastq, map or None
    @rtype: string
    """
    name = name.upper()
    if name.endswith(".GZ"):
        name = name[:-3]
    if name.endswith(".FASTA") or name.endswith("FA"):
        return "fasta"
    elif name.endswith(".FASTQ") or name.endswith("FQ"):  # fixes issue #5
        return "fastq"
    elif name.endswith(".MAP"):
        return "map"
    elif name.endswith(".SAM"):
        return "sam"
    elif name.endswith(".BAM"):
        return "bam"
    return None


def open_gzip(file_name):
    """Uses a separate zcat process to open a gzip
    stream.

    @param file_name: the name of the input file
    @return stream: gzip stream
    @rtype stream
    """
    if not isinstance(file_name, basestring):
        raise ValueError("The provided file name is not a string : %s" % (file_name))
    parameter = [__zcat(), file_name]
    zcat = subprocess.Popen(parameter, stdout=subprocess.PIPE)
    return zcat.stdout


def __zcat():
    """
    get the path to the zcat/gzcat executable
    or raise an exception
    """
    global __zcat_path
    if __zcat_path is not None:
        return __zcat_path

    __zcat_path = which("gzcat")
    if not __zcat_path:
        __zcat_path = which("zcat")
    if not __zcat_path:
        raise ValueError("Unable to find a zcat|gzcat executable in PATH!")
    return __zcat_path
