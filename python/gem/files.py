#!/usr/bin/env python
"""
gem.files handles opening files and streams
"""
from itertools import islice
import os
import logging
import subprocess
import __builtin__
import gem
import utils
from utils import which


__author__ = 'Thasso Griebel'
__zcat_path = None

class Parser(object):
    def __init__(self):
        self.read = gem.Read()

    def next(self, stream):
        """Implement this method
        and parse the next entry from the
        stream

        @param stream: the input stream
        @type stream: stream
        """
        pass


class parse_fasta(Parser):
    """Parse fasta entries from a stream
    """

    def next(self, stream):
        fasta_lines = list(islice(stream, 2))  # read in chunks of 2
        if not fasta_lines:
            return None
        self.read.type = "fasta"
        self.read.id = fasta_lines[0].rstrip()[1:]
        self.read.sequence = fasta_lines[1].rstrip()
        self.read.qualities = None
        self.read.summary = None
        self.read.mappings = None
        self.read.line = "".join(fasta_lines).strip()

        if len(self.read.sequence) <= 0:
            return self.next(stream)
        return self.read


class parse_fastq(Parser):
    """Parse fastq entries from a stream"""

    def next(self, stream):
        fastq_lines = list(islice(stream, 4))  # read in chunks of 2
        if not fastq_lines:
            return None
        self.read.type = "fastq"
        self.read.id = fastq_lines[0].rstrip()[1:]
        self.read.sequence = fastq_lines[1].rstrip()
        self.read.qualities = fastq_lines[3].rstrip()
        self.read.summary = None
        self.read.mappings = None
        self.read.line = "".join(fastq_lines).strip()

        if len(self.read.sequence) <= 0:
            return self.next(stream)
        return self.read


class parse_map(Parser):
    """Parse gem map entries from a stream"""

    def next(self, stream):
        line = stream.readline()
        if not line:
            return None
        self.line2read(line)
        return self.read

    def line2read(self, line):
        line = line.rstrip()
        split = line.split("\t")
        self.read.type = "map"
        self.read.id = split[0]
        self.read.sequence = split[1]
        self.read.line = line
        if len(split) == 5:
            # read with qualities
            self.read.qualities = split[2]
            self.read.summary = split[3]
            self.read.mappings = split[4]
        else:
            # read with no qualities
            self.read.summary = split[2]
            self.read.mappings = split[3]
            self.read.qualities = None

        if self.read.summary == None or self.read.summary in ['-', '*']:
            self.read.summary = "0"
        elif self.read.summary in ['+', '!']:
            self.read.summary = "%d" % gem._max_mappings
        return self.read

class parse_sam(Parser):
    """Parse gem map entries from a stream"""

    def next(self, stream):
        while True:
            line = stream.readline()
            if not line:
                return None
            if line.startswith("@"):
                continue
            line = line.rstrip()
            split = line.split("\t")
            self.read.type = "sam"
            self.read.line = line
            self.read.id = split[0]
            self.read.sequence = split[9]
            self.read.qualities = split[10]
            self.read.mappings = None
            self.read.summary = None
            self.read.line = line
            if self.read.qualities in ["", "*"]:
                self.read.qualities = None


            ## update id and add /1 /2 if its a paired read
            flag = int(split[1])
            if 0x1 & flag == 0x1 and not self.read.id.endswith("/1") and not self.read.id.endswith("/2"):
                ## multiple segments
                if 0x40 & flag == 0x40:
                    self.read.id += "/1"
                else:
                    self.read.id += "/2"

            chr = split[2]
            if chr != "*":
                if flag & 0x10 == 0x10:
                    ## seq is reverser complement
                    self.read.sequence = utils.reverseComplement(self.read.sequence)
                    if self.read.qualities:
                        self.read.qualities = self.read.qualities[::-1]

            return self.read


class ReadIterator(object):
    def __init__(self, stream, parser, filename=None, process=None, remove_after_iteration=False, quality=None, raw=False):
        """
        Create a ReadIterator from a stream with a given parser.
        If the filename is given, the iterator can be cloned to re-read
        the file. Cloning the iterator will not change the current state
        of this instance. The returned clone will start again from the
        beginning of the file.

        @param stream: the input stream
        @type stream: stream
        @param parser: the read parser
        @type type: Parser
        @param filename: optional name of the input file
        @type filename: string
        @param process: optinal process associated with this reader
        @type process: process
        @param remove_after_iteration: if set to True, the input file is removed after a completed iteration of the content
        @type remove_after_iteration: boolean
        """
        self.stream = stream
        self.parser = parser
        self.filename = filename
        self.process = process
        self.remove_after_iteration = remove_after_iteration
        self.quality = quality
        self.raw = raw

    def __iter__(self):
        return self

    def next(self):
        """ Delegs ates to the parser
        to get the next entry
        """
        ret = None
        if self.raw:
            ret = self.stream.readline()
        else:
            ret = self.parser.next(self.stream)
        if not ret:
            logging.debug("Read iterator closing stream on %s" % (self.filename))
            self.stream.close()
            if self.remove_after_iteration and self.filename:
                os.remove(self.filename)
            raise StopIteration()
        return ret

    def close(self):
        logging.debug("Read iterator closing stream on %s" % (self.filename))
        self.stream.close()
        if self.remove_after_iteration and self.filename:
            os.remove(self.filename)

    def clone(self):
        if not self.filename or self.remove_after_iteration:
            raise ValueError("No filename given or file is marked for deletion, this reader can not be cloned!")
        return ReadIterator(open_file(self.filename), self.parser.__class__(), self.filename, quality=self.quality)


## type to parser map
supported_types = {
    "fasta": parse_fasta,
    "fastq": parse_fastq,
    "map": parse_map,
    "sam": parse_sam,
    "bam": parse_sam
}


def open(input, type=None, process=None, remove_after_iteration=False, quality=None, raw=False):
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
    if type is None and is_string:
        type = _guess_type(input)
        if type is None:
            type = "map"
    if type not in supported_types.keys():
        raise ValueError("Unknown file type %s, call this method with a specific type (%s)" % (input, supported_types))

    stream = input
    if is_string:
        if type == "bam":
            process = subprocess.Popen(["samtools", "view", input],
                stdout=subprocess.PIPE,
                stderr=__builtin__.open("/dev/null", 'w'),
                close_fds=True)
            stream = process.stdout
        else:
            stream = open_file(input)
    else:
        input = None  ## reset filename

    return ReadIterator(stream, supported_types[type](), input, process=process,
                        remove_after_iteration=remove_after_iteration, quality=quality, raw=raw)


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
    elif name.endswith(".FASTQ") or name.endswith("FQ"): ## fixes issue #5
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
    __zcat_path = which("zcat")
    if not __zcat_path:
        __zcat_path = which("gzcat")
    if not __zcat_path:
        raise ValueError("Unable to find a zcat|gzcat executable in PATH!")
    return __zcat_path
