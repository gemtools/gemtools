#!/usr/bin/env python
"""Gem tools utilities
"""
import re

import subprocess
import logging
from threading import Thread

def read_to_sequence(read):
    """
    Transform the Read to a fasta/fastq seqeucen
    @param read: the read
    @type read: Read
    @return: fasta/q sequence representation
    @rtype: string
    """
    return "%s\n" % read.to_sequence()


def __parse_error_output(stream):
    """Parse the GEM error output stream and
    raise an exception if 'error' occurs"""
    for line in stream:
        if re.search("error", line):
            raise ValueError("GEM run error : %s " % (line))


def run_tool(params, input=None, output=None, name="", transform_fun=read_to_sequence, logfile=None):
    """
    Run the tool defined in the params array using a new process.
    The input is a ReadIterator and the method checks
    if the input comes from a file and can be used directly. If that
    is the case, the input file is passed to the
    parameters using '-i' parameter. If its not the case, a new Thread
    is used to pipe the input to the tool.


    @param params: the parameter list
    @type params: list
    @param input: optional input sequence that is passed through stdin
    @type input: sequence
    @param output: a string that is used as a log file for stdout of the process or None
    @type output: string
    @param name: optional name of the executied tool
    @type name: string
    @param name: optional transformation function that is used if
    @type name: string
    @param logfile: optional path to the log file
    @type logfile: string
    @raise: ValueError in case the execution failed

    """

    process_in = None
    process_out = None
    process_err = None

    ## handle input stream
    if input is not None:
        process_in = subprocess.PIPE

    ## handle output stream
    if output is not None:
        if isinstance(output, basestring):
            process_out = open(output, 'w')
        else:
            process_out = output
    else:
        process_out = subprocess.PIPE

    ## handle error stream
    if logfile is not None:
        process_err = logfile
        if isinstance(logfile, basestring):
            process_err = open(logfile, 'w')
    else:
        process_err = subprocess.PIPE

    logging.info("Starting %s :\n\t%s" % (name, " ".join(params)))
    process = subprocess.Popen(params, stdin=process_in, stdout=process_out, stderr=process_err, close_fds=True)


    if input is not None:
        input_thread = Thread(target=__write_input, args=(input, process.stdin, transform_fun))
        input_thread.start()

    if logfile is None:
        err_thread = Thread(target=__parse_error_output, args=(process.stderr,))
        err_thread.start()


    return process



def __write_input(sequence, stream, transformer=None):
    """
    Write the content form the given sequence to the stream,
    passing ith through the transformer function if specified

    @param sequence: the input sequence
    @type sequence: sequence
    @param stream: the output stream
    @type stream: stream
    @param transformer: optional transformer function
    @type transformer: function
    """
    for e in sequence:
        if transformer is not None:
            e = transformer(e)
        stream.write(str(e))
    stream.close()


def multisplit(s, seps):
    """Split a the string s using multiple separators"""
    res = [s]
    for sep in seps:
        s, res = res, []
        for seq in s:
            res += seq.split(sep)
    return res
