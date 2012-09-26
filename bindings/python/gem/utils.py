#!/usr/bin/env python
"""Gem tools utilities
"""
import os
import re
import string

import subprocess
import logging
from threading import Thread

import gem

def read_to_sequence(read):
    """
    Transform the Read to a fasta/fastq seqeucen
    @param read: the read
    @type read: Read
    @return: fasta/q sequence representation
    @rtype: string
    """
    return "%s\n" % read.to_sequence()


def read_to_map(read):
    """
    Transform the Read to a gem map line
    @param read: the read
    @type read: Read
    @return: fasta/q sequence representation
    @rtype: string
    """
    return "%s\n" % str(read)


def __parse_error_output(stream, name="Tool"):
    """Parse the GEM error output stream and
    raise an exception if 'error' occurs"""
    for line in stream:
        line = line.rstrip()
        if re.search("error", line):
            logging.error("%s raised an error !\n%s\n" % (name, line))



complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')

def reverseComplement(sequence):
    """Returns the reverse complement of the given sequence"""
    return sequence.translate(complement)[::-1]


def run_tools(tools, input=None, output=None, name="", transform_fun=read_to_sequence, logfile=None, raw_stream=False):
    """
    Run the tools defined in the tools list using a new process per tools.
    The input is a ReadIterator and the method checks
    if the input comes from a file and can be used directly. If that
    is the case, the input file is passed to the
    parameters using '-i' parameter. If its not the case, a new Thread
    is used to pipe the input to the tool.
    IF multiple tools are specified the stdout and stdin streams are connected
    automatically


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

    ## helper function to append a logger thread per process
    def append_logger(process, logfile):
        if logfile is None:
            tools_name = None
            if name:
                tools_name = name
            err_thread = Thread(target=__parse_error_output, args=(process.stderr,tools_name,))
            err_thread.start()

    process_in = None
    process_out = None
    process_err = None

    ## handle input stream
    if input is not None:
        if raw_stream:
            process_in = input.stream
        else:
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
    if logfile is not None and gem.log_output != gem.LOG_STDERR:
        process_err = logfile
        if isinstance(logfile, basestring):
            process_err = open(logfile, 'w')
    elif gem.log_output != gem.LOG_STDERR:
        process_err = subprocess.PIPE

    num_tools = len(tools)
    first_process = None
    current_process = None
    last_process = None

    for i, params in enumerate(tools):
        logging.info("Starting %s :\n\t%s" % (name, " ".join(params)))
        #print "Starting %s :\n\t%s" % (name, " ".join(params))
        p_in = process_in
        p_out = process_out

        if first_process is None:
            ## start the first process
            ## and set stdout to PIPE if there are more
            if num_tools > 1:
                p_out = subprocess.PIPE
            first_process = subprocess.Popen(params, stdin=p_in, stdout=p_out, stderr=process_err, close_fds=True)
            current_process = first_process
        else:
            ## add the next process
            ## and set process out if it is the last one
            p_out = process_out
            if i < num_tools - 1:
                p_out = subprocess.PIPE
            current_process = subprocess.Popen(params, stdin=current_process.stdout, stdout=p_out, stderr=process_err, close_fds=True)
        if gem.log_output != gem.LOG_STDERR:
            append_logger(current_process, logfile)
        last_process = current_process


    ## start the input thread
    if input is not None and not raw_stream:
        input_thread = Thread(target=__write_input, args=(input, first_process.stdin, transform_fun))
        input_thread.start()

    return last_process


def run_tool(params, input=None, output=None, name="", transform_fun=read_to_sequence, logfile=None, raw_stream=False):
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
    return run_tools([params], input, output, name, transform_fun, logfile, raw_stream)

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
        try:
            stream.write(str(e))
        except Exception, ex:
            logging.error("Failed to write %s to stream: %s" % (str(e), ex.message))
            exit(1)
    stream.close()


def multisplit(s, seps):
    """Split a the string s using multiple separators"""
    res = [s]
    for sep in seps:
        s, res = res, []
        for seq in s:
            res += seq.split(sep)
    return res


def which(program):
    """
    Find programm in path
    """
    import subprocess

    ## use which command
    try:
        params = ["which", program]
        output = subprocess.check_output(params)
        path = output.split("\n")[0]
        if path is None or len(path) == 0:
            return None
        return path
    except Exception:
        ## ignore exceptions and try path search
        return None


def find_in_path(program):
    """
    Explicitly search the PATH for the given program

    @param programm: the program
    @type programm: string
    @return: absolute path to the program or None
    @rtype: string
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def gzip(file):
    """Helper to call gzip on the given file name
    and compress the file
    """
    if subprocess.Popen(['gzip', file]).wait() != 0:
        raise ValueError("Error wile executing gzip on %s" % file)
    return "%s.gz" % file
