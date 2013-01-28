#!/usr/bin/env python
"""Gem tools utilities
"""
from Queue import Queue
import os
import re
import string

import subprocess
import logging
from threading import Thread

import gem
import sys

import datetime
import time

class Timer(object):
    """Helper class to take runtimes"""
    def __init__(self):
        self.start_time = time.time()

    def stop(self, message):
        end = datetime.timedelta(seconds=int(time.time() - self.start_time))
        if message is not None:
            logging.info(message % (str(end)))

class CommandException(Exception):
    pass


class Command(object):
    """Command base class to be registered 
    with the gem tools main command
    """
    def register(self, parser):
        pass
    def run(self, args):
        pass


class ProcessWrapper(object):
    def __init__(self, _parent, _threads, stdout=None):
        self._parent = _parent
        self._threads = _threads
        self.stdout = _parent.stdout
        self.stderr = _parent.stderr
        self.stdin = _parent.stdin
        if stdout is not None:
            self.stdout = stdout

    def wait(self):
        if self._threads is not None:
            for t in self._threads:
                t.join()
        return self._parent.wait()


class StreamWrapper(object):
    def __init__(self):
        self.queue = Queue()
        self.__closed = False
        self.__written = 0
        self.__read = 0

    def write(self, e):
        self.__written += 1
        self.queue.put(e)


    def readline(self, timeout=None):
        if self.__closed and self.__read == self.__written:
            return None

        if timeout is None:
            timeout = 1
        line = None
        try:
            line = self.queue.get(timeout=timeout)
        except Exception, e:
            if not self.__closed:
                if timeout > 30:
                    raise e
                else:
                    return self.readline(timeout=timeout*2)
        self.__read += 1
        return line

    def __iter__(self):
        return self

    def next(self):
        line = self.readline()
        if line is None:
            raise StopIteration()
        return line

    def close(self):
        self.__closed = True


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
    if read.type == "map" and read.line is not None:
        return "%s\n" % read.line
    return "%s\n" % str(read)


def __parse_error_output(stream, name="Tool"):
    """Parse the GEM error output stream and
    if 'error' occurs exists"""
    logging.debug("Error thread started for %s %s" % (name, stream))
    for line in stream:
        line = line.rstrip()
        if re.search("error", line):
            logging.error("%s raised an error !\n%s\n" % (name, line))
            exit(1)
    logging.debug("Error thread method finished for %s %s" % (name, stream))


complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')

def reverseComplement(sequence):
    """Returns the reverse complement of the given sequence"""
    return sequence.translate(complement)[::-1]


def run_tools(tools, input=None, output=None, name="", transform_fun=read_to_sequence, post_transform=None, logfile=None
              , raw_stream=False, path=None):
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
    @type output: sequence
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
            logging.debug("Appending stderr read thread to watch for errors in %s" % process)
            err_thread = Thread(target=__parse_error_output, args=(process.stderr, tools_name,))
            err_thread.start()

    process_in = None
    process_out = None
    process_err = None

    ## handle input stream
    if input is not None:
        if raw_stream and isinstance(input, gem.files.ReadIterator):
            process_in = input.stream
        else:
            process_in = subprocess.PIPE

    ## handle output stream
    if output is not None or post_transform is not None:
        if isinstance(output, basestring):
            process_out = open(output, 'w')
        else:
            process_out = output
    else:
        process_out = subprocess.PIPE

    ## handle error stream
    logging.debug("Global logging is set to %d" % gem.log_output)
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

    env = None
    if path is not None:
        env = {'PATH': path}
    for i, params in enumerate(tools):
        logging.info("Starting %s :\n\t%s" % (name, " ".join(params).replace("\t", "\\t")))
        #print "Starting %s :\n\t%s" % (name, " ".join(params))
        p_in = process_in
        p_out = process_out

        if first_process is None:
            ## start the first process
            ## and set stdout to PIPE if there are more
            if num_tools > 1 or post_transform is not None:
                p_out = subprocess.PIPE
            first_process = subprocess.Popen(params, stdin=p_in, stdout=p_out, stderr=process_err, close_fds=True,
                env=env)
            current_process = first_process
        else:
            ## add the next process
            ## and set process out if it is the last one
            p_out = process_out
            if i < num_tools - 1 or post_transform is not None:
                p_out = subprocess.PIPE
            current_process = subprocess.Popen(params, stdin=current_process.stdout, stdout=p_out, stderr=process_err,
                close_fds=True, env=env)
        if gem.log_output != gem.LOG_STDERR:
            append_logger(current_process, logfile)
        last_process = current_process

    ## start the input thread
    threads = []
    stdout = last_process.stdout

    if input is not None and not raw_stream:
        logging.debug("Starting input stream thread for %s " % (name))
        input_thread = Thread(target=__write_input, args=(input, first_process.stdin, transform_fun))
        input_thread.start()
        threads.append(input_thread)

    if post_transform is not None and isinstance(post_transform, (list, tuple)):
        if len(post_transform) == 0:
            post_transform = None

    if post_transform is not None:
        ## post_transform output
        logging.debug("Starting post transform thread for %s " % (name))
        thread_out = None
        if output is not None and isinstance(output, basestring):
            thread_out = open(output, 'w')
        else:
            thread_out = StreamWrapper()
        stdout = thread_out

        output_thread = Thread(target=__write_output, args=(last_process.stdout, thread_out, post_transform))
        output_thread.start()
        threads.append(output_thread)

    return ProcessWrapper(last_process, threads, stdout=stdout)


def run_tool(params, input=None, output=None, name="", transform_fun=read_to_sequence, post_transform=None, logfile=None
             , raw_stream=False):
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
    return run_tools([params], input, output, name, transform_fun, post_transform, logfile, raw_stream)


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
            logging.error("Failed to write %s to stream: %s" % (str(e), str(ex)))
            exit(1)
    stream.close()


def __write_output(output, stream, transformer=None):
    """
    Put lines from output through transformer function and write them to
    stream
    @param sequence: the input sequence
    @type sequence: sequence
    @param stream: the output stream
    @type stream: stream
    @param transformer: optional transformer function
    @type transformer: function
    """
    logging.info("Output thread method started for %s" % (stream))
    transformer_list = []
    if isinstance(transformer, (list, tuple)):
        transformer_list.extend(transformer)
    else:
        transformer_list.append(transformer)
    while True:
        e = output.readline()
        if e is None or len(e) == 0:
            break
        if transformer is not None:
            for t in transformer_list:
                e = t(e)
        try:
            stream.write(str(e))
        except Exception, ex:
            logging.error("Failed to write %s to stream: %s" % (str(e), str(ex)))
            stream.close()
            exit(1)
    stream.close()
    logging.info("Output thread method finished for %s" % (stream))


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
        process = subprocess.Popen(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.wait() != 0:
            return None
        output = process.communicate()[0]
        if output is None or len(output) == 0:
            return None

        path = output.split("\n")[0]
        if path is None or len(path) == 0:
            return None
        return path
    except Exception:
        ## ignore exceptions and try path search
        return None

def find_pair(file):
    """find another pair file or return none if it could not be
    found or a tuple of the clean name and the name of the second pair.
    """
    pairs = {
        "0.fastq.gz" : "1.fastq.gz",
        "1.fastq.gz" : "2.fastq.gz",
        "1.fastq" : "2.fastq",
        "0.fastq" : "1.fastq",
        "0.txt.gz" : "1.txt.gz",
        "1.txt.gz" : "2.txt.gz",
        "1.txt" : "2.txt",
        "0.txt" : "1.txt",
    }
    for k, v in pairs.items():
        if file.endswith(k):
            name = file[:-len(k)]
            other = name + v
            if name[-1] in (".", "-", "_"):
                name = name[:-1]
            return (os.path.basename(name), other)
    return (None, None)



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


def gzip(file, threads=1):
    """Helper to call gzip on the given file name
    and compress the file

    If threads is > 1, pigz has to be in path and is used
    """
    logging.debug("Starting GZIP compression for %s" % (file))

    if threads > 1 and which("pigz") is not None:
        if subprocess.Popen(['pigz', '-q', '-f', '-p', str(threads), file]).wait() != 0:
            raise ValueError("Error wile executing pigz on %s" % file)
    else:
        if subprocess.Popen(['gzip', '-f', '-q', file]).wait() != 0:
            raise ValueError("Error wile executing gzip on %s" % file)
    return "%s.gz" % file


class retriever(object):
    """Wrap around an instance of the gem retriever.
    The retriever instance is open until you explicitly
    close it
    """
    def __init__(self, index_hash):
        """Create a new instance specifying the index hash
        that should be used.
        """
        self.index_hash = index_hash
        self.__process = None
        if not os.path.exists(self.index_hash):
            raise ValueError("Index hash %s not found", self.index_hash)

    def __initialize_process(self):
        """Initialize the retriever instance"""
        if self.__process is not None:
            raise ValueError("Retriever instance is running already")
        pa = [gem.executables['gem-retriever'], 'query', self.index_hash]
        self.__process = subprocess.Popen(pa,
                                          stdin=subprocess.PIPE,
                                          stdout=subprocess.PIPE,
                                          shell=False)


    def get_junction(self, chr, strand, start, length):
        """Get junction site. Return a tuple left,right
        with left being the sequence from start to start+4
        and right being start+length-5 to start+length-1
        """
        if self.__process is None:
            self.__initialize_process()

        left = "%s\t%s\t%d\t%d\n"%(chr, strand, start, start+4)
        right = "%s\t%s\t%d\t%d\n"%(chr, strand, start+length-5, start+length-1)
        self.__process.stdin.write(left+right)
        self.__process.stdin.flush()
        g_left = self.__process.stdout.readline().rstrip()
        g_right = self.__process.stdout.readline().rstrip()
        return g_left, g_right

    def close(self):
        """Close the retriever instance"""
        if self.__process is not None:
            self.__process.stdin.close()
            self.__process.stdout.close()
            self.__process = None

