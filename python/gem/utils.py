#!/usr/bin/env python
"""Gem tools utilities and methods
to start external gem processes

In addition, the utilities class currently hosts
the command environment. If you want
to create a new gemtools command, create a subclass
of gem.utils.Command.
"""

import os
import string
import subprocess
import logging
from threading import Thread
import gem
import gem.gemtools as gt
import datetime
import time


class Timer(object):
    """Helper class to take runtimes
    To use this, create a new instance. The timer is
    started at creation time. Call stop() to stop
    timing. The time is logged to info level by default.

    For example:

        timer = Timer()
        ...<long running process>
        time.stop("Process finished in %s")
    """
    def __init__(self):
        """Create a new time and initialie it
        with current time"""
        self.start_time = time.time()

    def stop(self, message, loglevel="info"):
        """Stop timing and print result to logger.
        NOTE you have to add a "%s" into your message string to
        print the time
        """
        end = datetime.timedelta(seconds=int(time.time() - self.start_time))
        if message is not None:
            if loglevel is not None:
                ll = loglevel.upper()
                if ll == "info":
                    logging.info(message % (str(end)))
                elif ll == "warning":
                    logging.warning(message % (str(end)))
                elif ll == "error":
                    logging.error(message % (str(end)))
                elif ll == "debug":
                    logging.error(message % (str(end)))


class CommandException(Exception):
    """Exception thrown by gemtools commands"""
    pass


class Command(object):
    """Command base class to be registered
    with the gemtools main command. The command
    implementation has to implement two methods.

    register() which is called with the argparse parser to
    register new command line options, and
    run() wich is called with the parsed arguments.
    """
    def register(self, parser):
        """Add new command line options to the passed
        argparse parser

        parser -- the argparse parser
        """
        pass

    def run(self, args):
        """Run the command

        args -- the parsed arguments"""
        pass


# complement translation
__complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')


def reverseComplement(sequence):
    """Returns the reverse complement of the given DNA/RNA sequence"""
    return sequence.translate(__complement)[::-1]


class ProcessWrapper(object):
    """Class returned by run_tools that wraps around a list of processes and
    is able to wait. The wrapper is aware of the process log files and
    will do the cleanup around the process when after waiting.

    If a process does not exit with 0, its log file is printed to logger error.

    After the wait, all log files are deleted by default.
    """
    def __init__(self, keep_logfiles=False):
        """Create an empty process wrapper

        keep_logfiles -- if true, log files are not deleted
        """
        self.processlist = []
        self.process_startup = []
        self.logfiles = []
        self.keep_logfiles = keep_logfiles


    def run(self, command, input=None, output=None, env=None ):
        process_in = subprocess.PIPE
        if input is not None:
            process_in = input

        process_out = subprocess.PIPE
        if output is not None:
            process_out = output


        logging.debug("Starting process : %s ", command[0])
        process = subprocess.Popen(command, stdin=process_in, stdout=process_out, stderr=process_err, env=env, close_fds=True)


    def wait(self):
        """Wait for all processes in the process list to
        finish. If a process is exiting with non 0, the process
        log file is printed to logger error.

        All log files are delete if keep_logfiles is False
        """
        try:
            for i, process in enumerate(self.processlist):
                exit_value = process.wait()
                if exit_value != 0:
                    logging.error("Process %s exited with %d!" % (self.process_startup[i][0]), exit_value)
                    if self.logfiles[i] is not None:
                        with open(self.logfiles[i]) as f:
                            for line in f:
                                logging.error("%s" % (line.strip()))
                    return 1
            return 0
        finally:
            if not self.keep_logfiles:
                for logfile in self.logfiles:
                    if logfile is not None and os.path.exists(logfile):
                        logging.debug("Removing log file: %s" % (logfile))
                        os.remove(logfile)


def run_tools(tools, input=None, output=None, write_map=False, clean_id=True, append_extra=True, name=""):
    """
    Run the tools defined in the tools list using a new process per tool.
    The input must be a gem.gemtools.TemplateIterator that is used to get
    the input templates and write their string representation to the
    first process.

    The parameters 'write_map', clean_id', and 'append_extra' can be used to configure the
    output stream. If write_map is True, the output will be written
    in gem format, otherwise, it is transformed to fasta/fastq.
    If clean_id is set to True, the read pair information is always encoded as /1 /2 and
    any casava 1.8 information is dropped. If append_extra is set to False, no additional
    information will be printed to the read tag

    If output is a string or an open file handle, the
    stdout of the final process is piped to that file.

    tools        -- the list of tools to run. This is a list of lists.
    input        -- the input TemplateIterator
    output       -- optional output file name or open, writable file handle
    write_map    -- if true, write input in gem format
    clean_id     -- it true, /1 /2 pair identifiers are enforced
    append_extra -- if false, no additional information is printed in tag
    name         -- optional name for this process group
    """

    process_in = None
    process_out = None
    process_err = None

    ## handle input stream
    if input is not None:
        if isinstance(input, gem.files.ReadIterator):
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
    elif raw_stream:
        logging.debug("Passing raw stream to process")

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


def run_tool(params, input=None, output=None, name="", transform_fun=read_to_sequence, post_transform=None, logfile=None, raw_stream=False):
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
            stream.close()
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

