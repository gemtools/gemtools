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
import tempfile
import multiprocessing as mp


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


class ProcessError(Exception):
    """Error thrown by the process wrapper when there is a problem
    with either starting of the process or at runtime"""
    def __init__(self, message, process=None):
        """Initialize the exception with a message and the process wrapper that caused the
        exception

        message -- the error message
        process -- the underlying process wrapper
        """
        Exception.__init__(self, message)
        self.process = process


class Process(object):
    """Single process in a pipeline of processes"""

    def __init__(self, wrapper, commands, input=subprocess.PIPE, output=subprocess.PIPE, parent=None, env=None, logfile=None):
        """"Internal process"""
        self.wrapper = wrapper
        self.commands = commands
        self.input = input
        self.output = output
        self.env = env
        self.process = None
        self.logfile = logfile
        self.parent = parent
        self.process_input = None

    def run(self):
        """Start the process and return it"""
        stdin = self.input
        stdout = self.output
        stderr = None
        if self.parent is not None:
            stdin = self.parent.process.stdout
        if self.logfile is not None:
            stderr = open(self.logfile, 'wb')
        self.process_input = None
        if isinstance(stdin, ProcessInput):
            self.process_input = stdin
            stdin = subprocess.PIPE

        self.process = subprocess.Popen(self.commands, stdin=stdin, stdout=stdout, stderr=stderr, env=self.env, close_fds=True)

        if self.process_input is not None:
            self.process_input.start(self.process)

        return self.process

    def __str__(self):
        if self.commands is None:
            return "<process>"
        else:
            if isinstance(self.commands, (list, tuple)):
                return str(self.commands[0])
            return str(self.commands)

    def wait(self):
        """Wait for the process and return its exit value. If it did not exit with
        0, print the log file to error"""
        if self.process is None:
            raise ProcessError("Process was not started!", self)

        if self.process_input is not None:
            self.process_input.join()

        exit_value = self.process.wait()
        if exit_value is not 0:
            logging.error("Process %s exited with %d!" % (str(self), exit_value))
            if self.logfile is not None:
                with open(self.logfile) as f:
                    for line in f:
                        logging.error("%s" % (line.strip()))
        else:
            logging.debug("Process %s finished" % (str(self)))
        return exit_value

    def to_bash(self):
        if isinstance(self.commands, (list, tuple)):
                return " ".join(self.commands)
        return str(self.commands)


class ProcessWrapper(object):
    """Class returned by run_tools that wraps around a list of processes and
    is able to wait. The wrapper is aware of the process log files and
    will do the cleanup around the process when after waiting.

    If a process does not exit with 0, its log file is printed to logger error.

    After the wait, all log files are deleted by default.
    """

    def __init__(self, keep_logfiles=False, name=None, force_debug=False):
        """Create an empty process wrapper

        keep_logfiles -- if true, log files are not deleted
        """
        self.processes = []
        self.keep_logfiles = keep_logfiles
        self.name = name
        self.stdout = None
        self.stdin = None
        self.force_debug = force_debug

    def submit(self, command, input=subprocess.PIPE, output=subprocess.PIPE, env=None):
        """Run a command. The command must be list of command and its parameters.
        If input is specified, it is passed to the stdin of the subprocess instance.
        If output is specified, it is connected to the stdout of the underlying subprocess.
        Environment is optional and will be passed to the process as well.

        This is indened to be used in pipes and specifying output will close the pipe
        """
        logfile = None
        parent = None
        if len(self.processes) > 0:
            parent = self.processes[-1]
        if logging.getLogger().level is not logging.DEBUG and not self.force_debug:
            # create a temporary log file
            (fifo, logfile) = tempfile.mkstemp(suffix='.log', prefix=self.__command_name(command), dir=".")
            # close the file handle
            os.close(fifo)

        p = Process(self, command, input=input, output=output, env=env, logfile=logfile, parent=parent)
        self.processes.append(p)
        return p

    def __command_name(self, command):
        """Create a name for the given command. The name
        is either based on the specified wrapper name or on
        the first part of the command.

        The process list index is always appended.

        command -- the command
        """
        name = None
        if self.name is not None:
            name = self.name
        else:
            if isinstance(command, (list, tuple)):
                name = command[0].split()[0]
            else:
                name = str(command.split()[0])
        return "%s.%d" % (name, len(self.processes))

    def start(self):
        """Start the process pipe"""
        logging.info("Starting:\n\t%s" % (self.to_bash_pipe()))
        for p in self.processes:
            p.run()
        self.stdin = self.processes[0].process.stdin
        self.stdout = self.processes[-1].process.stdout

    def wait(self):
        """Wait for all processes in the process list to
        finish. If a process is exiting with non 0, the process
        log file is printed to logger error.

        All log files are delete if keep_logfiles is False
        """
        try:
            exit_value = 0
            for i, process in enumerate(self.processes):
                ev = process.wait()
                if exit_value is not 0:
                    exit_value = ev
            return ev
        finally:
            if not self.keep_logfiles:
                for p in self.processes:
                    if p.logfile is not None and os.path.exists(p.logfile):
                        logging.debug("Removing log file: %s" % (p.logfile))
                        os.remove(p.logfile)

    def to_bash_pipe(self):
        return " | ".join([p.to_bash() for p in self.processes])


class ProcessInput(object):
    """Helper class to write templates to a target process"""
    def __init__(self, input, write_map=False, clean_id=True, append_extra=True):
        self.write_map = write_map
        self.clean_id = clean_id
        self.append_extra = append_extra
        self.process = None
        self.input = input
        self.target_process = None

    def __write_input(self, process, q):
        """Internal method that writes templates ot the
        input stream of the process and puts a message on
        the queue when done to close the stream
        """
        of = None
        try:
            iterator = self.input
            if isinstance(iterator, gt.InputFile):
                iterator = iterator.templates()
            of = gt.OutputFile(stream=process.stdin)
            if self.write_map:
                of.write_map(iterator, clean_id=self.clean_id, append_extra=self.append_extra)
            else:
                of.write_fastq(iterator, clean_id=self.clean_id, append_extra=self.append_extra)
        finally:
            ## push to the queue to close the sreams
            logging.debug("External process writer finished, sending signal")
            q.put("Done")

    def __wait_for_process(self, q, target_process):
        """Internal method that is starte in a separate thread
        and waits for a signal from the queue before it closes the
        target_process input stream"""
        q.get()
        logging.debug("Signal received, closing input stream")
        target_process.stdin.close()

    def start(self, target_process):
        self.target_process = target_process
        q = mp.Queue()
        Thread(target=ProcessInput.__wait_for_process, args=(self, q, target_process)).start()
        self.process = mp.Process(target=ProcessInput.__write_input, args=(self, self.target_process, q,))
        self.process.start()

    def join(self):
        if self.process is not None:
            self.process.join()


def _prepare_input(input, write_map=False, clean_id=True, append_extra=True):
    if isinstance(input, basestring):
        return open(input, 'rb')
    if isinstance(input, file):
        return input
    else:
        return ProcessInput(input, write_map=write_map, clean_id=clean_id, append_extra=append_extra)


def _prepare_output(output):
    if isinstance(output, basestring):
        return open(output, 'wb')
    if isinstance(output, file):
        return output
    return subprocess.PIPE


def run_tools(tools, input=None, output=None, write_map=False, clean_id=False, append_extra=True, name=None, keep_logfiles=False, force_debug=False, env=None):
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
    p = ProcessWrapper(keep_logfiles=keep_logfiles, name=name, force_debug=force_debug)

    for i, commands in enumerate(tools):
        process_in = subprocess.PIPE
        process_out = subprocess.PIPE
        if i == 0:
            # prepare first process input
            process_in = _prepare_input(input, write_map=write_map, clean_id=clean_id, append_extra=append_extra)
        if i == len(tools) - 1:
            # prepare last process output
            process_out = _prepare_output(output)
        p.submit(commands, input=process_in, output=process_out, env=env)

    # start the run
    p.start()
    return p


def run_tool(tool, **kwargs):
    """
    Delegates to run_tools() with just a single tool
    """
    return run_tools([tool], **kwargs)


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
        "0.fastq.gz": "1.fastq.gz",
        "1.fastq.gz": "2.fastq.gz",
        "1.fastq": "2.fastq",
        "0.fastq": "1.fastq",
        "0.txt.gz": "1.txt.gz",
        "1.txt.gz": "2.txt.gz",
        "1.txt": "2.txt",
        "0.txt": "1.txt",
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

