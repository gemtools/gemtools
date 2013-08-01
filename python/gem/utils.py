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
import gem
import gem.gemtools as gt
import datetime
import time
import tempfile
import multiprocessing as mp
import json
import subprocess
import signal


# clobal process registry
# which is used to save multiprocessing.Process instances
# to be able to terminate them in case
# an error occured somewhere
_process_registry = set([])


def register_process(process):
    """Register a Process with the registry"""
    _process_registry.add(process)


def terminate_processes():
    """Do a hard cleanup and terminata all registered
    processes.
    """
    for p in _process_registry:
        if p.is_alive() and p._popen is not None:
            p.terminate()
            if p.is_alive():  # terminate send a SIGTERM, now kill it properly
                os.kill(p._popen.pid, signal.SIGKILL)


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
        self.end = None

    def stop(self, message, loglevel="info"):
        """Stop timing and print result to logger.
        NOTE you have to add a "%s" into your message string to
        print the time
        """
        self.end = datetime.timedelta(seconds=int(time.time() - self.start_time))
        if message is not None:
            if loglevel is not None:
                ll = loglevel.lower()
                if ll == "info":
                    logging.info(message % (str(self.end)))
                elif ll == "warning":
                    logging.warning(message % (str(self.end)))
                elif ll == "error":
                    logging.error(message % (str(self.end)))
                elif ll == "debug":
                    logging.error(message % (str(self.end)))


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
    def add_options(self, tool, parser):
        """Reads the tools JSON options and adds them to the
        command line parser. The commands 'tool_options' dict is
        set after execution and the parsed args can be translated to
        a command argument with get_command()
        """
        ## required parameters
        # read the json

        output = subprocess.Popen([gem.executables[tool], '-J'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE).communicate()[1]

        args = json.loads(output)
        groups = {}
        self.tool = tool
        self.tool_opts = {}
        if "options" in args:
            options = args['options']
            opt_parser = parser
            for opt in options:
                if opt.get('shortOption', "") == "h":
                    continue
                # get the group parser or create it
                if "group" in opt:
                    group = opt['group']
                    if len(group) > 0:
                        if group not in groups:
                            opt_parser = parser.add_argument_group(group)
                            groups[group] = opt_parser
                        else:
                            opt_parser = groups[group]

                shortOption = opt.get('shortOption', None)
                longOption = opt.get('longOption', None)
                optionType = opt.get('optionType', "")
                metaVar = opt.get("commandInfo", None)
                description = opt.get("description", None)
                all_options = []

                if shortOption is not None:
                    all_options.append("-" + shortOption)
                if longOption is not None:
                    all_options.append("--" + longOption)

                if metaVar is not None:
                    metaVar = metaVar.split(" ")[0]
                kwargs = {"dest": longOption,
                          "metavar": metaVar,
                          "help": description
                          }

                if opt.get('argumentType', None) is not None:
                    t = opt['argumentType']
                    if t == "int":
                        kwargs['type'] = int
                    elif t == "float":
                        kwargs['type'] = float

                if optionType == "noArgument":
                    kwargs['action'] = "store_true"
                    del kwargs['metavar']
                opt_parser.add_argument(*all_options, **kwargs)
                self.tool_opts[longOption] = kwargs

    def get_command(self, args):
        cmd = [gem.executables[self.tool]]
        dargs = vars(args)
        for k, v in dargs.items():
            if k in self.tool_opts and v is not None:
                if not isinstance(v, bool):
                    cmd.append("--" + k)
                    cmd.append(str(v))
                elif v:
                    cmd.append("--" + k)
        return cmd



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
        """"Initialize a single process. the process takes the outer wrapper, which is basically the pipeline
        integrating this process and the commands. Additionally, the process can have defined input and output.
        If a parent process is given, any input setting is overwritten and the stdout of
        the parent is used.

        wrapper  -- the outer pipeline
        commands -- the commands to be executed
        input    -- the input
        output   -- the output
        parent   -- optional parent process that defines the input
        env      -- optional environment definition
        logfile  -- optional path to the log file
        """
        self.wrapper = wrapper
        self.commands = commands
        self.input = input
        self.output = output
        self.env = env
        self.process = None
        self.logfile = logfile
        self.parent = parent
        self.input_writer = None

    def run(self):
        """Start the process and return it. If the input is a ProcessInput,
        the input of this process is set to pipe and the ProcessInput is written
        there.
        """
        stdin = self.input
        stdout = self.output
        stderr = None
        process_input = None

        # check the logfile
        if self.logfile is not None:
            if isinstance(self.logfile, basestring):
                logging.debug("Setting process log file to %s", self.logfile)
                stderr = open(self.logfile, 'wb')
            else:
                stderr = self.logfile

        #check outout
        if stdout is not None and isinstance(stdout, basestring):
            logging.debug("File output detected, opening output stream to %s", stdout)
            stdout = open(stdout, "wb")
        elif stdout is None:
            stdout = subprocess.PIPE

        # check input
        if self.parent is not None:
            logging.debug("Setting process input to parent output")
            stdin = self.parent.process.stdout
        else:
            if isinstance(stdin, ProcessInput):
                logging.debug("Process-Input detected, switching to PIPE input")
                process_input = stdin
                stdin = subprocess.PIPE

        logging.debug("Starting subprocess")
        self.process = subprocess.Popen(self.commands, stdin=stdin, stdout=stdout, stderr=stderr, env=self.env, close_fds=False)

        if process_input is not None:
            logging.debug("Starting process input writer")
            process_input.write(self.process)
            self.input_writer = process_input

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
        0, print the log file to error
        """
        if self.process is None:
            logging.error("Process was not started")
            raise ProcessError("Process was not started!", self)

        if self.input_writer is not None:
            logging.debug("Waiting for process input writer to finish")
            self.input_writer.wait()

        # wait for the process
        exit_value = self.process.wait()
        logging.debug("Process '%s' finished with %d", str(self), exit_value)
        if exit_value is not 0:
            logging.error("Process '%s' finished with %d", str(self), exit_value)
            if self.logfile is not None and isinstance(self.logfile, basestring):
                with open(self.logfile) as f:
                    for line in f:
                        logging.error("%s" % (line.strip()))
            raise ProcessError("Process '%s' finished with %d" % (str(self), exit_value))
        return exit_value

    def to_bash(self):
        """Returns the bash command representation
        """
        if isinstance(self.commands, (list, tuple)):
                return " ".join(self.commands)
        return str(self.commands)


class ProcessInput(object):
    """Takes a gemtools InputFile or filter and uses the write_stream() method
    to write the content to the given process stdin.
    """

    def __init__(self, input, write_map=False, clean_id=True, append_extra=True):
        """Initialize a new ProcessInput from a gemtools.InputFile or filter

        input         -- the input
        write_map     -- write in map format, otherwise writes in fasta/q
        clean_id      -- clean ids and ensure /1 /2 pairing ids
        append_extra  -- include any additional infomration stored in the read ids
        """
        self.write_map = write_map
        self.clean_id = clean_id
        self.append_extra = append_extra
        self.process = None
        self.input = input
        self.thread = None

    def __write_input(self):
        """Internal method that writes templates ot the
        input stream of the process and puts a message on
        the queue when done to close the stream
        """
        outfile = gt.OutputFile(self.process.stdin, clean_id=self.clean_id, append_extra=self.append_extra)
        self.input.write_stream(outfile, write_map=self.write_map, threads=2)

        try:
            if self.input.process is not None:
                self.input.process.wait()
        except Exception:
            pass

        outfile.close()
        self.process.stdin.close()

    def write(self, process):
        self.process = process
        logging.debug("Preparing process input stream -- clean_id: %s, append_extra: %s, write_map:%s" % (str(self.clean_id), str(self.append_extra), str(self.write_map)))
        self.thread = mp.Process(target=ProcessInput.__write_input, args=(self,))
        register_process(self.thread)
        #self.thread = Thread(target=ProcessInput.__write_input, args=(self,))
        self.thread.start()
        self.process.stdin.close()

    def wait(self):
        if self.thread is not None:
            self.thread.join()


class ProcessWrapper(object):
    """Class returned by run_tools that wraps around a list of processes and
    is able to wait. The wrapper is aware of the process log files and
    will do the cleanup around the process when after waiting.

    If a process does not exit with 0, its log file is printed to logger error.

    After the wait, all log files are deleted by default.
    """

    def __init__(self, keep_logfiles=False, name=None, force_debug=False, raw=None):
        """Create an empty process wrapper

        keep_logfiles -- if true, log files are not deleted
        """
        self.processes = []
        self.keep_logfiles = keep_logfiles
        self.name = name
        self.stdout = None
        self.stdin = None
        self.force_debug = force_debug
        self.raw = raw
        self.exit_value = None

    def submit(self, command, input=subprocess.PIPE, output=None, env=None, logfile=None):
        """Run a command. The command must be list of command and its parameters.
        If input is specified, it is passed to the stdin of the subprocess instance.
        If output is specified, it is connected to the stdout of the underlying subprocess.
        Environment is optional and will be passed to the process as well.

        This is indened to be used in pipes and specifying output will close the pipe
        """
        logfile = logfile
        parent = None
        if len(self.processes) > 0:
            parent = self.processes[-1]
        if logfile is None and logging.getLogger().level is not logging.DEBUG and not self.force_debug:
            # create a temporary log file
            tmpfile = tempfile.NamedTemporaryFile(suffix='.log', prefix=self.__command_name(command) + ".", delete=(not self.keep_logfiles))
            logfile = tmpfile.name
            tmpfile.close()

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
            name = name.split("/")[-1]
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
        if self.exit_value is not None:
            return self.exit_value
        try:
            if self.raw:
                for r in self.raw:
                    if r is not None:
                        r.wait()
            exit_value = 0
            for i, process in enumerate(self.processes):
                ev = process.wait()
                if exit_value is not 0:
                    exit_value = ev
            self.exit_value = exit_value
            return ev
        except Exception, e:
            if isinstance(e, OSError) and e.errno == 10:
                pass
            else:
                print "An error occured while waiting for one the child processes:", e
                self.exit_value = 1
        finally:
            if not self.keep_logfiles:
                for p in self.processes:
                    if p.logfile is not None and isinstance(p.logfile, basestring) and os.path.exists(p.logfile):
                        logging.debug("Removing log file: %s" % (p.logfile))
                        os.remove(p.logfile)

    def to_bash_pipe(self):
        return " | ".join([p.to_bash() for p in self.processes])


def _prepare_input(input, write_map=False, clean_id=True, append_extra=True):
    if isinstance(input, basestring):
        return open(input, 'rb')
    if isinstance(input, file):
        return input
    elif input is not None:
        return ProcessInput(input, write_map=write_map, clean_id=clean_id, append_extra=append_extra)
    return None


def _prepare_output(output):
    if isinstance(output, basestring):
        return output
    if isinstance(output, file):
        if output.name is not None and output.name is not "<fdopen>":
            output.close()
            return output.name
        raise ProcessError("Can not pass raw file descriptors")
    return None


def run_tools(tools, input=None, output=None, write_map=False, clean_id=False,
              append_extra=True, name=None, keep_logfiles=False,
              force_debug=False, env=None, raw=False, logfile=None):
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
    logfile      -- specify a filename or a string that is used as stderr
    """
    parent_process = None
    if raw:
        if isinstance(input, gt.InputFile):
            parent_process = [input.process]
    p = ProcessWrapper(keep_logfiles=keep_logfiles, name=name, force_debug=force_debug, raw=parent_process)

    for i, commands in enumerate(tools):
        process_in = subprocess.PIPE
        process_out = subprocess.PIPE
        if i == 0:
            # prepare first process input
            if raw:
                process_in = _prepare_input(input.raw_stream(), write_map=write_map, clean_id=clean_id, append_extra=append_extra)
            else:
                process_in = _prepare_input(input, write_map=write_map, clean_id=clean_id, append_extra=append_extra)

        if i == len(tools) - 1:
            # prepare last process output
            process_out = _prepare_output(output)
        p.submit(commands, input=process_in, output=process_out, env=env, logfile=logfile)

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
        "0.fq": "1.fq",
        "1.fq": "2.fq",
        "0.fq.gz": "1.fq.gz",
        "1.fq.gz": "2.fq.gz",
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
    def __init__(self, index):
        """Create a new instance specifying the index
        that should be used.
        """
        self.index = index
        self.__process = None
        if not os.path.exists(self.index):
            raise ValueError("Index %s not found", self.index)

    def __initialize_process(self):
        """Initialize the retriever instance"""
        if self.__process is not None:
            raise ValueError("Retriever instance is running already")
        pa = [gem.executables['gem-retriever'], self.index]
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


def get_max_read_length(input, threads=1, paired=False):
    p = [gem.executables['gt.stats'], "--threads", str(threads)]
    if paired:
        p.append("-p")

    process = run_tool(p, input=input, output=None, name="stats")
    import re
    pattern = re.compile(r'.*Reads.Length \(min,avg,max\) \(\d+,\d+,(\d+)\)')
    max_len = -1
    for line in process.processes[0].process.stdout:
        m = pattern.match(line)
        if m:
            max_len = int(m.groups()[0])
    if process.wait() != 0:
        raise ValueError("Error while getting max read length from stats!")
    return max_len
