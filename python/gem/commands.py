#!/usr/bin/env python
"""Gemtools commands"""

import imp
import sys
from sys import exit
import types
from textwrap import dedent
import os
import logging
from datetime import datetime, timedelta

import gem.utils
from gem.utils import CommandException

import jip
from jip.cli import colorize, YELLOW, GREEN, RED, BLUE

__VERSION__ = "1.8"

log = logging.getLogger("gemtools.command")


class cli(object):
    """Command decorator"""
    _cmd_registry = {}

    def __init__(self, name, inputs=None, outputs=None,
                 title=None, description=None, add_outputs=None,
                 pipeline=None):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.title = title
        self.description = description
        self.add_outputs = add_outputs
        self.pipeline = pipeline

    def __call__(self, cls):
        wrapped_function = False
        if isinstance(cls, types.FunctionType):
            wrapped_function = True
            ## wrap the function in a new class
            ref_fun = cls

            class __command_wrapper(gem.utils.Command):
                title = self.title
                description = self.description

                def run(self, args):
                    ref_fun(args)

                def register(self, parser):
                    pass

                @property
                def __doc__(self):
                    return ref_fun.__doc__

            cls = __command_wrapper

        cls.name = self.name
        cli._cmd_registry[self.name] = cls

        def _help(inst):
            from jip.tools import Tool
            return "%s\n%s" % (cls.description, Tool.help(inst))

        def _validate(inst):
            args = inst.options.to_dict()
            inst.instance.validate(args)
            for k, v in args.iteritems():
                opt = inst.options[k]
                if opt is not None:
                    opt.set(v)

        # register as jip tool
        jip_name = "gemtools_%s" % (self.name.replace("-", "_"))
        jip.tool(jip_name, inputs=self.inputs,
                 outputs=self.outputs,
                 validate="validate",
                 help=_help,
                 add_outputs=self.add_outputs,
                 pipeline=self.pipeline,
                 run='jip_run',
                 argparse='register' if not wrapped_function else None,
                 get_command='jip_command')(cls)
        return cls


def _load_plugins():
    """Load plugin modules form PYTHONPATH and dedicated folders

    Plugin modules are loaded in three flavors:

        1) we try to load a gemtools_plugins module from python path
        2) we check the GEMTOOLS_PLUGINS environment variable for
           folders that exists and load all modules that end in _plugin
           from those folders or all files directly
        3) we check the $HOME/.gemtools/plugins folder and load all files
           there as modules
    """
    # 1) try to load gemtools_plugins
    try:
        import gemtools_plugins
        log.info("gemtools_plugins module loaded")
    except ImportError:
        pass

    # 2) check GEMTOOLS_PLUGINS envorinment
    plugin_paths = os.getenv("GEMTOOLS_PLUGINS", None)
    if plugin_paths:
        for path in plugin_paths.split(":"):
            _load_plugins_from_folder(path)
    # 3) load $HOME/.gemtools/plugins folder
    _load_plugins_from_folder(
        os.path.join(os.getenv("HOME", ""), ".gemtools/plugins")
    )


def _load_plugins_from_folder(path):
    """Helper function that puts path in the python search path and
    loads modules that end in _plugin

    :param path: the root path
    """
    if not os.path.exists(path):
        return
    if os.path.isfile(path):
        ## load file as module
        module_name = os.path.basename(path)[:-3]
        imp.load_source(module_name, path)
    else:
        # we found a directory. put the directory in python path
        # and iterate all files, loading them as modules if
        # they end in _plugin.py
        sys.path.append(os.path.abspath(path))
        for root, dirs, files in os.walk(path):
            for f in [x for x in files if x.endswith("_plugin.py")]:
                module_dir = root.replace(path, "")
                if len(module_dir) > 0 and module_dir[0] == "/":
                    module_dir = module_dir[1:]
                module_dir = module_dir.replace("/", ".")
                module_name = f.replace(".py", "")
                if module_dir:
                    if "__init__.py" in files:
                        __import__(".".join([module_dir, module_name]))
                    else:
                        imp.load_source(
                            module_name,
                            os.path.join(root, f)
                        )
                else:
                    __import__(module_name)


def gemtools():
    """\
    The gemtools driver executes different sub-commands and pipelines.

    Usage:
        gemtools [--loglevel <log>] [--dry] [--show] <cmd> [<args>...]
        gemtools --help
        gemtools --version

    Options:
        <cmd>               The sub-command
        <args>              The sub-commands arguments
        --dry               Show an overview of the execution
        --show              Show the command line that will be executed
        -h, --help          Show this help message
        -v, --version       Show version information
        --loglevel <log>    Set the log level to one of
                            error|warn|info|debug
    """

    ## import all modules that contain commands
    import gem.production
    try:
        # load modules from plugin folders
        _load_plugins()

        from jip.vendor.docopt import docopt
        commands_string = """
        The following sub-commands are available:
        """
        cmds = []
        for k, v in cli._cmd_registry.iteritems():
            cmds.append("    %30s  %s" % (k, v.title))

        doc_string = dedent(gemtools.__doc__) + \
            dedent("\n".join([commands_string] + cmds))
        args = docopt(doc_string, options_first=True)

        if args["--version"]:
            print "GEMTools version %s" % __VERSION__
            exit(0)
        # set loglevel from environment
        gem.loglevel(os.getenv("_GT_LOGLEVEL", "ERROR"))

        if args["--loglevel"]:
            gem.loglevel(args["--loglevel"])

        cmd = args["<cmd>"]
        if not cmd in cli._cmd_registry:
            print >>sys.stderr, "gemtools command '%s' not found!" % cmd
            print >>sys.stderr, "See gemtools --help for a list of commands"
            sys.exit(1)

        jip_name = "gemtools_%s" % (cmd.replace("-", "_"))
        try:
            if os.getenv("_GT_EXEC"):
                # reset the options so the parser is triggered
                # again and the tool instance is properly initialized
                tool = jip.find(jip_name)
                tool._options = None
                tool.parse_args(args['<args>'])
                tool.run()
            else:
                if args['--dry'] or args['--show']:
                    tool = jip.find(jip_name)
                    # we handle --dry and --show separatly,
                    # create the jobs and call the show commands
                    try:
                        jobs = jip.create_jobs(tool, args=args['<args>'])
                    except jip.tools.ValidationError as va:
                        sys.stderr.write(str(va))
                        sys.stderr.write('\n')
                        sys.exit(1)

                    if args['--dry']:
                        jip.cli.show_dry(jobs, options=tool.options)
                    if args['--show']:
                        jip.cli.show_commands(jobs)
                    try:
                        jip.jobs.check_output_files(jobs)
                    except Exception as err:
                        print >>sys.stderr, "%s\n" % \
                            (jip.cli.colorize("Validation error!",
                                              jip.cli.RED))
                        print >>sys.stderr, str(err)
                        sys.exit(1)
                    return

                try:
                    _run_tool(jip_name, args['<args>'])
                except jip.ValidationError as va:
                    sys.stderr.write(str(va))
                    sys.stderr.write("\n")
                    sys.exit(1)
                except jip.ParserException as va:
                    sys.stderr.write(str(va))
                    sys.stderr.write("\n")
                    sys.exit(1)
                except Exception as va:
                    raise
        except CommandException as e:
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)
        except jip.ParserException as e:
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)
        except Exception as e:
            raise
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)
    except KeyboardInterrupt:
        exit(1)
    finally:
        pass


def _run_tool(jip_name, args):
    profile = jip.Profile()
    script = jip.find(jip_name)
    script.parse_args(args)
    try:
        threads = script.options['threads'].get(int)
        if threads > 1:
            profile.threads = threads
    except:
        pass
    jobs = jip.jobs.create_jobs(script, args=args, profile=profile)
    # assign job ids
    for i, j in enumerate(jobs):
        j.id = i + 1

    for exe in jip.jobs.create_executions(jobs):
        if exe.completed:
            print >>sys.stderr, colorize("Skipping", YELLOW), exe.name
        else:
            sys.stderr.write(colorize("Running", YELLOW) +
                             " {name:30} ".format(
                                 name=colorize(exe.name, BLUE)
                             ))
            sys.stderr.flush()
            start = datetime.now()
            success = jip.run_job(exe.job)
            end = timedelta(seconds=(datetime.now() - start).seconds)
            if success:
                print >>sys.stderr, colorize(exe.job.state, GREEN),\
                    "[%s]" % (end)
            else:
                print >>sys.stderr, colorize(exe.job.state, RED)
                sys.exit(1)


if __name__ == "__main__":
    gemtools()
