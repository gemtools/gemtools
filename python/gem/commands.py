#!/usr/bin/env python
"""Gemtools commands"""

import sys
from sys import exit
import types
from textwrap import dedent
import os

import gem.utils
from .utils import CommandException

import jip

__VERSION__ = "1.7"


_cmd_registry = None


class cli(object):
    """Command decorator"""
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
        global _cmd_registry
        if _cmd_registry is None:
            _cmd_registry = {}

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
        _cmd_registry[self.name] = cls

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
                 argparse='register' if not wrapped_function else None,
                 get_command='jip_command')(cls)
        return cls


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
        from jip.vendor.docopt import docopt
        commands_string = """
        The following sub-commands are available:
        """
        cmds = []
        for k, v in _cmd_registry.iteritems():
            cmds.append("    %30s  %s" % (k, v.title))

        doc_string = dedent(gemtools.__doc__) + \
            dedent("\n".join([commands_string] + cmds))
        args = docopt(doc_string, options_first=True)

        if args["--version"]:
            print "GEMTools version %s" % __VERSION__
            exit(0)

        if args["--loglevel"]:
            gem.loglevel(args["--loglevel"])

        cmd = args["<cmd>"]
        if not cmd in _cmd_registry:
            print >>sys.stderr, "gemtools command '%s' not found!" % cmd
            print >>sys.stderr, "See gemtools --help for a list of commands"
            sys.exit(1)

        jip_name = "gemtools_%s" % (cmd.replace("-", "_"))
        try:
            tool = jip.find(jip_name)
            if os.getenv("_GT_EXEC"):
                tool.parse_args(args['<args>'])
                tool.run()
            else:
                jip.run(tool, args['<args>'],
                        dry=args['--dry'], show=args['--show'])
        except CommandException as e:
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

if __name__ == "__main__":
    gemtools()
