#!/usr/bin/env python
"""Gemtools commands"""
from argparse import ArgumentParser
import gem.utils

import sys
from sys import exit
import types
from textwrap import dedent

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
            tool.parse_args(args["<args>"])
            tool.validate()
            if args['--dry'] or args['--show']:
                jip.run(tool, dry=args['--dry'], show=args['--show'])
            else:
                tool.run()
        except Exception as e:
            raise
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)

        #instance = _cmd_registry[cmd]()
        #try:
            #cmd_args = instance.parse_args([instance.name] + args["<args>"])
            #instance.run(cmd_args)
        #except gem.utils.CommandException, e:
            #sys.stderr.write("%s\n" % (str(e)))
            #exit(1)


        #parser = ArgumentParser(prog="gemtools",
                                #description="Gemtools driver to execute "
                                            #"different gemtools command and "
                                            #"pipelines")
        #parser.add_argument('--loglevel',
                            #default=None,
                            #help="Log level (error, warn, info, debug)")
        #parser.add_argument('-v', '--version',
                            #action='version',
                            #version='%(prog)s ' + __VERSION__)

        #commands = {
            #"index": gem.production.Index,
            #"rna-pipeline": gem.production.RnaPipeline,
            #"t-index": gem.production.TranscriptIndex,
            #"gtfcount": gem.production.GtfCount,
            #"merge": gem.production.Merge,
            #"convert": gem.production.Convert,
            #"gtf-junctions": gem.production.Junctions,
            #"denovo-junctions": gem.production.JunctionExtraction,
            #"stats": gem.production.Stats,
            #"filter": gem.production.Filter,
            #"report": gem.production.StatsReport,
            #"prepare": gem.production.PrepareInput,
            #"map": gem.production.GemMapper
        #}
        #instances = {}

        #if _cmd_registry is not None:
            #subparsers = parser.add_subparsers(
                #title="commands",
                #metavar="<command>",
                #description="Available commands",
                #dest="command"
            #)
            #for name, cmdClass in _cmd_registry.iteritems():
                #p = subparsers.add_parser(
                    #name,
                    #help=cmdClass.title,
                    #description=cmdClass.description
                #)
                #instances[name] = cmdClass()
                #instances[name].register(p)

        #args = parser.parse_args()
        #if args.loglevel is not None:
            #gem.loglevel(args.loglevel)

        #try:
            #instances[args.command].run(args)
        #except gem.utils.CommandException, e:
            #sys.stderr.write("%s\n" % (str(e)))
            #exit(1)
    except KeyboardInterrupt:
        exit(1)
    finally:
        pass

if __name__ == "__main__":
    gemtools()
