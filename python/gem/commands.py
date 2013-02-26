#!/usr/bin/env python
"""Gemtools commands"""
import argparse
import gem
import gem.production
import sys

__VERSION__ = "1.6"


def gemtools():
    try:
        parser = argparse.ArgumentParser(prog="gemtools",
                description="Gemtools driver to execute different gemtools command and pipelines"
                )
        parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __VERSION__)

        commands = {
            "index": gem.production.Index,
            "rna-pipeline": gem.production.RnaPipeline,
            "t-index": gem.production.TranscriptIndex,
            "merge": gem.production.Merge,
            "gtf-junctions": gem.production.Junctions,
            "denovo-junctions": gem.production.JunctionExtraction,
            "stats": gem.production.Stats
        }
        instances = {}

        subparsers = parser.add_subparsers(title="commands", metavar="<command>", description="Available commands", dest="command")
        for name, cmdClass in commands.items():
            p = subparsers.add_parser(name, help=cmdClass.title, description=cmdClass.description)
            instances[name] = cmdClass()
            instances[name].register(p)

        args = parser.parse_args()
        if args.loglevel is not None:
            gem.loglevel(args.loglevel)
        try:
            instances[args.command].run(args)
        except gem.utils.CommandException, e:
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)
    except KeyboardInterrupt:
        exit(1)
    finally:
        # cleanup
        gem.files._cleanup()
if __name__ == "__main__":
    gemtools()
