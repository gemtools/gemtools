#!/usr/bin/env python
"""Gem tools utilities
"""

import subprocess
import sys


def run_tool(input, output, params, name=""):
    """Run a tool with the given configuration. This expects
    the input to be a generator function that returns lines
    passed to the toos stdin. Output must be a file name."""

    if not isinstance(output, basestring):
        raise ValueError("Output must be a string")

    print >> sys.stderr, "Starting %s :\n\t%s" % (name, " ".join(params))

    process = subprocess.Popen(params,
                stdin=subprocess.PIPE,
                bufsize=0)

    for l in input:
        try:
            process.stdin.write(l.to_sequence())
            process.stdin.write('\n')
        except:
            raise
    process.stdin.close()
    ret = process.wait()
    if ret != 0:
        raise ValueError("%s execution failed" % name)


def get_mismatches(strata):
    """Parse the mismatch string and return the number
    of missmatches of the first mapping found
    or return -1 if no mapping was found
    """
    ## get the mappings
    mismatches = -1
    if strata is None:
        return -1
    if strata in ["-", "*", "+"]:
        return -1
    for idx, s in enumerate(multisplit(strata, [':', '+'])):
        if int(s) > 0:
            mismatches = idx
            break
    return mismatches


def multisplit(s, seps):
    """Split a the string s using multiple separators"""
    res = [s]
    for sep in seps:
        s, res = res, []
        for seq in s:
            res += seq.split(sep)
    return res
