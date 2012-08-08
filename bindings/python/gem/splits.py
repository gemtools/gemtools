#!/usr/bin/env python
#-t 200 --min-split-size 15 --refinement-step-size 2

import os
import subprocess
import types
import tempfile
import re

from . import filter as gemfilters
from threading import Thread
import gem

from gem.junctions import Exon, JunctionSite


def extract_denovo_junctions(gemoutput, index_hash, minsplit=4, maxsplit=2500000):
    splits2junctions_p = [
        'splits-2-junctions',
        str(minsplit),
        str(maxsplit)
    ]
    p = subprocess.Popen(splits2junctions_p, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True, bufsize=0)
    ## start the retriever
    retriever = subprocess.Popen(['gem-retriever', 'query', index_hash], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True, bufsize=0)

    ## start pipe thread
    input_thread = Thread(target=_pipe_geminput, args=(gemoutput, p))
    input_thread.start()

    ## read from process stdout and get junctions
    delta = 5
    sites = set([])
    for line in p.stdout:
        site = JunctionSite(line = line)
        que_don = __extract(delta, 1, site.descriptor[0], site.descriptor[1], site.descriptor[2])
        que_acc = __extract(delta, 0, site.descriptor[3], site.descriptor[4], site.descriptor[5])
        seq_don = __retrieve(retriever, que_don)
        seq_acc = __retrieve(retriever, que_acc)
        seq = seq_don+seq_acc
        if (re.search("GT......AG|GC......AG|ATATC...A.|GTATC...AT", seq)) or\
            (re.search("CT......AC|CT......GC|.T...GATAT|AT...GATAC", seq)):
            sites.add(site)


    ## wait for thread and process to finish
    input_thread.join()
    exit_value = p.wait()
    retriever.stdin.close()
    retriever.kill()
    if exit_value != 0:
        raise ValueError("Error while executing junction extraction")
    return sites


def __retrieve(retriever, query):
    retriever.stdin.write(query)
    retriever.stdin.write("\n")
    retriever.stdin.flush()
    result = retriever.stdout.readline().rstrip()
    return result

def __extract(delta, is_donor, chr, strand, pos):
    is_forw = strand == "+"
    #       return ((is_donor&&is_forw)||!is_forw&&!is_donor) ? chr"\t"str"\t"(pos+1)"\t"(pos+delta) :
    #                                                           chr"\t"str"\t"(pos-delta)"\t"(pos-1)

    if (is_donor and is_forw) or (not is_forw and not is_donor):
        return "%s\t%s\t%d\t%d"%(chr, strand, pos+1, pos+delta)
    else:
        return "%s\t%s\t%d\t%d"%(chr, strand, pos-delta, pos-1)



def _pipe_geminput(input, process):
    for read in input:
        process.stdin.write(str(read))
        process.stdin.write("\n")
    process.stdin.close()
