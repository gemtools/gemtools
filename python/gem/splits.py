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


def extract_denovo_junctions(gemoutput, minsplit=4, maxsplit=2500000, sites=None):
    splits2junctions_p = [
        gem.executables['splits-2-junctions'],
        str(minsplit),
        str(maxsplit)
    ]
    p = subprocess.Popen(splits2junctions_p, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True, bufsize=0)

    ## start pipe thread
    input_thread = Thread(target=_pipe_geminput, args=(gemoutput, p))
    input_thread.start()

    ## read from process stdout and get junctions
    if sites is None:
        sites = set([])

    for line in p.stdout:
        sites.add(JunctionSite(line = line))

    exit_value = p.wait()
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
        ## avoid printing max reads
        process.stdin.write(read.to_map(no_max_mappings=True))
        process.stdin.write("\n")
    process.stdin.close()



class append_xs_filter(object):
    ## for now, statically implemented junction sites
    # forward strand junctions
    forward = ["GT","GC","ATATC","GTATC"]
    # reverse strand junctions
    reverse = ["AG","AG","A.","AT"]

    # forward strand reverse complements
    forwardc = ["CT","CT",".T","AT"]
    # reverse strand reverse complements
    reversec = ["AC","GC","GATAT","GATAC"]

    # pattern to parse SAM cigar
    pat = re.compile("(\d+[A-Z=])")
    # pattern to identify split maps
    split_re = re.compile(".*N.*")

    def __init__(self, index_hash):
        """Create a new filter providing the .hash file
        of the genome reference
        """
        self.retriever = gem.utils.retriever(index_hash)

    def filter(self, input):
        if input[0] == "@" :
            return input

        line = input.rstrip()
        split = line.split("\t")
        start = int(split[3])
        sum = 0
        strandplus = 0
        strandminus = 0
        ## check for split map
        if not append_xs_filter.split_re.search(split[5]):
            return input

        # found split, check junction site
        for i in append_xs_filter.pat.findall(split[5]):
            if i[-1] == "N":
                ## get junction sequence
                ## todo : why only + strand ? because sam contains reverse complement in case of - strand ?
                left, right = self.retriever.get_junction(split[2], "+", start+sum, int(i[:-1]))

                strand = None

                ## we check both forward and
                ## reverse complement junctions here
                ## to make sure this is unique

                #check normal forward
                for j,f in enumerate(append_xs_filter.forward):
                    if re.search("^"+f, left):
                        if re.search(append_xs_filter.reverse[j]+"$", right ):
                            strand = "+"
                            strandplus += 1

                #check reverse complement forward
                for j,f in enumerate(append_xs_filter.forwardc):
                    if re.search("^"+f, left):
                        if re.search(append_xs_filter.reversec[j]+"$", right ):
                            strand = "-"
                            strandminus += 1
                # set XS to 0 for unknown/non unique site
                #if (strandplus > 0 and strandminus == 0) or (strandplus == 0 and strandminus > 0):
                if (strandplus == 0 and strandminus == 0) or (strandplus > 0 and strandminus > 0):
                    strand = "0"
                split.append("XS:A:"+strand)
                return "\t".join(split)+"\n"
