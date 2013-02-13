#!/usr/bin/env python
import re
import logging
import gem
from gem.junctions import JunctionSite


def extract_denovo_junctions(input, minsplit=4, maxsplit=2500000, sites=None, coverage=0, max_junction_matches=1, process=None, threads=1):
    """Extract denovo junctions from a split map run.

    gemoutput - a read iterator over gem splitmapper Output
    minsplit  - minimum length of a split to be considered, default is 4
    maxsplit  - maximum length of a split to be considered, defauilt is 2500000
    sites     - target set where the found junctions will be
                merged into, default is None
    coverage  - if > 0, a junction must be found > coverage times
                to be considered, default is 0 and therefore disabled
    """
    splits2junctions_p = [
        gem.executables['splits-2-junctions'],
        str(minsplit),
        str(maxsplit),
        str(max_junction_matches),
        str(threads)
    ]
    p = gem.utils.run_tool(splits2junctions_p, input=input, write_map=True)

    ## read from process stdout and get junctions
    if sites is None:
        sites = set([])
    initial_size = len(sites)
    local_sites = {}
    for line in p.stdout:
        js = JunctionSite(line=line)
        if coverage > 0:
            if js.hash not in local_sites:
                local_sites[js.hash] = js
                local_sites[js.hash].coverage = 1
            else:
                local_sites[js.hash].coverage += 1
        else:
            sites.add(js)

    if coverage > 0:
        for i, e in local_sites.items():
            if e.coverage >= coverage:
                sites.add(e)

    exit_value = p.wait()
    if process is not None:
        process.wait()
    if exit_value != 0:
        logging.error("Error while executing junction extraction")
        exit(1)
    logging.info("Junction extraction: Initial %d, Output %d (%d new)" % (initial_size, len(sites), (len(sites) - initial_size)))
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
        return "%s\t%s\t%d\t%d" % (chr, strand, pos + 1, pos + delta)
    else:
        return "%s\t%s\t%d\t%d" % (chr, strand, pos - delta, pos - 1)


def _pipe_geminput(input, process):
    for read in input:
        ## avoid printing max reads
        process.stdin.write(read.to_map(no_max_mappings=True))
        process.stdin.write("\n")
    process.stdin.close()


class append_xs_filter(object):
    ## for now, statically implemented junction sites
    # forward strand junctions
    forward = ["GT", "GC", "ATATC", "GTATC"]
    # reverse strand junctions
    reverse = ["AG", "AG", "A.", "AT"]

    # forward strand reverse complements
    forwardc = ["CT", "CT", ".T", "AT"]
    # reverse strand reverse complements
    reversec = ["AC", "GC", "GATAT", "GATAC"]

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
        if input[0] == "@":
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
                left, right = self.retriever.get_junction(split[2], "+", start + sum, int(i[:-1]))

                strand = None

                ## we check both forward and
                ## reverse complement junctions here
                ## to make sure this is unique

                #check normal forward
                for j, f in enumerate(append_xs_filter.forward):
                    if re.search("^" + f, left):
                        if re.search(append_xs_filter.reverse[j] + "$", right):
                            strand = "+"
                            strandplus += 1

                #check reverse complement forward
                for j, f in enumerate(append_xs_filter.forwardc):
                    if re.search("^" + f, left):
                        if re.search(append_xs_filter.reversec[j] + "$", right):
                            strand = "-"
                            strandminus += 1
                # set XS to 0 for unknown/non unique site
                #if (strandplus > 0 and strandminus == 0) or (strandplus == 0 and strandminus > 0):
                if (strandplus == 0 and strandminus == 0) or (strandplus > 0 and strandminus > 0):
                    strand = "0"
                split.append("XS:A:" + strand)
                return "\t".join(split) + "\n"
