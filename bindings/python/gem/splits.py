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


default_splice_consensus = [("GT","AG"),("CT","AC")]
extended_splice_consensus = [("GT","AG"),("CT","AC"),
                             ("GC","AG"),("CT","GC"),
                             ("ATATC","A."),(".T","GATAT"),
                             ("GTATC", "AT"),("AT","GATAC")
                            ]

default_filter = "same-chromosome,same-strand"

#$GEM_RNA_MAPPER
# -t 200 <-- matches theshold
# --min-split-size 15
# --refinement-step-size 2
# -s \"GT\"+\"AG\",\"CT\"+\"AC\",\"GC\"+\"AG\",\"CT\"+\"GC\",\"ATATC\"+\"A.\",\".T\"+\"GATAT\",\"GTATC\"+\"AT\",\"AT\"+\"GATAC\"
# -f ordered
# -I $GEM_IDX
# -i $GEM_PREFIX.1.$GEM_EXT
# -q $GEM_OFFSET
# -o $GEM_PREFIX.1.de-novo
# -T 8 > $GEM_PREFIX.1.de-novo.log 2>&1

def junction_detection(input,
                       output,
                       index,
                       index_hash,
                       junctions=0.02,
                       junctions_file=None,
                       filter="ordered",
                       refinement_step_size=2,
                       min_split_size=15,
                       matches_threshold=200,
                       splice_consensus=extended_splice_consensus,
                       quality=None,
                       threads=1,
                       tmpdir=None):
    splitmap = splitmapper(input,
                output,
                index,
                junctions=junctions,
                junctions_file=junctions_file,
                filter=filter,
                refinement_step_size=refinement_step_size,
                min_split_size=min_split_size,
                matches_threshold=matches_threshold,
                splice_consensus=splice_consensus,
                quality=quality,
                threads=threads,
                tmpdir=tmpdir)

    ## extract the junctions
    denovo_junctions = extract_denovo_junctions(splitmap, index_hash)
    return (gemfilters.gemoutput(output), denovo_junctions)


def splitmapper(input,
                output,
                index,
                junctions=0.02,
                junctions_file=None,
                filter=default_filter,
                refinement_step_size=2,
                min_split_size=15,
                matches_threshold=100,
                splice_consensus=None,
                quality=None,
                threads=1,
                tmpdir=None):
    """Start the GEM split mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. In case of a
    .map file, the input is converted from gem format
    to fastq and passed to the mapper.

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- string with the input file or a file handle or a generator
    output -- output file name or file handle
    index -- valid GEM2 index
    """

    ## check the index
    if index is None:
      raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
      raise ValueError("GEM index must be a string")
    if index.endswith(".gem"):
      index = index[:-4]
    # check output
    if output is None:
      raise ValueError("No output specified!")
    if not isinstance(output, basestring):
      raise ValueError("GEM output must be a string")
    if output.endswith(".map"):
      output = output[:-4]

    ## set default values
    if quality is None:
        quality = "offset-33"
    if splice_consensus is None:
        splice_consensus = default_splice_consensus

    splice_cons = None
    if isinstance(splice_consensus, basestring):
        splice_cons = splice_consensus
    else:
        ## translate the splice consensus tupel structure
        splice_cons = ",".join(['"%s"+"%s"'%(x[0],x[1]) for x in splice_consensus])

    inputfile = input
    fifo = None
    if isinstance(input, basestring):
      inputfile = input
    elif isinstance(input, types.GeneratorType):
      (fifo, inputfile) = tempfile.mkstemp(suffix=".fastq", prefix="splitmap-input", dir=tmpdir)
      fifo = open(inputfile, 'w')
      ## the splitmapper does not support fifo or piping from stdin
      ## so we have to write a tmp file
      try:
        for l in input:
          fifo.write(l)
        fifo.close()
      except:
        ## kill the process
        os.unlink(inputfile)
        raise

    pa = ['gem-rna-mapper',
         '-I', index,
         '-i', inputfile,
         '-o', output,
         '-q', quality,
         '--min-split-size',str(min_split_size),
         '--refinement-step-size', str(refinement_step_size),
         '--matches-threshold', str(matches_threshold),
         '-T', str(threads)
    ]

    if junctions_file is not None:
        pa.append("-J")
        pa.append(os.path.abspath(junctions_file))
        pa.append("-j")
        pa.append(str(junctions))
    if filter is not None:
        pa.append("-f")
        pa.append(filter)
    if splice_cons is not None and junctions_file is None:
        pa.append("-s")
        pa.append(splice_cons)


    print " ".join(pa)
    p = subprocess.Popen(pa)
    ret = p.wait()

    ## delete tmp file
    if fifo is not None:
      os.unlink(inputfile)
    if ret != 0:
      raise ValueError("GEM Splitmapper execution failed")
    return gem.gem_open(output+".map")


def extract_denovo_junctions(gemoutput, index_hash, minsplit=4, maxsplit=2500000):
    splits2junctions_p = [
        'splits-2-junctions',
        str(minsplit),
        str(maxsplit)
    ]
    p = subprocess.Popen(splits2junctions_p, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    ## start pipe thread
    input_thread = Thread(target=_pipe_geminput, args=(gemoutput, p))
    input_thread.start()

    ## read from process stdout and get junctions
    sites = []
    for line in p.stdout:
        site = JunctionSite(line = line)
        sites.append(site)

    ## wait for thread and process to finish
    input_thread.join()
    if p.wait() != 0:
        raise ValueError("Error while executing junction extraction")

    sorted_sites = set(sites)

    sites = []
    ## start the retriever
    retriever = subprocess.Popen(['gem-retriever', 'query', index_hash], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    delta = 5

    for j in sorted_sites:
        que_don = __extract(delta, 1, j.descriptor[0], j.descriptor[1], j.descriptor[2])
        que_acc = __extract(delta, 0, j.descriptor[3], j.descriptor[4], j.descriptor[5])
        seq_don = __retrieve(retriever, que_don)
        seq_acc = __retrieve(retriever, que_acc)
        seq = seq_don+seq_acc
        if (re.search("GT......AG|GC......AG|ATATC...A.|GTATC...AT", seq)) or\
           (re.search("CT......AC|CT......GC|.T...GATAT|AT...GATAC", seq)):
            sites.append(j)
    ## stop retriever
    retriever.stdin.close()
    retriever.kill()

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



    # cat $GEM_PREFIX.1.de-novo.map |
#  $GEM_SPLITS_2_JUNCTIONS 4 2500000 |
# awk '
# function invert(str){
#   if (str=="+"||str=="F") return "-";
#   else if (str=="-"||str=="R") return "+";
#   else exit
# }
# {
#   if (int($3)>int($6)) print $4"/"invert($5)"/"$6"/"$1"/"invert(\$2)"/"$3;
#   else print $1"/"$2"/"$3"/"$4"/"$5"/"$6}
# ' | LC_ALL=C sort | uniq -c |
# awk '{print $2"\t"$1}' |
# awk '
#   function extract(delta,is_donor,chr,str,pos){
#       is_forw=(str=="+");
#       return ((is_donor&&is_forw)||!is_forw&&!is_donor) ? chr"\t"str"\t"(pos+1)"\t"(pos+delta) :
#                                                           chr"\t"str"\t"(pos-delta)"\t"(pos-1)
#   }
#   BEGIN{
#       retr="gem-retriever query '$GEM_IDX'.hash"
#   }{
#       split($1,s,"/");
#       delta=5;
#       que_don=extract(delta,1,s[1],s[2],s[3]);
#       que_acc=extract(delta,0,s[4],s[5],s[6]);
#       print que_don |& retr;
#       retr |& getline seq_don;
#       print que_acc |& retr;
#       retr |& getline seq_acc;
#       seq=seq_don seq_acc;
#       if (seq~/GT......AG|GC......AG|ATATC...A.|GTATC...AT/ || seq~/CT......AC|CT......GC|.T...GATAT|AT...GATAC/)
#           print s[1]"/"s[2]"/"s[3]"/"s[4]"/"s[5]"/"s[6]"\t"$2}' > $GEM_PREFIX.1.de-novo.junctions.compare