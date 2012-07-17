#!/usr/bin/env python
#-t 200 --min-split-size 15 --refinement-step-size 2

import os
import sys
import subprocess
import types
import tempfile
import re

from . import filter as gemfilters


default_splice_consensus=[("GT","AG"),("CT","AC")]

def splitmapper(input,
                output,
                index,
                junctions=0.02,
                junctions_file=None,
                filter="same-chromosome,same-strand",
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
    output -- output file name or fiel handle
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
        ## teranslate the splcie consensus tupel structure
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


    ## joint splice site consensnsus
    #",".join(['"%s"+"%s"'%(x[0],x[1]) for x in a])

    print " ".join(pa)
    p = subprocess.Popen(pa)
    ret = p.wait()
    ## delete tmp file
    if fifo is not None:
      os.unlink(inputfile)
    if ret != 0:
      raise ValueError("GEM Splitmapper execution failed")
    return gemfilters.gemoutput(output+".map")


