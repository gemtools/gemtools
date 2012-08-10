#!/usr/bin/env python
"""Python wrapper around the GEM2 mapper that provides
ability to feed data into GEM and retreive the mappings"""
import os
import shutil
import subprocess
import sys
import logging


# set this to true to cpli qualities to
# read length and print a waring instead of raising an
# exception
import tempfile
import files
from . import utils
import junctions as gemjunctions
import splits
import filter as gemfilter

_trim_qualities = False
default_splice_consensus = [("GT","AG"),("CT","AC")]
extended_splice_consensus = [("GT","AG"),("CT","AC"),
    ("GC","AG"),("CT","GC"),
    ("ATATC","A."),(".T","GATAT"),
    ("GTATC", "AT"),("AT","GATAC")
]

default_filter = "same-chromosome,same-strand"

## paths to the executables
executables = {
    "gem-mapper": "gem-mapper",
    "gem-rna-mapper": "gem-rna-mapper",
    "gem-map-2-map": "gem-map-2-map",
    "gem-2-sam": "gem-2-sam",
    "samtools": "samtools",
    "gem-info": "gem-info",
    "splits-2-junctions": "splits-2-junctions",
    "gem-retriever": "gem-retriever",
}

class Read(object):
    """A single read. The read info covers
    its id, the raw sequence, qualities, the mapping summary
    and the actual mappings. Qualities, the summary, and the mappings
    are optional.

    The read can be transformed to GEM input using the to_sequence
    method. The __str__ implementation returns the GEM representation.
    """
    def __init__(self):
        """Create a new empty read the is supposed to be used as a
        container
        """
        self.id = None
        self.sequence = None
        self.qualities = None
        self.summary = None
        self.mappings = None

    def min_mismatches(self):
        """Parse the mismatch string and return the minimum number
        of mismatches of the first mapping found
        or return -1 if no mapping was found
        """
        ## get the mappings
        mismatches = -1
        if self.summary is None:
            return -1
        if self.summary in ["-", "*", "+"]:
            return -1

        for idx, s in enumerate(utils.multisplit(self.summary, [':', '+'])):
            if int(s) > 0:
                mismatches = idx
                break
        return mismatches


    def get_maps(self):
        if self.summary == None or self.summary in ['-','+','*']:
            return ([0], ["-"])
        sums = [int(x) for x in utils.multisplit(self.summary, [':', '+'])]
        maps = self.mappings.split(',')
        return (sums, maps)


    def __str__(self):
        """Returns GEM string representaton of the reads"""
        if self.qualities is None:
            return "%s\t%s\t%s\t%s" % (self.id,
                                           self.sequence,
                                           self.summary,
                                           self.mappings)
        else:
            return "%s\t%s\t%s\t%s\t%s" % (self.id,
                                           self.sequence,
                                           self.qualities,
                                           self.summary,
                                           self.mappings)

    def to_sequence(self):
        """Convert to sequence format. FastQ or FastA depending
        on qualities"""
        if self.qualities is None:
            ## print fasta
            return ">%s\n%s" % (self.id, self.sequence)
        else:
            return self.to_fastq()

    def to_fastq(self):
        """Convert to Fastq and fakes qualities if they are not present"""
        if len(self.sequence) <= 0:
            raise ValueError("Sequence length is < 0 for : \n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))

        ## do a sanity check for sequence and quality lengths
        qualities = self.qualities
        if qualities is None or len(qualities) == 0:
            ## fake the qualities
            qualities = '['*len(self.sequence)

        if len(self.sequence) != len(qualities) and self.qualities is not None:
            if _trim_qualities:
                logging.warn("Different sequence and quality sizes for : %s !! Trimming qualities to read length !" % (self.id))
                sizes = [len(self.qualities), len(self.sequence)]
                sizes.sort()
                self.qualities = self.qualities[:(sizes[0]-sizes[1])]
            else:
                raise ValueError("Different sequence and quality sizes for :\n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))

        return "@%s\n%s\n+\n%s" % (self.id, self.sequence, qualities)

    def length(self):
        return len(self.sequence)


def _prepare_index_parameter(index, gem_suffix=True):
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if gem_suffix:
        if not index.endswith(".gem"):
            index = index + ".gem"
    else:
        if index.endswith(".gem"):
            index = index[:-4]
    return index


def _prepare_splice_consensus_parameter(splice_consensus):
    """
    Convert the splice consensus tuple to
    valid gem parameter input.
    If the given splice_consensus is None, the
    default splice consensus is used
    """
    if splice_consensus is None:
        splice_consensus = default_splice_consensus
    if isinstance(splice_consensus, basestring):
        splice_cons = splice_consensus
    else:
        ## translate the splice consensus tupel structure
        splice_cons = ",".join(['"%s"+"%s"' % (x[0], x[1]) for x in splice_consensus])
    return splice_cons


def _prepare_quality_parameter(quality):
    ## check quality
    if quality is not None:
        quality = "offset-%d" % quality
    else:
        quality = 'ignore'

    return quality


def _write_sequence_file(input, tmpdir = None):
    """Takes a Read sequence and writes it to a tempfile
    in fastq or fasta format
    """
    (fifo, inputfile) = tempfile.mkstemp(suffix=".fastq", prefix="mapper_input", dir=tmpdir)
    fifo = open(inputfile, 'w')
    ## the splitmapper does not support fifo or piping from stdin
    ## so we have to write a tmp file
    try:
        for l in input:
            fifo.write(l.to_fastq())
            fifo.write("\n")
        fifo.close()
    except:
        fifo.close()
        ## kill the process
        os.unlink(inputfile)
        raise
    return inputfile


def validate():
    """Validate the gem executables
    """
    for exe, path in executables.items():
        exe_path = utils.which(path)
        found = exe_path is not None
        if found:
            print >> sys.stderr, "Executable '%s' : %s" % (exe, exe_path)
        else:
            print >> sys.stderr, "Executable '%s' : Unknown" % (exe)


def mapper(input, index, output=None,
           mismatches=0.04,
           delta=0,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=20,
           min_decoded_strata=0,
           min_matched_bases=0.80,
           threads=1):
    """Start the GEM mapper on the given input.
    If input is a file handle, it is assumed to
    provide fastq entries. If input is a string,
    it is checked for its extension. In case of a
    .map file, the input is converted from gem format
    to fastq and passed to the mapper.

    Output can be a string, which will be translated to
    the output file. In case output is a file handle,
    the GEM output is written there.

    input -- A ReadIterator with the input
    output -- output file name
    index -- valid GEM2 index
    mismatches - number or % mismatches, default=0.04
    delta -- strata after best <number> (default=0)
    quality -- one of 'ignore'|'offset-33'|'offset-64' defaults to offset-33
    quality_threshold <number> -- (default=26, that is e<=2e-3)
    max_decoded_matches -- maximum decoded matches, defaults to 20
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 0
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    """

    ## check the index
    index = _prepare_index_parameter(index)

    quality = _prepare_quality_parameter(quality)

    ## prepare the input
    pa = [executables['gem-mapper'], '-I', index,
          '-q', quality,
          '-m', str(mismatches),
          '-s', str(delta),
          '--max-decoded-matches', str(max_decoded_matches),
          '--min-decoded-strata', str(min_decoded_strata),
          '--min-matched-bases', str(min_matched_bases),
          '--gem-quality-threshold', str(quality_threshold),
          '-T', str(threads)
    ]

    ## run the mapper
    process = utils.run_tool(pa, input, output, name="GEM-Mapper")
    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process.wait() != 0:
            raise ValueError("GEM-Mapper execution failed!")
        return files.open(output, type="map", process=process)
    else:
        ## running in async mode, return iterator on
        ## the output stream
        return files.open(process.stdout, type="map", process=process)



def splitmapper(input,
                index,
                output=None,
                junctions=0.02,
                junctions_file=None,
                splice_consensus=None,
                filter=default_filter,
                refinement_step_size=2,
                min_split_size=15,
                matches_threshold=100,
                quality=33,
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
    index = _prepare_index_parameter(index, False)
    quality = _prepare_quality_parameter(quality)
    splice_cons = _prepare_splice_consensus_parameter(splice_consensus)

    input_file = _write_sequence_file(input, tmpdir=tmpdir)
    (fifo, output_file) = tempfile.mkstemp(suffix=".map", prefix="splitmap_output", dir=tmpdir)

    pa = [executables['gem-rna-mapper'],
          '-I', index,
          '-i', input_file,
          '-o', output_file[:-4],
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


    ## run the split-mapper
    process = utils.run_tool(pa, None, None, name="GEM-Split-Mapper")

    exit_value = process.wait()
    ## cleanup
    os.remove(input_file)

    if exit_value != 0:
        raise ValueError("GEM-Mapper execution failed, output file is : %s!" % (output_file))

    if output is not None:
        ## move temp file to specified output
        shutil.move(output_file, output)
        return files.open(output, type="map", process=process)
    else:
        return files.open(output_file, type="map", process=process, remove_after_iteration=True)


def extract_junctions(input,
                       index,
                       index_hash,
                       output=None,
                       junctions=0.02,
                       filter="ordered",
                       refinement_step_size=2,
                       min_split_size=15,
                       matches_threshold=200,
                       splice_consensus=extended_splice_consensus,
                       quality=33,
                       threads=1,
                       tmpdir=None,
                       merge_with=None,
                       keep_short_indels=True):
    ## run the splitmapper
    splitmap = splitmapper(input,
                           index,
                           output=None,
                           junctions=junctions,
                           filter=filter,
                           refinement_step_size=refinement_step_size,
                           min_split_size=min_split_size,
                           matches_threshold=matches_threshold,
                           splice_consensus=splice_consensus,
                           quality=quality,
                           threads=threads,
                           tmpdir=tmpdir)
    ## make sure we have an output file
    ## for the splitmap results
    output_file = output
    if output is None:
        (fifo, output_file) = tempfile.mkstemp(suffix=".map", prefix="splitmap_output", dir=tmpdir)

    ## helper filter to pass the read on to the junction extractor,
    ## kill the splitmap as long as it is no short indel and
    ## write the result to the output file
    def write_helper(reads):
        of = open(output_file, 'w')
        for read in reads:
            yield read
            if keep_short_indels:
                gemfilter.keep_short_indel(read)
            ## and write the read
            of.write(str(read))
            of.write("\n")
        of.close()

    denovo_junctions = splits.extract_denovo_junctions(write_helper(splitmap), index_hash)
    if merge_with is not None:
        to_merge = [x for x in merge_with]
        to_merge.append(denovo_junctions)
        denovo_junctions = gemjunctions.merge_junctions(to_merge)

    if output is None:
        ## return delete iterator
        return ( files.open(output_file, type="map", process=splitmapper, remove_after_iteration=True), denovo_junctions)
    else:
        return ( files.open(output_file, type="map", process=splitmapper), denovo_junctions)


def pairalign(input, index, output=None,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=20,
           min_decoded_strata=0,
           min_insert_size=0,
           max_insert_size=1000,
           max_edit_distance=0.08,
           min_matched_bases=0.80,
           max_extendable_matches="all",
           max_matches_per_extension=1,
           unique_pairing=False,
           map_both_ends=False,
           threads=1):

    ## check the index
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if not index.endswith(".gem"):
        index = index + ".gem"

    ## set default values
    quality = _prepare_quality_parameter(quality)

    pa = [executables['gem-mapper'],
         '-I', index,
         '-q', quality,
         '--gem-quality-threshold', str(quality_threshold),
         '--max-decoded-matches', str(max_decoded_matches),
         '--min-decoded-strata', str(min_decoded_strata),
         '--min-insert-size', str(min_insert_size),
         '--max-insert-size', str(max_insert_size),
         '-E', str(max_edit_distance),
         '--min-matched-bases', str(min_matched_bases),
         '--max-extendable-matches', str(max_extendable_matches),
         '--max-matches-per-extension', str(max_matches_per_extension),
         '-T', str(threads)
    ]
    if unique_pairing:
        pa.append("--unique-pairing")
    if map_both_ends:
        pa.append("--map-both-ends")

        ## run the mapper
    process = utils.run_tool(pa, input, output, name="GEM-Pair-align")
    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process.wait() != 0:
            raise ValueError("GEM-Pair-align execution failed!")
        return files.open(output, type="map", process=process)
    else:
        ## running in async mode, return iterator on
        ## the output stream
        return files.open(process.stdout, type="map", process=process)



def score(input,
          index,
          output=None,
          threads=1,
          scoring="+U,+u,-t,-s,-i,-a",
          validate_score="-s,-b,-i",
          validate_filter="2,25"):
    index = _prepare_index_parameter(index, gem_suffix=False)
    validate_p = [executables['gem-map-2-map'],
                '-I', index,
                '-v', '-r',
                '-s', validate_score,
                '-f', validate_filter,
                '-T', str(max(threads - 2, 1))
    ]
    score_p = [executables['gem-map-2-map'],
               '-I', index,
               '-s', scoring
    ]

    process = utils.run_tools([validate_p, score_p], input, output, "GEM-Score", utils.read_to_map)
    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process.wait() != 0:
            raise ValueError("GEM-Score execution failed!")
        return files.open(output, type="map", process=process)
    else:
        ## running in async mode, return iterator on
        ## the output stream
        return files.open(process.stdout, type="map", process=process)



def sam(input, index, output, single_end=False, compact=False, bam=False, sort_bam=True, threads=1):
    ## check the index
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if index.endswith(".gem"):
        index = index[:-4]
    outname = output
    if outname.endswith(".sam") or outname.endswith(".bam"):
        outname = outname[:-4]

    gem_2_sam_threads = max(threads - 1, 1)
    if bam:
        gem_2_sam_threads = max(threads - 2, 1)
        if sort_bam:
            gem_2_sam_threads = max(threads - 3, 1)


    gem_2_sam_p = [executables['gem-2-sam'],
                '-I', index,
                '-T', str(gem_2_sam_threads)
    ]
    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    sam_2_bam_p = [executables['samtools'],
               'view',
               '-S', '-b', '-'
    ]

    bam_sort_p = [executables['samtools'], 'sort', '-', outname]

    print >> sys.stderr, " ".join([str(x) for x in gem_2_sam_p])
    if bam:
        print >> sys.stderr, " ".join([str(x) for x in sam_2_bam_p])
        if sort_bam:
            print >> sys.stderr, " ".join([str(x) for x in bam_sort_p])


    # prepare outputs
    gem_2_sam_out = None
    if bam:
        gem_2_sam_out = subprocess.PIPE
    else:
        gem_2_sam_out = open('w', output)


    # start gem 2 stam conversion
    gem_2_sam = subprocess.Popen(gem_2_sam_p, stdin=subprocess.PIPE,
                           stdout=gem_2_sam_out, bufsize=-1)
    gem_2_sam_out = gem_2_sam.stdout
    sam_2_bam = None
    bam_sort = None
    if bam:
        sam_2_bam_out = subprocess.PIPE
        if not sort_bam:
            sam_2_bam_out = open(output,  "w")

        sam_2_bam = subprocess.Popen(sam_2_bam_p, stdin=gem_2_sam_out,
                             stdout=sam_2_bam_out, bufsize=-1, close_fds=True)
        if sort_bam:
            bam_sort = subprocess.Popen(bam_sort_p, stdin=sam_2_bam.stdout, close_fds=True)


    ## read from input and pipe to process
    for l in input:
        print >>gem_2_sam.stdin, "\t".join(l)

    gem_2_sam.stdin.close()
    if gem_2_sam.wait() != 0:
        gem_2_sam.stdout.close()
        raise ValueError("GEM 2 sam failed")
    gem_2_sam.stdout.close()

    if sam_2_bam is not None:
        if sam_2_bam.wait() != 0:
            sam_2_bam.stdout.close()
            raise ValueError("Sam to Bam conversion failed")
        sam_2_bam.stdout.close()
        if bam_sort is not None:
            if bam_sort.wait() != 0:
                raise ValueError("Bam sorting failed")


def merge(target, source, output):
    """Merge all mappings from the source files into the target file.
    The target file must contain all available mappings and the sort order
    of all files must be the same.

    target -- either an open file descriptor or a file name.
    source -- an open file descriptor to a single source file, a file name or a
              list of file names
    output -- either an open file descriptor that is used to write the output
              or the target file name
    """
    if target is None:
        raise ValueError("No target file specified")
    if source is None:
        raise ValueError("No source file specified")
    if output is None:
        raise ValueError("No output file specified")
    if isinstance(target, basestring):
        target = open(target, 'r')
    if isinstance(output, basestring):
        output = open(output, 'w')

    handles = [source]
    if isinstance(source, basestring):
        handles = [open(source, 'r')]
    elif isinstance(source, (list, tuple)):
        handles = [open(s, 'r') for s in source]

    lines = []
    for x in handles:
        lines.append(None)
    ## merge the files
    for mapping in target:
        id = mapping.split("\t")[0]
        ## read one line from each handle
        for i, h in enumerate(handles):
            source_line = lines[i]
            if source_line is None:
                source_line = h.readline()
                lines[i] = source_line
            if not source_line:
                continue

            ss = source_line.split("\t")
            sid = ss[0]
            if id == sid:
                mis = utils.get_mismatches(ss[3])
                if mis >= 0:
                    mapping = source_line  # mattchin mapping, replace
                ## read next into cache
                lines[i] = h.readline()
        output.write(mapping)

    ### close everything
    target.close()
    map(lambda x: x.close(), handles)
    output.close()
