#!/usr/bin/env python
"""Python wrapper around the GEM2 mapper that provides
ability to feed data into GEM and retreive the mappings"""
import subprocess
import sys
import logging


# set this to true to cpli qualities to
# read length and print a waring instead of raising an
# exception
import files
from gem import utils

_trim_qualities = False

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
        print self.summary, self.mappings
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
            ## do a sanity check for sequence and quality lengths
            if len(self.sequence) != len(self.qualities):
                if _trim_qualities:
                    logging.warn("Different sequence and quality sizes for : %s !! Trimming qualities to read length !" % (self.id))
                    sizes = [len(self.qualities), len(self.sequence)]
                    sizes.sort()
                    self.qualities = self.qualities[:(sizes[0]-sizes[1])]
                else:
                    raise ValueError("Different sequence and quality sizes for :\n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))
            if len(self.sequence) <= 0:
                raise ValueError("Sequence length is < 0 for : \n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))
            return "@%s\n%s\n+\n%s" % (self.id, self.sequence, self.qualities)

    def length(self):
        return len(self.sequence)


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
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if not index.endswith(".gem"):
        index = index + ".gem"

    ## set default values
    if quality is not None:
        quality = "offset-%d" % quality


    ## prepare the input
    pa = ['gem-mapper', '-I', index,
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
                output,
                index,
                junctions=0.02,
                junctions_file=None,
                filter=None,
                refinement_step_size=2,
                min_split_size=15,
                matches_threshold=100,
                splice_consensus=None,
                quality=None,
                threads=1,
                tmpdir=None):
    return splits.splitmapper(input,
                output,
                index,
                junctions,
                junctions_file,
                filter,
                refinement_step_size,
                min_split_size,
                matches_threshold,
                splice_consensus,
                quality,
                threads,
                tmpdir)


def stats(input, output=None):
    """Run the gemtools stats in the mapping input

    input -- mapping input file or file descriptor or generator function
             that produces mappings
    output -- output prefix
    """
    ## prepare the input
    input_generator = filter.prepare_input(input)
    pa = ['gemtools', '-s', '-p', output, '-']
    ## run the mapper
    utils.run_tool(input_generator, output, pa, "GEM-Stats")


#
#-b|--map-both-ends                       (default=false)
#     --min-insert-size <number>               (default=0)
#     --max-insert-size <number>               (default=1000)
#     -E <max_edit_distance>|<%_differences>   (default=0.08)
#     --min-matched-bases <number>|<%>         (default=0.80)
#     --extension-triggering-mismatches <number>|<%>
#                                              (default=0.00)
#     --max-extendable-matches <number>|'all'  (default=all)
#     --max-matches-per-extension <number>     (default=1)
#     --unique-pairing                         (default=false)

def pairalign(input, output, index,
           quality=None, quality_threshold=26,
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
    if quality is None:
        quality = "offset-33"

    ## prepare the input
    input_generator = filter.prepare_input(input)

    pa = ['gem-mapper',
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
    utils.run_tool(input_generator, output, pa, "GEM-pair-aligner")

    ## return a gen generator on the output
    return filter.gemoutput(output)





def score(input, index, output, threads=1):
    ## check the index
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if index.endswith(".gem"):
        index = index[:-4]

    output_fd = open(output, "w")
    validate_p = ['gem-map-2-map',
                '-I', index,
                '-v', '-r',
                '-s', '-s,-b,-i',
                '-f', '2,25',
                '-S',
                '-T', str(max(threads - 2, 1))
    ]
    score_p = ['gem-map-2-map',
               '-I', index,
               '-s', '+U,+u,-t,-s,-i,-a'
    ]

    print >> sys.stderr, " ".join(validate_p)
    print >> sys.stderr, " ".join(score_p)

    validate = subprocess.Popen(validate_p, stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE, bufsize=-1)
    score = subprocess.Popen(score_p, stdin=validate.stdout,
                             stdout=output_fd, bufsize=-1, close_fds=True)

    ## read from input and pipe to process
    for l in input:
        print >> validate.stdin, "\t".join(l)
    validate.stdin.close()

    if validate.wait() != 0 or score.wait() != 0:
        output_fd.close()
        raise ValueError("GEM map 2 map scoring failed")
    output_fd.close()
    return filter.gemoutput(output)


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


    gem_2_sam_p = ['gem-2-sam',
                '-I', index,
                '-T', str(gem_2_sam_threads)
    ]
    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    sam_2_bam_p = ['samtools',
               'view',
               '-S', '-b', '-'
    ]

    bam_sort_p = ['samtools', 'sort', '-', outname]

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
