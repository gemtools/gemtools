#!/usr/bin/env python
"""Python wrapper around gemtools
In addition to gemtoosl support, the gem module
provides ways to start the GEM mapper directly.
"""
import os
import subprocess
import sys
import types
from  itertools import islice

import gemtools as gt

from . import filter
from . import junctions as gemjunctions
from . import splits
from . import utils


class Template(object):
    """One line in a gem map file that
    contains the tag, alignments, counters, map_relations
    """
    def __init__(self, gt_template=None):
        self._gt_template = gt_template
        self.tag = None
        self.alignments = []
        self.counters = []
        self.multimaps = None

        if(gt_template):
            self.tag = gt_template.tag
            # fill alignments ?
            self.alignments = [Alignment(self, gt.gt_template_get_block(gt_template, i))
                    for i in xrange(0,
                                    gt.gt_template_get_num_blocks(gt_template))]
            self.counters = [gt.gt_template_get_counter(gt_template, i + 1)
                    for i in xrange(0, gt.gt_template_get_num_counters(gt_template))]

    def __str__(self):
        if(self.tag):
            return "\n".join([ str(a) for a in self.alignments])
        else:
            return "Empty Template"


class Alignment(object):
    """Alignment"""
    def __init__(self, tmpl=None, gt_alignment=None):
        self.tag = None
        self.qualities = None
        self.read = None
        self.counters = []
        self.maps = []
        if gt_alignment is not None:
            ## initialize from alignment
            tag = gt.gt_alignment_get_tag(gt_alignment)
            if tag:
                self.tag = tag
            else:
                self.tag = tmpl.tag
            self.read = gt.gt_alignment_get_read(gt_alignment)
            self.qualities = gt.gt_alignment_get_qualities(gt_alignment)
            self.counters = [gt.gt_alignment_get_counter(gt_alignment, i + 1)
                    for i in xrange(0, gt.gt_alignment_get_num_counters(gt_alignment))]
            self.maps = [Mapping(self, gt.gt_alignment_get_map(gt_alignment, i))
                    for i in xrange(0, gt.gt_alignment_get_num_maps(gt_alignment))]

    def __str__(self):
        return "%s\t%s\t%s\t%s" % (self.tag, self.read, self.qualities, ":".join([str(s) for s in self.counters]))


class Mapping(object):
    def __init__(self, alignment=None, gt_map=None):
        self.name = None
        self.position = 0
        self.length = 0
        self.strand = "+"
        self.distance = 0
        self.score = 0
        self.mismatches = []
        self.next_map = None
        if gt_map:
            self.name = gt.gt_map_get_seq_name(gt_map)
            self.position = gt.gt_map_get_position(gt_map)
            strand = gt.gt_map_get_direction(gt_map)
            if strand == 1:
                self.strand = "-"
            self.length = gt.gt_map_get_base_length(gt_map)
            self.score = gt.gt_map_get_score(gt_map)

    def __str__(self):
        return "%s:%s%d:%d" % (self.name, self.strand, self.position, self.length)


class Mismatch(object):
    def __init__(self):
        self.type = None
        self.position = 0
        self.character = None
        self.size = 0


class Read(object):
    """Represets a single read. The read info covers
    its id, the raw sequence, qualities, the mapping summary
    and the actual mappings. Qualitites, the summary, and the mappings
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
        self.template = None
        self.alignments = None

    def from_template(self, template, alignment):
        """Initialize a Read from the given gemtools
        template"""
        self.id = gt.gt_template_get_tag(template)
        self.sequence = gt.gt_alignment_get_read(alignment)
        self.qualities = gt.gt_alignment_get_qualities(alignment)
        self.template = template
        self.alignment = alignment


    def get_mismatches(self):
        """Parse the mismatch string and return the number
        of missmatches of the first mapping found
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

    def __str__(self):
        """Returns GEM string representaton of the reads"""
        if self.qualities is None:
            pass
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
            return "@%s\n%s\n+\n%s" % (self.id, self.sequence, self.qualities)


def _from_fastq(fastq):
    """Generator function that return a Read
    created from a fastq file

    Arguments
    ---
    fastq - a string to a file or an open file descriptor or an iterable
    """
    if fastq is None:
        raise ValueError("No fastq input specified")
    if isinstance(fastq, basestring):
        ## open a file
        fastq = open(fastq, "r")

    read = Read()
    while True:
        fastq_lines = list(islice(fastq, 4))  # read in chunks of 4
        if not fastq_lines:
            break
        read.id = fastq_lines[0].rstrip()[1:]
        read.sequence = fastq_lines[1].rstrip()
        read.qualities = fastq_lines[3].rstrip()
        yield read
    try:
        fastq.close()
    except:
        pass  # ignore, exception in case fastq is not file descriptor


def _from_fasta(fastq):
    """Generator function that return a Read
    created from a fasta file

    Arguments
    ---
    fastq - a string to a file or an open file descriptor or an iterable
    """
    if fastq is None:
        raise ValueError("No fasta input specified")
    if isinstance(fastq, basestring):
        ## open a file
        fastq = open(fastq, "r")

    read = Read()
    while True:
        fastq_lines = list(islice(fastq, 2))  # read in chunks of 4
        if not fastq_lines:
            break
        read.id = fastq_lines[0].rstrip()[1:]
        read.sequence = fastq_lines[1].rstrip()
        yield read
    try:
        fastq.close()
    except:
        pass  # ignore, exception in case fastq is not file descriptor


def _from_map(input):
    """Generator function that return a Read
    created from a map file

    Arguments
    ---
    input - a string to a file or an open file descriptor
    """
    if input is None:
        raise ValueError("No map input specified")

    infile = None
    if isinstance(input, basestring):
        ## open a file
        infile = gt.gt_input_file_open(input, False)
    else:
        ## open stream
        infile = gt.gt_input_stream_open(input)

    py_tmpl = Template()
    template = gt.gt_template_new()
    map_input = gt.gt_buffered_map_input_new(infile)
    while True:
        value = gt.gtpy_fill_template(py_tmpl, template, map_input)
        if not value:
            break
        print "Yielding", py_tmpl, py_tmpl.tag
        yield py_tmpl

    gt.gt_buffered_map_input_close(map_input)
    gt.gt_input_file_close(infile)


def _open_file(input):
    """Generator that opens the given file and yields
    the lines. The function can handle gzipped files
    and uses zcat to read the content (extra process)
    """
    if not isinstance(input, basestring):
        raise ValueError("Zcat input is not a stirng")

    fd = None
    if input.endswith(".gz"):
        zcat = subprocess.Popen(["zcat", "-f", input],
                                    stdout=subprocess.PIPE)
        fd = zcat.stdout
    else:
        fd = open(input, 'r')
    return fd


def _create_input(input):
    """Check the input and return it if its a generator,
    otherwise open as file
    """
    if isinstance(input, types.GeneratorType):
        return input
    else:
        return gem_open(input)


def gem_open(input):
    """Generator function that opens the input and
    yields Read objects. Automatic type detection is
    done using the file extension. .map .fastq .fa and
    .fastq are supported file formats. If a generator function
    is passed, the content of the generator is
    yield

    ---
    Arguments:
    input - single string or a list of strings with the file names
    """
    inputs = input
    if isinstance(input, basestring):
        inputs = [input]

    for file in inputs:
        content = _open_file(file)
        gen = None
        # remove .gz from name
        if file.endswith(".gz"):
            file = file[:-3]
        if file.endswith(".map"):
            gen = _from_map(content)
        elif file.endswith(".fastq"):
            gen = _from_fastq(content)
        elif file.endswith(".fasta") or file.endswith(".fa"):
            gen = _from_fasta(content)
        else:
            raise ValueError("Unsupported file %s" % file)
        for read in gen:
            yield read


def trim(reads, left_trim=0, right_trim=0):
    """Generator function that trims reads"""
    for read in reads:
        if right_trim > 0:
            read.sequence = read.sequence[left_trim:-right_trim]
            if read.qualities is not None:
                read.qualities = read.qualities[left_trim:-right_trim]
        elif left_trim > 0:
            read.sequence = read.sequence[left_trim:]
            if  read.qualities is not None:
                read.qualities = read.qualities[left_trim:]
        yield read


def unmapped(reads, exclude=-1):
    """Yield only unmapped reads and reads
    that have only mappings with mismatches >= exclude
    """
    for read in reads:
        mis = read.get_mismatches()
        if mis < 0:
            yield read
        elif exclude > 0 and mis >= exclude:
            yield read


def mapper(input, output, index,
           mismatches=0.04, delta=0,
           quality=None, quality_threshold=26,
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

    input -- string with the input file or a file handle or a generator
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
    if quality is None:
        quality = "offset-33"

    out_prefix = output
    if out_prefix.endswith(".map"):
        out_prefix = out_prefix[:-4]

    ## create output directories
    base = os.path.dirname(os.path.abspath(output))
    if not os.path.exists(base):
        os.makedirs(base)

    ## prepare the input
    pa = ['gem-mapper', '-I', index,
         '-o', out_prefix,
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
    utils.run_tool(_create_input(input), output, pa, "GEM-Mapper")
    return _from_map(output)


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


def extract_junctions(annotation, output=None):
    """Extract junctions from given gtf annotation. This checks for
    and existing junctions file and returns it if it exists. Otherwise
    the annotation is parsed and junctions are extracted before the name of the
    file is returned
    """
    if output == None:
        output = os.path.basename(annotation) + ".junctions"
    if os.path.exists(output):
        return os.path.abspath(output)

    out_fd = open(output, 'wb')
    for j in gemjunctions.from_gtf(annotation):
        print >>out_fd, str(j)
    out_fd.close()
    return os.path.abspath(output)


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


def sam(input, index, output, single_end=False, compact=False, threads=1):
    ## check the index
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    if index.endswith(".gem"):
        index = index[:-4]

    output_fd = open(output, 'w')

    gem_2_sam_p = ['gem-2-sam',
                '-I', index,
                '-T', str(threads)
    ]
    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    sam_2_bam_p = ['samtools',
               'view', index,
               '-s', '+U,+u,-t,-s,-i,-a'
    ]
    print >> sys.stderr, " ".join(gem_2_sam_p)
    print >> sys.stderr, " ".join(sam_2_bam_p)

    gem_2_sam = subprocess.Popen(gem_2_sam_p, stdin=subprocess.PIPE,
                           stdout=output_fd, bufsize=-1)
    #sam_2_bam = subprocess.Popen(score_p, stdin=validate.stdout,
    #                         stdout=output_fd, bufsize=-1, close_fds=True)
    ## read from input and pipe to process
    for l in input:
        print >>gem_2_sam.stdin, "\t".join(l)
    gem_2_sam.stdin.close()

    if gem_2_sam.wait() != 0:
        output_fd.close()
        raise ValueError("GEM 2 sam failed")
    output_fd.close()
    return filter.gemoutput(output)


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
            while True:
                source_line = lines[i]
                if source_line is None:
                    source_line = h.readline()
                lines[i] = source_line
                if not source_line:
                    break
                ss = source_line.split("\t")
                sid = ss[0]
                if id == sid:
                    mis = utils.get_mismatches(ss[3])
                if mis >= 0:
                    mapping = source_line  # mattchin mapping, replace
                ## read next into cache
                lines[i] = h.readline()
                break
        output.write(mapping)

    ### close everything
    target.close()
    map(lambda x: x.close(), handles)
    output.close()
