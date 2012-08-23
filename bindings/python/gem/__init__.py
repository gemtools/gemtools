#!/usr/bin/env python
"""Python wrapper around the GEM2 mapper that provides
ability to feed data into GEM and retreive the mappings"""
import os
import shutil
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
        self.line = None

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
                #logging.warn("Different sequence and quality sizes for : %s !! Trimming qualities to read length !" % (self.id))
                sizes = [len(qualities), len(self.sequence)]
                sizes.sort()
                qualities = self.qualities[:(sizes[0]-sizes[1])]
            else:
                raise ValueError("Different sequence and quality sizes for :\n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))

        return "@%s\n%s\n+\n%s" % (self.id, self.sequence, qualities)


    def sam_2_gem(self):
        """Extract minimal information from the sam line
        and create a GEM mapping that can be passed to the
        realigner. NOTE: Without realigning the
        mapping will be invalid !!"""
        s = self.line.split("\t")
        flag = int(s[1])
        chr = s[2]

        ## update read id if paired end
        if 0x1 & flag == 0x1:
            ## multiple segments
            if 0x40 & flag == 0x40:
                self.id += "/1"
            else:
                self.id += "/2"

        if chr == "*":
            ## unmapped read
            self.summary = 0
            self.mappings = '-'
        else:
            pos = int(s[3])
            strand = "+"
            if flag & 0x10 == 0x10:
                ## seq is reverser complement
                self.sequence = utils.reverseComplement(self.sequence)
                if self.qualities:
                    self.qualities = self.qualities[::-1]
                strand = '-'

            self.summary = "1"
            self.mappings = "%s:%s:%d:%d" %(chr, strand, pos, len(self.sequence))


    def length(self):
        return len(self.sequence)


def _prepare_index_parameter(index, gem_suffix=True):
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    file_name = index
    if not file_name.endswith(".gem"):
        file_name = file_name+".gem"
    if not os.path.exists(file_name):
        raise ValueError("Index file not found : %s"%file_name)

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
    count = 0
    try:
        for l in input:
            fifo.write(l.to_fastq())
            fifo.write("\n")
            count += 1
        fifo.close()
    except:
        fifo.close()
        ## kill the process
        os.unlink(inputfile)
        raise
    return inputfile, count


def _prepare_output(process, output=None, type="map", name="GEM", remove_after_iteration=False):
    """If output is not None, this waits for the process to
    finish and opens a ReadIterator
    on the output file using the given type.
    If output is not given, this opens a ReadIterator on the
    stdout stream of the process.
    """
    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process.wait() != 0:
            raise ValueError("%s execution failed!" % name)
        return files.open(output, type=type, process=process, remove_after_iteration=remove_after_iteration)
    else:
        ## running in async mode, return iterator on
        ## the output stream
        return files.open(process.stdout, type=type, process=process, remove_after_iteration=remove_after_iteration)


def validate_executables():
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
           max_big_indel_length=15,
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
          '--max-big-indel-length', str(max_big_indel_length),
          '-T', str(threads)
    ]

    ## run the mapper
    process = utils.run_tool(pa, input, output, name="GEM-Mapper")
    return _prepare_output(process, output, type="map", name="GEM-Mapper")


def splitmapper(input,
                index,
                output=None,
                mismatches=0.04,
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

    input_file, read_count = _write_sequence_file(input, tmpdir=tmpdir)
    read_count = max(1, read_count)
    (fifo, output_file) = tempfile.mkstemp(suffix=".map", prefix="splitmap_output", dir=tmpdir)

    pa = [executables['gem-rna-mapper'],
          '-I', index,
          '-i', input_file,
          '-o', output_file[:-4],
          '-q', quality,
          '-m', str(mismatches),
          '--min-split-size',str(min_split_size),
          '--refinement-step-size', str(refinement_step_size),
          '--matches-threshold', str(matches_threshold),
          '-T', str(min(threads, read_count))
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
        raise ValueError("GEM-Mapper execution failed, output file is : %s, tmp input was : %s" % (output_file, input_file))

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
                       mismatches=0.04,
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
                           mismatches=mismatches,
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
    return _prepare_output(process, output, type="map", name="GEM-Pair-aling")


def realign(input,
             index,
             output=None,
             threads=1,):
    index = _prepare_index_parameter(index, gem_suffix=False)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-r',
                  '-T', str(threads)
    ]
    process = utils.run_tool(validate_p, input, output, "GEM-Validate", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Validate")


def validate(input,
             index,
             output=None,
             validate_score="-s,-b,-i",
             validate_filter="2,25",
             threads=1,):
    index = _prepare_index_parameter(index, gem_suffix=False)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-v', '-r',
                  '-s', validate_score,
                  '-f', validate_filter,
                  '-T', str(max(threads - 2, 1))
    ]
    process = utils.run_tool(validate_p, input, output, "GEM-Validate", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Validate")



def score(input,
          index,
          output=None,
          scoring="+U,+u,-t,-s,-i,-a",
          threads=1):
    index = _prepare_index_parameter(index, gem_suffix=False)
    score_p = [executables['gem-map-2-map'],
               '-I', index,
               '-s', scoring,
               '-T', str(threads)
    ]
    process = utils.run_tool(score_p, input, output, "GEM-Score", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Score")


def validate_and_score(input,
          index,
          output=None,
          scoring="+U,+u,-t,-s,-i,-a",
          validate_score="-s,-b,-i",
          validate_filter="2,25",
          threads=1):
    validator = validate(input, index, None, validate_score, validate_filter, max(threads-2, 1))
    return score(validator, index, output, scoring, 1)


def gem2sam(input, index, output=None, single_end=False, compact=False, threads=1):
    index = _prepare_index_parameter(index, False)
    gem_2_sam_p = [executables['gem-2-sam'],
                   '-I', index,
                   '-T', str(threads)
    ]
    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")
    process = utils.run_tool(gem_2_sam_p, input, output, "GEM-2-SAM", utils.read_to_map)
    return _prepare_output(process, output, "sam", name="GEM-2-SAM")


def sam2bam(input, output=None, sorted=False, tmpdir=None):
    sam2bam_p = ['samtools', 'view', '-S', '-b', '-']
    tools = [sam2bam_p]
    out_name = output
    if sorted:
        # we can not pipe samtools sort :(
        if out_name is not None:
            if out_name.endswith('.bam'):
                out_name = out_name[:-4]
        else:
            (fifo, out_name) = tempfile.mkstemp(suffix="", prefix="bam_sort", dir=tmpdir)
        bam_sort = ['samtools', 'sort', '-', out_name]
        out_name = out_name + ".bam"
        tools.append(bam_sort)
    process = utils.run_tools(tools, input, output, "SAM-2-BAM", raw_stream=True)
    return _prepare_output(process, out_name, type="bam", name="SAM-2-BAM", remove_after_iteration=(output is None))


class merger(object):
    """Merge all mappings from the source files into the target file.
    The target file must contain all available mappings and the sort order
    of all files must be the same.

    target -- either an open file descriptor or a file name.
    source -- an open file descriptor to a single source file, a file name or a
              list of file names
    """
    def __init__(self, target, source):
        if target is None:
            raise ValueError("No target file specified")
        if source is None:
            raise ValueError("No source file specified")

        self.target = target
        self.source = source
        self.reads = []
        for x in self.source:
            self.reads.append(None)

    def __iter__(self):
        return self

    def __get_next_read(self, stream):
        try:
            return stream.next()
        except :
            return None

    def next(self):
        target_read = self.__get_next_read(self.target)
        if not target_read:
            raise StopIteration()
        ## read one line from each handle
        source_read = None
        result_read = target_read
        for i, h in enumerate(self.source):
            source_read = self.reads[i]
            if source_read is None:
                source_read = self.__get_next_read(h)
                self.reads[i] = source_read
            if not source_read:
                continue

            if target_read.id == source_read.id:
                mis = source_read.min_mismatches()
                if mis >= 0:
                    result_read = source_read  # matching mapping, replace
                    ## read next into cache
                self.reads[i] = self.__get_next_read(h)
        return result_read

    def merge(self, output):
        of = open(output, 'w')
        for read in self:
            of.write(str(read))
            of.write("\n")
        of.close()
        return files.open(output, type="map")

