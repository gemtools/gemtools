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
import pkg_resources
import splits
import filter as gemfilter
import gem.gemtools as gt

LOG_NOTHING = 1
LOG_STDERR = 2

log_output = LOG_NOTHING
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.WARN)

_trim_qualities = False
default_splice_consensus = [("GT", "AG"), ("CT", "AC")]
extended_splice_consensus = [("GT", "AG"), ("CT", "AC"),
    ("GC", "AG"), ("CT", "GC"),
    ("ATATC", "A."), (".T", "GATAT"),
    ("GTATC", "AT"), ("AT", "GATAC")
]

default_filter = "same-chromosome,same-strand"

## use the bundled executables
use_bundled_executables = True
## max mappings to replace mapping counts for + and ! summaries
__max_mappings = 999999999
## filter to work around GT-32 and #006 in gem-map-2-map 
__awk_filter = ["awk", "-F", "\t", '{if($4 == "*" || $4 == "-"){print $1"\t"$2"\t"$3"\t0\t"$5}else{if($4 == "!" || $4 == "+"){print $1"\t"$2"\t"$3"\t'+str(__max_mappings)+'\t"$5}else{print}}}']

class execs_dict(dict):
    """Helper dict that resolves bundled binaries"""
    def __getitem__(self, item):
        if use_bundled_executables and pkg_resources.resource_exists("gem", "gembinaries/%s"%item):
            f = pkg_resources.resource_filename("gem", "gembinaries/%s"%item)
            return f
        return dict.__getitem__(self, item)

## paths to the executables
executables = execs_dict({
    "gem-indexer": "gem-indexer",
    "gem-mapper": "gem-mapper",
    "gem-rna-mapper": "gem-rna-mapper",
    "gem-map-2-map": "gem-map-2-map",
    "gem-2-sam": "gem-2-sam",
    "samtools": "samtools",
    "gem-info": "gem-info",
    "splits-2-junctions": "splits-2-junctions",
    "gem-retriever": "gem-retriever",
    })


def loglevel(level):
    """Simple way to set the current log level globally for the root logger.
    Accepts either 'debug','info','warning', 'error'

    Log levels debug also ensures executable output is written to stderr
    """
    global log_output
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)
    logging.getLogger().setLevel(numeric_level)
    if level.upper() in ["DEBUG"]:
        log_output = LOG_STDERR


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
        self.type = None
        self.template = gt.Template()

    def fill(self, other):
        """Fill this read with the content of another read"""
        self.id = other.id
        self.sequence = other.sequence
        self.qualities = other.qualities
        self.summary = other.summary
        self.mappings = other.mappings
        self.line = other.line
        self.type = other.type

    def min_mismatches(self):
        """Parse the mismatch string and return the minimum number
        of mismatches of the first mapping found
        or return -1 if no mapping was found
        """
        ## get the mappings
        mismatches = -1
        if self.summary is None:
            return -1
        if self.summary in ["-", "*", "+", "!"]:
            return -1

        for idx, s in enumerate(utils.multisplit(self.summary, [':', '+'])):
            if int(s) > 0:
                mismatches = idx
                break
        return mismatches


    def get_maps(self):
        if self.summary == None or self.summary in ['-', '*']:
            return ([0], ["-"])
        elif self.summary in ['+', '!']:
            return ([__max_mappings], ["-"])

        sums = [int(x) for x in utils.multisplit(self.summary, [':', '+'])]
        maps = self.mappings.split(',')
        return (sums, maps)


    def _fill_template(self):
        return self.template.fill(self.line)


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
        ## add support for paired reads
        isPaired = self.sequence.find(" ") >=0
        sequence_split = None
        append1 = ""
        append2 = ""
        if isPaired:
            append1 = "/1"
            append2 = "/2"
            sequence_split = self.sequence.split(" ")

        if self.qualities is None:
            ## print fasta
            if isPaired:
                return ">%s%s\n%s\n>%s%s\n%s" % (self.id, append1, sequence_split[0], self.id, append2, sequence_split[1])
            return ">%s\n%s" % (self.id, self.sequence)
        else:
            return self.to_fastq()


    def to_fastq(self):
        """Convert to sequence format. FastQ or FastA depending
        on qualities"""
        ## add support for paired reads
        isPaired = self.sequence.find(" ") >=0
        sequence_split = None
        quality_split = None
        append1 = ""
        append2 = ""
        if isPaired:
            append1 = "/1"
            append2 = "/2"
            sequence_split = self.sequence.split(" ")
            if self.qualities is not None:
                quality_split = self.qualities.split(" ")
            else:
                quality_split = [None, None]
        if isPaired:
            r1 = self._to_fastq(self.id+append1, sequence_split[0], quality_split[0], _trim_qualities)
            r2 = self._to_fastq(self.id+append2, sequence_split[1], quality_split[1], _trim_qualities)
            return "%s\n%s" % (r1, r2)
        return "%s" % self._to_fastq(self.id, self.sequence, self.qualities, _trim_qualities)


    def _to_fastq(self, id, sequence, qualities, trim_qualities=False):
        """Convert to Fastq and fakes qualities if they are not present"""
        if len(sequence) <= 0:
            raise ValueError("Sequence length is < 0 for : \n%s\n%s\n%s" % (self.id, self.sequence, self.qualities))

        ## do a sanity check for sequence and quality lengths
        if qualities is None or len(qualities) == 0:
            ## fake the qualities
            qualities = '[' * len(sequence)

        if len(sequence) != len(qualities) and qualities is not None:
            if trim_qualities:
                ql = len(qualities)
                sl = len(sequence)
                if ql < sl:
                    # extend qualities
                    qualities = "%s%s" % (qualities, ('[' * (sl - ql)))
                else:
                    # cut qualities
                    qualities = qualities[:(sl - ql)]
            else:
                raise ValueError(
                    "Different sequence and quality sizes for :\n%s\n%s\n%s" % (id, sequence, qualities))
        return "@%s\n%s\n+\n%s" % (id, sequence, qualities)


    def length(self):
        if self.sequence.find(" ") >=0:
            return len(self.sequence.split(" ")[0])
        return len(self.sequence)


def _prepare_index_parameter(index, gem_suffix=True):
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    file_name = index
    if not file_name.endswith(".gem"):
        file_name = file_name + ".gem"
    if not os.path.exists(file_name):
        raise ValueError("Index file not found : %s" % file_name)

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
    if quality is not None and quality in ["offset-33", "offset-64", "ignore"]:
        return quality
    if quality is not None and quality not in ["none", "ignore"]:
        quality = "offset-%d" % quality
    else:
        quality = 'ignore'

    return quality


def _write_sequence_file(input, tmpdir=None):
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


def _prepare_output(process, output=None, type="map", name="GEM", remove_after_iteration=False, quality=None):
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
        return files.open(output, type=type, process=process, remove_after_iteration=remove_after_iteration, quality=quality)
    else:
        ## running in async mode, return iterator on
        ## the output stream
        return files.open(process.stdout, type=type, process=process, remove_after_iteration=remove_after_iteration, quality=quality)


def validate_executables():
    """Validate the gem executables
    """
    for exe, path in executables.items():
        path = executables[exe]
        exe_path = utils.which(executables[exe])
        found = exe_path is not None
        if found:
            print >> sys.stderr, "Executable '%s' (%s) : %s" % (exe, path, exe_path)
        else:
            print >> sys.stderr, "Executable '%s' (%s) : Unknown" % (exe, path)


def mapper(input, index, output=None,
           mismatches=0.04,
           delta=0,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=20,
           min_decoded_strata=2,
           min_matched_bases=0.80,
           max_big_indel_length=15,
           max_edit_distance=0.20,
           mismatch_alphabet="ACGT",
           trim=None,
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
    max_edit_distance -- max edit distance, 0.20 per default
    max_decoded_matches -- maximum decoded matches, defaults to 20
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 2 2
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    trim -- tuple or list that specifies left and right trimmings
    """

    ## check the index
    index = _prepare_index_parameter(index)

    if quality is None and isinstance(input, files.ReadIterator):
        quality = input.quality
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
          '--mismatch-alphabet', mismatch_alphabet,
          '-T', str(threads)
    ]

    if max_edit_distance > 0:
        pa.append("-e")
        pa.append("%s"%str(max_edit_distance))

    trim_c = [executables['gem-map-2-map'], '-c']
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    # workaround for GT-32 - filter away the !
    # buld list of tools
    tools = [pa, __awk_filter]
    if trim is not None:
        tools.append(trim_c)

    ## run the mapper
    process = None
    if len(tools) == 1:
        process = utils.run_tool(tools[0], input, output, "GEM-Mapper", utils.read_to_sequence)
    else:
        process = utils.run_tools(tools, input, output, "GEM-Mapper", utils.read_to_sequence)

    return _prepare_output(process, output, type="map", name="GEM-Mapper", quality=quality)


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
                mismatch_strata_delta=1,
                quality=33,
                trim=None,
                post_validate=True,
                mismatch_alphabet="ACGT",
                threads=1):
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
    if quality is None and isinstance(input, files.ReadIterator):
        quality = input.quality
    quality = _prepare_quality_parameter(quality)
    splice_cons = _prepare_splice_consensus_parameter(splice_consensus)

    pa = [executables['gem-rna-mapper'],
          '-I', index,
          '-q', quality,
          '-m', str(mismatches),
          '--min-split-size', str(min_split_size),
          '--refinement-step-size', str(refinement_step_size),
          '--matches-threshold', str(matches_threshold),
          '--strata-after-first', str(mismatch_strata_delta),
          '--mismatch-alphabet', mismatch_alphabet,
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
        pa.append("-c")
        pa.append(splice_cons)

    trim_c = [executables['gem-map-2-map'], '-c']
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    tools = [pa, __awk_filter]
    if trim is not None:
        tools.append(trim_c)

    ## run the mapper
    process = None
    original_output = output
    if post_validate:
        output = None
    if len(tools) == 1:
        process = utils.run_tool(tools[0], input, output, "GEM-Split-Mapper", utils.read_to_sequence)
    else:
        process = utils.run_tools(tools, input, output, "GEM-Split-Mapper", utils.read_to_sequence)

    splitmap_out = _prepare_output(process, output, type="map", name="GEM-Split-Mapper", quality=quality)

    if post_validate:
        return validate(splitmap_out, index, original_output, threads=threads)

    return splitmap_out


def extract_junctions(input,
                      index,
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
                      merge_with=None,
                      min_split=4,
                      max_split=2500000,
                      keep_short_indels=True,
                      tmpdir=None):
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
        post_validate=False,
        threads=threads)
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

    denovo_junctions = splits.extract_denovo_junctions(
        write_helper(splitmap),
        minsplit=min_split,
        maxsplit=max_split,
        sites=merge_with)

    if output is None:
        ## return delete iterator
        return (
        files.open(output_file, type="map", process=splitmapper, remove_after_iteration=True), denovo_junctions)
    else:
        return files.open(output_file, type="map", process=splitmapper), denovo_junctions


def pairalign(input, index, output=None,
              quality=33,
              quality_threshold=26,
              max_decoded_matches=20,
              min_decoded_strata=1,
              min_insert_size=0,
              max_insert_size=1000,
              max_edit_distance=0.30,
              min_matched_bases=0.80,
              max_extendable_matches=0,
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
    if quality is None and isinstance(input, files.ReadIterator):
        quality = input.quality

    quality = _prepare_quality_parameter(quality)

    pa = [executables['gem-mapper'],
          '-p',
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
          '--max-extensions-per-match', str(max_matches_per_extension),
          '-T', str(threads)
    ]
    if unique_pairing:
        pa.append("--unique-pairing")
    if map_both_ends:
        pa.append("--map-both-ends")

        ## run the mapper
    process = utils.run_tool(pa, input, output, "GEM-Pair-align", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Pair-align", quality=quality)


def realign(input,
            index,
            output=None,
            threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=False)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-r',
                  '-T', str(threads)
    ]
    quality = None
    if isinstance(input, files.ReadIterator):
        quality = input.quality
    process = utils.run_tool(validate_p, input, output, "GEM-Validate", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Validate", quality=quality)


def validate(input,
             index,
             output=None,
             validate_score=None, # "-s,-b,-i"
             validate_filter=None, # "2,25"
             threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=False)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-v', '-r',
                  '-T', str(max(threads, 1))
    ]

    if validate_score is not None:
        validate_p.extend(["-s", validate_score])

    if validate_filter is not None:
        validate_p.extend(['-f', validate_filter])

    quality = None
    if isinstance(input, files.ReadIterator):
        quality = input.quality

    process = utils.run_tool(validate_p, input, output, "GEM-Validate", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Validate", quality=quality)


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
    quality = None
    if isinstance(input, files.ReadIterator):
        quality = input.quality
    process = utils.run_tool(score_p, input, output, "GEM-Score", utils.read_to_map)
    return _prepare_output(process, output, type="map", name="GEM-Score", quality=quality)


def validate_and_score(input,
                       index,
                       output=None,
                       scoring="+U,+u,-t,-s,-i,-a",
                       validate_score="-s,-b,-i",
                       validate_filter="2,25",
                       threads=1):
    validator = validate(input, index, None, validate_score, validate_filter, max(threads / 2, 1))
    return score(validator, index, output, scoring, max(threads / 2, 1))


def gem2sam(input, index=None, output=None, single_end=False, compact=False, threads=1, quality=None, check_ids=True):
    if index is not None:
        index = _prepare_index_parameter(index, True)
    gem_2_sam_p = [executables['gem-2-sam'],
                   '-T', str(threads)
    ]
    if index is not None:
        gem_2_sam_p.extend(['-I', index])

    if quality is None and isinstance(input, files.ReadIterator):
        quality = input.quality

    quality = _prepare_quality_parameter(quality)
    if quality is not None:
        gem_2_sam_p.extend(["-q", quality])

    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    # GT-25 transform id's
    transform = utils.read_to_map
    if check_ids and not single_end:
        def t(read):
            split = read.id.split()
            if len(split) > 1:
                if split[1][0] in ("1", "2"):
                    read.id = "%s/%s" % (split[0], split[1][0])
                else:
                    raise ValueError("Unable to identify read id pair counter from %s" % read.id)
            return utils.read_to_map(read)
        transform = t

    process = utils.run_tool(gem_2_sam_p, input, output, "GEM-2-SAM", transform)
    return _prepare_output(process, output, "sam", name="GEM-2-SAM", quality=quality)


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
    quality = None
    if isinstance(input, files.ReadIterator):
        quality = input.quality

    return _prepare_output(process, out_name, type="bam", name="SAM-2-BAM", remove_after_iteration=(output is None), quality=quality)


def index(input, output, content="dna", threads=1):
    """Run teh gem-indexer on the given input. Input has to be the path
    to a single fasta file that contains the genome to be indexed.
    Output should be the path to the target index file. Note that
    the gem index has to end in .gem and the prefix is added if necessary and
    the returned path will always be the correct path to the index.

    The method checks for the existence of the target index file
    and will NOT override but exit silently without recreating the index!

    Returns the absolute path to the resulting index file
    """
    indexer_p = [
        executables['gem-indexer'],
        '-t', str(threads),
        '--content-type', content.lower()
    ]

    if isinstance(input, basestring):
        if not os.path.exists(input):
            raise ValueError("Indexer input file %s not found" % input)
        indexer_p.extend(["-i", input])
    else:
        raise ValueError("The indexer wrapper can not handle the input %s, pass a file or a list of files" % input )

    existing = output
    if existing[-4:] != ".gem": existing = "%s.gem" % existing
    if os.path.exists(existing):
        logging.warning("Index %s already exists, skipping indexing" % existing)
        return os.path.abspath(existing)


    # indexer takes the prefix
    if output[-4:] == ".gem":
        output = output[:-4]
    indexer_p.extend(['-o', output])

    # the indexer need the other indexer tools in PATH
    path="%s:%s" % (os.path.dirname(executables['gem-indexer']), os.getenv("PATH"))

    process = utils.run_tools([indexer_p], input=None, output=None, name="gem-indexer", raw_stream=True, path=path)
    if process.wait() != 0:
        raise ValueError("Error while executing the gem-indexer")
    return os.path.abspath("%s.gem" % output)


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
        self.result_read = Read()
        for x in self.source:
            self.reads.append(None)

    def __iter__(self):
        return self

    def __get_next_read(self, stream):
        try:
            return stream.next()
        except:
            return None

    def next(self):
        target_read = self.__get_next_read(self.target)
        if not target_read:
            raise StopIteration()
            ## read one line from each handle
        source_read = None
        self.result_read.fill(target_read)
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
                    self.result_read.fill(source_read)
                ## read next into cache
                self.reads[i] = self.__get_next_read(h)
        return self.result_read

    def merge(self, output):
        of = open(output, 'w')
        for read in self:
            of.write(str(read))
            of.write("\n")
        of.close()
        return files.open(output, type="map")
