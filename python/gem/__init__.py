#!/usr/bin/env python
"""Python wrapper around the GEM2 mapper that provides
ability to feed data into GEM and retrieve the mappings"""
import os
import sys
import logging

import tempfile
import files
from . import utils
import pkg_resources
import splits
import filter as gemfilter
import gem.gemtools as gt
import gem

LOG_NOTHING = 1
LOG_STDERR = 2

# default logger configuration
log_output = LOG_NOTHING
logging.basicConfig(format='%(asctime)-15s %(levelname)s: %(message)s', level=logging.WARN)


default_splice_consensus = [("GT", "AG"), ("CT", "AC")]
extended_splice_consensus = [("GT", "AG"), ("CT", "AC"),
    ("GC", "AG"), ("CT", "GC"),
    ("ATATC", "A."), (".T", "GATAT"),
    ("GTATC", "AT"), ("AT", "GATAC")
]

#default filter
default_filter = "same-chromosome,same-strand"

## use the bundled executables
use_bundled_executables = True

## max mappings to replace mapping counts for + and ! summaries
_max_mappings = 999999999

## filter to work around GT-32 and #006 in gem-map-2-map
__awk_filter = ["awk", "-F", "\t", '{if($4 == "*" || $4 == "-"){print $1"\t"$2"\t"$3"\t0\t"$5}else{if($4 == "!" || $4 == "+"){print $1"\t"$2"\t"$3"\t' + str(_max_mappings) + '\t"$5}else{print}}}']


class execs_dict(dict):
    """Helper dictionary that resolves bundled binaries
    based on the configuration. We check first for GEM_PATH
    environment variable. If its set, and points to a directory with
    the executable, the path to that executable is returned.
    Next, we check the use_bundled_executable flag. If that is true(default)
    the path to the bundled executable is returned.
    If nothing is found, the plain executable name is returned and we
    assume it can be found in PATH

    """
    def __getitem__(self, item):
        # check if there is an environment variable set
        # to specify the path to the GEM executables
        base_dir = os.getenv("GEM_PATH", None)
        if base_dir is not None:
            file = "%s/%s" % (base_dir, item)
            if os.path.exists(file):
                logging.debug("Using binary from GEM_PATH : %s" % file)
                return file

        if use_bundled_executables and pkg_resources.resource_exists("gem", "gembinaries/%s" % item):
            f = pkg_resources.resource_filename("gem", "gembinaries/%s" % item)
            logging.debug("Using bundled binary : %s" % f)
            return f
        logging.debug("Using binary from PATH: %s" % item)
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
    "compute-transcriptome": "compute-transcriptome",
    "transcriptome-2-genome": "transcriptome-2-genome",
    })


def loglevel(level):
    """Simple way to set the current log level globally for the root logger.
    Accepts either 'debug','info','warning', 'error'

    Log levels debug also ensures executable output is written to stderr

    level -- one of debug, info, warn, error
    """
    global log_output
    numeric_level = level
    if isinstance(level, basestring):
        numeric_level = getattr(logging, level.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(level=numeric_level)
    logging.getLogger().setLevel(numeric_level)


def _prepare_index_parameter(index, gem_suffix=True):
    """Prepares the index file and checks that the index
    exists. The function throws a IOError if the index file
    can not be found.

    index      -- the path to the index file
    gem_suffix -- if true, the function ensures that the index ends in .gem,
                  otherwise, it ensures that the .gem suffix is removed.

    """
    if index is None:
        raise ValueError("No valid GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM index must be a string")
    file_name = index

    if not file_name.endswith(".gem"):
        file_name = file_name + ".gem"

    if not os.path.exists(file_name):
        raise IOError("Index file not found : %s" % file_name)

    if gem_suffix:
        if not index.endswith(".gem"):
            index = index + ".gem"
    else:
        if index.endswith(".gem"):
            index = index[:-4]
    return index


def _prepare_splice_consensus_parameter(splice_consensus):
    """Convert the splice consensus tuple to
    valid gem parameter input.

    If the given splice_consensus is None, the
    default splice consensus is used

    splice_consensus -- if string, it is just passed on, if its a list of tuples like ("A", "G") it
                        is translated into gem representation
    """
    if splice_consensus is None:
        splice_consensus = default_splice_consensus
    if isinstance(splice_consensus, basestring):
        splice_cons = splice_consensus
    else:
        ## translate the splice consensus tupel structure
        splice_cons = ",".join(['"%s"+"%s"' % (x[0], x[1]) for x in splice_consensus])
    return splice_cons


def _prepare_quality_parameter(quality, input=None):
    """Prepare and returnn the quality parameter for gem runs. If the
    input is None, qualities are disabled, if it is a string and
    valid parameter it is returned as is. Otherwise it as to be 33
    or 64 and will be translated to a valid gem parameter

    quality -- the quality offset, 33|64, None to ignore or a string that
                is either offset-33|offset-64|ignore
    input   -- optional input that can be checkd for a quality parameter. This
               currently works only for gemtools.TemplateIterator or gemtools.InputFile
               instances
    """
    if quality is None and isinstance(input, (gt.TemplateIterator, gt.InputFile)):
        quality = input.quality

    ## check quality
    if quality is not None and quality in ["offset-33", "offset-64", "ignore"]:
        return quality
    if quality is not None and quality not in ["none", "ignore"]:
        i = int(quality)
        if i not in [33, 64]:
            raise ValueError("%s is not a valid quality value, try None, 33 or 64" % (str(quality)))
        quality = "offset-%d" % int(quality)
    else:
        quality = 'ignore'

    return quality


def _prepare_output(process, output=None, quality=None, delete_after_iterate=False):
    """Creates a new gem.gemtools.Inputfile from the given process.
    If output is specivied, the function blocks and waits for the process to finish
    successfully before the InputFile is created on the specified output.

    Otherwise, a stream based InputFile is created using the process stdout
    stream.

    If quality is specified it is passed on to the input file.

    process -- the process
    output  -- output file name or None when process stdout should be used
    quality -- optional quality that is passed to the InputFile
    """
    if output is not None:
        # we are writing to a file
        # wait for the process to finish
        if process is not None and process.wait() != 0:
            raise ValueError("Execution failed!")
        logging.debug("Opening output file %s" % (output))
        if output.endswith(".bam"):
            return gt.InputFile(stream=gem.files.open_bam(output), file_name=output, quality=quality, process=process, delete_after_iterate=delete_after_iterate)
        return gt.InputFile(file_name=output, quality=quality, process=process, delete_after_iterate=delete_after_iterate)
    else:
        logging.debug("Opening output stream")
        ## running in async mode, return iterator on
        ## the output stream
        return gt.InputFile(stream=process.stdout, quality=quality, process=process, delete_after_iterate=delete_after_iterate)


def validate_executables():
    """Validate the gem executables and
    print the paths to the executables in use
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
           min_decoded_strata=1,
           min_matched_bases=0.80,
           max_big_indel_length=15,
           max_edit_distance=0.20,
           mismatch_alphabet="ACGT",
           trim=None,
           unique_mapping=False,
           threads=1,
           extra=None,
           key_file=None,
           force_min_decoded_strata=False
           ):
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
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 1
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    trim -- tuple or list that specifies left and right trimmings
    extra -- list of additional parameters added to gem mapper call
    """

    ## prepare inputs
    index = _prepare_index_parameter(index)
    quality = _prepare_quality_parameter(quality, input)

    if delta >= min_decoded_strata and not force_min_decoded_strata:
        logging.warning("Changing min-decoded-strata from %s to %s to cope with delta of %s" % (
            str(min_decoded_strata), str(delta + 1), str(delta)))
        min_decoded_strata = delta + 1

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

    if unique_mapping:
        pa.append("--unique-mapping")

    if max_edit_distance > 0:
        pa.append("-e")
        pa.append("%s" % str(max_edit_distance))

    ## extend with additional parameters
    _extend_parameters(pa, extra)

    if key_file is not None:
        threads = max(1, threads / 2)

    trim_c = [executables['gem-map-2-map'], '-c', '-T', str(threads)]
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    # workaround for GT-32 - filter away the !
    # build list of tools
    tools = [pa]
    if unique_mapping:
        tools.append(__awk_filter)

    if trim is not None:
        tools.append(trim_c)

    # convert to genome coordinates if mapping to transcriptome
    if key_file is not None:
        convert_to_genome = [executables['transcriptome-2-genome'], key_file, str(threads)]
        tools.append(convert_to_genome)

    ## run the mapper
    process = utils.run_tools(tools, input=input, output=output, name="GEM-Mapper")
    return _prepare_output(process, output=output, quality=quality)


def transcript_mapper(input, indices, key_files, output=None,
           mismatches=0.04,
           delta=0,
           quality=33,
           quality_threshold=26,
           max_decoded_matches=100,
           min_decoded_strata=1,
           min_matched_bases=0.80,
           max_big_indel_length=15,
           max_edit_distance=0.20,
           mismatch_alphabet="ACGT",
           trim=None,
           threads=1,
           extra=None,
           ):
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
    min_decoded_strata -- strata that are decoded fully (ignoring max decoded matches), defaults to 1
    min_matched_bases -- minimum number (or %) of matched bases, defaults to 0.80
    trim -- tuple or list that specifies left and right trimmings
    extra -- list of additional parameters added to gem mapper call
    """

    if not isinstance(indices, (list, tuple)):
        indices = [indices]

    if not isinstance(key_files, (list, tuple)):
        key_files = [key_files]

    outputs = []
    output_files = []
    for i, index in enumerate(indices):
        output_file = output
        if len(indices) > 1:
            (fifo, output_file) = tempfile.mkstemp(suffix=".map", prefix="transcript_mapping_output", dir=".")
            os.close(fifo)

        output_files.append(output_file)
        outputs.append(mapper(input.clone(), index,
           key_file=key_files[i],
           output=output_file,
           mismatches=mismatches,
           delta=delta,
           quality=quality,
           quality_threshold=quality_threshold,
           max_decoded_matches=max_decoded_matches,
           min_decoded_strata=min_decoded_strata,
           min_matched_bases=min_matched_bases,
           max_big_indel_length=max_big_indel_length,
           max_edit_distance=max_edit_distance,
           mismatch_alphabet=mismatch_alphabet,
           trim=trim,
           threads=threads,
           extra=extra,
           force_min_decoded_strata=True
           )
        )
    if len(indices) > 1:
        merged = merger(outputs[0], outputs[1:])
        if(output is not None):
            merged.merge(output, threads=threads, same_content=True)
            for f in output_files:
                os.remove(f)
            return _prepare_output(None, output=output, quality=quality)
        return merged
    else:
        return _prepare_output(None, output=output, quality=quality)


def splitmapper(input,
                index,
                output=None,
                mismatches=0.04,
                splice_consensus=extended_splice_consensus,
                filter=default_filter,
                refinement_step_size=2,
                min_split_size=15,
                matches_threshold=100,
                strata_after_first=1,
                mismatch_alphabet="ACGT",
                quality=33,
                trim=None,
                filter_splitmaps=True,
                post_validate=True,
                threads=1,
                extra=None):
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
    index = _prepare_index_parameter(index, gem_suffix=True)
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
          '-s', str(strata_after_first),
          '--mismatch-alphabet', mismatch_alphabet,
          '-T', str(threads)
    ]
    min_threads = int(round(max(1, threads / 2)))

    if filter is not None:
        pa.append("-f")
        pa.append(filter)
    if splice_cons is not None:
        pa.append("-c")
        pa.append(splice_cons)

    ## extend with additional parameters
    _extend_parameters(pa, extra)

    trim_c = [executables['gem-map-2-map'], '-c', '-T', str(min_threads)]
    if trim is not None:
        ## check type
        if not isinstance(trim, (list, tuple)) or len(trim) != 2:
            raise ValueError("Trim parameter has to be a list or a tuple of size 2")
        input = gemfilter.trim(input, trim[0], trim[1], append_label=True)

    tools = [pa]
    if filter_splitmaps:
        tools.append(__awk_filter)
    if trim is not None:
        tools.append(trim_c)

    ## run the mapper
    process = None
    original_output = output
    if post_validate:
        output = None

    process = utils.run_tools(tools, input=input, output=output, name="GEM-Split-Mapper")
    splitmap_out = _prepare_output(process, output=output, quality=quality)

    if post_validate:
        return validate(splitmap_out, index, original_output, threads=threads)

    return splitmap_out


def extract_junctions(input,
                      index,
                      filter="ordered,non-zero-distance",
                      mismatches=0.04,
                      refinement_step_size=2,
                      min_split_size=15,
                      matches_threshold=75,
                      splice_consensus=extended_splice_consensus,
                      strata_after_first=1,
                      quality=33,
                      threads=1,
                      merge_with=None,
                      min_split=4,
                      max_split=2500000,
                      coverage=0,
                      max_junction_matches=5,
                      tmpdir=None,
                      extra=None):
    ## run the splitmapper
    splitmap = splitmapper(input,
        index,
        output=None,
        filter=filter,
        mismatches=mismatches,
        refinement_step_size=refinement_step_size,
        min_split_size=min_split_size,
        matches_threshold=matches_threshold,
        splice_consensus=splice_consensus,
        quality=quality,
        strata_after_first=strata_after_first,
        filter_splitmaps=False,
        post_validate=False,
        threads=threads,
        extra=extra)

    denovo_junctions = splits.extract_denovo_junctions(
        splitmap.raw_stream(),  # pass the raw stream
        minsplit=min_split,
        maxsplit=max_split,
        coverage=coverage,
        sites=merge_with,
        max_junction_matches=max_junction_matches,
        process=splitmap.process)
    return denovo_junctions


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
              threads=1,
              extra=None):
    ## check the index
    index = _prepare_index_parameter(index)
    quality = _prepare_quality_parameter(quality, input)

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

    ## extend with additional parameters
    _extend_parameters(pa, extra)

    if unique_pairing:
        pa.append("--unique-pairing")
    if map_both_ends:
        pa.append("--map-both-ends")

    ## run the mapper and trim away all the unused stuff from the ids
    process = utils.run_tool(pa, input=input, output=output, name="GEM-Pair-align", write_map=True, clean_id=True, append_extra=False)
    return _prepare_output(process, output=output, quality=quality)


def realign(input,
            index,
            output=None,
            threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=True)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-r',
                  '-T', str(threads)
    ]
    process = utils.run_tool(validate_p, input=input, output=output, name="GEM-Realign", write_map=True)
    return _prepare_output(process, output=output)


def validate(input,
             index,
             output=None,
             validate_score=None,  # "-s,-b,-i"
             validate_filter=None,  # "1,2,25"
             threads=1, ):
    index = _prepare_index_parameter(index, gem_suffix=True)
    validate_p = [executables['gem-map-2-map'],
                  '-I', index,
                  '-v', '-r',
                  '-T', str(max(threads, 1))
    ]
    if validate_score is not None:
        validate_p.extend(["-s", validate_score])
    if validate_filter is not None:
        validate_p.extend(['-f', validate_filter])

    process = utils.run_tool(validate_p, input=input, output=output, name="GEM-Validate", write_map=True)
    return _prepare_output(process, output=output)


def score(input,
          index,
          output=None,
          scoring="+U,+u,-s,-t,+1,-i,-a",
          filter=None,  # "1,2,25"
          quality=None,
          threads=1):
    """Score the input. In addition, you can specify a tuple with (<score_strata_to_keep>,<max_strata_distance>,<max_alignments>) to
    filter the result further.
    """

    quality = _prepare_quality_parameter(quality)
    index = _prepare_index_parameter(index, gem_suffix=True)
    score_p = [executables['gem-map-2-map'],
               '-I', index,
               '-s', scoring,
               '-T', str(threads)
    ]

    if filter is not None:
        score_p.append("-f")
        ff = filter
        if not isinstance(filter, basestring):
            ff = ",".join(filter)
        score_p.append(ff)

    raw = False
    if isinstance(input, gt.InputFile):
        raw = True

    process = utils.run_tool(score_p, input=input, output=output, name="GEM-Score", write_map=True, raw=raw)
    return _prepare_output(process, output=output)


def gem2sam(input, index=None, output=None,
    single_end=False, compact=False, threads=1,
    quality=None, check_ids=True):
    # if output is None:
    #     raise ValueError("You have to specify a sam output file! We are working on the piping support!")

    if index is not None:
        index = _prepare_index_parameter(index, gem_suffix=True)

    gem_2_sam_p = [executables['gem-2-sam'],
                   '-T', str(threads)
    ]
    if index is not None:
        gem_2_sam_p.extend(['-I', index])

    quality = _prepare_quality_parameter(quality, input)
    if quality is not None:
        gem_2_sam_p.extend(["-q", quality])

    if single_end:
        gem_2_sam_p.append("--expect-single-end-reads")
    if compact:
        gem_2_sam_p.append("-c")

    # GT-25 transform id's
    process = utils.run_tool(gem_2_sam_p, input=input, output=output, name="GEM-2-sam", write_map=True, clean_id=True, append_extra=False)
    # if process.wait() != 0:
    #     raise ValueError("GEM-2-SAM execution failed!")
    # return output
    return _prepare_output(process, output=output, quality=quality)


def sam2bam(input, output=None, sorted=False, tmpdir=None, mapq=None):
    sam2bam_p = ['samtools', 'view', '-S', '-b']
    if mapq is not None:
        sam2bam_p.append("-q")
        sam2bam_p.append(str(mapq))
    sam2bam_p.append('-')

    tools = [sam2bam_p]
    out_name = output
    delete_after_iterate = False
    if sorted:
        if out_name is not None:
            if out_name.endswith('.bam'):
                out_name = out_name[:-4]
        else:
            tmpfile = tempfile.NamedTemporaryFile(prefix="sorting", suffix=".bam")
            tmpfile.close()
            out_name = tmpfile.name[:-4]
            delete_after_iterate = True


        bam_sort = ['samtools', 'sort', '-']
        bam_sort.append(out_name)
        out_name = out_name + ".bam"
        tools.append(bam_sort)

    process = utils.run_tools(tools, input=input, output=output, name="SAM-2-BAM", raw=True)
    return _prepare_output(process, output=out_name, quality=33, delete_after_iterate=delete_after_iterate)


def compute_transcriptome(max_read_length, index, junctions, substract=None):
    """Compute the transcriptome based on a set of junctions. You can optionally specify
    a *substract* junction set. In that case only junctions not in substract are passed
    to compute the transcriptome.
    The function returns a tuple of a .fa file with the transcriptome genome and a
    .keys file with the translation table.

    max_read_length -- the maximum read length
    index -- path to the gem index
    junctions -- path to the junctions file
    substract -- additional juntions that are not taken into account and substracted from the main junctions
    """
    transcriptome_p = [
        executables['compute-transcriptome'],
        str(max_read_length),
        index,
        junctions
    ]
    if substract is not None:
        transcriptome_p.append(substract)

    process = utils.run_tools([transcriptome_p], input=None, output=None, name="compute-transcriptome")
    if process.wait() != 0:
        raise ValueError("Error while computing transcriptome")

    return (os.path.abspath("%s.fa" % junctions), os.path.abspath("%s.keys" % junctions))


def index(input, output, content="dna", threads=1):
    """Run the gem-indexer on the given input. Input has to be the path
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
        '-T', str(threads),
        '--content-type', content.lower()
    ]

    if isinstance(input, basestring):
        if not os.path.exists(input):
            raise ValueError("Indexer input file %s not found" % input)
        indexer_p.extend(["-i", input])
    else:
        raise ValueError("The indexer wrapper can not handle the input %s, pass a file or a list of files" % input)

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
    path = "%s:%s" % (os.path.dirname(executables['gem-indexer']), os.getenv("PATH"))

    process = utils.run_tools([indexer_p], name="gem-indexer", env={"PATH": path})
    if process.wait() != 0:
        raise ValueError("Error while executing the gem-indexer")
    return os.path.abspath("%s.gem" % output)


def hash(input, output):
    """Run the gem-retriever on the given input and create a hash
    version of a genome that can be used by the retriever to query
    the reference by chromosome and coordinates
    """
    p = [
        executables['gem-retriever'],
        'hash', input, output
    ]
    process = utils.run_tools([p], name="gem-retriever")
    if process.wait() != 0:
        raise ValueError("Error while executing the gem-retriever")
    return os.path.abspath(output)


class merger(object):
    """Merge all mappings from the source files into the target file.
    The target file must contain all available mappings and the sort order
    of all files must be the same.

    target -- either an open file descriptor or a file name.
    source -- an open file descriptor to a single source file, a file name or a
              list of file names
    exclusive - if set the True, the next mappings are take
                exlusively if the reads is not mapped at all. The default
                is False, where all mappings are merged
    paired -- merging paried reads
    """

    def __init__(self, target, source, paired=False):
        if target is None:
            raise ValueError("No target file specified")
        if source is None:
            raise ValueError("No source file specified")

        self.target = target
        self.source = source
        self.paired = paired

    def __iter__(self):
        return gt.merge(self.target, self.source)

    def merge(self, output, threads=1, paired=False, same_content=False):
        if same_content:

            files = [self.target.file_name]
            for f in self.source:
                files.append(f.file_name)

            compressed = False
            for f in files:
                if f.endswith(".gz") or f.endswith(".bz"):
                    compressed = True
                    break

            if not compressed:
                logging.debug("Using merger with %d threads and same content" % (threads))
                logging.debug("Files to merge: %s" % (" ".join(files)))
                params = [executables['gem-map-2-map'], '-T', str(threads), '-M', ",".join(files)]
                process = utils.run_tool(params, input=None, output=output, name="GEM-Merge")
                if process.wait() != 0:
                    logging.error("Merging failed !")
                    raise ValueError("Mergin failed!")
                return gt.InputFile(file_name=output)

        if len(self.source) == 1:
            logging.debug("Using paired merger with %d threads and same content %s" % (threads, str(same_content)))
            logging.debug("Files to merge: %s and %s" % (self.target.file_name, self.source[0].file_name))

            merger = gt.merge(self.target, self.source, init=False)
            merger.merge_pairs(self.target.file_name, self.source[0].file_name, output, same_content, threads)
        else:
            logging.debug("Using unpaired merger")
            merger = gt.merge(self.target, self.source, init=True)
            of = gt.OutputFile(file_name=output)
            of.write_map(merger, clean_id=False, append_extra=True)
            of.close()
        return gt.InputFile(file_name=output)


def _is_i3_compliant(stream):
    """Reads lines from the input stream and scans "flags" lines
    and returns true if the flags are compatible with the GEM
    i3 bundle. This is usually filled with the content of /proc/cpuinfo
    to determine the current systems capabilities

    The input is a stream that must provide a readline method"""
    i3_flags = set(["popcnt", "ssse3", "sse4_1", "sse4_2"])
    cpu_flags = set([])
    for line in iter(stream.readline, ''):
        line = line.rstrip()
        if line.startswith("flags"):
            for e in line.split(":")[1].strip().split(" "):
                cpu_flags.add(e)
    return i3_flags.issubset(cpu_flags)


def _extend_parameters(pa, extra=None):
    """Extend parameter array pa
    with extra parameters. If extra is a string, it is
    split by space and parameters are added separatly,
    otherwise it is assumed that extra is a list and
    it is appended to the parameter array as is.
    """
    if extra is not None:
        if isinstance(extra, (basestring,)):
            pa.extend(extra.split(" "))
        else:
            pa.extend(extra)
