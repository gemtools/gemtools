#!/usr/bin/env
"""Pipeline utilities"""
import os
import logging
import gem

from gem.utils import Timer


class PipelineError(Exception):
    """Exception thrown by the mapping pipeline"""
    pass


class MappingPipeline(object):
    """General mapping pipeline class."""

    def __init__(self):
        self.files = []
        self.mappings = []

        self.name = None  # target name
        self.index = None  # genome index
        self.output_dir = None  # Output directory
        self.annotation = None  # GTF annotation to use
        self.threads = 1  # number of threads
        self.junctioncoverage = 2  # junction coverage
        self.maxlength = 150  # max read length
        self.transcript_index = None  # transcriptome index
        self.transcript_keys = None  # transcriptome keys file
        self.denovo_index = None  # the denovo index to use
        self.denovo_keys = None  # the denovo keys to use
        self.quality = 33  # quality offset
        self.junctions_file = None  # file with both denovo and GTF junctions
        self.junction_length = 500000

        self.genome_mismatches = 0.06
        self.genome_quality_threshold = 26
        self.genome_max_decoded_matches = 20
        self.genome_min_decoded_strata = 1
        self.genome_min_matched_bases = 0.80
        self.genome_max_big_indel_length = 15
        self.genome_max_edit_distance = 0.20
        self.genome_mismatch_alphabet = "ACGT"
        self.genome_strata_after_best = 1

        self.transcript_mismatches = self.genome_mismatches
        self.transcript_quality_threshold = 26
        self.transcript_max_decoded_matches = 50
        self.transcript_min_decoded_strata = 1
        self.transcript_min_matched_bases = 0.80
        self.transcript_max_big_indel_length = 15
        self.transcript_max_edit_distance = 0.20
        self.transcript_mismatch_alphabet = "ACGT"
        self.transcript_strata_after_best = 1

        self.junction_mismatches = 0.04
        self.junctions_max_junction_matches = 5
        self.junctions_min_split_length = 4
        self.junctions_max_split_length = 500000
        self.junctions_refinement_step_size = 2
        self.junctions_min_split_size = 15
        self.junctions_matches_threshold = 75

        self.pairing_quality_threshold = 26
        self.pairing_max_decoded_matches = 20
        self.pairing_min_decoded_strata = 1
        self.pairing_min_insert_size = 0
        self.pairing_max_insert_size = self.junctions_max_split_length
        self.pairing_max_edit_distance = 0.30
        self.pairing_min_matched_bases = 0.80
        self.pairing_max_extendable_matches = 0
        self.pairing_max_matches_per_extension = 1

        self.scoring_scheme = "+U,+u,-s,-t,+1,-i,-a"
        self.filter = (1, 2, 25)

        self.compress = True
        self.remove_temp = True
        self.bam_mapq = 0
        self.bam_create = True
        self.bam_sort = True
        self.bam_index = True

        self.pipeline = []

    def add_step(self, name, function, args=None, kwargs=None):
        if name is None:
            raise PipelineError("Please specify a name for the pipeline step")
        if function is None:
            raise PipelineError("Please specify a function to call")
        try:
            if isinstance(function, basestring):
                getattr(self, function)
            else:
                function = function.__name__
        except:
            raise PipelineError("Unable to resolve function %s" % (function))

        self.pipeline.append({
            "name": name,
            "function": function,
            "args": args,
            "kwargs": kwargs
        })

    def run(self):
        logging.info("Running pipeline with %d steps" % (len(self.pipeline)))
        for i, step in enumerate(self.pipeline):
            logging.info("Step %d - %s" % (i + 1, step["name"]))
            function = self.__get_function(step["function"])
            args = step["args"]
            kwargs = step["kwargs"]
            function(*args, **kwargs)

    def __get_function(self, function):
        try:
            return getattr(self, function)
        except:
            raise PipelineError("Unable to resolve function %s" % (function))

    def initialize(self):
        ## initialize defaults
        if self.output_dir is None:
            logging.debug("Switching output directory to %s" % os.getcwd())
            self.output_dir = os.getcwd()

        ## check if output is specified and create folder in case it does not exists
        if not os.path.exists(self.output_dir):
            logging.debug("Creating output directory  %s" % self.output_dir)
            os.makedirs(self.output_dir)

        if self.annotation is not None:
            if self.transcript_index is None:
                self.transcript_index = self.annotation + ".gem"
                logging.debug("Setting trascriptome index to %s" % (self.transcript_index))

            if self.transcript_keys is None:
                self.transcript_keys = self.transcript_index[:-4] + ".junctions.keys"
                logging.debug("Setting trascriptome keys to %s" % (self.transcript_keys))

        if self.transcript_mismatches is None:
            self.transcript_mismatches = self.genome_mismatches

    def cleanup(self):
        """Delete all remaining temporary and intermediate files
        """
        if self.remove_temp:
            for f in self.files:
                if os.path.exists(f):
                    logging.info("Removing intermediate file : %s" % (f))
                    os.remove(f)

    def create_file_name(self, suffix, file_suffix="map", final=False):
        """Create a result file name"""
        file = ""
        if suffix is not None and len(suffix) > 0:
            file = "%s/%s_%s.%s" % (self.output_dir, self.name, suffix, file_suffix)
        else:
            file = "%s/%s.%s" % (self.output_dir, self.name, file_suffix)
        if not final:
            self.files.append(file)
        return file

    def _guess_input_name(self, input):
        try:
            return input.file_name
        except Exception:
            return "STREAM"

    def merge(self, suffix, same_content=False):
        """Merge current set of mappings and delete last ones"""
        out = self.create_file_name(suffix)
        if os.path.exists(out):
            logging.warning("Merge target exists, skipping merge : %s" % (out))
            return gem.files.open(out, quality=self.quality)

        logging.info("Merging %d mappings into %s" % (len(self.mappings), out))
        timer = Timer()
        merged = gem.merger(
            self.mappings[0].clone(),
            [m.clone() for m in self.mappings[1:]]
        ).merge(out, self.threads, same_content=same_content)

        if self.remove_temp:
            for m in self.mappings:
                logging.info("Removing temporary mapping %s ", m.filename)
                os.remove(m.filename)

        self.mappgins = []
        self.mappings.append(merged)
        timer.stop("Merging finished in %s")
        return merged

    def mapping_step(self, input, suffix, trim=None):
        """Single mapping step where the input is passed on
        as is to the mapper. The output name is created based on
        the dataset name in the parameters the suffix.

        The function returns an opened ReadIterator on the resulting
        mapping

        input -- mapper input
        suffix -- output name suffix
        trim -- optional trimming options, i.e. '0,20'
        """
        timer = Timer()
        mapping_out = self.create_file_name(suffix)
        if os.path.exists(mapping_out):
            logging.warning("Mapping step target exists, skip mapping set : %s" % (mapping_out))
            mapping = gem.files.open(mapping_out, quality=self.quality)
            self.mappings.append(mapping)
            return mapping
        input_name = self._guess_input_name(input)
        logging.debug("Mapping from %s to %s" % (input_name, mapping_out))
        mapping = gem.mapper(input,
                self.index,
                mapping_out,
                mismatches=self.genome_mismatches,
                quality_threshold=self.genome_quality_threshold,
                max_decoded_matches=self.genome_max_decoded_matches,
                min_decoded_strata=self.genome_min_decoded_strata,
                min_matched_bases=self.genome_min_matched_bases,
                max_big_indel_lengt=self.genome_max_big_indel_length,
                max_edit_distance=self.genome_max_edit_distance,
                mismatch_alphabet=self.genome_mismatch_alphabet,
                delta=self.genome_strata_after_best,
                trim=trim,
                quality=self.quality,
                threads=self.threads)
        self.mappings.append(mapping)
        timer.stop("Mapping step finished in %s")
        return mapping

    def transcript_mapping_step(self, input, suffix, index=None, key=None, trim=None):
        """Single transcript mapping step where the input is passed on
        as is to the mapper. The output name is created based on
        the dataset name in the parameters the suffix.

        The function returns an opened ReadIterator on the resulting
        mapping

        input -- mapper input
        suffix -- output name suffix
        trim -- optional trimming options, i.e. '0,20'
        """
        timer = Timer()
        mapping_out = self.create_file_name(suffix + "_transcript")
        if os.path.exists(mapping_out):
            logging.warning("Transcript mapping step target exists, skip mapping set : %s" % (mapping_out))
            mapping = gem.files.open(mapping_out, quality=self.quality)
            self.mappings.append(mapping)
            return mapping
        input_name = self._guess_input_name(input)
        logging.debug("Mapping transcripts from %s to %s" % (input_name, mapping_out))

        indices = [self.transcript_index, self.denovo_index]
        keys = [self.transcript_keys, self.denovo_keys]

        if index is not None:
            indices = [index]
            if key is None:
                key = index[:-4] + ".junctions.keys"
            keys = [key]

        mapping = gem.transcript_mapper(
                                input,
                                indices,
                                keys,
                                mapping_out,
                                mismatches=self.transcript_mismatches,
                                quality_threshold=self.transcript_quality_threshold,
                                max_decoded_matches=self.transcript_max_decoded_matches,
                                min_decoded_strata=self.transcript_min_decoded_strata,
                                min_matched_bases=self.transcript_min_matched_bases,
                                max_big_indel_lengt=self.transcript_max_big_indel_length,
                                max_edit_distance=self.transcript_max_edit_distance,
                                mismatch_alphabet=self.transcript_mismatch_alphabet,
                                delta=self.transcript_strata_after_best,
                                trim=trim,
                                quality=self.quality,
                                threads=self.threads,
                                )
        self.mappings.append(mapping)
        timer.stop("Transcript-Mapping step finished in %s")
        return mapping

    def gtf_junctions(self):
        """check if there is a .junctions file for the given annotation, if not,
        create it. Returns a tuple of the set of junctions and the output file.
        """
        timer = Timer()
        gtf_junctions = self.annotation + ".junctions"
        out = None
        junctions = None
        if os.path.exists(gtf_junctions):
            logging.info("Loading existing junctions from %s" % (gtf_junctions))
            out = gtf_junctions
            junctions = set(gem.junctions.from_junctions(gtf_junctions))
        else:
            out = self.create_file_name("gtf", file_suffix="junctions")
            if os.path.exists(out):
                logging.info("Loading existing junctions from %s" % (out))
                junctions = set(gem.junctions.from_junctions(out))
            else:
                logging.info("Extracting junctions from %s" % (self.annotation))
                junctions = set(gem.junctions.from_gtf(self.annotation))
                gem.junctions.write_junctions(junctions, out, self.index)

        logging.info("%d Junctions from GTF" % (len(junctions)))
        timer.stop("GTF-Junctions prepared in %s")
        return (junctions, out)

    def create_denovo_transcriptome(self, input):
        """Squeeze input throw the split mapper to extract denovo junctions
        and create a transcriptome out of the junctions"""
        index_denovo_out = self.create_file_name("denovo_transcripts", file_suffix="gem")
        junctions_out = self.create_file_name("all", file_suffix="junctions")
        denovo_keys = junctions_out + ".keys"
        if os.path.exists(index_denovo_out) and os.path.exists(denovo_keys):
            logging.warning("Transcriptome index and keys found, skip creating : %s" % (index_denovo_out))
            ## add .fa and .keys to list of files to delete
            self.files.append(junctions_out + ".fa")
            self.files.append(junctions_out + ".keys")
            self.files.append(index_denovo_out[:-4] + ".log")
            self.denovo_keys = denovo_keys
            self.denovo_index = index_denovo_out
            return index_denovo_out

        (gtf_junctions, junctions_gtf_out) = self.gtf_junctions()

        timer = Timer()
        ## get de-novo junctions
        logging.info("Getting de-novo junctions")

        ## add .fa and .keys to list of files to delete
        self.files.append(junctions_out + ".fa")
        self.files.append(junctions_out + ".keys")

        denovo_junctions = gem.extract_junctions(
            input,
            self.index,
            mismatches=self.junction_mismatches,
            threads=self.threads,
            strata_after_first=0,
            coverage=self.junctioncoverage,
            min_split=self.junctions_min_split_length,
            max_split=self.junctions_max_split_length,
            refinement_step_size=self.junctions_refinement_step_size,
            min_split_size=self.junctions_min_split_size,
            matches_threshold=self.junctions_matches_threshold,
            max_junction_matches=self.junctions_max_junction_matches,
            )

        logging.info("Denovo Junctions %d" % len(denovo_junctions))
        filtered_denovo_junctions = set(gem.junctions.filter_by_distance(denovo_junctions, self.junction_length))
        logging.info("Denovo Junction passing distance (%d) filter %d (%d removed)" % (self.junction_length,
            len(filtered_denovo_junctions), (len(denovo_junctions) - len(filtered_denovo_junctions))))

        junctions = gtf_junctions.union(filtered_denovo_junctions)
        logging.info("Total Junctions %d" % (len(junctions)))
        timer.stop("Denovo Transcripts extracted in %s")
        timer = Timer()
        gem.junctions.write_junctions(junctions, junctions_out, self.index)

        logging.info("Computing denovo transcriptome")
        (denovo_transcriptome, denovo_keys) = gem.compute_transcriptome(self.maxlength, self.index, junctions_out, junctions_gtf_out)
        logging.info("Indexing denovo transcriptome")
        timer.stop("Transcriptome generated in %s")
        timer = Timer()
        idx = gem.index(denovo_transcriptome, index_denovo_out, threads=self.threads)
        timer.stop("Transcriptome indexed in %s")
        # add indexer log file to list of files
        self.files.append(idx[:-4] + ".log")
        self.denovo_keys = denovo_keys
        self.denovo_index = idx
        # do not add the denovo mapping but make sure file goes to temp files
        #self.mappings.append(denovo_mapping)
        #self.files.append(denovo_out)
        return idx

    def pair_align(self, input, final=False):
        n = "paired"
        if final:
            n = ""
        paired_out = self.create_file_name(n, final=final)
        if os.path.exists(paired_out):
            logging.warning("Paired alignment found, skip pairing : %s" % (paired_out))
            mapping = gem.files.open(paired_out, quality=self.quality)
            self.mappings.append(mapping)
            return mapping

        logging.info("Running pair aligner")
        timer = Timer()

        paired_mapping = gem.pairalign(input, self.index, None,
          quality_threshold=self.pairing_quality_threshold,
          max_decoded_matches=self.pairing_max_decoded_matches,
          min_decoded_strata=self.pairing_min_decoded_strata,
          min_insert_size=self.pairing_min_insert_size,
          max_insert_size=self.junctions_max_split_length,
          max_edit_distance=self.pairing_max_edit_distance,
          min_matched_bases=self.pairing_min_matched_bases,
          max_extendable_matches=self.pairing_max_extendable_matches,
          max_matches_per_extension=self.pairing_max_matches_per_extension,
          threads=max(self.threads - 2, 1),
          quality=self.quality)
        scored = gem.score(paired_mapping, self.index, paired_out, threads=min(2, self.threads), quality=self.quality)
        timer.stop("Pair-Align and scoring finished in %s")
        return scored

    def compress(self):
        """Compress the final alignment"""
        if len(self.files) > 0:
            file_name = self.files[-1]
            if os.path.exists(file_name + ".gz"):
                logging.warning("Compressed file found, skipping compression: %s" % file_name)
            else:
                logging.info("Compressing mapping")
                timer = Timer()
                gem.utils.gzip(self.files[-1], threads=self.threads)
                timer.stop("Results compressed in %s")
        else:
            logging.warning("No files found, skip compressions")

    def create_bam(self, input, sort=True):
        logging.info("Converting to sam/bam")
        sam_out = self.create_file_name("", file_suffix="bam", final=True)
        if os.path.exists(sam_out):
            logging.warning("SAM/BAM exists, skip creating : %s" % (sam_out))
        else:
            timer = Timer()
            sam = gem.gem2sam(input, self.index, threads=max(1, int(self.threads / 2)))
            gem.sam2bam(sam, sam_out, sorted=sort)
            timer.stop("BAM file created in %s")
