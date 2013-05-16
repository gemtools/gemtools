from gemapi cimport *

import os
import sys
import multiprocessing
import string


from cython.parallel import parallel, prange

cdef class TemplateFilter(object):
    """Filter templates. Override the filter
    method to modify the filter.
    """

    cpdef bool filter(self, Template template):
        """Default implementation always returns true,
        but you can modify this to exclude tempaltes completely.

        Implementations are also allowed to modify the template
        """
        return True


cdef class filter_unmapped(TemplateFilter):
    """Filter for unmapped reads or reads with mismatches >= max_mismatches"""
    # number of mismatches allowed to be not marked as unmapped
    cdef int64_t max_mismatches

    def __init__(self, int64_t max_mismatches=GT_ALL):
        """Initialize a new unmapped filter. If max_mismatches is
        set, this defines the maximum number of mismatches a template
        can have to be still treated as mapped.

        max_mismatches -- number of allowed max mismatches for a template
        """
        self.max_mismatches = max_mismatches

    cpdef bool filter(self, Template template):
        return template.is_unmapped(self.max_mismatches)

cdef class filter_unique(TemplateFilter):
    """Filter for unique mappings up to given uniqueness level
    """
    # the max level
    cdef uint64_t level

    def __init__(self, uint64_t level=0):
        self.level = level

    cpdef bool filter(self, Template template):
        template_level = template.level(self.level)
        return template_level > 0 and template_level >= self.level


cdef class filter_trim(TemplateFilter):
    """Trimming filter that always returns true, but trims the
    template
    """
    cdef uint64_t left
    cdef uint64_t right
    cdef uint64_t min_length
    cdef bool set_extra

    def __init__(self, uint64_t left=0, uint64_t right=0, uint64_t min_length=10, bool set_extra=True):
        self.left = left
        self.right = right
        self.min_length = min_length
        self.set_extra = set_extra

    cpdef bool filter(self, Template template):
        gt_template_trim(template.template, self.left, self.right, self.min_length, self.set_extra)
        return True

cpdef unmapped(source, int64_t max_mismatches=GT_ALL):
    """Wrapper function that creates a filtered iterator"""
    return filter(source, filter_unmapped(max_mismatches))

cpdef unique(source, uint64_t level=0):
    """Wrapper function that creates a filtered iterator"""
    return filter(source, filter_unique(level))

cpdef trim(source, uint64_t left=0, uint64_t right=0, uint64_t min_length=10, bool set_extra=True):
    """Wrapper function that creates a filtered iterator"""
    return filter(source, filter_trim(left=left, right=right, min_length=min_length, set_extra=set_extra))

cdef class filter(object):
    """Filtering iterator that takes a source
    and can apply a filter
    """
    cdef object source
    cdef object template_filter

    def __cinit__(self, source, template_filter):
        self.source = source
        if isinstance(template_filter, (list, tuple)):
            self.template_filter = template_filter
        else:
            self.template_filter = [template_filter]

    def __iter__(self):
        self.source = self.source.__iter__()
        return self

    def __next__(self):
        cdef bool filtered = False
        while True:
            filtered = False
            n = self.source.__next__()
            for f in self.template_filter:
                if not f.filter(n):
                    filtered = True
                    break
            if not filtered:
                return n


    cpdef write_stream(self, OutputFile output, bool write_map=False, uint64_t threads=1):
        """Write the content of this filter to the output file.

        output   -- the output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        threads       -- number of threads to use (if supported by the iterator)
        """
        for t in self:
            output.write(t, write_map=write_map)



cdef class interleave(object):
    """Interleaving iterator over templates from a list of input files
    """
    # the input files
    cdef object files
    # counter
    cdef int64_t i
    # number of iterators
    cdef int64_t length
    # number of threads
    cdef int64_t threads
    # interleave
    cdef bool interleave

    def __init__(self, files, interleave=True, uint64_t threads=1):
        self.files = files
        self.i = 0
        self.length = len(files)
        self.interleave = interleave
        self.threads = threads

    def __iter__(self):
        # initialize iterators
        self.files = [f.__iter__() for f in self.files]
        return self

    def __next__(self):
        cdef int64_t mises = 0
        if self.interleave:
            while True:
                try:
                    return self.files[self.i].__next__()
                except StopIteration:
                    mises += 1
                    if mises >= self.length:
                        raise StopIteration()
                finally:
                    self.i += 1
                    if self.i >= self.length:
                        self.i = 0
        else:
            while True:
                try:
                    return self.files[self.i].__next__()
                except StopIteration:
                    self.i += 1
                    mises += 1
                    if mises >= self.length or self.i >= self.length:
                        raise StopIteration()

    cpdef write_stream(self, OutputFile output, bool write_map=False, uint64_t threads=1):
        """Write the content interleaved to the output file

        output_file   -- the output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        threads       -- number of threads to use (if supported by the iterator)
        """
        __run_write_stream(self.files, output, write_map, max(threads, self.threads), self.interleave, None)

    cpdef close(self):
        try:
            for f in self.files:
                f.close()
        except Exception:
            pass



cdef class cat(interleave):

    def __init__(self, files, uint64_t threads=1):
        interleave.__init__(self, files, interleave=False, threads=threads)


cdef class merge(object):
    """Merge a master stream with a list of child streams"""
    cdef object inputs
    cdef object buffers

    def __init__(self, inputs):
        self.inputs = inputs
        self.buffers = []

    def __iter__(self):
        self.inputs = [s.__iter__() for s in self.inputs]
        for i, t in enumerate(self.inputs[1:]):
            self.buffers.append(t.__next__())
        return self

    def __next__(self):
        cdef Template m = self.inputs[0].__next__()
        cdef Template o
        cdef buffers = self.buffers
        for i, t in enumerate(self.inputs[1:]):
            o = buffers[i]
            if o is None:
                try:
                    buffers[i] = t.__next__()
                except:
                    pass
            if o is not None and o.same(m):
                m.merge(o)
                try:
                    buffers[i] = t.__next__()
                except:
                    pass
        return m

    cpdef write_stream(self, OutputFile output, bool write_map=False, uint64_t threads=1, bool async=False):
        """Write the content of this input stream to the output file
        file.

        output   -- the output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        interleave    -- interleave muliple inputs
        threads       -- number of threads to use (if supported by the iterator)
        """
        # cdef gt_output_file* output_file = output.output_file
        # cdef uint64_t num_inputs = len(self.inputs)
        # cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        # cdef bool clean_id = output.clean_id
        # cdef bool append_extra = output.append_extra
        # for i in range(num_inputs):
        #     inputs[i] = (<InputFile> self.inputs[i])._open()

        # with nogil:
        #     gt_merge_files_synch(output_file, threads, num_inputs, inputs);

        # output.close()
        # for i in range(num_inputs):
        #     gt_input_file_close(inputs[i])
        # free(inputs)
        __run_write_stream(self.inputs, output, write_map, threads, False, None, __write_merge_stream, async)




cdef class OutputFile:
    """The OutputFile can write content to a file or
    or a stream.
    """
    # the target output file
    cdef gt_output_file* output_file
    # the target
    cdef readonly object target
    # the map attributes
    cdef gt_output_map_attributes* map_attributes
    # the fasta attributes
    cdef gt_output_fasta_attributes* fasta_attributes
    # clean read ids
    cdef bool clean_id
    # append extra information
    cdef bool append_extra
    # list of filters
    cdef object filters

    def __init__(self, target, bool clean_id=False, bool append_extra=True):
        """Initialize the output file from the given target. The
        target can be either a string a stream. If init_buffer is
        true, the output buffer is initialized. The output can be configured
        to clean the ids and append extra information. Clean ids encode the
        pair information as /1 and /2. Extra information is everything after
        the first space in the id or the everything after the casava id in case
        of tempalte tags encoded as casava >= 1.8+.

        target       -- the target file or stream
        clean_id     -- ensure /1 /2 read pair encoding
        append_extra -- append additional infomration to the id
        init_buffer  -- if true, the buffered output will be initialized, default True
        """
        self.target = target
        self.map_attributes = gt_output_map_attributes_new()
        self.fasta_attributes = gt_output_fasta_attributes_new()
        self.clean_id = clean_id
        self.append_extra = append_extra
        self.filters = None
        # init attributes
        gt_output_map_attributes_set_print_extra(self.map_attributes, append_extra)
        gt_output_map_attributes_set_print_casava(self.map_attributes, not clean_id)
        gt_output_fasta_attributes_set_print_extra(self.fasta_attributes, append_extra)
        gt_output_fasta_attributes_set_print_casava(self.fasta_attributes, not clean_id)

        # open the outout file or stream
        if isinstance(target, basestring):
            self._open_file(<char*> target)
        else:
            self._open_stream(<file> target)

    cpdef bool is_stream(self):
        return isinstance(self.target, basestring)

    def add_filter(self, filter):
        if self.filters is None:
            self.filters = []
        self.filters.append(filter)

    def __dealloc__(self):
        self.close()
        gt_output_map_attributes_delete(self.map_attributes)
        gt_output_fasta_attributes_delete(self.fasta_attributes)

    cpdef _open_stream(self, file stream):
        """Initialize this instance from a stream

        stream      -- the output stream
        init_buffer -- if true, initialize the output buffer
        """
        self.output_file = gt_output_stream_new(PyFile_AsFile(stream), SORTED_FILE)

    cpdef _open_file(self, char* file_name):
        """Initialize this instance from a file

        file_name   -- the output file name
        init_buffer -- if true, initialize the output buffer
        """
        self.output_file = gt_output_file_new(file_name, SORTED_FILE)

    cpdef close(self):
        """Close the output file"""
        if self.output_file is not NULL:
            gt_output_file_close(self.output_file)
            self.output_file = NULL
        if not isinstance(self.target, basestring):
            self.target.close()

    cpdef write(self, Template template, write_map=True):
        """Write a single template to this otuout file.

        template  -- the source tempalte
        write_map -- writes map or fastq/a format
        """
        if self.filters is not None:
            for f in self.filters:
                if not f.filter(template):
                    return
        # write a single template
        if write_map:
            gt_output_map_ofprint_template(self.output_file, template.template, self.map_attributes)
        else:
            gt_output_fasta_ofprint_template(self.output_file, template.template, self.fasta_attributes)



cdef class InputFile(object):
    """GEMTools input file. The input file extends TemplateIterator
    and can be used directly to iterate templates. To access the underlying
    alignment, use the alignments() function to get an alignment iterator.
    """
    cdef readonly object source
    # the source file name if specified
    cdef readonly object filename
    # force reading paired reads
    cdef readonly bool force_paired_reads
    # memory map the file
    cdef readonly bool mmap_file
    # the process that creates the content of the file or stream
    cdef readonly object process
    # the quality offset
    cdef readonly object quality

    # parsing attributes
    # the buffered input file
    # the source input file
    cdef gt_input_file* input_file
    #
    cdef gt_buffered_input_file* buffered_input
    # parser attributes
    cdef gt_generic_parser_attr* parser_attr
    # the template instance that is used to iterate templates
    cdef readonly Template template
    # remove scores when printing
    cdef public bool remove_scores

    def __init__(self, source, bool mmap_file=False, bool force_paired_reads=False, object quality=None, object process=None):
        """Initialize a new input file on the source. The source can be either a file name
        or a stream.

        """
        self.source = source
        self.force_paired_reads = force_paired_reads
        self.mmap_file = mmap_file
        self.process = process
        self.quality = quality
        self.template = Template()
        self.remove_scores = False
        if isinstance(source, basestring):
            self.filename = source
        # make sure memory mapping is disabled for compressed files and
        # # streams
        if self.filename is None or self.filename.endswith(".gz") or self.filename.endswith(".bz2"):
                self.mmap_file = False

    def __dealloc__(self):
        if self.buffered_input is not NULL:
            gt_buffered_input_file_close(self.buffered_input)
            self.buffered_input = NULL
        if self.parser_attr is not NULL:
            free(self.parser_attr)
            self.parser_attr = NULL
        if self.input_file is not NULL:
            gt_input_file_close(self.input_file)
            self.input_file = NULL

    cdef gt_input_file* _open(self):
        """Open a new gt_input_file"""
        if self.filename is not None:
            return gt_input_file_open(<char*>self.filename, self.mmap_file)
        else:
            return gt_input_stream_open(PyFile_AsFile(self.source))

    def raw_sequence_stream(self):
        """Return true if this is a file based
        input file, uncompressed and fastq/q format
        """
        if self.filename is None:
            return False
        cdef gt_input_file* infile = self._open()
        valid = False
        if infile.file_type == REGULAR_FILE and infile.file_format == FASTA:
            valid = True
        gt_input_file_close(infile)
        return valid

    def clone(self):
        """If this is a file based, this returns a new
        instance of the input file, otherwise an
        exception is thrown as we can not
        clone stream based inputs.
        """
        if self.filename is None:
            raise ValueError("Can not clone a stream based input file")
        else:
            return InputFile(self.source, mmap_file=self.mmap_file, force_paired_reads=self.force_paired_reads, quality=self.quality, process=self.process)

    def raw_stream(self):
        """Return the raw stream on this input file.
        In case this is stream based, the stream is returned,
        otherwise a new file handle is opened on
        the input file.
        """
        if self.filename is None:
            return self.source
        else:
            return open(self.filename, "rb")

    def __iter__(self):
        """Initialize buffers and prepare for iterating"""
        ##
        self.input_file = self._open()
        self.buffered_input = gt_buffered_input_file_new(self.input_file)
        self.parser_attr = gt_input_generic_parser_attributes_new(self.force_paired_reads)
        return self

    def __next__(self):
        if self._next() == GT_STATUS_OK:
            return self.template
        else:
            raise StopIteration()

    cpdef gt_status _next(self):
        """Internal iterator method"""
        cdef gt_status s = gt_input_generic_parser_get_template(self.buffered_input, self.template.template, self.parser_attr)
        if s != GT_STATUS_OK:
            if self.process is not None:
                # if this is a stream based process, make sure we clean up
                self.process.wait()
        return s

    cpdef write_stream(self, OutputFile output, bool write_map=False, uint64_t threads=1):
        """Write the content of this input stream to the output file
        file.

        output   -- the output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        interleave    -- interleave muliple inputs
        threads       -- number of threads to use (if supported by the iterator)
        """
        __run_write_stream([self], output, write_map, threads, True, self.process, remove_scores=self.remove_scores)

    cpdef close(self):
        if self.buffered_input is not NULL:
            gt_buffered_input_file_close(self.buffered_input)
            self.buffered_input = NULL
        if self.parser_attr is not NULL:
            free(self.parser_attr)
            self.parser_attr = NULL
        if self.input_file is not NULL:
            gt_input_file_close(self.input_file)
            self.input_file = NULL


cpdef __run_write_stream(source, OutputFile output, bool write_map=False, uint64_t threads=1, bool interleave=True, parent=None, function=__write_stream, bool async=False, bool remove_scores=False):
    import gem.utils
    process = multiprocessing.Process(target=function, args=(source, output, write_map, threads, interleave, remove_scores))
    gem.utils.register_process(process)
    process.start()

    if not async:
        process.join()
        if parent is not None:
            parent.wait()
    return process

cpdef __write_stream(source, OutputFile output, bool write_map=False, uint64_t threads=1, bool interleave=True, bool remove_scores=False):
    cdef gt_output_file* output_file = output.output_file
    cdef uint64_t num_inputs = len(source)
    cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
    cdef bool clean_id = output.clean_id
    cdef bool append_extra = output.append_extra
    cdef uint64_t use_threads = threads

    for i in range(num_inputs):
        inputs[i] = (<InputFile> source[i])._open()

    with nogil:
        gt_write_stream(output_file, inputs, num_inputs, append_extra, clean_id, interleave, use_threads, write_map, remove_scores)

    output.close()
    for i in range(num_inputs):
        gt_input_file_close(inputs[i])
    free(inputs)

cpdef __write_merge_stream(source, OutputFile output, bool write_map=False, uint64_t threads=1, bool interleave=True):
    cdef gt_output_file* output_file = output.output_file
    cdef uint64_t num_inputs = len(source)
    cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
    cdef uint64_t use_threads = threads

    for i in range(num_inputs):
        inputs[i] = (<InputFile> source[i])._open()

    with nogil:
        gt_merge_files_synch(output_file, threads, num_inputs, inputs);
    output.close()
    for i in range(num_inputs):
        gt_input_file_close(inputs[i])
    free(inputs)


cdef _create_alignment(gt_alignment* ali):
    a = Alignment(initialize=False)
    a.alignment = ali
    return a

cdef class Alignment:
    """Wrapper class around gt_alignment"""
    cdef gt_alignment* alignment
    cdef bool initialize

    def __cinit__(self, initialize=True):
        self.initialize = initialize
        if initialize:
            self.alignment = gt_alignment_new()

    def __dealloc__(self):
        if self.initialize:
            gt_alignment_delete(self.alignment)


    property tag:
        def __get__(self):
            return gt_alignment_get_tag(self.alignment)
        def __set__(self, value):
            gt_alignment_set_tag(self.alignment, value, len(value))

    property pair:
        def __get__(self):
            return gt_alignment_get_pair(self.alignment)

    property counters:
        def __get__(self):
            return gt_alignment_get_num_counters(self.alignment)

    property num_maps:
        def __get__(self):
            return gt_alignment_get_num_maps(self.alignment)

    property mcs:
        def __get__(self):
            return gt_alignment_get_mcs(self.alignment)
        def __set__(self, value):
            gt_alignment_set_mcs(self.alignment, value)

    property has_qualities:
        def __get__(self):
            return gt_alignment_has_qualities(self.alignment)

    property not_unique_flag:
        def __get__(self):
            return gt_alignment_get_not_unique_flag(self.alignment)
        def __set__(self, value):
            gt_alignment_set_not_unique_flag(self.alignment, value)

    property length:
        def __get__(self):
            return self._length()

    property read:
        def __get__(self):
            return self._read()

    property qualities:
        def __get__(self):
            return self._qualities()

    cdef uint64_t _length(self):
        return gt_alignment_get_read_length(self.alignment)

    cdef char* _read(self):
        return gt_alignment_get_read(self.alignment)

    cdef char* _qualities(self):
        return gt_alignment_get_qualities(self.alignment)


    cpdef uint64_t get_counter(self, uint64_t stratum):
        return gt_alignment_get_counter(self.alignment, stratum)

    cpdef set_counter(self, uint64_t stratum, uint64_t value):
        gt_alignment_set_counter(self.alignment, stratum, value)

    cpdef uint64_t get_num_maps(self):
        return gt_alignment_get_num_maps(self.alignment)


    def to_map(self):
        cdef gt_string* gt = self._to_map()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_fasta(self):
        cdef gt_string* gt = self._to_fasta()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_fastq(self):
        cdef gt_string* gt = self._to_fastq()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_sequence(self):
        cdef gt_string* gt = self._to_sequence()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s


    cdef gt_string* _to_map(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_map_attributes* attr = gt_output_map_attributes_new()
        gt_output_map_sprint_alignment(s,self.alignment, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_map_attributes_delete(attr)
        return s

    cdef gt_string* _to_fasta(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_fasta_attributes* attr = gt_output_fasta_attributes_new()
        gt_output_fasta_attributes_set_format(attr, F_FASTA)
        gt_output_fasta_sprint_alignment(s,self.alignment, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_fasta_attributes_delete(attr)
        return s

    cdef gt_string* _to_fastq(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_fasta_attributes* attr = gt_output_fasta_attributes_new()
        gt_output_fasta_sprint_alignment(s,self.alignment, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_fasta_attributes_delete(attr)
        return s

    cdef gt_string* _to_sequence(self):
        if gt_alignment_has_qualities(self.alignment):
            return self._to_fastq()
        else:
            return self._to_fasta()

    cpdef int64_t get_min_mismatches(self):
        """Return minimum number of mismatches or -1 if
        the alignment is unmapped"""
        cdef gt_alignment* alignment = self.alignment
        cdef uint64_t counter = 0
        cdef uint64_t c = gt_alignment_get_num_counters(alignment)
        cdef int64_t i = 0
        for i in range(c):
            counter = gt_alignment_get_counter(alignment, i)
            if counter > 0:
                return i
        return -1

    cpdef get_pair(self):
        """Return 0 for unpaired or 1 or 2"""
        return gt_alignment_get_pair(self.alignment)

    cpdef int64_t level(self, uint64_t max_level=GT_ALL):
        cdef gt_alignment* alignment = self.alignment
        cdef uint64_t counter = 0
        cdef uint64_t c = gt_alignment_get_num_counters(alignment)
        cdef int64_t i, j = 0
        cdef int64_t level = 0
        for i in range(c):
            counter = gt_alignment_get_counter(alignment, i)
            if counter == 1:
                for j in range(i+1, c):
                    if level >= max_level:
                        return level
                    counter = gt_alignment_get_counter(alignment, j)
                    if counter > 0:
                        return <int64_t> (j-(i+1))
                    else:
                        level += 1
                return <int64_t> (c - (i+1))
            elif counter > 1:
                return -1
        return -1

    cpdef maps(self):
        return [_create_map(gt_alignment_get_map(self.alignment, i)) for i in range(self.num_maps)]

cdef _create_map(gt_map* _map):
    m = Map(initialize=False)
    m._map = _map
    return m

cdef class Map:
    cdef gt_map* _map
    cdef bool initialize

    def __cinit__(self, initialize=True):
        self.initialize = initialize
        if initialize:
            self._map = gt_map_new()

    def __dealloc__(self):
        if self.initialize:
            gt_map_delete(self._map)

    property seqname:
        """The sequence name of the genomic sequence (i.e. Chromosome)"""
        def __get__(self):
            return gt_map_get_seq_name(self._map);

    property strand:
        """The strand, either FORWARD or REVERSE"""
        def __get__(self):
            return gt_map_get_strand(self._map);



cdef class Template:
    """Wrapper class around gt_tempalte"""
    cdef gt_template* template

    def __cinit__(self):
        self.template = gt_template_new()

    def __dealloc__(self):
        gt_template_delete(self.template)

    property tag:
        def __get__(self):
            return gt_template_get_tag(self.template)
        def __set__(self, value):
            gt_template_set_tag(self.template, value, len(value))

    property pair:
        def __get__(self):
            return gt_template_get_pair(self.template)

    property blocks:
        def __get__(self):
            return gt_template_get_num_blocks(self.template)

    property num_alignments:
        def __get__(self):
            return gt_template_get_num_blocks(self.template)

    property counters:
        def __get__(self):
            return gt_template_get_num_counters(self.template)

    property num_maps:
        def __get__(self):
            return gt_template_get_num_mmaps(self.template)

    property mcs:
        def __get__(self):
            return gt_template_get_mcs(self.template)
        def __set__(self, value):
            gt_template_set_mcs(self.template, value)

    property has_qualities:
        def __get__(self):
            return gt_template_has_qualities(self.template)

    property not_unique_flag:
        def __get__(self):
            return gt_template_get_not_unique_flag(self.template)
        def __set__(self, value):
            gt_template_set_not_unique_flag(self.template, value)

    property length:
        def __get__(self):
            return self._length()

    property read:
        def __get__(self):
            return self._read()

    property qualities:
        def __get__(self):
            return self._qualities()

    cdef uint64_t _length(self):
        cdef gt_alignment* alignment = self._get_alignment(0)
        return gt_alignment_get_read_length(alignment)

    cdef char* _read(self):
        cdef gt_alignment* alignment = self._get_alignment(0)
        return gt_alignment_get_read(alignment)

    cdef char* _qualities(self):
        cdef gt_alignment* alignment = self._get_alignment(0)
        return gt_alignment_get_qualities(alignment)

    cdef gt_alignment* _get_alignment(self, uint64_t posistion):
        return gt_template_get_block(self.template, posistion)

    cpdef uint64_t get_counter(self, uint64_t stratum):
        return gt_template_get_counter(self.template, stratum)

    cpdef set_counter(self, uint64_t stratum, uint64_t value):
        gt_template_set_counter(self.template, stratum, value)

    cpdef uint64_t get_num_maps(self):
        return gt_template_get_num_mmaps(self.template)

    cpdef bool is_unmapped(self, uint64_t max_mismatches=GT_ALL):
        return not gt_template_is_thresholded_mapped(self.template, max_mismatches)

    def to_map(self):
        cdef gt_string* gt = self._to_map()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_fasta(self):
        cdef gt_string* gt = self._to_fasta()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_fastq(self):
        cdef gt_string* gt = self._to_fastq()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s

    def to_sequence(self):
        cdef gt_string* gt = self._to_sequence()
        s = ""+gt_string_get_string(gt)
        gt_string_delete(gt)
        return s


    cdef gt_string* _to_map(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_map_attributes* attr = gt_output_map_attributes_new()
        gt_output_map_sprint_template(s,self.template, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_map_attributes_delete(attr)
        return s

    cdef gt_string* _to_fasta(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_fasta_attributes* attr = gt_output_fasta_attributes_new()
        gt_output_fasta_attributes_set_format(attr, F_FASTA)
        gt_output_fasta_sprint_template(s,self.template, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_fasta_attributes_delete(attr)
        return s

    cdef gt_string* _to_fastq(self):
        cdef gt_string* s = gt_string_new(512)
        cdef gt_output_fasta_attributes* attr = gt_output_fasta_attributes_new()
        gt_output_fasta_sprint_template(s,self.template, attr)
        gt_string_set_length(s, gt_string_get_length(s)-1)
        gt_string_append_eos(s)
        gt_output_fasta_attributes_delete(attr)
        return s

    cdef gt_string* _to_sequence(self):
        if gt_template_has_qualities(self.template):
            return self._to_fastq()
        else:
            return self._to_fasta()

    cpdef int64_t get_min_mismatches(self):
        """Return minimum number of mismatches or -1 if
        the template is unmapped"""
        cdef gt_template* template = self.template
        cdef uint64_t counter = 0
        cdef uint64_t c = gt_template_get_num_counters(template)
        cdef int64_t i = 0
        for i in range(c):
            counter = gt_template_get_counter(template, i)
            if counter > 0:
                return i
        return -1

    cpdef merge(self, Template other):
        """Merge the other template into this one"""
        cdef gt_template* tmpl = gt_template_union_template_mmaps(self.template, other.template)
        gt_template_delete(self.template)
        self.template = tmpl

    cpdef get_pair(self):
        """Return 0 for unpaired or 1 or 2"""
        return gt_template_get_pair(self.template)

    cpdef same(self, Template other):
        return strcmp(gt_template_get_tag(self.template), gt_template_get_tag(other.template)) == 0 and self.get_pair() == other.get_pair()

    cpdef int64_t level(self, uint64_t max_level=GT_ALL):
        cdef gt_template* template = self.template
        cdef uint64_t counter = 0
        cdef uint64_t c = gt_template_get_num_counters(template)
        cdef int64_t i, j = 0
        cdef int64_t level = 0
        for i in range(c):
            counter = gt_template_get_counter(template, i)
            if counter == 1:
                for j in range(i+1, c):
                    if level >= max_level:
                        return level
                    counter = gt_template_get_counter(template, j)
                    if counter > 0:
                        return <int64_t> (j-(i+1))
                    else:
                        level += 1
                return <int64_t> (c - (i+1))
            elif counter > 1:
                return -1
        return -1

    cpdef parse(self, char* string):
        gt_input_map_parse_template(string, self.template)
        return self

    cpdef alignments(self):
        return [_create_alignment(gt_template_get_block(self.template, i)) for i in range(self.num_alignments)]

cdef class Stats(object):
    cdef gt_stats* stats
    cdef bool best_map
    cdef bool paired

    def __init__(self, bool best_map, bool paired):
        self.best_map =best_map
        self.paired = paired
        self.stats = gt_stats_new()

    def __dealloc__(self):
        gt_stats_delete(self.stats)

    cpdef read(self, input, uint64_t threads=1):
        if self.best_map:
            __calculate_stats(None, self, input, threads)
        else:
            __calculate_stats(self, None, input, threads)

    cpdef write(self, output):
        gt_stats_print_stats(PyFile_AsFile(output), self.stats, self.paired)

    property min_length:
        def __get__(self):
            return self.stats.min_length
    property max_length:
        def __get__(self):
            return self.stats.max_length
    property total_bases:
        def __get__(self):
            return self.stats.total_bases
    property total_bases_aligned:
        def __get__(self):
            return self.stats.total_bases_aligned
    property mapped_min_length:
        def __get__(self):
            return self.stats.mapped_min_length
    property mapped_max_length:
        def __get__(self):
            return self.stats.mapped_max_length
    property num_blocks:
        def __get__(self):
            return self.stats.num_blocks
    property num_reads:
        def __get__(self):
            if self.paired:
                return 2*self.num_blocks
            else:
                return self.num_blocks
    property num_alignments:
        def __get__(self):
            return self.stats.num_alignments
    property num_maps:
        def __get__(self):
            return self.stats.num_maps
    property num_mapped:
        def __get__(self):
            return self.stats.num_mapped
    property num_mapped_total:
        def __get__(self):
            if self.paired:
                return self.num_mapped * 2
            else:
                return self.num_mapped
    property nt_counting:
        def __get__(self):
            cdef uint64_t * c = self.stats.nt_counting
            return {"A": c[0], "C":c[1], "G":c[2], "T":c[3], "N":c[4]}
    property mmap:
        def __get__(self):
            cdef uint64_t * c = self.stats.mmap
            return [c[i] for i in range(8)]
    property mmap_description:
        def __get__(self):
            return ["0", "5", "10", "50", "100", "500", "1000", "more"]
    property uniq:
        def __get__(self):
            cdef uint64_t * c = self.stats.uniq
            return [c[i] for i in range(10)]
    property uniq_description:
        def __get__(self):
            return ["0", "1", "2", "3", "10", "50", "100", "500", "more", "rest"]
    property maps_profile:
        def __get__(self):
            return StatsMapProfile(self)
    property splits_profile:
        def __get__(self):
            return StatsSplitsProfile(self)

    property __dict__:
        def __get__(self):
            return {
                "min_length": self.min_length,
                "max_length": self.max_length,
                "total_bases": self.total_bases,
                "total_bases_aligned": self.total_bases_aligned,
                "mapped_min_length": self.mapped_min_length,
                "mapped_max_length": self.mapped_max_length,
                "num_blocks": self.num_blocks,
                "num_reads": self.num_reads,
                "num_alignments": self.num_alignments,
                "num_maps": self.num_maps,
                "num_mapped": self.num_mapped,
                "nt_counting": self.nt_counting,
                "mmap": self.mmap,
                "mmap_description": self.mmap_description,
                "uniq": self.uniq,
                "uniq_description": self.uniq_description,
                "maps_profile": self.maps_profile.__dict__,
                "splits_profile": self.splits_profile.__dict__,
            }


cdef class StatsSplitsProfile(object):
    cdef gt_splitmaps_profile* profile

    def __init__(self, Stats stats):
        self.profile = stats.stats.splitmaps_profile

    property num_mapped_with_splitmaps:
        def __get__(self):
            return self.profile.num_mapped_with_splitmaps
    property num_mapped_only_splitmaps:
        def __get__(self):
            return self.profile.num_mapped_only_splitmaps
    property total_splitmaps:
        def __get__(self):
            return self.profile.total_splitmaps
    property total_junctions:
        def __get__(self):
            return self.profile.total_junctions
    property num_junctions:
        def __get__(self):
            return [self.profile.num_junctions[i] for i in range(4)]
    property num_junctions_description:
        def __get__(self):
            return ["1", "2", "3", "more"]
    property length_junctions:
        def __get__(self):
            return [self.profile.length_junctions[i] for i in range(6)]
    property length_junctions_description:
        def __get__(self):
            return ["100", "1000", "5000", "10000", "50000", "more"]
    property junction_position:
        def __get__(self):
            return [self.profile.junction_position[i] for i in range(100)]
    property pe_sm_sm:
        def __get__(self):
            return self.profile.pe_sm_sm
    property pe_sm_rm:
        def __get__(self):
            return self.profile.pe_sm_rm
    property pe_rm_rm:
        def __get__(self):
            return self.profile.pe_rm_rm
    property __dict__:
        def __get__(self):
            return {
                "num_mapped_with_splitmaps": self.num_mapped_with_splitmaps,
                "num_mapped_only_splitmaps": self.num_mapped_only_splitmaps,
                "total_splitmaps": self.total_splitmaps,
                "total_junctions": self.total_junctions,
                "num_junctions": self.num_junctions,
                "num_junctions_description": self.num_junctions_description,
                "length_junctions": self.length_junctions,
                "length_junctions_description": self.length_junctions_description,
                "junction_position": self.junction_position,
                "pe_sm_sm": self.pe_sm_sm,
                "pe_sm_rm": self.pe_sm_rm,
                "pe_rm_rm": self.pe_rm_rm,
            }



cdef class StatsMapProfile(object):
    cdef gt_maps_profile* profile

    def __init__(self, Stats stats):
        self.profile = stats.stats.maps_profile

    property mismatches_description:
        def __get__(self):
            return ["1", "2","3", "4", "5", "6", "7", "8", "9", "10", "20", "50", "more"]
    property mismatches:
        def __get__(self):
            return [self.profile.mismatches[i] for i in range(14)]
    property levenshtein:
        def __get__(self):
            return [self.profile.levenshtein[i] for i in range(14)]
    property insertion_length:
        def __get__(self):
            return [self.profile.insertion_length[i] for i in range(14)]
    property deletion_length:
        def __get__(self):
            return [self.profile.deletion_length[i] for i in range(14)]
    property errors_events:
        def __get__(self):
            return [self.profile.errors_events[i] for i in range(14)]
    property total_mismatches:
        def __get__(self):
            return self.profile.total_mismatches
    property total_levenshtein:
        def __get__(self):
            return self.profile.total_levenshtein
    property total_indel_length:
        def __get__(self):
            return self.profile.total_indel_length
    property total_errors_events:
        def __get__(self):
            return self.profile.total_errors_events
    property error_position:
        def __get__(self):
            return [self.profile.error_position[i] for i in range(1000)]
    property total_bases:
        def __get__(self):
            return self.profile.total_bases
    property total_bases_matching:
        def __get__(self):
            return self.profile.total_bases_matching
    property total_bases_trimmed:
        def __get__(self):
            return self.profile.total_bases_trimmed
    property inss_fine_grain:
        def __get__(self):
            return [self.profile.inss_fine_grain[i] for i in range(GT_STATS_INSS_FG_RANGE)]
    property inss:
        def __get__(self):
            return [self.profile.inss[i] for i in range(16)]
    property inss_description:
        def __get__(self):
            return ["(-inf, 0)", "(-100, 0)", "(0, 100]", "(100, 200]", "(200, 300]", "(300, 400]", "(400, 500]", "(500, 600]", "(600, 700]", "(700, 800]", "(800, 900]", "(900, 1000]", "(1000, 2000]", "(2000, 5000]", "(5000, 10000]", "(1000, inf]"]
    property misms_transition:
        def __get__(self):
            return [self.profile.misms_transition[i] for i in range(5*14)] # GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_RANGE */
    property qual_score_misms:
        def __get__(self):
            return [self.profile.qual_score_misms[i] for i in range(256)]
    property qual_score_errors:
        def __get__(self):
            return [self.profile.qual_score_errors[i] for i in range(256)]
    property misms_1context:
        def __get__(self):
            return [self.profile.misms_1context[i] for i in range(GT_STATS_MISMS_1_CONTEXT_RANGE)]
    property misms_2context:
        def __get__(self):
            return [self.profile.misms_2context[i] for i in range(GT_STATS_MISMS_2_CONTEXT_RANGE)]
    property indel_transition_1:
        def __get__(self):
            return [self.profile.indel_transition_1[i] for i in range(GT_STATS_INDEL_TRANSITION_1_RANGE)]
    property indel_transition_2:
        def __get__(self):
            return [self.profile.indel_transition_2[i] for i in range(GT_STATS_INDEL_TRANSITION_2_RANGE)]
    property indel_transition_3:
        def __get__(self):
            return [self.profile.indel_transition_3[i] for i in range(GT_STATS_INDEL_TRANSITION_3_RANGE)]
    property indel_transition_4:
        def __get__(self):
            return [self.profile.indel_transition_4[i] for i in range(GT_STATS_INDEL_TRANSITION_4_RANGE)]
    property indel_1context:
        def __get__(self):
            return [self.profile.indel_1context[i] for i in range(GT_STATS_INDEL_1_CONTEXT)]
    property indel_2context:
        def __get__(self):
            return [self.profile.indel_2context[i] for i in range(GT_STATS_INDEL_2_CONTEXT)]
    property __dict__:
        def __get__(self):
            return {
                "mismatches_description": self.mismatches_description,
                "mismatches": self.mismatches,
                "levenshtein": self.levenshtein,
                "insertion_length": self.insertion_length,
                "deletion_length": self.deletion_length,
                "errors_events": self.errors_events,
                "total_mismatches": self.total_mismatches,
                "total_levenshtein": self.total_levenshtein,
                "total_indel_length": self.total_indel_length,
                "total_errors_events": self.total_errors_events,
                "error_position": self.error_position,
                "total_bases": self.total_bases,
                "total_bases_matching": self.total_bases_matching,
                "total_bases_trimmed": self.total_bases_trimmed,
                "inss": self.inss,
                "inss_fine_grain": self.inss_fine_grain,
                "inss_description": self.inss_description,
                "misms_transition": self.misms_transition,
                "qual_score_misms": self.qual_score_misms,
                "qual_score_errors": self.qual_score_errors,
                "misms_1context": self.misms_1context,
                "misms_2context": self.misms_2context,
                "indel_transition_1": self.indel_transition_1,
                "indel_transition_2": self.indel_transition_2,
                "indel_transition_3": self.indel_transition_3,
                "indel_transition_4": self.indel_transition_4,
                "indel_1context": self.indel_1context,
                "indel_2context": self.indel_2context,
            }


cpdef read_stats(source, Stats all=None, Stats best=None, uint64_t threads=1):
    """Calculate stats from input and do this optionally for two stats
    at once, one for all mappings, and one for only the best mappings.
    """
    #__run_stats(all, best, source, threads)
    __calculate_stats(all, best, source, threads)


cpdef __calculate_stats(Stats all_stats, Stats best_stats, source, uint64_t threads=1):
    # import gem.utils
    # process = multiprocessing.Process(target=__calculate_stats_process, args=(all_stats, best_stats, source, threads))
    # gem.utils.register_process(process)
    # process.start()
    # process.join()
    __calculate_stats_process(all_stats, best_stats, source, threads)


cpdef __calculate_stats_process(Stats all_stats, Stats best_stats, source, uint64_t threads=1):
    cdef gt_input_file* input = <gt_input_file*> (<InputFile>source)._open()
    cdef uint64_t use_threads = threads
    cdef gt_stats* target_all = NULL
    cdef gt_stats* target_best = NULL
    cdef bool paired = True
    if all_stats is not None:
        target_all = all_stats.stats
        paired = all_stats.paired
    if best_stats is not None:
        target_best = best_stats.stats
        paired = best_stats.paired

    with nogil:
        gt_stats_fill(input, target_all, target_best, threads, paired)



cpdef __write_filter(source, OutputFile output, uint64_t threads=1, params=None):
    cdef gt_output_file* output_file = output.output_file
    cdef gt_input_file* input_file = (<InputFile> source)._open()
    cdef uint64_t use_threads = threads
    if params is None:
        params = {}

    cdef gt_filter_params p
    p.max_matches = params.get("max_matches", 0)
    p.min_event_distance = params.get("min_event_distance", 0)
    p.max_event_distance = params.get("max_event_distance", UINT64_MAX)
    p.min_levenshtein_distance = params.get("min_levenshtein_distance", 0)
    p.max_levenshtein_distance = params.get("max_levenshtein_distance", UINT64_MAX)
    p.max_inss = params.get("max_inss", INT64_MAX)
    p.min_inss = params.get("min_inss", INT64_MIN);
    p.filter_by_strand = params.get("filter_strand", False)
    p.keep_unique = params.get("keep_unique", False)
    p.min_score = params.get("min_score", 0)
    p.filter_groups = params.get("filter_groups", False)
    p.group_1 = params.get("group_1", False)
    p.group_2 = params.get("group_2", False)
    p.group_3 = params.get("group_3", False)
    p.group_4 = params.get("group_4", False)
    p.close_output = params.get("close_output", True)
    if "annotation" in params and params["annotation"] is not None:
        p.annotation = params["annotation"]
    else:
        p.annotation = ""
    with nogil:
        gt_filter_stream(input_file, output_file, use_threads, &p)
    gt_input_file_close(input_file)


cpdef filter_map(input, OutputFile output, params, uint64_t threads=1,
                 background_process=True):
    if background_process:
        import gem.utils
        process = multiprocessing.Process(target=__write_filter, args=(input, output, threads, params))
        gem.utils.register_process(process)
        process.start()
        process.join()
        return process
    else:
        params["close_output"] = False
        __write_filter(input, output, threads, params)
        return None

