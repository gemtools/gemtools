from gemapi cimport *
import os
import sys


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
    # interleave
    cdef bool interleave

    def __init__(self, files, interleave=True):
        self.files = files
        self.i = 0
        self.length = len(files)
        self.interleave = interleave

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

        cdef gt_output_file* output_file = output.output_file
        cdef uint64_t num_inputs = len(self.files)
        cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        cdef bool clean_id = output.clean_id
        cdef bool append_extra = output.append_extra
        cdef bool interleave = self.interleave
        for i in range(num_inputs):
            inputs[i] = (<InputFile> self.files[i])._open()
        with nogil:
            gt_write_stream(output_file, inputs, num_inputs, append_extra, clean_id, interleave, threads, write_map)
        output.close()
        for i in range(num_inputs):
            gt_input_file_close(inputs[i])
        free(inputs)



cdef class cat(interleave):
    def __init__(self, files):
        interleave.__init__(self, files, interleave=False)


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

    cpdef write_stream(self, OutputFile output, bool write_map=False, uint64_t threads=1):
        """Write the content of this input stream to the output file
        file.

        output   -- the output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        interleave    -- interleave muliple inputs
        threads       -- number of threads to use (if supported by the iterator)
        """
        cdef gt_output_file* output_file = output.output_file
        cdef uint64_t num_inputs = len(self.inputs)
        cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        cdef bool clean_id = output.clean_id
        cdef bool append_extra = output.append_extra
        for i in range(num_inputs):
            inputs[i] = (<InputFile> self.inputs[i])._open()

        with nogil:
            gt_merge_files_synch(output_file, threads, num_inputs, inputs);
        output.close()


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
        cdef gt_output_file* output_file = output.output_file
        cdef uint64_t num_inputs = 1
        cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        cdef bool clean_id = output.clean_id
        cdef bool append_extra = output.append_extra
        inputs[0] = self._open()

        with nogil:
            gt_write_stream(output_file, inputs, num_inputs, append_extra, clean_id, True, threads, write_map)
        output.close()
        gt_input_file_close(inputs[0])
        free(inputs)




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
