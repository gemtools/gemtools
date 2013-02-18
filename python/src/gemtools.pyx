from gemapi cimport *
import os
import sys

cdef class TemplateIterator:
    """Base class for iterating template streams.
    The base implementation can also write to
    stream, but rather slow. If a specific implementation
    can speed up the writing, override the write() method.

    NOTE that the __next__ implementation reuses the
    iterators Template instance.
    """
    # the buffered input file
    cdef gt_buffered_input_file* buffered_input
    # parser attributes
    cdef gt_generic_parser_attr* parser_attr
    # the source input file
    cdef gt_input_file* input_file
    # the template instance that is used to iterate templates
    cdef readonly Template template
    # optional source iterator
    cdef TemplateIterator source
    # the quality offset of the input
    cdef readonly object quality
    # the source file name
    cdef readonly object filename

    def __iter__(self):
        """Returns self"""
        return self

    def __next__(self):
        """Fill the template instance and return it"""
        if self._next() == GT_STATUS_OK:
            return self.template
        else:
            raise StopIteration()

    cdef _init(self, TemplateIterator source):
        """Initialize from another template iterator"""
        self.input_file = source.input_file
        self.buffered_input = source.buffered_input
        self.parser_attr = source.parser_attr
        self.template = source.template
        self.source = source
        self.quality = source.quality
        self.filename = source.filename

    cpdef gt_status _next(self):
        """Implement this to read the next template from a
        stream
        """
        pass

    cpdef write_stream(self, output_file, write_map=False, clean_id=False, append_extra=True, threads=1):
        """Write the content of the this template iterator to the output file.

        output_file   -- the target output file
        write_map     -- if true, write map, otherwise write fasta/q sequence
        clean_id      -- clean id's and ensure /1 /2 endings for paired reads
        append_extra  -- append any additional information to the read id
        threads       -- number of threads to use (if supported by the iterator)
        """

    cpdef write(self, output, bool clean_id=False, bool append_extra=True, bool write_map=False, uint64_t threads=1):
        """Override this if there is a faster implementation"""
        cdef OutputFile of = OutputFile(output)

        if write_map:
            of.write_map(self, clean_id=clean_id, append_extra=append_extra)
        else:
            of.write_fastq(self, clean_id=clean_id, append_extra=append_extra)
        if isinstance(output, basestring):
            of.close()

    def close(self):
        """Close the template iterator. If a source is
        set, this closes the source iterator.
        """
        if self.source is not None:
            self.source.close()


cdef class SimpleTemplateIterator(TemplateIterator):
    """Simple iterator that iterates a list of templates"""
    cdef object templates
    cdef int64_t index
    cdef int64_t length

    def __init__(self, templates):
        self.templates = templates
        self.index = 0
        self.length = len(templates)

    cpdef gt_status _next(self):
        if self.index < self.length:
            self.index += 1
            self.template = self.templates[self.index]
            return GT_STATUS_OK
        return GT_STATUS_FAIL


cdef class interleave(TemplateIterator):
    """Pairwise interleave of two streams of equal size"""
    cdef object iterators
    cdef uint64_t current
    cdef uint64_t length

    def __init__(self, iterators):
        self.iterators = [i.__iter__() for i in iterators]
        self.length = len(iterators)
        self.current = 0
        cdef TemplateIterator main = <TemplateIterator> self.iterators[0]
        self._init(main)

    cpdef gt_status _next(self):
        cdef TemplateIterator ci = <TemplateIterator> self.iterators[self.current]
        self.template = ci.template
        cdef gt_status s = ci._next()
        self.current += 1
        if self.current >= self.length:
            self.current = 0
        return s

    cpdef write_fastq(self, output, bool clean_id, bool append_extra, uint64_t threads):
        cdef OutputFile of
        if isinstance(output, basestring):
            of = OutputFile(file_name=output)
        else:
            of = OutputFile(stream=output)

        cdef gt_output_file* output_file = of.output_file
        cdef uint64_t num_inputs = len(self.iterators)
        cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        for i in range(num_inputs):
            inputs[i] = (<TemplateIterator> self.iterators[i]).input_file
        with nogil:
            gt_write_sequence_stream(output_file, inputs, num_inputs, append_extra, clean_id, True, threads)
        #if isinstance(output, basestring):
        of.close()


cdef class cat(TemplateIterator):
    """Cat iterators"""
    cdef object iterators
    cdef uint64_t current
    cdef uint64_t length
    cdef TemplateIterator current_itereator

    def __init__(self, iterators):
        self.iterators = [i.__iter__() for i in iterators]
        self.length = len(iterators)
        self.current = 0
        self.current_itereator = <TemplateIterator> self.iterators[0]
        self._init(self.current_itereator)
        self.template = self.current_itereator.template

    cpdef gt_status _next(self):
        cdef gt_status s = self.current_itereator._next()
        if s != GT_STATUS_OK:
            if self.current < self.length - 1:
                self.current += 1
                self.current_itereator = <TemplateIterator> self.iterators[self.current]
                self.template = self.current_itereator.template
                return self._next()
            else:
                return GT_STATUS_FAIL
        return s


cdef class unmapped(TemplateIterator):
    """Filter for unmapped reads or reads with mismatches >= max_mismatches"""
    cdef TemplateIterator iterator
    cdef int64_t max_mismatches


    def __init__(self, iterator, int64_t max_mismatches=GT_ALL):
        self.iterator = iterator.__iter__()
        self.max_mismatches = max_mismatches
        self._init(iterator)
        if max_mismatches < 0:
            max_mismatches = GT_ALL

    cpdef gt_status _next(self):
        cdef gt_status s = self.iterator._next()
        while( s == GT_STATUS_OK ):
            if self.template.is_unmapped(self.max_mismatches):
                return s
            s = self.iterator._next()
        return GT_STATUS_FAIL

cdef class unique(TemplateIterator):
    """Yield unique reads"""
    cdef TemplateIterator iterator
    cdef uint64_t level


    def __init__(self, iterator, uint64_t level=0):
        self.iterator = iterator.__iter__()
        self.level = level
        self._init(iterator)

    cpdef gt_status _next(self):
        cdef gt_status s = self.iterator._next()
        cdef int64_t template_level = 0
        while( s == GT_STATUS_OK ):
            template_level = self.template.level(self.level)
            if template_level > 0 and template_level >= self.level:
                return s
            s = self.iterator._next()
        return GT_STATUS_FAIL


cdef class trim(TemplateIterator):
    """Trim reads"""
    cdef TemplateIterator iterator
    cdef uint64_t left
    cdef uint64_t right
    cdef uint64_t min_length
    cdef bool set_extra

    def __init__(self, iterator, uint64_t left=0, uint64_t right=0, uint64_t min_length=10, bool set_extra=True):
        self.iterator = iterator.__iter__()
        self.left = left
        self.right = right
        self.min_length = min_length
        self.set_extra = set_extra
        self._init(iterator)

    cpdef gt_status _next(self):
        if self.iterator._next() == GT_STATUS_OK:
            if self.left > 0 or self.right > 0:
                # trim the read
                gt_template_trim(self.template.template, self.left, self.right, self.min_length, self.set_extra)
            return GT_STATUS_OK
        return GT_STATUS_FAIL


cdef class merge(TemplateIterator):
    """Merge a master stream with a list of child streams"""
    cdef object children
    cdef object states
    cdef object merge_master
    cdef object merge_children

    def __init__(self, master, children, bool init=True):
        if init:
            self._init(master.__iter__())
            self.children = [c.__iter__() for c in children]
            self.states = []
            for c in self.children:
                s = c._next()
                self.states.append(s)
        else:
            self.merge_master = master
            self.merge_children = children


    cpdef merge_synch(self, output, uint64_t threads):
        cdef gt_output_file* output_file

        if isinstance(output, basestring):
            output_file = gt_output_file_new(<char*>output, SORTED_FILE)
        else:
            output_file = gt_output_stream_new(PyFile_AsFile(output), SORTED_FILE)

        cdef gt_input_file* master = gt_input_file_open(self.merge_master.file_name, False)
        cdef uint64_t num_slaves = len(self.merge_children) + 1
        cdef gt_input_file** slaves = <gt_input_file**>malloc( num_slaves *sizeof(gt_input_file*))
        slaves[0] = master
        for i in range(len(self.merge_children)):
            slaves[i+1] = gt_input_file_open(self.merge_children[i].file_name, False)

        with nogil:
            gt_merge_files_synch(output_file, threads, num_slaves, slaves)

        for i in range(num_slaves):
            gt_input_file_close(slaves[i])
        free(slaves)
        gt_output_file_close(output_file)



    cpdef gt_status _next(self):
        if self.source._next() == GT_STATUS_OK:
            for i, s in enumerate(self.states):
                other = self.children[i]
                if s == GT_STATUS_OK and self.template.same(other.template):
                    self.template.merge(other.template)
                    self.states[i] = other._next()
            return GT_STATUS_OK
        return GT_STATUS_FAIL


cdef class OutputFile:
    """The OutputFile can write content to a file or
    or a stream.
    """
    # the target output file
    cdef gt_output_file* output_file
    # buffered output
    cdef gt_buffered_output_file* buffered_output
    # the target
    cdef readonly object target

    def __init__(self, target):
        """Initialize the output file from the given target. The
        target can be either a string a stream. If init_buffer is
        true, the output buffer is initialized

        target      -- the target file or stream
        init_buffer -- if true, the buffered output will be initialized, default True
        """
        self.target = target
        if isinstance(target, basestring):
            self._open_file(<char*> target)
        else:
            self._open_stream(stream)

        if init_buffer:
            self.buffered_output = gt_buffered_output_file_new(self.output_file)

    def __dealloc__(self):
        self.close()

    cpdef init_buffer(self):
        """Ensure that the output buffer is initialized"""
        if self.buffered_output is NULL:
            self.buffered_output = gt_buffered_output_file_new(self.output_file)

    cpdef _open_stream(self, file stream, init_buffer=True):
        """Initialize this instance from a stream

        stream      -- the output stream
        init_buffer -- if true, initialize the output buffer
        """
        self.output_file = gt_output_stream_new(PyFile_AsFile(stream), SORTED_FILE)

    cpdef def _open_file(self, char* file_name, init_buffer=True):
        """Initialize this instance from a file

        file_name   -- the output file name
        init_buffer -- if true, initialize the output buffer
        """
        self.output_file = gt_output_file_new(file_name, SORTED_FILE)

    cpdef close(self):
        """Close the output file"""
        if self.buffered_output is not NULL:
            gt_buffered_output_file_close(self.buffered_output)
            gt_buffered_output_file_delete(self.buffered_output)
            self.buffered_output = NULL
        if self.output_file is not NULL:
            gt_output_file_close(self.output_file)
            gt_output_file_delete(self.output_file)
            self.output_file = NULL

    cpdef write_stream(self, iterator, write_map=False, clean_id=False, append_extra=True, threads=1):
        """Write the content of the given template iterator to this output
        file.

        iterator      -- the template iterator providing the source templates
        write_map     -- if true, write map, otherwise write fasta/q sequence
        clean_id      -- clean id's and ensure /1 /2 endings for paired reads
        append_extra  -- append any additional information to the read id
        threads       -- number of threads to use (if supported by the iterator)
        """
        with nogil:
            iterator.write(self, clean_id=clean_id, append_extra=append_extra, write_map=write_map, threads=threads)



cdef class InputFile(object):
    """GEMTools input file. The input file extends TemplateIterator
    and can be used directly to iterate templates. To access the underlying
    alignment, use the alignments() function to get an alignment iterator.
    """
    cdef readonly object source
    # the source file name if specified
    cdef readonly char* filename
    # force reading paired reads
    cdef readonly bool force_paired_reads
    # memory map the file
    cdef readonly bool mmap_file
    # delete the file after iteraation
    cdef readonly bool delete_after_iterate
    # the process that creates the content of the file or stream
    cdef readonly object process
    # the quality offset
    cdef readonly object quality

    def __init__(self, source, bool mmap_file=False, bool force_paired_reads=False, object quality=None, object process=None, bool delete_after_iterate=False):
        """Initialize a new input file on the source. The source can be either a file name
        or a stream.

        """
        self.source = source
        self.force_paired_reads = force_paired_reads
        self.mmap_file = mmap_file
        self.process = process
        self.quality = quality
        self.delete_after_iterate = delete_after_iterate

        if isinstance(source, basestring):
            self.filename = <char*> source

        # make sure memory mapping is disabled for compressed files and
        # # streams
        if self.filename is None or self.filename.endswith(".gz") or self.filename.endswith(".bz2"):
                self.mmap_file = False

    cdef gt_input_file* open(self):
        """Open a new gt_input_file"""
        if self.filename is not None:
            return gt_input_file_open(self.filename, self.mmap_file)
        else:
            return gt_input_stream_open(PyFile_AsFile(self.source))

    def raw_sequence_stream(self):
        """Return true if this is a file based
        input file, uncompressed and fastq/q format
        """
        if self.filename is None:
            return False
        cdef gt_input_file* infile = self.open()
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
            return InputFile(self.source, mmap_file=self.mmap_file, force_paired_reads=self.force_paired_reads, quality=self.quality, process=self.process, delete_after_iterate=self.delete_after_iterate)

    def raw_stream(self):
        """Return the raw stream on this input file.
        In case this is stream based, the stream is returned,
        otherwise a new file handle is opened on
        the input file.
        """
        if self.filename is None:
            return PyFile_AsFile(self.source)
        else:
            return open(self.file_name, "rb")

    def __iter__(self):
        """Returns a new template iterator in the input"""
        return self._templates()

    cpdef InputFileTemplateIterator templates(self):
        """Returns a new template iterator in the input"""
        return InputFileTemplateIterator(self)


cdef class InputFileTemplateIterator(TemplateIterator):
    # delete the source file after iteration
    cdef bool delete_after_iterate
    # underlying process
    cdef object process

    def __cinit__(self, InputFile input_file):
        self.input_file = input_file._input_file()
        if input_file.file_name is not None:
            self.file_name = input_file.file_name

        cdef gt_generic_parser_attr* gp = gt_input_generic_parser_attributes_new(input_file.force_paired_reads)
        self.parser_attr = gp
        self.template = Template()
        self.quality = input_file.quality
        self.delete_after_iterate = input_file.delete_after_iterate
        self.process = input_file.process
        #
        #initialize the buffer
        #
        self.buffered_input = gt_buffered_input_file_new(self.input_file)


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

    cpdef gt_status _next(self):
        """Internal iterator method"""
        cdef gt_status s = gt_input_generic_parser_get_template(self.buffered_input, self.template.template, self.parser_attr)

        if s != GT_STATUS_OK :
            if self.delete_after_iterate:
                os.remove(self.file_name)
            if self.process is not None:
                self.process.wait()
        return s

    cpdef write_fastq(self, output, bool clean_id, bool append_extra, uint64_t threads):
        cdef OutputFile of
        if isinstance(output, basestring):
            of = OutputFile(file_name=output)
        else:
            of = OutputFile(stream=output)

        cdef gt_output_file* output_file = of.output_file
        cdef uint64_t num_inputs = 1
        cdef gt_input_file** inputs = <gt_input_file**>malloc( num_inputs *sizeof(gt_input_file*))
        inputs[0] = self.input_file
        with nogil:
            gt_write_sequence_stream(output_file, inputs, num_inputs, append_extra, clean_id, True, threads)
        #if isinstance(output, basestring):
        of.close()


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
