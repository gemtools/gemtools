from gemapi cimport *


cdef class TemplateIterator:
    """Base class for iterating template streams"""
    cdef gt_buffered_input_file* buffered_input
    cdef gt_generic_parser_attr* parser_attr
    cdef gt_input_file* input_file
    cdef Template template
    cdef TemplateIterator source
    cdef readonly object quality

    def __iter__(self):
        return self

    def __next__(self):
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

    cdef gt_status _next(self):
        pass

    def close(self):
        pass

cdef class AlignmentIterator:
    """Base class to iterate alignments"""
    cdef gt_buffered_input_file* buffered_input
    cdef gt_generic_parser_attr* parser_attr
    cdef gt_input_file* input_file
    cdef Alignment alignment
    cdef AlignmentIterator source

    def __iter__(self):
        return self

    cdef _init(self, AlignmentIterator source):
        """Initialize from another alignment iterator"""
        self.input_file = source.input_file
        self.buffered_input = source.buffered_input
        self.parser_attr = source.parser_attr
        self.template = source.template
        self.source = source

    def __next__(self):
        if self._next() == GT_STATUS_OK:
            return self.alignment
        else:
            raise StopIteration()

    cdef gt_status _next(self):
        pass

cdef class SimpleTemplateIterator(TemplateIterator):
    """Simple iterator that iterates a list of templates"""
    cdef object templates
    cdef int64_t index
    cdef int64_t length

    def __init__(self, templates):
        self.templates = templates
        self.index = 0
        self.length = len(templates)

    cdef gt_status _next(self):
        if self.index < self.length:
            self.index += 1
            self.template = self.templates[self.index]
            print "TEMPLATE SET RETURN OKEY"
            return GT_STATUS_OK
        print "FINISH ITERATING"
        return GT_STATUS_FAIL


cdef class interleave(TemplateIterator):
    """Pairwise interleave of two streams of equal size"""
    cdef object iterators
    cdef uint64_t current
    cdef uint64_t length

    def __init__(self, iterators):
        self.iterators = iterators
        self.length = len(iterators)
        self.current = 0
        cdef TemplateIterator main = <TemplateIterator> self.iterators[0]
        self._init(main)

    cdef gt_status _next(self):
        cdef TemplateIterator ci = <TemplateIterator> self.iterators[self.current]
        self.template = ci.template
        cdef gt_status s = ci._next()
        self.current += 1
        if self.current >= self.length:
            self.current = 0
        return s

cdef class unmapped(TemplateIterator):
    """Filter for unmapped reads or reads with mismatches >= max_mismatches"""
    cdef TemplateIterator iterator
    cdef int64_t max_mismatches


    def __init__(self, iterator, int64_t max_mismatches=GT_ALL):
        self.iterator = iterator
        self.max_mismatches = max_mismatches
        self._init(iterator)

    cdef gt_status _next(self):
        cdef gt_status s = self.iterator._next()
        while( s == GT_STATUS_OK ):
            if self.template.is_unmapped(self.max_mismatches):
                return s
            s = self.iterator._next()
        return GT_STATUS_FAIL

cdef class unique(TemplateIterator):
    """Pairwise interleave of two streams of equal size"""
    cdef TemplateIterator iterator
    cdef uint64_t level


    def __init__(self, iterator, uint64_t level=0):
        self.iterator = iterator
        self.level = level
        self._init(iterator)

    cdef gt_status _next(self):
        cdef gt_status s = self.iterator._next()
        cdef int64_t template_level = 0
        while( s == GT_STATUS_OK ):
            template_level = self.template.level(self.level)
            if template_level > 0 and template_level >= self.level:
                return s
            s = self.iterator._next()
        return GT_STATUS_FAIL


cdef class merge(TemplateIterator):
    """Merge a master stream with a list of child streams"""
    cdef TemplateIterator master
    cdef object children

    def __init__(self, TemplateIterator master, children):
        self._init(master)
        self.children = []
        # initialize children
        for c in children:
            i = c.__iter__()
            children.append(i)

    cdef gt_status _next(self):
        cdef gt_status s = self.source._next()
        while( s == GT_STATUS_OK ):
            for c in self.children:
                if self.template.same(c.template):
                    self.template.merge(c.template)
                    c._next()
            return s


cdef class OutputFile:
    """Wrapper around gt_output_file"""
    cdef gt_output_file* output_file
    cdef gt_buffered_output_file* buffered_output
    cdef bool do_close

    def __init__(self,file_name=None, file stream=None):
        self.do_close = True
        if file_name is not None:
            self.open_file(<char*> file_name)
        elif stream is not None:
            self.open_stream(stream)

    def open_stream(self, file stream):
        self.output_file = gt_output_stream_new(PyFile_AsFile(stream), SORTED_FILE)
        self.do_close = False
        self.buffered_output = gt_buffered_output_file_new(self.output_file)
        return self

    def open_file(self, char* file_name):
        self.output_file = gt_output_file_new(file_name, SORTED_FILE)
        self.buffered_output = gt_buffered_output_file_new(self.output_file)
        return self

    def close(self):
        if self.buffered_output is not NULL:
            gt_buffered_output_file_close(self.buffered_output)
        if self.do_close:
            gt_output_file_close(self.output_file)

    def write_fastq(self, iterator, clean_id=False, append_extra=True):
        self._write_fastq(iterator, clean_id, append_extra)

    def write_map(self, iterator, print_scores=True, clean_id=False, append_extra=True):
        self._write_map(iterator, print_scores, clean_id, append_extra)

    cdef void _write_fastq(self, TemplateIterator iterator, bool clean_id, bool append_extra):
        """Write templated from iterator as fastq"""
        gt_buffered_input_file_attach_buffered_output(iterator.buffered_input, self.buffered_output)
        cdef gt_output_fasta_attributes* attributes = gt_output_fasta_attributes_new()
        gt_output_fasta_attributes_set_print_extra(attributes, append_extra)
        gt_output_fasta_attributes_set_print_casava(attributes, not clean_id)
        if(iterator._next() == GT_STATUS_OK):
            ## check fastq / fasta by qualities
            if not iterator.template.has_qualities:
                gt_output_fasta_attributes_set_format(attributes, F_FASTA)

            ## print the first one
            gt_output_fasta_bofprint_template(self.buffered_output, iterator.template.template, attributes)
            # print the rest
            while(iterator._next() == GT_STATUS_OK):
                gt_output_fasta_bofprint_template(self.buffered_output, iterator.template.template, attributes)
        self.close()

    cdef void _write_map(self, TemplateIterator iterator, bool print_scores, bool clean_id, bool append_extra):
        """Write templated from iterator as fastq"""
        gt_buffered_input_file_attach_buffered_output(iterator.buffered_input, self.buffered_output)
        cdef gt_output_map_attributes* attributes = gt_output_map_attributes_new()
        gt_output_map_attributes_set_print_scores(attributes, print_scores)
        gt_output_map_attributes_set_print_casava(attributes, not clean_id)
        gt_output_map_attributes_set_print_extra(attributes, append_extra)
        while(iterator._next() == GT_STATUS_OK):
            gt_output_map_bofprint_template(self.buffered_output, iterator.template.template, attributes)
        self.close()


cdef class InputFile(TemplateIterator):
    """GEMTools input file. The input file extends TemplateIterator
    and can be used directly to iterate templates. To access the underlying
    alignment, use the alignments() function to get an alignment iterator.
    """
    cdef readonly object file_name
    cdef readonly bool force_paired_reads
    cdef readonly bool mmap_file
    cdef readonly object process
    cdef readonly object stream
    cdef InputFileTemplateIterator iterator

    def __init__(self,
        file_name=None,
        object stream=None,
        bool mmap_file=True,
        bool force_paired_reads=False,
        object quality=None,
        object process=None
        ):
        if file_name is not None:
            self.file_name = file_name
        else:
            self.file_name = None

        self.force_paired_reads = force_paired_reads
        self.mmap_file = mmap_file
        self.process = process
        self.stream = stream
        self.quality = quality
        if file_name is None or file_name.endswith(".gz") or file_name.endswith(".bz2"):
            self.mmap_file = False
        self.iterator = self._templates()

    cdef gt_input_file* _input_file(self):
        if self.stream is not None:
                return gt_input_stream_open(PyFile_AsFile(self.stream))
        else:
            if self.file_name is not None:
                return gt_input_file_open(<char*>self.file_name, self.mmap_file)


    def raw_stream(self):
        if self.stream is not None:
            return self.stream
        else:
            return open(self.file_name, "rb")

    def __iter__(self):
        return self._templates()

    cdef gt_status _next(self):
        return self.iterator._next()

    def templates(self):
        return self._templates()

    cdef InputFileTemplateIterator _templates(self):
        return InputFileTemplateIterator(self)

    def alignments(self):
        return self._alignments()

    cdef InputFileAlignmentIterator _alignments(self):
        return InputFileAlignmentIterator(self)


cdef class InputFileTemplateIterator(TemplateIterator):

    def __cinit__(self, InputFile input_file):
        self.input_file = input_file._input_file()
        cdef gt_generic_parser_attr* gp = gt_input_generic_parser_attributes_new(input_file.force_paired_reads)
        self.parser_attr = gp
        self.template = Template()
        self.quality = input_file.quality
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

    cdef gt_status _next(self):
        """Internal iterator method"""
        return gt_input_generic_parser_get_template(self.buffered_input, self.template.template, self.parser_attr)


cdef class InputFileAlignmentIterator(AlignmentIterator):

    def __cinit__(self, InputFile input_file):
        self.input_file = input_file._input_file()
        cdef gt_generic_parser_attr* gp = gt_input_generic_parser_attributes_new(False)
        self.parser_attr = gp
        self.alignment = Alignment()
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


    cdef gt_status _next(self):
        """Internal iterator method"""
        return gt_input_generic_parser_get_alignment(self.buffered_input, self.alignment.alignment, self.parser_attr)


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
        return gt_template_is_thresholded_mapped(self.template, max_mismatches)

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

    cpdef int64_t get_pair(self):
        """Return 0 for unpaired or 1 or 2"""
        return gt_template_get_pair(self.template)

    cpdef same(self, Template other):
        return gt_template_get_tag(self.template) == gt_template_get_tag(other.template) and self.get_pair() == other.get_pair()

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


cdef class Alignment:
    """Wrapper class around alignments"""
    cdef gt_alignment* alignment

    def __cinit__(self):
        self.alignment = gt_alignment_new()

    def __dealloc__(self):
        gt_alignment_delete(self.alignment)

    property tag:
        def __get__(self):
            return gt_alignment_get_tag(self.alignment)
        def __set__(self, value):
            gt_alignment_set_tag(self.alignment, value, len(value))

