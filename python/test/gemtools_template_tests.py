#!/usr/bin/env python

import gem.gemtools as gt
from testfiles import testfiles

test_mapping = testfiles["test.map"]
test_zipped_mapping = testfiles["test.map.gz"]
test_fastq = testfiles["test.fastq"]


def test_template_attribute_reading_first_mapping():
    infile = gt.InputFile(testfiles["test.map"])
    for tmpl in infile:
        print "Template"
        assert tmpl.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13924:124292", tmpl.tag
        assert tmpl.blocks == 1, tmpl.blocks
        assert tmpl.counters == 3, tmpl.counters
        assert tmpl.mcs == 2, tmpl.mcs
        assert tmpl.level() == -1, tmpl.level()
        break


def test_template_uniqness_level():
    infile = gt.InputFile(testfiles["test.map"])
    levels = [t.level() for t in infile]
    assert levels == [-1, 37, 0, -1, 0, 0, 0, 0, 0, 0], levels


def test_template_interleave():
    infile1 = gt.InputFile(testfiles["test.map"])
    infile2 = gt.InputFile(testfiles["test.map"])
    interleave = gt.interleave([infile1, infile2])
    levels = [t.level() for t in interleave]
    assert levels == [-1, -1, 37, 37, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], levels


def test_template_unmapped_filter_length():
    infile = gt.InputFile(testfiles["test.map.gz"])
    assert 1 == len([x for x in gt.unmapped(infile, 1)]), len([x for x in gt.unmapped(infile, 1)])
    assert 5 == len([x for x in gt.unmapped(infile, 0)]), len([x for x in gt.unmapped(infile, 0)])
    assert 1 == len([x for x in gt.unmapped(infile, 2)]), len([x for x in gt.unmapped(infile, 2)])
    assert 0 == len([x for x in gt.unmapped(infile, 3)]), len([x for x in gt.unmapped(infile, 3)])
    assert 0 == len([x for x in gt.unmapped(infile, 4)]), len([x for x in gt.unmapped(infile, 4)])


def test_template_unique_filter_level():
    infile = gt.InputFile(testfiles["test.map.gz"])
    assert 1 == len([x for x in gt.unique(infile, 20)]), "Should be length 1 buyt is: " + str([x.level() for x in gt.unique(infile, 20)])

# def test_fastq_writer_process():
#     infile = gt.InputFile(testfiles["test_merge_source_2.map"])
#     outfile = gt.OutputFile()


def test_template_map_parsing():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.tag == "A"


def test_template_length():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.length == 3


def test_template_read():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.read == "AAA"


def test_template_qualities():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.qualities == "###"


def test_template_no_qualities():
    template = gt.Template()
    template.parse("A/1\tAAA\t\t0\t-\n")
    assert template.qualities == ""


def test_template_to_map():
    template = gt.Template()
    template.parse("A/1\tAAA\t\t0\t-\n")
    assert template.to_map() == "A/1\tAAA\t\t0\t-", "Not '%s'" % template.to_map()


def test_template_to_fasta():
    template = gt.Template()
    template.parse("A/1\tAAA\t\t0\t-\n")
    assert template.to_fasta() == ">A/1\nAAA", "Not '%s'" % template.to_fasta()


def test_template_to_fastq():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.to_fastq() == "@A/1\nAAA\n+\n###", "Not '%s'" % template.to_fastq()


def test_template_to_sequence():
    template = gt.Template()
    template.parse("A/1\tAAA\t###\t0\t-\n")
    assert template.to_sequence() == "@A/1\nAAA\n+\n###", "Not '%s'" % template.to_sequence()
    template.parse("A/1\tAAA\t\t0\t-\n")
    assert template.to_fasta() == ">A/1\nAAA", "Not '%s'" % template.to_sequence()


def test_template_merge():
    template_1 = gt.Template()
    template_1.parse("A/1\tAAA\t###\t0\t-\n")
    template_2 = gt.Template()
    template_2.parse("A/1\tAAA\t###\t1\tchr1:+:50:3")
    template_1.merge(template_2)
    assert template_1.to_map() == "A/1\tAAA\t###\t1\tchr1:+:50:3", "Not '%s'" % template_1.to_map()

#
#def test_template_counters_list():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    for tmpl in infile:
#        counters = []
#        for c in tmpl.counters():
#            counters.append(c)
#        assert counters == [0, 4]
#        # test list direct list conversion
#        assert list(tmpl.counters()) == [0, 4]
#        assert len(tmpl.counters()) == 2
#
#
#def test_template_block_list():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    for tmpl in infile:
#        blocks = []
#        for c in tmpl.blocks():
#            blocks.append(c)
#        assert len(blocks) == 2
#        ##assert c[0] != c[1]
#        ## test that we have two different alignments
#        assert blocks[0].tag != blocks[1].tag
#        assert blocks[0].tag.endswith("/1")
#        assert blocks[1].tag.endswith("/2")
#        assert blocks[0].read == "GCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG"
#        assert blocks[1].read == "TTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT"
#        assert len(tmpl.blocks()) == 2
#
#
### test without iterating the mismatch blocks
#def test_template_mapping_iteration_with_paired_splitmap_number_of_blocks():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    assert infile != None
#    block_count = 0
#    for template in infile:
#        for block in template.mappings():
#            block_count += 1
#    assert block_count == 4
#
#
## test with iterating the mismatch blocks
#def test_template_mapping_iteration_with_paired_splitmap():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    assert infile != None
#    block_count = 0
#    for template in infile:
#        for block in template.mappings():
#            block_count += 1
#            for m, junction, distance in block:
#                assert m.seq_name == "chr7"
#    assert block_count == 4
#
#
#def test_template_mappings_iterator():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    c = 0
#    for tmpl in infile:
#        maps = []
#        for mappings in tmpl.mappings():
#            maps.append(mappings)
#            if c == 0:
#                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
#                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
#                assert [m[2] for m in mps] == [334, 65, -1]
#            elif c == 1:
#                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
#                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
#                assert [m[2] for m in mps] == [334, 399, -1]
#            elif c == 2:
#                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
#                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
#                assert [m[2] for m in mps] == [1554014, 399, -1]
#            c += 1
#    assert c == 4
#
#
## def test_template_mappings_iterator_outside_loop_throws_exception():
##     print("Error iteration")
##     infile = gt.open_file(testfiles.paired_w_splitmap)
##     c = 0
##     for tmpl in infile:
##         maps = []
##         for mappings in tmpl.mappings():
##             maps.append(mappings)
##     for m in maps:
##         for (mis, skip_type, skip_distance) in m:
##             print mis
##     assert c == 4
###
