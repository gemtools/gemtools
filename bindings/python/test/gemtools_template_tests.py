#!/usr/bin/env python

import gem.gemtools as gt
import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


def test_template_attribute_reading():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    for tmpl in infile:
        assert tmpl.tag == "HWI-ST661:131:C051TACXX:2:1101:1653:2244"
        assert tmpl.num_blocks == 2
        assert tmpl.num_counters == 2
        assert tmpl.max_complete_strata == 0


def test_template_counters_list():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    for tmpl in infile:
        counters = []
        for c in tmpl.counters():
            counters.append(c)
        assert counters == [0, 4]
        # test list direct list conversion
        assert list(tmpl.counters()) == [0, 4]
        assert len(tmpl.counters()) == 2


def test_template_block_list():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    for tmpl in infile:
        blocks = []
        for c in tmpl.blocks():
            blocks.append(c)
        assert len(blocks) == 2
        ##assert c[0] != c[1]
        ## test that we have two different alignments
        assert blocks[0].tag != blocks[1].tag
        assert blocks[0].tag.endswith("/1")
        assert blocks[1].tag.endswith("/2")
        assert blocks[0].read == "GCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG"
        assert blocks[1].read == "TTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT"
        assert len(tmpl.blocks()) == 2


## test without iterating the mismatch blocks
def test_template_mapping_iteration_with_paired_splitmap_number_of_blocks():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    assert infile != None
    block_count = 0
    for template in infile:
        for block in template.mappings():
            block_count += 1
    assert block_count == 4


## test with iterating the mismatch blocks
def test_template_mapping_iteration_with_paired_splitmap():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    assert infile != None
    block_count = 0
    for template in infile:
        for block in template.mappings():
            block_count += 1
            for m, junction, distance in block:
                assert m.seq_name == "chr7"
    assert block_count == 4


def test_template_mappings_iterator():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    c = 0
    for tmpl in infile:
        maps = []
        for mappings in tmpl.mappings():
            maps.append(mappings)
            if c == 0:
                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
                assert [m[2] for m in mps] == [334, 65, -1]
            elif c == 1:
                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
                assert [m[2] for m in mps] == [334, 399, -1]
            elif c == 2:
                mps = [(mis, skip_type, skip_distance) for mis, skip_type, skip_distance in mappings]
                assert [m[1] for m in mps] == [gt.SKIP_SPLICE, gt.SKIP_INSERT, gt.SKIP_NO_JUNCTION]
                assert [m[2] for m in mps] == [1554014, 399, -1]
            c += 1
    assert c == 4


