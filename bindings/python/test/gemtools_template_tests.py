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
            print "Original tag", c.tag
        assert len(blocks) == 2
        ##assert c[0] != c[1]
        ## test that we have two different alignments
        assert blocks[0].tag != blocks[1].tag
        assert blocks[0].tag.endswith("/1")
        assert blocks[1].tag.endswith("/2")
        assert blocks[0].read == "GCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG"
        assert blocks[1].read == "TTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT"
        assert len(tmpl.blocks()) == 2

