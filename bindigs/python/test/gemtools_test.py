#!/usr/bin/env python

import gem.gemtools as gt
import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


def test_map_2_fastq_conversion():
    infile = gt.open_file(test_mapping)
    assert infile != None
    lines = []
    counter = 0
    for template in infile:
        for alignment in template.blocks:
            fq = alignment.to_sequence()
            assert fq != None
            counter += 1
            for l in fq.split("\n"):
                lines.append(l+"\n");
    assert counter == 10
    of = open(test_fastq, 'r')
    original_lines = of.readlines()
    assert original_lines == lines

