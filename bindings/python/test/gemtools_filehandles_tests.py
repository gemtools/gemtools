#!/usr/bin/env python

import gem.gemtools as gt
import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


def test_open_map_file():
    infile = gt.open_file(testfiles.bedconvert_map)
    assert infile != null
