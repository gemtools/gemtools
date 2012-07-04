#!/usr/bin/env python

import os

dir = os.path.dirname(os.path.abspath(__file__))
print dir
testfiles_dir = os.path.realpath(dir + "/../../../testdata")
print testfiles_dir
test_fastq = testfiles_dir + "/test.fastq"
assert os.path.exists(test_fastq)
test_map = testfiles_dir + "/test.map"
assert os.path.exists(test_map)
paired_w_splitmap = testfiles_dir + "/paired_w_splitmap.map"
assert os.path.exists(test_map)
test_map_gz = testfiles_dir + "/test.map.gz"
assert os.path.exists(test_map_gz)
index = testfiles_dir + "/genome.gem"
assert os.path.exists(index)
