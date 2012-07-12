#!/usr/bin/env python
import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


#class ConversionTestCase(unittest.TestCase):
#    def test_map_2_fastq_conversion_lib(self):
#        of = gem.gem_open(test_mapping)
#        counter = 0
#       lines = []
#        for read in of:
#            fastq_lines = read.to_sequence()
#            for l in fastq_lines.split("\n"):
#                lines.append(l + "\n")
#            counter += 1
#        assert counter == 10
#
#        of = open(test_fastq, 'r')
#        original_lines = of.readlines()
#        assert original_lines == lines
#
#    def test_zipped_map_2_fastq_conversion_lib(self):
#        of = gem.gem_open(test_zipped_mapping)
#        counter = 0
#        lines = []
#        for read in of:
#            fastq_lines = read.to_sequence()
#            for l in fastq_lines.split("\n"):
#                lines.append(l + "\n")
#            counter += 1
 #       assert counter == 10
#
#        of = open(test_fastq, 'r')
#        original_lines = of.readlines()
#       assert original_lines == lines
