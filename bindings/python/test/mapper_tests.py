#!/usr/bin/env python

import testfiles
import unittest
import os

import gem


#class MapperTests(unittest.TestCase):
#    def test_simple_mapper_run(self):
#        if os.path.exists("testresults/output.map"):
#            os.remove("testresults/output.map")
#        output = gem.mapper(testfiles.test_fastq, "testresults/output.map", testfiles.index)
#        assert os.path.exists("testresults/output.map")
#        assert output is not None
#        counter = 0
#        for read in output:
#           counter += 1
#        assert counter == 10
#
#    def test_simple_mapper_run_with_trimming(self):
#        if os.path.exists("testresults/output.map"):
#            os.remove("testresults/output.map")
#        output = gem.mapper(testfiles.test_fastq, "testresults/output.map", testfiles.index)

#        for read in gem.unmapped(output):
#            print str(read)
        #output = gem.mapper(gem.trim(gem.unmapped(output, 3), 0, 10), "testresults/output.map", testfiles.index)
        #assert os.path.exists("testresults/output.map")
        #assert output is not None
        #counter = 0
        #for read in output:
        #    counter += 1
        #    assert len(read.sequence) == 65
        #assert counter == 2
