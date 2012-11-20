#!/usr/bin/env python

import gem.gemtools as gt
from testfiles import testfiles

test_mapping = testfiles["test.map"]
test_zipped_mapping = testfiles["test.map.gz"]
test_fastq = testfiles["test.fastq"]

#def test_maps_paired_alignment_iteratoion():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    maps = []
#    for tmpl in infile:
#        for alignment_block in tmpl.blocks():
#            for mappings_iterator in alignment_block.mappings():
#                for m, t, d in mappings_iterator:
#                    maps.append(m)
#    assert maps[0].seq_name == "chr7"
#    assert maps[0].length == 37
#    assert maps[0].global_length == 409    
#    assert maps[0].base_length == 37
#    assert maps[0].num_mismatches == 0
#    assert len(maps[0].mismatches()) == 0
#
#def test_map_mismatches():
#    infile = gt.open_file(testfiles["paired_sm_mm.map"])
#    maps = []
#    for tmpl in infile:
#        for alignment_block in tmpl.blocks():
#            for mappings_iterator in alignment_block.mappings():
#                for m, t, d in mappings_iterator:
#                    maps.append(m)
#    assert maps[0].seq_name == "chr6"
#    assert maps[0].length == 75
#    assert maps[0].global_length == 75    
#    assert maps[0].base_length == 75
#    assert maps[0].num_mismatches == 1
#    assert len(maps[0].mismatches()) == 1
#    mm = list(maps[0].mismatches())
#    assert mm[0].type == gt.MISMS
#    assert mm[0].position == 73
#    assert mm[0].base == "T"
#    assert mm[0].size == -1
    
    
