#!/usr/bin/env python

import gem.gemtools as gt
import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


def test_alignment_attribute_reading_in_loop():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = []
    for tmpl in infile:
        for block in tmpl.blocks():
            alis.append({
                "tag":block.tag,
                'read':block.read,
                'qualities':block.qualities,
                'mcs':block.max_complete_strata
            })
    assert alis[0]["tag"] == "HWI-ST661:131:C051TACXX:2:1101:1653:2244/1"  
    assert alis[0]["read"] == "GCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG"  
    assert alis[0]["qualities"] == "CCCFFFFFHHHHHJJIEHJIIIJJGIH>4=ADDDDDDBDD@DDBDBDDDBDBBDDDDDDDBBDDCC>CBD@BDDD"  
    assert alis[0]["mcs"] == 0  


    assert alis[1]["tag"] == "HWI-ST661:131:C051TACXX:2:1101:1653:2244/2"
    assert alis[1]["read"] == "TTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT"
    assert alis[1]["qualities"] == "CCCFFFFFHHHHHJJGIIJJJIGIJJJFHIJIICHIIGIHIIJIJHFFDD?BDDDDDDDD@BDDDDDDBDDDDDD"
    assert alis[1]["mcs"] == 0


def test_alignment_attribute_reading_in_list():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = [ali for tmpl in infile for ali in tmpl.blocks()]
    assert alis[0].tag == "HWI-ST661:131:C051TACXX:2:1101:1653:2244/1"  
    assert alis[0].read == "GCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG"  
    assert alis[0].qualities == "CCCFFFFFHHHHHJJIEHJIIIJJGIH>4=ADDDDDDBDD@DDBDBDDDBDBBDDDDDDDBBDDCC>CBD@BDDD"  
    assert alis[0].max_complete_strata == 0  


    assert alis[1].tag == "HWI-ST661:131:C051TACXX:2:1101:1653:2244/2"
    assert alis[1].read == "TTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT"
    assert alis[1].qualities == "CCCFFFFFHHHHHJJGIIJJJIGIJJJFHIJIICHIIGIHIIJIJHFFDD?BDDDDDDDD@BDDDDDDBDDDDDD"
    assert alis[1].max_complete_strata == 0


def test_alignment_to_sequencet():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = [ali for tmpl in infile for ali in tmpl.blocks()]
    assert alis[0].to_sequence() == "@HWI-ST661:131:C051TACXX:2:1101:1653:2244/1\nGCGCGGCCGGGACCGCAGAGCCCCGGGAGCCCGCTCGAGGAGGAGCGGCAGACGCAGCGCTCTAAACCGCAGCCG\n+\nCCCFFFFFHHHHHJJIEHJIIIJJGIH>4=ADDDDDDBDD@DDBDBDDDBDBBDDDDDDDBBDDCC>CBD@BDDD"
    assert alis[1].to_sequence() == "@HWI-ST661:131:C051TACXX:2:1101:1653:2244/2\nTTCCGCTTGGTGCTCTCGCTGCAGCGGTTCAGGATGAGGTCGGCGCTCGGCCGCGGGGGCACCGCCGGCTGCGGT\n+\nCCCFFFFFHHHHHJJGIIJJJIGIJJJFHIJIICHIIGIHIIJIJHFFDD?BDDDDDDDD@BDDDDDDBDDDDDD"

def test_alignment_counters():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = [ali for tmpl in infile for ali in tmpl.blocks()]
    print [list(c.counters()) for c in alis ]
    assert len(alis[0].counters()) == 2    
    assert len(alis[1].counters()) == 1
    assert list(alis[0].counters()) == [0, 4]
    assert list(alis[1].counters()) == [4]


def test_alignment_mappings_counts():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = [ali for tmpl in infile for ali in tmpl.blocks()]
    mappings = [ali.mappings()  for ali in alis]
    assert len(mappings) == 2
    maps = [list(m) for m in mappings]
    assert len(maps[0]) == 4
    assert len(maps[1]) == 4

def test_alignment_mappings_content():
    infile = gt.open_file(testfiles.paired_w_splitmap)
    alis = [ali for tmpl in infile for ali in tmpl.blocks()]
    mappings = [ali.mappings()  for ali in alis]
    assert len(mappings) == 2
    maps = [list(m) for m in mappings]

    #chr7:R74572684<*334>38::chr7:F74572619,
    #chr7:F74203012<*334>38::chr7:R74203411,
    #chr7:F72649332<*1554014>38::chr7:R72649731,
    #chr7:F72649332<*334>38::chr7:R72649731
    m1 = list(maps[0][0])
    assert m1[0][0].seq_name == "chr7"
    assert m1[0][0].position == 74572684
    assert m1[0][0].length == 37
    assert m1[0][0].direction == gt.REVERSE
    assert m1[0][2] == 334
    assert m1[0][1] == gt.SKIP_SPLICE
    assert m1[1][0].seq_name == "chr7"
    assert m1[1][0].position == 74572684 + m1[0][0].length + 334
    assert m1[1][0].direction == gt.REVERSE
    assert m1[1][2] == -1
    assert m1[1][1] == gt.SKIP_NO_JUNCTION

    m1 = list(maps[0][1])
    assert m1[0][0].seq_name == "chr7"
    assert m1[0][0].position == 74203012
    assert m1[0][0].direction == gt.FORWARD
    assert m1[0][2] == 334
    assert m1[0][1] == gt.SKIP_SPLICE
    assert m1[1][0].position == m1[0][0].position + m1[0][0].length + m1[0][2]
    
