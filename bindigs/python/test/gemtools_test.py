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
                lines.append(l + "\n")
    assert counter == 10
    of = open(test_fastq, 'r')
    original_lines = of.readlines()
    assert original_lines == lines


def test_template_mapping_iteration():    
    infile = gt.open_file(test_mapping)
    assert infile != None
    map_count = 0
    map_mismatches = [1, 1, 2, 3, 41, 0, 1, 1, 1, 2, 2, 0, 1, 1, 0, 0, 0]
    map_chrs = ["chr11", "chrX", "chr1", "chrM", 
                "chr5", "chrM", "chr1", "chr6", 
                "chr16", "chr19", "chr1", "chr1",
                "chr2", "chr9", "chr15", "chr12", 
                "chr1"]
    map_positions = [77597507,  108297381, 237766340, 11809,
                     134997659, 7109, 567659, 113902972, 61089473, 7731168,
                     192685585, 160654905, 162359617, 38397301, 72492866,
                     123774072, 91852848]
    read_mismatches = []
    read_chrs = []
    read_positions = []
    for template in infile:            
        for mappings in template.mappings:
            for m, junction, distance in mappings:
                assert junction == 0
                # assert distance == -1
                map_count += 1
                read_mismatches.append(m.num_mismatches)
                read_chrs.append(m.seq_name)
                read_positions.append(m.position)

    assert map_count == 17
    assert map_mismatches == read_mismatches
    assert map_chrs == read_chrs
    assert map_positions == read_positions


def test_template_alignment_mapping_iteration():    
    infile = gt.open_file(test_mapping)
    assert infile != None
    map_count = 0
    map_mismatches = [1, 1, 2, 3, 41, 0, 1, 1, 1, 2, 2, 0, 1, 1, 0, 0, 0]
    map_chrs = ["chr11", "chrX", "chr1", "chrM", 
                "chr5", "chrM", "chr1", "chr6", 
                "chr16", "chr19", "chr1", "chr1",
                "chr2", "chr9", "chr15", "chr12", 
                "chr1"]
    map_positions = [77597507,  108297381, 237766340, 11809,
                     134997659, 7109, 567659, 113902972, 61089473, 7731168,
                     192685585, 160654905, 162359617, 38397301, 72492866,
                     123774072, 91852848]
    read_mismatches = []
    read_chrs = []
    read_positions = []
    for template in infile: 
        for block in template.blocks:                       
            for mappings in block.mappings: 
                for m, junction, distance in mappings:                    
                    assert junction == 0
                    # assert distance == -1
                    map_count += 1
                    read_mismatches.append(m.num_mismatches)
                    read_chrs.append(m.seq_name)
                    read_positions.append(m.position)

    assert map_count == 17
    assert map_mismatches == read_mismatches
    assert map_chrs == read_chrs
    assert map_positions == read_positions


def test_template_mapping_iteration_with_paired_splitmap():    
    infile = gt.open_file(testfiles.paired_w_splitmap)
    assert infile != None
    block_count = 0
    for template in infile:         
        for block in template.mappings:            
            block_count += 1            
            for m, junction, distance in block:                
                assert m.seq_name == "chr7"                
                print "%s:%d %d %d" % (m.seq_name, m.position, junction, distance)
                # assert distance == -1
    print block_count
    assert block_count == 4
