#!/usr/bin/env python

#import gem
#import gem.gemtools as gt
from testfiles import testfiles

test_mapping = testfiles["test.map"]
test_zipped_mapping = testfiles["test.map.gz"]
test_fastq = testfiles["test.fastq"]



#def test_get_template():
#    infile = gem.files.open(test_mapping)
#    for read in infile:
#        assert read._get_template() is not None
#        #print read.template.max_complete_strata
#
#def test_gt_merge_templates():
#    read_1 = gem.Read()
#    read_1.id = "ID"
#    read_1.line = "ID\tACGT\t####\t1\tchr1:-:20:4"
#    read_2 = gem.Read()
#    read_2.id = "ID"
#    read_2.line = "ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1"
#    print gt.merge_templates(read_1._get_template(), read_2._get_template())

#def test_map_2_fastq_conversion():
#    infile = gt.open_file(test_mapping)
#    assert infile != None
#    lines = []
#    counter = 0
#    for template in infile:
#        for alignment in template.blocks():
#            fq = alignment.to_sequence()
#            assert fq != None
#            counter += 1
#            for l in fq.split("\n"):
#                lines.append(l + "\n")
#    assert counter == 10
#    ## commented for now as the library now cuts tags and ids ?
##    of = open(test_fastq, 'r')
##    original_lines = of.readlines()
##    assert original_lines == lines
#
#
#def test_template_mapping_iteration():
#    infile = gt.open_file(test_mapping)
#    assert infile != None
#    map_count = 0
#    map_mismatches = [1, 1, 2, 3, 41, 0, 1, 1, 1, 2, 2, 0, 1, 1, 0, 0, 0]
#    map_chrs = ["chr11", "chrX", "chr1", "chrM",
#                "chr5", "chrM", "chr1", "chr6",
#                "chr16", "chr19", "chr1", "chr1",
#                "chr2", "chr9", "chr15", "chr12",
#                "chr1"]
#    map_positions = [77597507,  108297381, 237766340, 11809,
#                     134997659, 7109, 567659, 113902972, 61089473, 7731168,
#                     192685585, 160654905, 162359617, 38397301, 72492866,
#                     123774072, 91852848]
#    read_mismatches = []
#    read_chrs = []
#    read_positions = []
#    for template in infile:
#        for mappings in template.mappings():
#            for m, junction, distance in mappings:
#                assert junction == 0
#                # assert distance == -1
#                map_count += 1
#                read_mismatches.append(m.num_mismatches)
#                read_chrs.append(m.seq_name)
#                read_positions.append(m.position)
#
#    assert map_count == 17
#    assert map_mismatches == read_mismatches
#    assert map_chrs == read_chrs
#    assert map_positions == read_positions
#
#
#def test_template_alignment_mapping_iteration():
#    infile = gt.open_file(test_mapping)
#    assert infile != None
#    map_count = 0
#    map_mismatches = [1, 1, 2, 3, 41, 0, 1, 1, 1, 2, 2, 0, 1, 1, 0, 0, 0]
#    map_chrs = ["chr11", "chrX", "chr1", "chrM",
#                "chr5", "chrM", "chr1", "chr6",
#                "chr16", "chr19", "chr1", "chr1",
#                "chr2", "chr9", "chr15", "chr12",
#                "chr1"]
#    map_positions = [77597507,  108297381, 237766340, 11809,
#                     134997659, 7109, 567659, 113902972, 61089473, 7731168,
#                     192685585, 160654905, 162359617, 38397301, 72492866,
#                     123774072, 91852848]
#    read_mismatches = []
#    read_chrs = []
#    read_positions = []
#    for template in infile:
#        for block in template.blocks():
#            for mappings in block.mappings():
#                for m, junction, distance in mappings:
#                    assert junction == 0
#                    # assert distance == -1
#                    map_count += 1
#                    read_mismatches.append(m.num_mismatches)
#                    read_chrs.append(m.seq_name)
#                    read_positions.append(m.position)
#
#    assert map_count == 17
#    assert map_mismatches == read_mismatches
#    assert map_chrs == read_chrs
#    assert map_positions == read_positions
#
#
#def test_template_mapping_iteration_with_paired_splitmap():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    assert infile != None
#    block_count = 0
#    for template in infile:
#        for block in template.mappings():
#            block_count += 1
#            for m, junction, distance in block:
#                assert m.seq_name == "chr7"
#                #print "%s:%d %d %d" % (m.seq_name, m.position, junction, distance)
#                # assert distance == -1
#    assert block_count == 4
#
#
#def test_template_alignment_mapping_iteration_with_paired_splitmap():
#    infile = gt.open_file(testfiles["paired_w_splitmap.map"])
#    assert infile != None
#    block_count = 0
#    for template in infile:
#        for block in template.blocks():
#            for mapping in block.mappings():
#                block_count += 1
#                for m, junction, distance in mapping:
#                    assert m.seq_name == "chr7"
#                    #print "%s:%d %d %d" % (m.seq_name, m.position, junction, distance)
#                # assert distance == -1
#    assert block_count == 8
#
#def test_conversion_to_bed():
#    infile = gt.open_file(testfiles["bedconvert.map"])
#    for tempalte in infile:
#        read = "/1"
#        beds = []
#        for alignment in tempalte.blocks():
#            rname = alignment.tag
#            if not rname.endswith("/1"):
#                    rname += read
#                    read = "/2"
#
#            for mapping in alignment.mappings():
#                block_count = 0
#                block_sizes = []
#                block_starts = []
#                start = -1
#                end = -1
#                name = None
#                strand = "+"
#                for m, junction, distance in mapping:
#                    if start == -1:
#                        name = m.seq_name
#                        start = m.position - 1
#                        end = m.position + m.global_length - 1
#                        if m.direction != 0:
#                            strand = "-"
#                    block_count += 1
#                    block_starts.append(m.position - 1 - start)
#                    block_sizes.append(m.length)
#                if not (start + block_starts[-1] + block_sizes[-1]) == end:
#                    assert (start + block_starts[-1] + block_sizes[-1]) == end
#                beds.append("%s\t%d\t%d\t%s\t0\t%s\t.\t.\t0,0,0\t%d\t%s\t%s" % (name, start, end, rname, strand,
#                                                                            block_count,
#                                                                            ",".join([str(c) for c in block_sizes]),
#                                                                            ",".join([str(c) for c in block_starts])))
##        for line in set(beds):
##            print line
