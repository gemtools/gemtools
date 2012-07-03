#!/usr/bin/env python
import unittest
import gemtools as gt
import gem

import testfiles

test_mapping = testfiles.test_map
test_zipped_mapping = testfiles.test_map_gz
test_fastq = testfiles.test_fastq


class ConversionTestCase(unittest.TestCase):
    def test_map_2_fastq_conversion(self):
        infile = gt.gt_input_file_open(test_mapping, False)
        template = gt.gt_template_new()
        map_input = gt.gt_buffered_map_input_new(infile)

        error_code = 1
        counter = 0
        lines = []
#        while True:
#
#            error_code = gt.gt_buffered_map_input_get_template(map_input, template)
#            if error_code == gt.GT_BMI_FAIL:
#                continue
#            if error_code == 0:
#                break
#            counter += 1
#            num_blocks = gt.gt_template_get_num_blocks(template)
#            for i in range(0, num_blocks):
#                alignment = gt.gt_template_get_block(template, 0)
#                tag = gt.gt_template_get_tag(template)
#                read = gt.gt_alignment_get_read(alignment)
#                quals = gt.gt_alignment_get_qualities(alignment)
#                lines.append("@" + tag + "\n")
#                lines.append(read + "\n")
#                lines.append("+\n")
#                lines.append(quals + "\n")
#        assert counter == 10
#        of = open(test_fastq, 'r')
#        original_lines = of.readlines()
#        assert original_lines == lines


if __name__ == "__main__":
    import sys
    import gempy
    import gem
    import gemtools

#    infile = gt.gt_input_stream_open(sys.stdin)
#    gt_template = gt.gt_template_new()
#    map_input = gt.gt_buffered_map_input_new(infile)
#    while True:
#        status = gt.gt_buffered_map_input_get_template(map_input, gt_template)
#        if status == 0:
#            break
#        t = gem.Template(gt_template)
#        for a in t.alignments:
#            print a.tag

    for tmpl in gempy.myiter(sys.stdin):
         print "Tag", tmpl.tag
         print "MCS", tmpl.max_complete_strata
         print ":".join([str(x) for x in tmpl.counters])
         for block in tmpl.blocks:
             print "Block tag", block.tag

