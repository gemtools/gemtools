#!/usr/bin/env python
import sys
import unittest
import gemtools as gt
import gem

test_mapping = "../../testdata/test.map"
test_zipped_mapping = "../../testdata/test_zipped.map.gz"
test_fastq = "../../testdata/test.fastq"


class ConversionTestCase(unittest.TestCase):
    def test_map_2_fastq_conversion(self):
        infile = gt.gt_input_file_open(test_mapping, False)
        template = gt.gt_template_new()
        map_input = gt.gt_buffered_map_input_new(infile)

        error_code = 1
        counter = 0
        lines = []
        while True:

            error_code = gt.gt_buffered_map_input_get_template(map_input, template)
            if error_code == gt.GT_BMI_FAIL:
                continue
            if error_code == 0:
                break
            counter += 1
            num_blocks = gt.gt_template_get_num_blocks(template)
            for i in range(0, num_blocks):
                alignment = gt.gt_template_get_block(template, 0)
                tag = gt.gt_template_get_tag(template)
                read = gt.gt_alignment_get_read(alignment)
                quals = gt.gt_alignment_get_qualities(alignment)
                lines.append("@" + tag + "\n")
                lines.append(read + "\n")
                lines.append("+\n")
                lines.append(quals + "\n")
        assert counter == 10
        of = open(test_fastq, 'r')
        original_lines = of.readlines()
        assert original_lines == lines

    def test_map_2_fastq_conversion_lib(self):
        of = gem.gem_open(test_mapping)
        counter = 0
        lines = []
        for read in of:
            fastq_lines = read.to_sequence()
            for l in fastq_lines.split("\n"):
                lines.append(l + "\n")
            counter += 1
        assert counter == 10

        of = open(test_fastq, 'r')
        original_lines = of.readlines()
        assert original_lines == lines

    def test_zipped_map_2_fastq_conversion_lib(self):
        of = gem.gem_open(test_zipped_mapping)
        counter = 0
        lines = []
        for read in of:
            fastq_lines = read.to_sequence()
            for l in fastq_lines.split("\n"):
                lines.append(l + "\n")
            counter += 1
        assert counter == 10

        of = open(test_fastq, 'r')
        original_lines = of.readlines()
        assert original_lines == lines
