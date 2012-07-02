#!/usr/bin/python

import sys
import gem
import gemtools as gt


def convert_direct():
    infile = gt.gt_input_stream_open(sys.stdin)
    template = gt.gt_template_new()
    map_input = gt.gt_buffered_map_input_new(infile)
    error_code = 1
    while True:
        error_code = gt.gt_buffered_map_input_get_template(map_input, template)
        if error_code == gt.GT_BMI_FAIL:
            continue
        if error_code == 0:
            break
        num_blocks = gt.gt_template_get_num_blocks(template)
        for i in xrange(0, num_blocks):
            alignment = gt.gt_template_get_block(template, 0)
            print "@%s\n%s\n+\n%s" % (gt.gt_template_get_tag(template),
                                     gt.gt_alignment_get_read(alignment),
                                     gt.gt_alignment_get_qualities(alignment))


def convert_lib(file):
    for read in gem.gem_open(file):
        print read.to_sequence()


if __name__ == "__main__":

    #convert_direct()

    if len(sys.argv) != 2:
        print >> sys.stderr, "Usage: gem2fastq.py <gemfile>"
        exit(1)
    convert_lib(sys.argv[1])
