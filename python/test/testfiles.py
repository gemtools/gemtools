#!/usr/bin/env python

import os

## setup testfiles for easier access
__base_dir = os.path.dirname(os.path.abspath(__file__))
__testfiles_dir = os.path.realpath(__base_dir + "/../testdata")

testfiles = {}

for file in os.listdir(__testfiles_dir):
    testfiles[file] = "%s/%s" % (__testfiles_dir, file)

