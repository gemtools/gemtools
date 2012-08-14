#!/usr/bin/env python
import os
import shutil
from nose.tools import with_setup
import gem
from gem import files
from gem import filter
from gem import junctions
from testfiles import testfiles

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'

index = testfiles["genome.gem"]
results_dir = None

def setup_func():
    global results_dir
    results_dir = "test_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    results_dir = os.path.abspath(results_dir)


def cleanup():
    shutil.rmtree(results_dir, ignore_errors=True)


@with_setup(setup_func, cleanup)
def test_realigning_BWA_input():
    input = files.open(testfiles["subsetBWA.sam"])
    realigned = gem.realign(filter.sam_2_map(input), index)
    sam = gem.gem2sam(realigned, index, single_end=True)
    count = 0
    for read in sam:
        count += 1
    assert count == 100
