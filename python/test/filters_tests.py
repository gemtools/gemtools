import os
import shutil
from nose.tools import with_setup
import sys

import gem.gemtools as gt
from gem import files
from gem import filter
from testfiles import testfiles

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'

results_dir = None


def setup_func():
    global results_dir
    results_dir = "test_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    results_dir = os.path.abspath(results_dir)


def cleanup():
    shutil.rmtree(results_dir, ignore_errors=True)


def test_cat():
    reads_1 = files.open(testfiles["reads_1.fastq"])
    reads_2 = files.open(testfiles["reads_2.fastq"])
    num_reads = sum(1 for r in filter.cat([reads_1, reads_2]))
    assert num_reads == 20000


@with_setup(setup_func, cleanup)
def test_cat_to_file():
    target = results_dir + "/catted.fastq"
    reads_1 = files.open(testfiles["reads_1.fastq"])
    reads_2 = files.open(testfiles["reads_2.fastq"])
    out = gt.OutputFile(target)
    gt.cat([reads_1, reads_2]).write_stream(out)
    c = 0
    ids = []
    with open(target) as f:
        for l in f:
            if c % 4 == 0:
                ids.append(l.strip())
            c += 1
    assert len(set(ids)) == 20000
