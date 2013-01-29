import gem
from gem import files
from testfiles import testfiles
import logging
import shutil
import os
from nose.tools import with_setup

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'

results_dir = None


def setup_func():
    logging.basicConfig(format='%(asctime)-15s %(levelname)s: %(message)s', level=logging.INFO)
    global results_dir
    results_dir = "test_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    results_dir = os.path.abspath(results_dir)


def cleanup():
    shutil.rmtree(results_dir, ignore_errors=True)


@with_setup(setup_func, cleanup)
def test_file_merge():
    reads_1 = files.open(testfiles["test.map"])
    reads_2 = files.open(testfiles["test.map"])
    merged = gem.merger(reads_1, [reads_2]).merge(results_dir + "/merge_result.map")
    num_reads = sum(1 for r in merged)
    assert num_reads == 10
