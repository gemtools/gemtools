import gem
from gem import files
from testfiles import testfiles
import logging
import shutil
import os
from nose.tools import with_setup
import subprocess

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
def test_file_merge_async():
    reads_1 = files.open(testfiles["test.map"])
    reads_2 = files.open(testfiles["test.map"])
    merged = gem.merge(reads_1, [reads_2])
    num_reads = sum(1 for r in merged)
    assert num_reads == 10


@with_setup(setup_func, cleanup)
def test_file_merge_pairwise_same():
    reads_1 = files.open(testfiles["test.map"])
    reads_2 = files.open(testfiles["test.map"])
    merged = gem.merge(reads_1, [reads_2],
                       output=results_dir + "/merge_result.map",
                       threads=8, paired=True)
    num_reads = sum(1 for r in merged)
    assert num_reads == 10


@with_setup(setup_func, cleanup)
def test_file_merge_pairwise_same_content_big():
    reads_1 = files.open(testfiles["20t.map.gz"])
    reads_2 = files.open(testfiles["20t.map.gz"])
    merged = gem.merge(reads_1, [reads_2],
                       output=results_dir + "/merge_result.map",
                       threads=8, same_content=True)
    num_reads = sum(1 for r in merged)
    assert num_reads == 20000


@with_setup(setup_func, cleanup)
def test_file_merge_pairwise_same_content_big_uncompressed():
    subprocess.call("cp %s %s; gunzip %s;" % (testfiles["20t.map.gz"], results_dir, results_dir + "/20t.map.gz"), shell=True)
    reads_1 = files.open(results_dir + "/20t.map")
    reads_2 = files.open(results_dir + "/20t.map")
    merged = gem.merge(reads_1, [reads_2],
                       output=results_dir + "/merge_result.map",
                       threads=8, same_content=True)
    num_reads = sum(1 for r in merged)
    assert num_reads == 20000


@with_setup(setup_func, cleanup)
def test_file_merge_pairwise_same_content_big_sub():
    reads_1 = files.open(testfiles["20t.map.gz"])
    reads_2 = files.open(testfiles["20t_sub.map.gz"])
    merged = gem.merge(reads_1, [reads_2],
                       output=results_dir + "/merge_result.map",
                       threads=8, paired=True)
    num_reads = sum(1 for r in merged)
    assert num_reads == 20000
