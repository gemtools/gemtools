#!/usr/bin/env python
import os
import shutil
from nose.tools import with_setup
import gem
from gem import files
from gem import filter
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
def test_sync_mapper_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index, results_dir + "/result.mapping")
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is not None
    assert mappings.filename == results_dir + "/result.mapping"
    assert sum(1 for x in mappings) == 10000


@with_setup(setup_func, cleanup)
def test_async_mapper_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is None
    assert sum(1 for x in mappings) == 10000


@with_setup(setup_func, cleanup)
def test_async_mapper_pipes_with_filter():
    def on_half_filter(reads):
        count = 0
        for read in reads:
            count += 1
            if count % 2 == 0: yield read


    input = files.open(testfiles["reads_1.fastq"])
    initial = gem.mapper(input, index)
    mappings = gem.mapper(on_half_filter(initial), index, results_dir + "/piped_out.mapping")

    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is not None
    assert mappings.filename == results_dir + "/piped_out.mapping"
    assert os.path.exists(results_dir + "/piped_out.mapping")
    assert sum(1 for x in mappings) == 5000

@with_setup(setup_func, cleanup)
def test_interleaved_mapper_run():
    input1 = files.open(testfiles["reads_1.fastq"])
    input2 = files.open(testfiles["reads_2.fastq"])
    mappings = gem.mapper(gem.filter.interleave([input1, input2]), index)
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is None
    assert sum(1 for x in mappings) == 20000

@with_setup(setup_func, cleanup)
def test_interleaved_pair_aligner_run():
    input1 = files.open(testfiles["reads_1.fastq"])
    input2 = files.open(testfiles["reads_2.fastq"])
    mappings = gem.mapper(gem.filter.interleave([input1, input2]), index)
    paired = gem.pairalign(mappings, index)
    assert paired is not None
    assert sum(1 for x in mappings) == 20000  ## test dataset does not pair at all

