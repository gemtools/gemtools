#!/usr/bin/env python

import subprocess
import gem.gemtools as gt
from testfiles import testfiles
import os
import shutil
from nose.tools import with_setup
from gem import files

results_dir = None

def setup_func():
    global results_dir
    results_dir = "test_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    results_dir = os.path.abspath(results_dir)


def cleanup():
    shutil.rmtree(results_dir, ignore_errors=True)


def test_open_file():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    count = 0
    for tmpl in infile:
        count += 1
        assert tmpl.tag == "WI-ST472_0084:6:1101:1179:2208#CTTGTA", tmpl.tag
        assert tmpl.blocks == 2, tmpl.blocks
    assert count == 1


def test_open_stream():
    p = subprocess.Popen(["cat", testfiles["bedconvert.map"]], stdout=subprocess.PIPE)
    infile = gt.InputFile(p.stdout)
    p.wait()
    count = 0
    for tmpl in infile:
        count += 1
        assert tmpl.tag == "WI-ST472_0084:6:1101:1179:2208#CTTGTA", tmpl.tag
        assert tmpl.blocks == 2, tmpl.blocks
    assert count == 1


def test_open_raw_stream_from_stream():
    p = subprocess.Popen(["cat", testfiles["bedconvert.map"]], stdout=subprocess.PIPE)
    infile = gt.InputFile(p.stdout)
    p.wait()
    count = 0
    for line in infile.raw_stream():
        count += 1
        assert line.startswith("WI-ST472_0084:6:1101:1179:2208#CTTGTA")
    assert count == 1


def test_open_raw_stream_from_file():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    count = 0
    for line in infile.raw_stream():
        count += 1
        assert line.startswith("WI-ST472_0084:6:1101:1179:2208#CTTGTA")
    assert count == 1


def test_iterating_input_file():
    infile = gt.InputFile(testfiles["reads_1.fastq"])
    count = 0
    for line in infile:
        count += 1
    assert count == 10000, count


def test_iterating_interleaved():
    infile_1 = gt.InputFile(testfiles["reads_1.fastq"])
    infile_2 = gt.InputFile(testfiles["reads_2.fastq"])
    count = 0
    for line in gt.interleave([infile_1, infile_2]):
        count += 1
    assert count == 20000, count


def test_one_level_filter_chain():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    filtered = gt.filter(infile, gt.filter_unique(2))
    count = 0
    for t in filtered:
        count += 1
    assert count == 1


def test_unique_filter():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    filtered = gt.unique(infile, 2)
    count = 0
    for t in filtered:
        count += 1
    assert count == 1


def test_chain_filter_indirect():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    filtered = gt.trim(gt.unique(infile, 2), left=10, right=10)
    count = 0
    for t in filtered:
        count += 1
        assert len(t.read) == 81
    assert count == 1


@with_setup(setup_func, cleanup)
def test_writing_input_file():
    source = files.open(testfiles["test_merge_target.map"])
    target = results_dir + "/write_input.fastq"
    out = gt.OutputFile(target)
    source.write_stream(out, write_map=False)
    with open(target) as f:
        lines = f.readlines()
        assert len(lines) == 40


@with_setup(setup_func, cleanup)
def test_writing_interleaved_file():
    source1 = files.open(testfiles["reads_1.fastq"])
    source2 = files.open(testfiles["reads_2.fastq"])
    target = results_dir + "/write_interleaved.fastq"
    out = gt.OutputFile(target)
    gt.interleave([source1, source2]).write_stream(out, write_map=False)
    with open(target) as f:
        lines = f.readlines()
        assert len(lines) == 80000

