#!/usr/bin/env python
import os
import shutil
from nose.tools.nontrivial import with_setup
from gem import Read
import gem
from testfiles import testfiles

index = testfiles["genome.gem"]

def setup_func():
    global results_dir
    results_dir = "test_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    results_dir = os.path.abspath(results_dir)


def cleanup():
    shutil.rmtree(results_dir, ignore_errors=True)

def test_read_to_sequence_wo_qualities():
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT"
    assert read.to_sequence() == ">myid\nACGT"


def test_read_to_sequence_w_qualities():
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT"
    read.qualities = "####"
    assert read.to_sequence() == "@myid\nACGT\n+\n####"

def test_read_to_fastq_wo_qualities():
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT"
    assert read.to_fastq() == "@myid\nACGT\n+\n[[[["

def test_read_to_fastq_w_qualities_trimmed():
    gem._trim_qualities=True
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT"
    read.qualities = "[[[[["
    assert read.to_fastq() == "@myid\nACGT\n+\n[[[[", read.to_fastq()
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT"
    read.qualities = "[["
    assert read.to_fastq() == "@myid\nACGT\n+\n[[[[", read.to_fastq()
    gem._trim_qualities=False

def test_read_paired_to_sequence_wo_qualities():
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT CGTA"
    assert read.to_sequence() == ">myid/1\nACGT\n>myid/2\nCGTA"

def test_read_paired_to_sequence_w_qualities():
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT CGTA"
    read.qualities = "[[[[ ####"
    assert read.to_sequence() == "@myid/1\nACGT\n+\n[[[[\n@myid/2\nCGTA\n+\n####"

def test_read_paired_to_sequence_w_qualities_trimmed():
    gem._trim_qualities=True
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT CGTA"
    read.qualities = "[[[ ###"
    assert read.to_sequence() == "@myid/1\nACGT\n+\n[[[[\n@myid/2\nCGTA\n+\n###["
    read = Read()
    read.id = "myid"
    read.sequence = "ACGT CGTA"
    read.qualities = "[[[[[ #####"
    assert read.to_sequence() == "@myid/1\nACGT\n+\n[[[[\n@myid/2\nCGTA\n+\n####"
    gem._trim_qualities=False


def test_without_NH_sam_field():
    p = gem.files.parse_map()
    read = p.line2read("ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1")
    sam = gem.gem2sam([read], compact=True, quality=33, append_nh=False)
    assert sam is not None
    for r in sam:
        f = r.line.split("\t")
        assert len(f) == 16


def test_NH_sam_field():
    p = gem.files.parse_map()
    read = p.line2read("ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1")
    sam = gem.gem2sam([read], compact=True, quality=33, append_nh=True)
    assert sam is not None
    for r in sam:
        f = r.line.split("\t")
        assert len(f) == 17


@with_setup(setup_func, cleanup)
def test_gem2sam_execution():
    input = gem.files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    sam = gem.gem2sam(mappings, index, compact=True, append_nh=True)
    assert sam is not None
    assert sam.process is not None
    assert sam.filename is None
    count = 0
    for read in sam:
        count += 1
    assert count == 10000
