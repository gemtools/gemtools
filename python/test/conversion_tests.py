#!/usr/bin/env python
import os
import shutil
from nose.tools.nontrivial import with_setup
from gem import Read
import gem
from testfiles import testfiles

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


def test_sam_fields_number():
    p = gem.files.parse_map()
    read = p.line2read("ID\tACGT\t####\t1:1\tchr1:-:20:4,chr9:+:50:2C1")
    sam = gem.gem2sam([read], compact=True, quality=33)
    assert sam is not None
    for r in sam:
        f = r.line.split("\t")
        assert len(f) == 17


@with_setup(setup_func, cleanup)
def test_gem2sam_execution():
    input = gem.files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    sam = gem.gem2sam(mappings, index, compact=True)
    assert sam is not None
    assert sam.process is not None
    assert sam.filename is None
    count = 0
    for read in sam:
        #print read.line.strip()
        count += 1
    assert count == 10000, "Count 10000!=%d" % count

@with_setup(setup_func, cleanup)
def test_gem2sam_execution_to_file():
    input = gem.files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    result = results_dir+"/test_sam.sam"
    sam = gem.gem2sam(mappings, index, output=result, compact=True)
    assert sam is not None
    assert sam.process is not None
    assert sam.filename == result
    assert os.path.exists(result)
    count = 0
    for read in sam:
        #print read.line.strip()
        count += 1
    assert count == 10000, "Count 10000!=%d" % count

@with_setup(setup_func, cleanup)
def test_gem2sam_sam2bam():
    input = gem.files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    sam = gem.gem2sam(mappings, index, compact=True)
    result = results_dir+"/test_sam.bam"
    bam = gem.sam2bam(sam, output=result)
    assert os.path.exists(result)
    count = 0
    for l in gem.files.open(result):
        count += 1
    assert count == 10000, "Count 10000!=%d" % count
