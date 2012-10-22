#!/usr/bin/env python
from gem import Read
import gem

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
