from gem import files
from gem import filter
from testfiles import testfiles

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'



def test_fastq_filter_interleave():
    reads_1 = files.open(testfiles["reads_1.fastq"])
    reads_2 = files.open(testfiles["reads_2.fastq"])
    num_reads = sum(1 for r in filter.interleave([reads_1, reads_2]))
    assert num_reads == 20000


def test_fastq_filter_unmapped():
    reads = files.open(testfiles["reads_1.fastq"])
    num_reads = sum(1 for r in reads.clone())
    sum_length = sum(1 for r in filter.unmapped(reads.clone()))
    assert num_reads == sum_length
    assert sum_length == 10000


def test_fastq_trim_right():
    reads = files.open(testfiles["reads_1.fastq"])
    sum_length = sum(r.length for r in filter.trim(reads, 0, 20))
    assert sum_length == 550000, sum_length


def test_fastq_trim_left():
    reads = files.open(testfiles["reads_1.fastq"])
    sum_length = sum(r.length for r in filter.trim(reads, 20, 0))
    assert sum_length == 550000, sum_length


def test_fastq_trim_both():
    reads = files.open(testfiles["reads_1.fastq"])
    sum_length = sum(r.length for r in filter.trim(reads, 10, 10))
    assert sum_length == 550000, sum_length


def test_fastq_filtering():
    reads = files.open(testfiles["reads_1.fastq"])
    num_reads = sum(1 for r in reads.clone())
    sum_length = sum(r.length for r in reads.clone())
    assert num_reads == 10000, num_reads
    assert sum_length == 750000, sum_length


def test_iterating_fastq():
    reader = files.open(testfiles["reads_1.fastq"])
    assert reader is not None
    assert len(list(reader)) == 10000
    reader = files.open(testfiles["reads_1.fastq.gz"])
    assert reader is not None
    assert len(list(reader)) == 10000


def test_reader_cloning():
    reader = files.open(testfiles["reads_1.fastq"])
    assert len(list(reader)) == 10000
    clone = reader.clone()
    assert len(list(clone)) == 10000


def test_open_file():
    reader = files.open_file(testfiles["reads_1.fastq"])
    assert reader is not None
    assert len(reader.readlines()) == 40000


def test_guess_type():
    assert files._guess_type("file.fa") == "fasta"
    assert files._guess_type("file.fasta") == "fasta"
    assert files._guess_type("file.fastq") == "fastq"
    assert files._guess_type("file.map") == "map"
    assert files._guess_type("file.abc") == None
    assert files._guess_type("") == None
    assert files._guess_type(".fa") == "fasta"


def test_open_gzip_file():
    reader = files.open_gzip(testfiles["reads_1.fastq.gz"])
    assert reader is not None
    lines = reader.readlines()
    assert lines is not None
    assert len(lines) == 40000, len(lines)

