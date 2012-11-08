from gem import files
from gem import filter
from testfiles import testfiles

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'


def test_cat():
    reads_1 = files.open(testfiles["reads_1.fastq"])
    reads_2 = files.open(testfiles["reads_2.fastq"])
    num_reads = sum(1 for r in filter.cat([reads_1, reads_2]))
    assert num_reads == 20000


