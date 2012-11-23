from gem import files
from gem import filter
from testfiles import testfiles
import filecmp
import os

__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'


def test_cat():
    reads_1 = files.open(testfiles["reads_1.fastq"])
    reads_2 = files.open(testfiles["reads_2.fastq"])
    num_reads = sum(1 for r in filter.cat([reads_1, reads_2]))
    assert num_reads == 20000

def test_basic_stats():
	# Prepare input files
	maps               = files.open(testfiles["test_unmapped.map"])
	correct_stats_file = testfiles["stats.txt"]
	stats_file         = "stats.txt"
	
	# Run test
	bs = filter.BasicStats(stats_file)
	for aln in bs.basic_stats(maps):
		pass
	bs.print_stats()
	
	# Check results of test
	assert filecmp.cmp(correct_stats_file, stats_file)
	
	# Remove generated test result
	os.remove(stats_file)
