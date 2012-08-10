#!/usr/bin/env python
import os
import shutil
from nose.tools import with_setup
import gem
from gem import files
from gem import filter
from gem import junctions
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
    mappings = gem.mapper(filter.interleave([input1, input2]), index)
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is None
    assert sum(1 for x in mappings) == 20000


@with_setup(setup_func, cleanup)
def test_interleaved_pair_aligner_run():
    input1 = files.open(testfiles["reads_1.fastq"])
    input2 = files.open(testfiles["reads_2.fastq"])
    mappings = gem.mapper(filter.interleave([input1, input2]), index)
    paired = gem.pairalign(mappings, index)
    assert paired is not None
    assert sum(1 for x in paired) == 20000  ## test dataset does not pair at all


@with_setup(setup_func, cleanup)
def test_sync_splitmapper_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.splitmapper(input, index, results_dir + "/splitmap_out.mapping")
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename == results_dir + "/splitmap_out.mapping"
    assert os.path.exists(results_dir + "/splitmap_out.mapping")
    assert sum(1 for x in mappings) == 10000


@with_setup(setup_func, cleanup)
def test_async_splitmapper_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.splitmapper(input, index)
    assert mappings is not None
    assert mappings.process is not None
    assert mappings.filename is not None
    print mappings.filename
    assert os.path.exists(mappings.filename)
    assert mappings.remove_after_iteration
    assert sum(1 for x in mappings) == 10000
    assert not os.path.exists(mappings.filename)


@with_setup(setup_func, cleanup)
def test_junction_extraction_from_gtf():
    gtf_junctions = list(junctions.from_gtf(testfiles["refseq.gtf"]))
    assert len(gtf_junctions) == 260
    junctions.write_junctions(gtf_junctions, results_dir + "/annotation.junctions")

    ## reread the junctions
    reread = junctions.from_junctions(results_dir + "/annotation.junctions")
    assert reread is not None
    assert len(reread) == 260
    assert len(set(gtf_junctions).intersection(set(reread))) == 260

    ## merge junctions
    merged = junctions.merge_junctions([gtf_junctions, reread])
    assert merged is not None
    assert len(merged) == 260


@with_setup(setup_func, cleanup)
def test_junction_extraction_from_splitmap():
    input = files.open(testfiles["reads_1.fastq"])
    index_hash = testfiles["genome.hash"]
    index = testfiles["genome.gem"]
    gtf_junctions = list(junctions.from_gtf(testfiles["refseq.gtf"]))
    (splitmap, jj) = gem.extract_junctions(input, index, index_hash, merge_with=[gtf_junctions])
    assert splitmap is not None
    assert junctions is not None
    assert len(jj) == 260
    assert sum(1 for x in splitmap) == 10000


@with_setup(setup_func, cleanup)
def test_sync_score_and_validate_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    scored = gem.validate_and_score(mappings, index, results_dir + "/scored.mapping")
    assert scored is not None
    assert scored.process is not None
    assert scored.filename is not None
    assert scored.filename == results_dir + "/scored.mapping"
    count = 0
    for read in scored:
        count += 1
    assert count == 10000


@with_setup(setup_func, cleanup)
def test_gem2sam_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index)
    sam = gem.gem2sam(mappings, index, compact=True)
    assert sam is not None
    assert sam.process is not None
    assert sam.filename is None
    count = 0
    for read in sam:
        count += 1
    assert count == 10000
