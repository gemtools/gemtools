#!/usr/bin/env python
import os
import shutil
from nose.tools import with_setup
import gem
from gem import files
from gem import filter
from gem import junctions
from testfiles import testfiles
import gem.gemtools as gt
import sys

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
def test_merging_maps_files():
    input = testfiles["test_merge_target.map"]
    source_1 = testfiles["test_merge_source_1.map"]
    source_2 = testfiles["test_merge_source_2.map"]
    #result = gem.mapper(input, index, results_dir + "/merged.mapping")
    merger = gem.merge(input, [source_1, source_2])
    count = 0
    for read in merger:
        count += 1
        if read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13866:124450":
            assert read.to_map().split("\t")[3] == "0+1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr2:-:162359617:74G"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13753:124452":
            assert read.to_map().split("\t")[3] == "0+1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr9:+:38397301:54G20"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:14211:124259":
            assert read.to_map().split("\t")[3] == "1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr15:+:72492866:75"
    assert count == 10


@with_setup(setup_func, cleanup)
def test_merging_maps():
    input = files.open(testfiles["test_merge_target.map"])
    source_1 = files.open(testfiles["test_merge_source_1.map"])
    source_2 = files.open(testfiles["test_merge_source_2.map"])
    #result = gem.mapper(input, index, results_dir + "/merged.mapping")
    merger = gem.merge(input, [source_1, source_2])
    count = 0
    for read in merger:
        count += 1
        if read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13866:124450":
            assert read.to_map().split("\t")[3] == "0+1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr2:-:162359617:74G"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13753:124452":
            assert read.to_map().split("\t")[3] == "0+1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr9:+:38397301:54G20"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:14211:124259":
            assert read.to_map().split("\t")[3] == "1", read.to_map()
            assert read.to_map().split("\t")[4] == "chr15:+:72492866:75"
    assert count == 10


@with_setup(setup_func, cleanup)
def test_merging_maps_chr21():
    s = ["chr21_mapping_initial.map",
          "chr21_mapping_denovo.map",
          "chr21_mapping_initial_split.map",
          "chr21_mapping_trim_20.map",
          "chr21_mapping_trim_20_split.map"
         ]
    fs = []
    for f in s:
        fs.append(files.open(testfiles[f]))
    m = gem.merge(fs[0], fs[1:])
    count = 0
    ms = 0
    for read in m:
        count += 1
        ms += read.num_maps
    assert count == 2000, count
    assert ms == 3590, ms

@with_setup(setup_func, cleanup)
def test_merging_maps_to_file():
    input = files.open(testfiles["test_merge_target.map"])
    source_1 = files.open(testfiles["test_merge_source_1.map"])
    source_2 = files.open(testfiles["test_merge_source_2.map"])
    result = results_dir + "/merged.mapping"
    merger = gem.merge(input, [source_1, source_2], output=result)
    count = 0
    for read in merger:
        count += 1
        if read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13866:124450 1:N:0:GCCAAT":
            assert read.to_map().split("\t")[3] == "0+1"
            assert read.to_map().split("\t")[4] == "chr2:-:162359617:74G"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:13753:124452 1:N:0:GCCAAT":
            assert read.to_map().split("\t")[3] == "0+1"
            assert read.to_map().split("\t")[4] == "chr9:+:38397301:54G20"
        elif read.tag == "HWI-ST661:153:D0FTJACXX:2:1102:14211:124259 1:N:0:GCCAAT":
            assert read.to_map().split("\t")[3] == "1"
            assert read.to_map().split("\t")[4] == "chr15:+:72492866:75"
    assert count == 10
    assert os.path.exists(result)

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


# @with_setup(setup_func, cleanup)
# def test_async_mapper_pipes_with_filter():
#     def on_half_filter(reads):
#         count = 0
#         for read in reads:
#             count += 1
#             if count % 2 == 0: yield read


#     input = files.open(testfiles["reads_1.fastq"])
#     initial = gem.mapper(input, index)
#     mappings = gem.mapper(on_half_filter(initial), index, results_dir + "/piped_out.mapping")

#     assert mappings is not None
#     assert mappings.process is not None
#     assert mappings.filename is not None
#     assert mappings.filename == results_dir + "/piped_out.mapping"
#     assert os.path.exists(results_dir + "/piped_out.mapping")
#     assert sum(1 for x in mappings) == 5000


@with_setup(setup_func, cleanup)
def test_writing_fastq():
    output = results_dir + "/print_fastq.out"
    input1 = files.open(testfiles["reads_1.fastq"])
    input1.write_stream(gt.OutputFile(output), write_map=False)
    s = 0
    with open(output) as f:
        for l in f:
            s += 1
    assert s == 40000


@with_setup(setup_func)
def test_writing_fastq_interleaved():
    output = results_dir + "/print_fastq.out"
    input1 = files.open(testfiles["reads_1.fastq"])
    input2 = files.open(testfiles["reads_2.fastq"])
    filter.interleave([input1, input2]).write_stream(gt.OutputFile(output), write_map=False)
    s = 0
    with open(output) as f:
        for l in f:
            s += 1
    assert s == 80000, s


@with_setup(setup_func, cleanup)
def test_interleaving():
    input1 = files.open(testfiles["reads_1.fastq"])
    input2 = files.open(testfiles["reads_2.fastq"])
    s = 0
    for t in filter.interleave([input1, input2]):
        s += 1
    assert s == 20000


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
    assert sum(1 for x in paired) == 20000  # test dataset does not pair at all


@with_setup(setup_func, cleanup)
def test_sync_splitmapper_execution():
    gem.loglevel("debug")
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
    assert mappings.filename is None
    assert sum(1 for x in mappings) == 10000


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
    index = testfiles["genome.gem"]
    gtf_junctions = set(junctions.from_gtf(testfiles["refseq.gtf"]))
    jj = gem.extract_junctions(input, index, merge_with=gtf_junctions)
    assert junctions is not None
    assert len(jj) == 260


@with_setup(setup_func, cleanup)
def test_quality_pass_on_execution():
    input = files.open(testfiles["reads_1.fastq"])
    mappings = gem.mapper(input, index, output=results_dir+"/quality_passon_mapping.map")
    assert mappings.quality == "offset-33", "Quality should be 'offset-33' but is %s" % (str(mappings.quality))


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


@with_setup(setup_func, cleanup)
def test_sam2bam_no_sort_sync_execution():
    result = results_dir + "/reads_1.bam"
    input = files.open(testfiles["reads_1.sam"])
    bam = gem.sam2bam(input, result, False)
    assert bam is not None
    assert os.path.exists(result)
    assert sum(1 for x in bam) == 10000


@with_setup(setup_func, cleanup)
def test_sam2bam_sort_sync_execution():
    result = results_dir + "/reads_1.bam"
    input = files.open(testfiles["reads_1.sam"])
    bam = gem.sam2bam(input, result, True)
    assert bam is not None
    assert os.path.exists(result)
    assert sum(1 for x in bam) == 10000


@with_setup(setup_func, cleanup)
def test_sam2bam_sort_async_execution():
    input = files.open(testfiles["reads_1.sam"])
    bam = gem.sam2bam(input, sorted=True)
    assert bam is not None
    assert sum(1 for x in bam) == 10000


def test_parameter_extension():
    pa = []
    gem._extend_parameters(pa, "--test me")
    assert pa == ["--test", "me"]
    pa = []
    gem._extend_parameters(pa, ["--test", "me"])
    assert pa == ["--test", "me"]
