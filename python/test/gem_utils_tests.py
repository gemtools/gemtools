#!/usr/bin/env python
import gem.utils as gu
import gem.gemtools as gt
from testfiles import testfiles
import tempfile
import os
import sys
test_mapping = testfiles["test.map"]
test_zipped_mapping = testfiles["test.map.gz"]
test_fastq = testfiles["test.fastq"]


def test_pipeing_simple_process_with_file_handle():
    of = open(test_mapping, "rb")
    p = gu.run_tools([["cat", "-"], ["wc", "-l"]], input=of)
    lines = 0
    while(True):
        s = p.stdout.readline()
        if s is None or len(s) == 0:
            break
        lines = int(s.strip())
    assert lines == 10
    assert p.wait() == 0
    of.close()


def test_pipeing_simple_process_with_file():
    p = gu.run_tools([["cat", "-"], ["wc", "-l"]], input=test_mapping)
    lines = 0
    while(True):
        s = p.stdout.readline()
        if s is None or len(s) == 0:
            break
        lines = int(s.strip())
    assert lines == 10
    assert p.wait() == 0


def test_pipeing_simple_process_with_file_output_file():
    (f, out) = tempfile.mkstemp()
    os.close(f)
    p = gu.run_tools([["cat", "-"], ["wc", "-l"]],
        input=test_mapping, output=out)
    lines = 0
    assert p.wait() == 0
    assert os.path.exists(out)
    with open(out) as f:
        lines = int(f.readline().strip())
    assert lines == 10
    os.remove(out)


def test_pipeline_fastaq_input():
    ff = gt.InputFile(test_mapping)
    p = gu.run_tools([["cat", "-"]], input=ff, force_debug=True)
    lines = 0
    while(True):
        s = p.stdout.readline()
        if s is None or len(s) == 0:
            break
        lines += 1
    assert p.wait() == 0
    assert lines == 40


def test_pipeline_map_input():
    ff = gt.InputFile(test_mapping)
    p = gu.run_tools([["cat", "-"]], input=ff, force_debug=True, write_map=True, clean_id=True)
    lines = 0
    while(True):
        s = p.stdout.readline()
        if s is None or len(s) == 0:
            break
        lines += 1
    assert p.wait() == 0
    assert lines == 10

