#!/usr/bin/env python

import subprocess
import gem.gemtools as gt
from testfiles import testfiles


def test_open_file():
    infile = gt.InputFile(testfiles["bedconvert.map"])
    count = 0
    for tmpl in infile.templates():
        count += 1
        assert tmpl.tag == "WI-ST472_0084:6:1101:1179:2208#CTTGTA", tmpl.tag
        assert tmpl.blocks == 2, tmpl.blocks
    assert count == 1


def test_open_stream():
    p = subprocess.Popen(["cat", testfiles["bedconvert.map"]], stdout=subprocess.PIPE)
    infile = gt.InputFile(stream=p.stdout)
    p.wait()
    count = 0
    for tmpl in infile.templates():
        count += 1
        assert tmpl.tag == "WI-ST472_0084:6:1101:1179:2208#CTTGTA", tmpl.tag
        assert tmpl.blocks == 2, tmpl.blocks
    assert count == 1


def test_open_raw_stream_from_stream():
    p = subprocess.Popen(["cat", testfiles["bedconvert.map"]], stdout=subprocess.PIPE)
    infile = gt.InputFile(stream=p.stdout)
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
