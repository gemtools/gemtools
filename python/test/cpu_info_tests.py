import subprocess

import gem
from testfiles import testfiles


def test_i3_compliance():
    file = testfiles["cpuinfo_i3.txt"]
    of = open(file, "r")
    assert gem._is_i3_compliant(of) == True
    of.close()

def test_core2_compliance():
    file = testfiles["cpuinfo_core2.txt"]
    of = open(file, "r")
    assert gem._is_i3_compliant(of) == False
    of.close()

