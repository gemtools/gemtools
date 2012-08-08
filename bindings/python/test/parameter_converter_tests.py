__author__ = 'Thasso Griebel <thasso.griebel@gmail.com>'

import gem

def test_splice_consensus_conversion():
    assert gem._prepare_splice_consensus_parameter(None) == '"GT"+"AG","CT"+"AC"'
    assert gem._prepare_splice_consensus_parameter([("AA", "BB"), ("CC", "DD")]) == '"AA"+"BB","CC"+"DD"'


def test_quality_parameter():
    assert gem._prepare_quality_parameter(None) == 'ignore'
    assert gem._prepare_quality_parameter(33) == 'offset-33'
    assert gem._prepare_quality_parameter(64) == 'offset-64'


def test_index_parameter():
    assert gem._prepare_index_parameter("my_index.gem") == "my_index.gem"
    assert gem._prepare_index_parameter("my_index") == "my_index.gem"
    assert gem._prepare_index_parameter("my_index", False) == "my_index"
    assert gem._prepare_index_parameter("my_index.gem", False) == "my_index"