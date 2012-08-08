#!/usr/bin/env python
import os
import re
import subprocess
import gem.filter

"""Module that handles junction sites and provides extraction
methods to prepare junction files used by the split mapper.

Junctions can be extracted from GTF file and split mappings.
"""

class Exon(object):
    """Exons representation"""
    def __init__(self, chr, id, start, stop, strand):
        """Create a new exon

        Arguments:
        chr -- the name of the chromosome
        id  -- the transcript id
        start -- the genomic start position
        stop  -- the genomic stop position
        strand -- the strand, +,-,F,R are allowed
        """
        self.chr = chr
        self.id = id
        self.start = start
        self.stop = stop
        self.strand = strand
        if self.strand == "F": self.strand = "+"
        elif self.strand == "R": self.strand = "-"
        assert self.strand in ["+", "-"]

        if self.start > self.stop:
            self.start, self.stop = self.stop, self.start
            if self.strand == "+": self.strand = "-"
            else: self.strand = "+"

class Junction(object):
    """A junction holds a list of exons and can
    create a list of compatible sites. The string representation
    is in the valid gem format and sorted by exon start positions
    """
    def __init__(self):
        """Create a new junction
        """
        self.exons = []

    def append(self, exon):
        self.exons.append(exon)

    def sites(self):
        self.exons.sort(key=lambda x: x.start)
        last_exon = None
        for i, exon in enumerate(self.exons):
            if i != 0:
                yield JunctionSite(last_exon, exon)
            last_exon = exon

class JunctionSite(object):
    """A junction site is created from two exons. The
    string representation is in the proper format and
    and junction sites are comparable with == and !=
    """
    def __init__(self, start_exon=None, stop_exon=None, line=None):
        """Create a junction site from the start and stop exon"""
        if start_exon is not None:
            self.descriptor = self.__descriptor(start_exon, stop_exon)
        else:
            ## parse from line
            self.descriptor = self.__descriptor_from_line(line)

        self.hash = hash(str(self))


    def __descriptor_from_line(self, line):
        s = line.split("\t")
        chr_1 = s[0]
        strand_1 = s[1]
        pos_1 = long(s[2])
        chr_2 = s[3]
        strand_2 = s[4]
        pos_2 = long(s[5])

        if pos_1 > pos_2:
            pos_1, pos_2 = pos_2, pos_1
            strand_1 = self.__invert(strand_1)
            strand_2 = self.__invert(strand_2)
        desc = [chr_1, strand_1, pos_1, chr_2, strand_2, pos_2]
        return desc

    def __invert(self, strand):
        if strand == "+": return "-"
        else: return "+"


    def __descriptor(self, start_exon, stop_exon):
        desc = []
        if start_exon.strand == "+":
            desc.extend([start_exon.chr, start_exon.strand, start_exon.stop])
        else:
            desc.extend([start_exon.chr, start_exon.strand, start_exon.start])
        if stop_exon.strand == "+":
            desc.extend([stop_exon.chr, stop_exon.strand, stop_exon.start])
        else:
            desc.extend([stop_exon.chr, stop_exon.strand, stop_exon.stop])
        return desc

    def __str__(self):
        return "\t".join([str(s) for s in self.descriptor])

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.descriptor == other.descriptor)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        s = self.descriptor
        o = other.descriptor
        if s[0] < o[0]:
            return True
        if s[0] == o[0] and s[2] < o[2]:
            return True
        ## compare the string representation
        return str(s) < str(o)

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return self > other or self == other

    def __hash__(self):
        return self.hash

def filter_junctions(index, junctions, output):
    chrs = _get_chromosomes(index)
    jf = open(output, 'w')
    with open(os.path.abspath(junctions), 'r') as f:
        for line in f:
            split = line.split("\t")
            if split[0] in chrs and split[3] in chrs:
                jf.write(line)
    jf.close()
    return os.path.abspath(output)


def write_junctions(junctions, index, output):
    chrs = _get_chromosomes(index)
    jf = open(output, 'w')
    s = None
    for js in junctions:
        if s is None:
            s = set(js)
        else:
            s = s.union(js)

    for junction in s:
        if junction.descriptor[0] in chrs and junction.descriptor[3] in chrs:
            jf.write(str(junction))
            jf.write("\n")
    jf.close()
    return os.path.abspath(output)


def _get_chromosomes(index):
    """Return the list of chromosome names contained in the given index file"""
    if not isinstance(index, basestring):
        raise ValueError("The index must be a string")
    if index.endswith(".gem"):
        index = index[:-4]
    pa = subprocess.Popen('gem-info '+index, stdout=subprocess.PIPE,shell=True)
    chrs = []
    for line in pa.stdout:
        line = line.strip()
        split = re.split("\s+",line)
        if len(split) == 5 and split[0][0] == "#":
            chrs.append(split[1][1:-1])
    if pa.wait() != 0:
        raise ValueError("gem-info failed")
    return sorted(set(chrs))


def from_gtf(annotation):
    """Extract junctions from given gtf annotation"""
    in_fd = gem.filter.zcat(annotation)  ##open(annotation, 'rb')
    trans_re = re.compile('.*transcript_id "(.*)";.*')
    junctions = {}
    for line in in_fd:
      split = line.split("\t")
      if split[2] != "exon":
        continue
      id = trans_re.match(split[8]).group(1)+split[0]
      exon = Exon(split[0], id, int(split[3]), int(split[4]), split[6])
      junctions.setdefault(id, Junction()).append(exon)

    for j in sorted(set((
                site for sites in
                (junction.sites() for k,junction in junctions.items())
                for site in sites))):
        yield j
