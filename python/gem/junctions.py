#!/usr/bin/env python
"""Module that handles junction sites and provides extraction
methods to prepare junction files used by the split mapper.

Junctions can be extracted from GTF file and split mappings.
"""
import os
import re
import subprocess
import sys
import gem.filter
import gem
import logging

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
                if exon.strand == "-":
                    yield JunctionSite(exon, last_exon)
                else:
                    yield JunctionSite(last_exon, exon)
            last_exon = exon


class JunctionSite(object):
    """A junction site is created from two exons. The
    string representation is in the proper format and
    and junction sites are comparable with == and !=
    """

    def __init__(self, start_exon=None, stop_exon=None, line=None, coverage=0):
        """Create a junction site from the start and stop exon"""
        if start_exon is not None:
            self.descriptor = self.__descriptor(start_exon, stop_exon)
        else:
            ## parse from line
            self.descriptor = self.__descriptor_from_line(line)

        self.coverage = coverage
        self.hash = hash(self._unique())  # was from str(self)

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
            if strand_1 == "+":
                strand_1 = "-"
            else:
                strand_1 = "+"
            if strand_2 == "+":
                strand_2 = "-"
            else:
                strand_2 = "+"
        desc = [chr_1, strand_1, pos_1, chr_2, strand_2, pos_2]
        return desc

    def _unique(self):
        d = self.descriptor
        return "%s %d %s %d" % (d[0], min(d[2], d[5]), d[3], max(d[2], d[5]))

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
    """
    Filter a junctions file by chromosomes
    contained in the given index. The result is written to
    output

    @param index: path to the gem index
    @type index: string
    @param junctions: path to junctions file
    @type junctions: string
    @param output: path to the output file
    @type output: output
    """
    chrs = _get_chromosomes(index)
    jf = open(output, 'w')
    with open(os.path.abspath(junctions), 'r') as f:
        for line in f:
            split = line.split("\t")
            if split[0] in chrs and split[3] in chrs:
                jf.write(line)
    jf.close()


def merge_junctions(junction_sets):
    """
    Merge the given set of junctions and return
    the union

    @param junction_sets: list of list of junctions
    @type junction_sets: list
    @return: merged junctions
    @rtype: set
    """
    s = set([])
    for js in junction_sets:
        if s is None:
            s = set(js)
        else:
            s = s.union(js)
    return s


def write_junctions(junctions, output=None, index=None):
    """Write the given list of JunctionSites to output file.
    If an index is specified, only junctions that are contained
    in the index are written

    @param junctions: the list of junction sites
    @type junctions: list
    @param output: filename or none for stdout
    @type output: string
    @param index: path to the index file
    @type index: string
    @return: void
    @rtype: void
    """
    chrs = None
    if index:
        chrs = _get_chromosomes(index)
    jf = sys.stdout
    if output:
        jf = open(output, 'w')

    ## write the junctions
    for junction in junctions:
        if (chrs is None) or (junction.descriptor[0] in chrs and junction.descriptor[3] in chrs):
            jf.write(str(junction))
            jf.write("\n")
    if output:
        jf.close()


def _get_chromosomes(index):
    """Return the list of chromosome names contained in the given index file"""
    if not isinstance(index, basestring):
        raise ValueError("The index must be a string")
#    if index.endswith(".gem"):
#        index = index[:-4]
    pa = subprocess.Popen(gem.executables['gem-info'] + ' ' + index, stdout=subprocess.PIPE, shell=True)
    chrs = []
    for line in pa.stdout:
        line = line.strip()
        split = re.split("\s+", line)
        if len(split) == 5 and split[0][0] == "#":
            chrs.append(split[1][1:-1])
    if pa.wait() != 0:
        raise ValueError("gem-info failed")
    return sorted(set(chrs))


def filter_by_distance(junctions, min_distance, max_distance):
    """Yields the junction sites that have a distance less than equal max_distance"""
    for j in junctions:
        d = abs(j.descriptor[2] - j.descriptor[5])
        if min_distance <= d and d <= max_distance:
            yield j


def from_junctions(junctions_file):
    """
    Read junctions file and return a list of JunctionSites
    @param junctions_file: path to the junctions file
    @type junctions_file: string
    @return: list of junction sites
    @rtype: list
    """
    junctions = []
    with open(os.path.abspath(junctions_file), 'r') as f:
        for line in f:
            junctions.append(JunctionSite(line=line))
    return junctions


def from_gtf(annotation):
    """Extract junctions from given gtf annotation"""
    def __get_id(tr):
        if tr.startswith("transcript_id"):
            t = tr.split(" ")[1]
            if t[-1] == ";": t = t[:-1]
            if t[0] == '"': t = t[1:]
            if t[-1] == '"': t = t[:-1]
            return t
        return None

    def __extract_transcript(cell, check_first=0):
        split = [s.strip() for s in cell.split(";")]
        t = None
        if len(split) < check_first:
            t = __get_id(split[check_first])
        if t is not None:
            return t
        for s in split:
            t = __get_id(s)
            if t is not None:
                return t
        return None

    in_fd = gem.files.open_file(annotation)  ##open(annotation, 'rb')
    #trans_re = re.compile('.*transcript_id "(.*)";.*')
    junctions = {}
    line_count = 0
    for line in in_fd:
        if len(line.strip()) == 0 or line.strip()[0] == '#':
            logging.debug("Skipping comment line")
            # skip comment
            continue
        line_count += 1
        split = line.split("\t")
        if len(split) < 8:
            logging.warning("GTF line %d has < 8 fields, skipping" % line_count)
            continue
        if split[2] != "exon":
            continue
        #id = trans_re.match(split[8]).group(1) + split[0]
        try:
            id = __extract_transcript(split[8].strip())
            if id is None:
                logging.error("Failed to extract a transcript id from line %d" % line_count)
                continue
            id = id + split[0]
        except:
            logging.error("Failed to extract a transcript id from line %d" % line_count)
            continue
        start = int(split[3])
        end = int(split[4])
        if start != end:
            exon = Exon(split[0], id, start, end, split[6])
            if exon.stop - exon.start > 1:
                junctions.setdefault(id, Junction()).append(exon)
    logging.debug("Found %d raw sites, extracting unique" % (len(junctions)))
    unique = sorted(set((
        site for sites in
        (junction.sites() for k, junction in junctions.items())
        for site in sites)))
    logging.debug("Found %d sites" % (len(unique)))
    for j in unique:
        yield j
