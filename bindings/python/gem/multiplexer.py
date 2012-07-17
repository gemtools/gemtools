#!/usr/bin/env python
"""Multiplex GEM input files and clean CASAVA 1.8 ID's"""
import sys
import subprocess


def _clean_id(id, label, filtered=False, read=0):
    """Take a read ID and clean it up.
    If label is specified, a GEM label is added.
    If filtered is true, also reads not passing
    casava filter are includes
    """
    s = id.split(" ")
    if read > 0 and not id.endswith("/1") and not id.endswith("/2"):
        id = id.rstrip() + "/" + str(read) + "\n"
    elif len(s) > 1 and (s[1].startswith("1:") or s[1].startswith("2:")):
        if not filtered and s[1][2] == "Y":
            return None
        id = s[0] + "/" + s[1][0] + "\n"
    if label:
        id = "@" + str(label) + " L " + id[1:]
    return id


def _multiplex_entry(input, target, label, filtered=False, read=0):
    """Read fastq format from input stream,
    cleanup the fastq header and + line and
    write to the output stream
    """
    line_count = 0
    inc = True
    for line in input:
        line_count += 1
        if inc:
            if line_count == 1:
                line = _clean_id(line, label, filtered, read)
                if not line:
                    inc = False
                    continue
            if line_count == 3:
                line = "+\n"
            ## write the line
            target.write(line)

        if line_count == 4:
            return True
    return False


def _multiplex_entry_generator(input, label, filtered=False):
    """Read fastq format from input stream,
    cleanup the fastq header and + line and
    write to the output stream
    """
    line_count = 0
    inc = True
    lines = []
    for line in input:
        line_count += 1
        if line_count == 1:
            line = _clean_id(line, label, filtered)
            if not line:
                inc = False
                continue
        if line_count == 3:
            line = "+\n"
        ## write the line
        if inc:
            lines.append(line)

        if line_count == 4:
            if inc:
                return "".join(lines)
            else:
                line_count = 0
                lines = []
                inc = True
    return False


def _open_fd(file):
    if(file.endswith(".gz")):
        return subprocess.Popen(["zcat", file], stdout=subprocess.PIPE).stdout
    else:
        return open(file, 'r')


def multiplex(files, target=None, labels=True, filtered=False, sort=True, add_ids=False):
    if target is None:
        target = sys.stdout
    clfd = False
    if isinstance(target, basestring):
        target = open(target, 'w')
        clfd = True
    if not files:
        raise ValueError("No files specified")
    if len(files) > 2:
        raise ValueError("Can not merge more than 2 files")
    if sort:
        files.sort()
    handles = map(_open_fd, files)

    ## start multiplexing
    label = 0
    if labels:
        label = 1
    second = len(handles) > 1
    count = 0
    read = 0
    second_read = 0
    if second and add_ids:
        read = 1
        second_read = 2

    while _multiplex_entry(handles[0], target, label, filtered, read):
        if label:
            label += 1
        if second:
            _multiplex_entry(handles[1], target, label, filtered, second_read)
        if label:
            label += 1
        count += 1
        if count % 100000 == 0:
            sys.stderr.write("#Processed %d reads\n" % count)

    for h in handles:
        h.close()
    if clfd:
        target.close()


def multiplex_generator(files, labels=True, filtered=False, sort=True):
    if not files:
        raise ValueError("No files specified")
    if len(files) > 2:
        raise ValueError("Can not merge more than 2 files")
    if sort:
        files.sort()
    handles = map(_open_fd, files)

    ## start multiplexing
    label = 0
    if labels:
        label = 1
    second = len(handles) > 1

    while True:
        l1 = _multiplex_entry_generator(handles[0], label, filtered)
        if not l1:
            break
        l2 = ""
        if label:
            label += 1
        if second:
            l2 = _multiplex_entry_generator(handles[1], label, filtered)
            if label:
                label += 1
        yield l1 + l2

    for h in handles:
        h.close()
