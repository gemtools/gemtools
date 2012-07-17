#!/usr/bin/env python
import sys
import gem.gemtools as gt

if __name__ == "__main__":
    infile = gt.open_stream(sys.stdin)
    for tempalte in infile:
        beds = []
        read = "/1"
        for alignment in tempalte.blocks():
            rname = alignment.tag
            if not rname.endswith("/1"):
                    rname += read
                    read = "/2"
            mm = alignment.mappings()
            for mapping in mm:
                block_count = 0
                block_sizes = []
                block_starts = []
                start = -1
                end = -1
                name = None
                strand = "+"
                for m, junction, distance in mapping:
                    if start == -1:
                        name = m.seq_name
                        start = m.position - 1
                        end = m.position + m.global_length - 1
                        if m.direction != 0:
                            strand = "-"
                    block_count += 1
                    block_starts.append(m.position - 1 - start)
                    block_sizes.append(m.length)
                assert start + block_starts[-1] + block_sizes[-1] == end
                beds.append("%s\t%d\t%d\t%s\t0\t%s\t.\t.\t0,0,0\t%d\t%s\t%s" % (name, start, end, rname, strand,
                                                                            block_count,
                                                                            ",".join([str(c) for c in block_sizes]),
                                                                            ",".join([str(c) for c in block_starts])))
        for line in set(beds):
            print line
