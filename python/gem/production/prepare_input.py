#!/usr/bin/env python
import os
import argparse

import gem


def main():
    parser = argparse.ArgumentParser()

    ## required parameters
    parser.add_argument('-i', '--index', dest="index", help='The path to the gem index', required=True)
    parser.add_argument('-a', '--annotation', dest="annotation", help='The annotation file', required=True)
    parser.add_argument('-m', '--max-read-length', dest="maxlength", help='The maximum read length', required=True)
    parser.add_argument('-t', '--threads', dest="threads", help='Number of threads', default=2)
    parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")

    ## parsing command line arguments
    args = parser.parse_args()

    if args.loglevel is not None:
        gem.loglevel(args.loglevel)

    junctions_out = os.path.basename(args.annotation) + ".junctions"
    index_out = os.path.basename(args.annotation) + ".gem"

    print "Loading Junctions"
    junctions = set(gem.junctions.from_gtf(args.annotation))
    print "%d Junctions loaded from annotation " % (len(junctions))
    gem.junctions.write_junctions(junctions, junctions_out)
    print "Junctions writen to %s " % (junctions_out)

    print "Computing transcriptome..."
    (transcriptome, keys) = gem.compute_transcriptome(args.maxlength, args.index, junctions_out)

    print "Indexing transcriptome"
    gem.index(transcriptome, index_out, threads=args.threads)
    print "Done"


if __name__ == "__main__":
    main()
