#!/usr/bin/env python3

import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Make genome file with contig lengths using info '
        'from SAM/BAM header')
parser.add_argument('aln_file', help='input file')

args = parser.parse_args()

samfile = pysam.AlignmentFile(args.aln_file)
for rname in samfile.references:
	rlen = samfile.get_reference_length(rname)
	sys.stdout.write("{}\t{}\n".format(rname, rlen))
