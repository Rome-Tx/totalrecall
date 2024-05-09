#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import HTSRead, IntervalWrapper, intersect_iter

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify nearby positive/negative breakpoints')
parser.add_argument('clip_read_file',
    help='TSV file with info for clipped reads')
parser.add_argument('locus_tsv', help='TSV file with loci of interest')
parser.add_argument('-e', '--ext-length', type=int, default=0,
    help='extension length for regions of interest')

args = parser.parse_args()

class Locus:
    def __init__(self, line):
        S = line.split()
        self.ch = S[0]
        self.pos = int(S[1])

with open(args.clip_read_file) as fid, open(args.locus_tsv) as locus_fid:
    locus_iter = (IntervalWrapper(loc.ch, loc.pos - args.ext_length, 
        loc.pos + args.ext_length, loc)
        for loc in (Locus(line) for line in locus_fid))
    cr_iter = (IntervalWrapper(cr.rname, cr.ref_locus, cr.ref_locus, 
        cr) for cr in (HTSRead.from_line(line) for line in fid))
    for cr, bp in intersect_iter(cr_iter, locus_iter, one_by_one=False):
        sys.stdout.write(f"{cr}\n")
