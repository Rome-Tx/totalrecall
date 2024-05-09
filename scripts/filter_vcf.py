#!/usr/bin/env python3

import sys
import argparse

from TRHelper.Aux import VCF
from TRHelper.Align import polya_score_full_length
from TRHelper.Constants import POLYA_R, POLYA_FULL_Q, POLYA_MINSCORE_FULL


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='filter VCF')
parser.add_argument('-l', '--transposon-length', type=int, default=30,
    help='minimal transposon length')
parser.add_argument('-L', '--line-only', action='store_true',
    help='only keep LINE elements')
parser.add_argument('-n', '--num-reads', type=int, default=3,
    help='minimal number of reads')


args = parser.parse_args()

for line in sys.stdin:
    if line.startswith("#"):
        sys.stdout.write(line)
        continue
    vc = VCF.from_line(line)
    if args.line_only and "LINE1" not in vc.alt:
        # only keep LINEs
        continue
    if "CLIPPED.CLUSTERS" in vc.filters or "LOW_FRAC_OF_DEPTH" in vc.filters:
        continue
    if vc.info.SVLEN < args.transposon_length:
        continue
    if vc.info.BP5NR < args.num_reads:
        continue
    if vc.info.BP3NR >= args.num_reads or polya_score_full_length(
    vc.info.BP3SEQ, POLYA_R, POLYA_FULL_Q) >= POLYA_MINSCORE_FULL:
        sys.stdout.write(f"{vc}\n")
