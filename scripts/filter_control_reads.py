#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import (NamedBreakpoint, HTSRead, IntervalWrapper,
    intersect_iter)
from TRHelper.Align import check_vs_normal

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify nearby positive/negative breakpoints')
parser.add_argument('breakpoint_file',
    help='TSV file with named breakpoint info')
parser.add_argument('control_read_file',
    help='TSV file with info for clipped reads from the control sample')
parser.add_argument('-d', '--max-breakpoint-dist', type=int, default=0,
    help='maximal distance to control clipped read')
parser.add_argument("-l", "--min-length", type=int, default=5,
    help="minimal length of the clipped tail of the read in the control "
    "sample")

args = parser.parse_args()

with open(args.breakpoint_file) as fid, open(
args.control_read_file) as control_fid:
    bp_iter = (IntervalWrapper(bp.ch, bp.pos - args.max_breakpoint_dist,
        bp.pos + args.max_breakpoint_dist, bp) for bp in
        (NamedBreakpoint(line) for line in fid))
    cr_iter = (IntervalWrapper(cr.rname, cr.ref_locus, cr.ref_locus, 
        cr) for cr in (HTSRead.from_line(line) for line in control_fid))
    for bp, control_reads in intersect_iter(bp_iter, cr_iter,
    one_by_one=False):
        # check that both are the same; i.e., R/R or L/L
        # and that the sequences match
        if check_vs_normal(bp.seq, (r.cseq for r in control_reads if
        bp.bkp==r.bkp and len(r.cseq)>=args.min_length)):
            sys.stdout.write(bp.name.strip("@RL") + "\n")
