#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import (NamedBreakpoint, HTSRead, IntervalWrapper,
    intersect_iter, read_iter)
from TRHelper.Align import calc_num_clusters
from TRHelper.Constants import CC_DIST_SET, CC_CLUSTERING_CUTOFF

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify nearby positive/negative breakpoints')
parser.add_argument('breakpoint_file',
    help='TSV file with named breakpoint info')
parser.add_argument('clip_file',
    help='TSV file with info for clipped reads from the case sample')
parser.add_argument("-l", "--min-length", type=int, default=10,
    help="minimal length of the clipped tail of the read in the case "
    "sample")
parser.add_argument("-m", "--masked-regions", help="file with masked regions")

args = parser.parse_args()

max_breakpoint_dist = max(CC_DIST_SET)

with open(args.breakpoint_file) as fid:
    bp_iter = (IntervalWrapper(bp.ch, bp.pos - max_breakpoint_dist,
        bp.pos + max_breakpoint_dist, bp) for bp in
        (NamedBreakpoint(line) for line in fid))
    cr_iter = (IntervalWrapper(cr.rname, cr.ref_locus, cr.ref_locus, 
        cr) for cr in read_iter(args.clip_file, args.masked_regions))
    for bp, other_reads in intersect_iter(bp_iter, cr_iter,
    one_by_one=False, yield_all=True):
        # check that both are the same; i.e., R/R or L/L
        # and that the sequences match
        relevant_reads = tuple(r for r in other_reads if bp.bkp==r.bkp 
            and len(r.cseq)>=args.min_length)
        N = [calc_num_clusters(bp.seq, (r.cseq for r in relevant_reads if
            abs(bp.pos - r.ref_locus) <= dist),
            num_clusters_cutoff=CC_CLUSTERING_CUTOFF)
            for dist in CC_DIST_SET]
        sys.stdout.write("{}\t{}\n".format(bp.name,
            "\t".join(str(n) for n in N)))
