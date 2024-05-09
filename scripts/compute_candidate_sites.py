#!/usr/bin/env python3

import sys
from math import sqrt, erfc
from itertools import chain
import argparse
from TRHelper.Aux import (Breakpoint, TargetSite, Evidence, IntervalWrapper,
    intersect_iter, read_gzipped_genome)

sys.setrecursionlimit(10000)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify nearby positive/negative breakpoints')
parser.add_argument('-d', '--max-breakpoint-dist', type=int, default=50,
    help='maximal breakpoint distance')
parser.add_argument('-n', '--no-orphan', action='store_true',
    help='do not call orphan non-LTR transductions')
parser.add_argument('-r', '--min-reads', type=int, default=1,
    help='minimal number of reads supporting breakpoint')
parser.add_argument('breakpoint_file', help='TSV file with breakpoint info')
parser.add_argument('genome_gz_file', help='gzipped genome fasta file')

args = parser.parse_args()

N = 1

gh = read_gzipped_genome(args.genome_gz_file)

LTR_ANNOT_SET = {Evidence.LTR3, Evidence.LTR5}
nonLTR_5SET = {Evidence.LINE, Evidence.LINE_INV, Evidence.SINE, Evidence.SVA}
TS_LEN_SIGMA = 15
TS_DIST_DENOM = 1 / (TS_LEN_SIGMA * sqrt(2))
TS_LEN_NO_PENALTY = 20

def ts_len_score(bpL, bpR):
    L = bpL.pos - bpR.pos + 1 if bpR.pos <= bpL.pos else bpL.pos -bpR.pos - 1
    if L <= TS_LEN_NO_PENALTY:
        return 1
    return erfc(TS_DIST_DENOM * (L - TS_LEN_NO_PENALTY))

def ts_score(bpL, bpR):
    """
    Odds ratio for the pairing between the L and R breakpoints
    """
    if bpL.annotation is not None and bpL.annotation.transposon == "LTR":
        if bpR.annotation is not None and bpR.annotation.transposon == "LTR":
            if (set((bpL.annotation.evidence, bpR.annotation.evidence))
            == LTR_ANNOT_SET):
                # one 5' end and one 3' end
                return 2**(bpL.annotation.bitscore + bpR.annotation.bitscore
                    ) *  ts_len_score(bpL, bpR)
            else:
                # both 5' end or both 3' end - not a valid retrovirus
                return 0
        else:
            # only one LTR end - not a valid retrovirus
            return 0
    elif bpL.annotation is not None and bpL.annotation.evidence == Evidence.polyA:
        if bpR.annotation is None:
            # orphan transduction - one poly(A) endsys.stderr.write("1\n")
            return 2**(bpL.annotation.bitscore) *  ts_len_score(bpL, bpR)
        elif bpR.annotation.evidence in nonLTR_5SET:
            # non-LTR transposon with a poly(A)
            return 2**(bpL.annotation.bitscore + bpR.annotation.bitscore
                ) *  ts_len_score(bpL, bpR)
        else:
            # invalid pairing
            return 0
    elif bpR.annotation is not None and bpR.annotation.evidence == Evidence.polyA:
        if bpL.annotation is None:
            # orphan transduction - one poly(A) end
            return 2**(bpR.annotation.bitscore) *  ts_len_score(bpL, bpR)
        elif bpL.annotation.evidence in nonLTR_5SET:
            # transposon
            return 2**(bpL.annotation.bitscore + bpR.annotation.bitscore
                ) *  ts_len_score(bpL, bpR)
        else:
            # invalid pairing
            return 0
    else:
        return 0


def bp_cluster_generator(fname, dist):
    """
    Generator of the pre-clusters -- breakpoints within some distance
    Parameters
    fname: TSV file name with breakpoint info
    dist: max distance between breakpoints in a pre-cluster
    """
    bpl = []
    with open(fname) as fid:
        for line in fid:
            bp = Breakpoint.from_line(line)
            if not bpl:
                bpl = [bp]
                continue
            prev = bpl[-1]
            if bp.ch != prev.ch or bp.pos - prev.pos > dist:
                yield bpl
                bpl = [bp]
            else:
                bpl.append(bp)
    # correctly process the last cluster
    if bpl:
        yield bpl


def ts_generator(bplist):
    bpL, bpR, max_likelihood = None, None, 0
    iter_L = (IntervalWrapper(bp.ch, bp.pos - args.max_breakpoint_dist,
        bp.pos + args.max_breakpoint_dist, bp) for bp in bplist
        if bp.bkp == "L")
    iter_R = (IntervalWrapper(bp.ch, bp.pos, bp.pos, bp) for bp in bplist
        if bp.bkp == "R")
    for bpl, bpr in intersect_iter(iter_L, iter_R):
        cur_likelihood = ts_score(bpl, bpr) 
        if cur_likelihood > max_likelihood:
            bpL, bpR, max_likelihood = bpl, bpr, cur_likelihood
    if max_likelihood > 0:
        # found some valid pairing(s), selected the best
        bpL.used = True
        bpR.used = True
        yield (bpL, bpR)
        start, end = sorted((bpL.pos, bpR.pos))
        bplist_L = [bp for bp in bplist if bp.pos < start]
        bplist_R = [bp for bp in bplist if bp.pos > end]
        yield from ts_generator(bplist_L)
        yield from ts_generator(bplist_R)


def ts_2step_generator(bplist):
    bplist_hc = [bp for bp in bplist if bp.nr > 1]  # at least two reads
    ts_hc = list(ts_generator(bplist_hc))
    ts_hc.sort(key = lambda x: min(bp.pos for bp in x))
    bplist_other = [bp for bp in bplist if not bp.used]
    ts_other = list(ts_generator(bplist_other))
    ts_other.sort(key = lambda x: min(bp.pos for bp in x))
    yield from ts_hc
    hc_iter = (IntervalWrapper("ch", min(bpl.pos, bpr.pos),
        max(bpl.pos, bpr.pos), None) for bpl, bpr in ts_hc)
    other_iter = (IntervalWrapper("ch", min(bpl.pos, bpr.pos),
        max(bpl.pos, bpr.pos), (bpl, bpr)) for bpl, bpr in ts_other)
    for (bpL, bpR), x in intersect_iter(other_iter, hc_iter, one_by_one=False,
    yield_all=True):
        if not x and (bpL.nr > 1 or bpR.nr > 1):
            yield (bpL, bpR)


for bpl_full in bp_cluster_generator(args.breakpoint_file,
args.max_breakpoint_dist):
    bpl = [bp for bp in bpl_full if bp.nr >= args.min_reads
        and (bp.annotation is not None or not args.no_orphan)]
    ltr_bpl = [bp for bp in bpl if bp.annotation is not None
        and bp.annotation.transposon == "LTR"]
    nonltr_bpl = [bp for bp in bpl if bp.annotation is None
        or bp.annotation.transposon != "LTR"]
    for bp1, bp2 in chain(ts_2step_generator(ltr_bpl),
    ts_2step_generator(nonltr_bpl)):
        name = f"IS.{N}"
        ts = TargetSite.from_breakpoints(name, bp1, bp2, gh)
        sys.stdout.write(f"{ts}\n")
        N += 1
