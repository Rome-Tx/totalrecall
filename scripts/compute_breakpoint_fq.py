#!/usr/bin/env python3

import sys
import argparse

from TRHelper.Aux import HTSRead
from TRHelper.Align import (cluster_clipped_tails, polya_score_qual,
    polya_score_full_length_qual, compute_qual_scores)
from TRHelper.Constants import (POLYA_R, POLYA_Q, POLYA_FULL_Q,
    POLYA_MINSCORE_1P, POLYA_MINSCORE_FULL, POLYA_5TRIM, MAX_Q,
    MAX_PRECLUSTER_LEN, MAX_CLUSTER_INPUT_LEN)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify putative breakpoints from clipped reads')
parser.add_argument('infile', help='TSV file with clipped reads')
parser.add_argument('-d', '--dist', type=int, default=5,
     help='maximal allowed distance between the clipping coordinates '
     'of reads within the pre-cluster')
parser.add_argument("-F", "--polya-fl-mismatch", type=int, default=POLYA_FULL_Q,
    help="mismatch penalty for the full length poly(A) search")
parser.add_argument("-l", "--min-length", type=int, default=5,
    help="minimal length of the clipped tail in each read")
parser.add_argument("-L", "--min-length-consensus", type=int, default=15,
    help="minimal length of the consensus clipped sequence of a cluster")
parser.add_argument("-m", "--polya-match", type=int, default=POLYA_R,
    help="match score for poly(A) search")
parser.add_argument("-M", "--polya-mismatch", type=int, default=POLYA_Q,
    help="mismatch penalty for the poly(A) search")
parser.add_argument("-n", "--min-num-reads", type=int, default=2,
    help="minimal number of reads within the cluster")
parser.add_argument('-o', '--outfile-polya',
    help='output file with coordinates of breakpoints having a poly(A) tail')
parser.add_argument('-O', '--outfile-mask',
    help='output file with coordinates of regions masked for having too large '
    'clusters of clipped reads')
parser.add_argument("-q", "--min-mapq", type=int, default=1,
    help="minimal MAPQ of the breakpoint")
parser.add_argument("-s", "--min-polya-score", type=int,
    default=POLYA_MINSCORE_1P,
    help="minimal alignment score for sequence to be considered poly(A)")
parser.add_argument("-S", "--min-polya-fl-score", type=int,
    default=POLYA_MINSCORE_FULL,
    help="minimal full length alignment score for sequence to be "
    "considered poly(A)")
parser.add_argument("-t", "--polya-trim", type=int, default=POLYA_5TRIM,
    help="number of bases allowed to be trimmed at the 5' end "
    "for the poly(A) search")

args = parser.parse_args()

ofid = open(args.outfile_polya, "w") if args.outfile_polya else None
ofid_mask = open(args.outfile_mask, "w") if args.outfile_mask else None

Smatch, Smism = compute_qual_scores(args.polya_match, args.polya_mismatch,
    MAX_Q)
Smatch_FL, Smism_FL = compute_qual_scores(args.polya_match,
    args.polya_fl_mismatch, MAX_Q)


def precluster_iter(fname, dist, L):
    """
    Iterator over pre-clusters
    Parameters
    fname: TSV file name
    dist: max distance between adjacent reads in the pre-cluster
    L: minimal length for the clipped tail
    """
    rl = {"R":[], "L":[]}
    huge_cluster = {"R": False, "L":False}
    r = None  # needed as a guard for empty input file
    with open(fname) as fid:
        for line in fid:
            r = HTSRead.from_line(line)
            if len(r.cseq) < L:
                continue
            if not rl[r.bkp]:
                rl[r.bkp] = [r]
                huge_cluster[r.bkp] = False
                continue
            prev = rl[r.bkp][-1]
            if (r.rname != prev.rname or r.ref_locus - prev.ref_locus > dist):
                if len(rl[r.bkp]) >= args.min_num_reads:
                    if huge_cluster[r.bkp]:
                        yield rl[r.bkp], True
                    else:
                        for read_list in cluster_clipped_tails(rl[r.bkp],
                        max_cluster_input_len=MAX_CLUSTER_INPUT_LEN):
                            yield read_list, False
                rl[r.bkp] = [r]
                huge_cluster[r.bkp] = False
            elif huge_cluster[r.bkp]:
                rl[r.bkp][-1] = r
            else:
                rl[r.bkp].append(r)
                if len(rl[r.bkp]) >= MAX_PRECLUSTER_LEN:
                    huge_cluster[r.bkp] = True
    if r is not None:
        # correctly process the last clusters
        for bkp in ("R", "L"):
            if len(rl[bkp]) >= args.min_num_reads:
                if huge_cluster[r.bkp]:
                    yield rl[r.bkp], True
                else:
                    for read_list in cluster_clipped_tails(rl[r.bkp],
                    max_cluster_input_len=MAX_CLUSTER_INPUT_LEN):
                        yield read_list, False

N = 1
# aggregate reads and check the tail sequences
for rl, huge_cluster_flag in precluster_iter(args.infile, args.dist,
args.min_length):
    if huge_cluster_flag:
        if ofid_mask is not None:
            ch = rl[0].rname
            start = rl[0].ref_locus
            end = rl[-1].ref_locus
            ofid_mask.write(f"{ch}\t{start}\t{end}\n")
        continue
    if len(rl) <- args.min_num_reads:
        continue
    MQ = max(r.mapq for r in rl)
    if MQ < args.min_mapq:
        continue
    for r in rl[:2]:
        polya_score_fl = polya_score_full_length_qual(r.cseq, r.cqual,
            Smatch_FL, Smism_FL)
        if len(r.cseq) < args.min_length_consensus:
            continue
        if ofid is not None and (polya_score_qual(r.cseq, r.cqual,
        Smatch, Smism, args.polya_trim) >= args.min_polya_score
        or polya_score_fl >= args.min_polya_fl_score):
            ofid.write(f"{r.rname}\t{r.ref_locus}\n")
        else:
            sys.stdout.write(f"@read{N}:{r.rname}:{r.ref_locus}\n{r.cseq}\n"
                f"+\n{r.cqual}\n")
            N += 1
