#!/usr/bin/env python3

import sys
import argparse

from TRHelper.Aux import HTSRead, Breakpoint, read_iter
from TRHelper.Align import cluster_clipped_tails, polya_score_qual
from TRHelper.Constants import MAX_CLUSTER_INPUT_LEN

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Identify putative breakpoints from clipped reads')
parser.add_argument('infile', help='TSV file with clipped reads')
parser.add_argument('bp_file', help='output file for TSV representation of '
    'breakpoints')
parser.add_argument('fa_file', help='output FASTA file for clipped sequences')
parser.add_argument('-d', '--dist', type=int, default=5,
     help='maximal allowed distance between the clipping coordinates '
     'of reads within the pre-cluster')
parser.add_argument("-l", "--min-length", type=int, default=2,
    help="minimal length of the clipped tail in each read")
parser.add_argument("-L", "--min-length-consensus", type=int, default=10,
    help="minimal length of the consensus clipped sequence of a cluster")
parser.add_argument("-n", "--min-num-reads", type=int, default=2,
    help="minimal number of reads within the cluster")
parser.add_argument("-m", "--masked-regions", help="file with masked regions")
parser.add_argument("-q", "--min-mapq", type=int, default=0,
    help="minimal MAPQ of the breakpoint")

args = parser.parse_args()
        

def precluster_iter(fname, dist, L):
    """
    Iterator over pre-clusters
    Parameters
    fname: TSV file name
    dist: max distance between adjacent reads in the pre-cluster
    L: minimal length for the clipped tail
    """
    rl = {"R":[], "L":[]}
    r = None  # needed as a guard for empty input file
    for r in read_iter(fname, args.masked_regions):
        if len(r.cseq) < L:
            continue
        if not rl[r.bkp]:
            rl[r.bkp] = [r]
            continue
        prev = rl[r.bkp][-1]
        if r.rname != prev.rname or r.ref_locus - prev.ref_locus > dist:
            if len(rl[r.bkp]) >= args.min_num_reads:
                for read_list in cluster_clipped_tails(rl[r.bkp],
                max_cluster_input_len=MAX_CLUSTER_INPUT_LEN):
                    yield read_list
            rl[r.bkp] = [r]
        else:
            rl[r.bkp].append(r)
    if r is not None:
        # correctly process the last clusters
        for bkp in ("R", "L"):
            if len(rl[bkp]) >= args.min_num_reads:
                for read_list in cluster_clipped_tails(rl[bkp],
                max_cluster_input_len=MAX_CLUSTER_INPUT_LEN):
                    yield read_list

N=0
# aggregate reads and check the tail sequences
with open(args.bp_file, "w") as fh_bp, open(args.fa_file, "w") as fh_fa:
    for rl in precluster_iter(args.infile, args.dist, args.min_length):
        if len(rl) <- args.min_num_reads:
            continue
        if len(rl[0].cseq) < args.min_length_consensus:
            continue
        MQ = max(r.mapq for r in rl)
        if MQ < args.min_mapq:
            continue
        N += 1
        bp = Breakpoint.from_cluster(rl, N)
        fh_bp.write(f"{bp}\n")
        fh_fa.write(f"{bp.clipped_fasta()}")
