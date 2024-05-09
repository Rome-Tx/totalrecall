#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement
from TRHelper import Align
from TRHelper.Aux import NamedBreakpoint

LC_SCORE_THRESHOLD = 5

SUBSEQ_LEN = 10
NUM_OFFSETS = 5

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='compute low complexity filters')
parser.add_argument('input_bn', help='basename ')
parser.add_argument('-s', '--lc-score', type=int, default=LC_SCORE_THRESHOLD,
    help='threshold for running low complexity score')

args = parser.parse_args()

nbp_file = args.input_bn + '.ins.nbp.tsv'
dust_file = args.input_bn + '.bps.dust'


def lc_max_score(seq):
    S = 0
    MS = 0  # maximal score
    for ch in seq:
        if ch in "ACGT":
            S -= 1
        else:
            S += 1
            MS = max(S, MS)
    return MS
        
    
with open(nbp_file) as fid:
    bps = {bp.name:bp.seq for bp in (NamedBreakpoint(line) for line in fid)}
        

for rec in SeqIO.parse(dust_file, "fasta"):
    seq = str(rec.seq)
    rid = rec.id
    L = len(seq) // 2
    seqL = reverse_complement(seq[:L])
    seqR = seq[L:]
    # check for masking of the genomic sequence
    if int(lc_max_score(seqL) >= args.lc_score) + int(lc_max_score(seqR
    ) >= args.lc_score) == 1:
    # exactly one side is low complexity
        sys.stdout.write("{}\n".format(rid.strip("RL@")))
        continue
    # check for identity between the clipped tail and the genomic sequence
    # (will be the case for repeat read with clipped end)
    seqClip = bps[rid]
    if len(seqClip) < SUBSEQ_LEN:
        continue
    seqClip = seqClip[:SUBSEQ_LEN]
    num_sequences = min(NUM_OFFSETS, L - len(seqClip) + 1)
    genSeqL = [seqL[L-len(seqClip)-k:L-k] for k in range(num_sequences)]
    if Align.check_vs_normal(seqClip, genSeqL):
        sys.stdout.write("{}\n".format(rid.strip("RL@")))
