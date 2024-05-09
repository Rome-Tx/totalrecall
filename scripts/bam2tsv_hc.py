#!/usr/bin/env python3

import sys
import argparse
from collections import deque
import pysam
from Bio.Seq import reverse_complement
from TRHelper.Aux import HTSRead, Breakpoint
from TRHelper.Align import compute_clipping_from_qual
from TRHelper.Constants import Q_OFFSET

BAM_SOFT_CLIP = pysam.CSOFT_CLIP.value
BAM_HARD_CLIP = pysam.CHARD_CLIP.value

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Convert BAM/SAM to TSV')
parser.add_argument('aln_file', help='input file')
parser.add_argument('-l', '--min-map-length', type=int, default=10,
    help='minimal distance between the start and end coordinates')
parser.add_argument('-n', '--sample-name', default="NA",
    help='sample name')
parser.add_argument('-t', '--threads', type=int, default=1,
    help='number of threads to use')
parser.add_argument('-q', '--min-end-qual', type=int, default=3,
    help='minimal quality of the end base of the clipped region')
parser.add_argument('-Q', '--min-mapq', type=int, default=0,
    help='minimal mapping quality for the read to be output')


args = parser.parse_args()


def aln_group_generator(aln_iterator):
    """
    Assume that aln_iterator is sorted/collated by read name
    and yield the group of alignments for each read
    """
    AG = deque()
    for aln in aln_iterator:
        if aln.is_unmapped:
            continue
        if aln.mapping_quality < args.min_mapq:
            continue
        if aln.reference_end - aln.reference_start < args.min_map_length:
            continue
        if not AG or AG[-1].query_name == aln.query_name:
            AG.append(aln)
        else:
            yield AG
            AG = deque()
            AG.append(aln)
    if AG:
        yield AG


def primary_alignment(AL):
    for aln in AL:
        if not aln.is_supplementary:
            return aln
    return None


def compute_clipping_info(r, qual, min_qual):
    """
    find out how many bases are soft clipped on each end
    """
    len_L, len_R = 0, 0
    if r.is_supplementary:
        if (r.cigartuples[0][0] == BAM_HARD_CLIP or
        r.cigartuples[0][0] == BAM_SOFT_CLIP):
            len_L = r.cigartuples[0][1]
        if (r.cigartuples[-1][0] == BAM_HARD_CLIP or
        r.cigartuples[-1][0] == BAM_SOFT_CLIP):
            len_R = r.cigartuples[-1][1]
    else:
        if r.cigartuples[0][0] == BAM_SOFT_CLIP:
            len_L = r.cigartuples[0][1]
        if r.cigartuples[-1][0] == BAM_SOFT_CLIP:
            len_R = r.cigartuples[-1][1]
    if r.is_reverse:
        len_L, len_R = len_R, len_L
    c_left, c_right = compute_clipping_from_qual(qual, len_L, len_R, min_qual)
    return len_L, len_R, c_left, c_right


def clipped_htsr_generator(AL):
    pa = primary_alignment(AL)
    if pa is None:
        return
    seq, qual = pa.get_forward_sequence(), bytes(pa.get_forward_qualities())
    for aln in AL:
        len_L, len_R, c_left, c_right = compute_clipping_info(aln, qual,
            args.min_end_qual)
        if c_left == c_right == 0:
            continue
        if c_left >= c_right:
            # tie -> assume the 5' end is the major clipped end
            p2 = len_L
            p1 = p2 - c_left
            clip_seq = seq[p1:p2]
            clip_seq = reverse_complement(clip_seq)
            clip_qual = "".join(chr(q + Q_OFFSET) for q in qual[p1:p2])
            clip_qual = clip_qual[::-1]
            bkp = "L" if aln.is_reverse else "R"
        else:
            p1 = len(seq) - c_right
            p2 = p1 + c_right
            clip_seq = seq[p1:p2]
            clip_qual = "".join(chr(q + Q_OFFSET) for q in qual[p1:p2])
            bkp = "R" if aln.is_reverse else "L"
        yield HTSRead.from_aligned_segment(aln, c_right if aln.is_reverse else
            c_left, c_left if aln.is_reverse else c_right, bkp,
            args.sample_name, clip_seq, clip_qual)

# process the input sam file
samfile = pysam.AlignmentFile(sys.stdin if args.aln_file == "-"
    else args.aln_file, threads=args.threads)
for AG in aln_group_generator(samfile):
    R1 = [aln for aln in AG if aln.is_read1]
    R2 = [aln for aln in AG if aln.is_read2]
    for htsr in clipped_htsr_generator(R1):
        sys.stdout.write(f"{htsr}\n")
    for htsr in clipped_htsr_generator(R2):
        sys.stdout.write(f"{htsr}\n")
