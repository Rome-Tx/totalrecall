#!/usr/bin/env python3

import sys
import argparse
import pysam
from Bio.Seq import reverse_complement
from TRHelper.Aux import HTSRead, Breakpoint

BAM_SOFT_CLIP = pysam.CSOFT_CLIP.value
from TRHelper.Constants import Q_OFFSET

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Convert BAM/SAM to TSV')
parser.add_argument('aln_file', help='input file')
parser.add_argument("-c", "--as-clipped", action='store_true',
                    help='compute mapping locus using clipping '
                    'info instead of discordant mate info')
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


# def compute_clipping_from_qual(Q, min_qual):
    # max_score_pos = 0
    # max_score = 0
    # score = 0
    # for pos, q in enumerate(Q, 1):
        # if q >= min_qual:
            # score += 1
            # if score > max_score:
                # max_score = score
                # max_score_pos = pos
        # else:
            # score -= 2
    # return max_score_pos


def soft_clip_length(r, min_qual):
    """
    find out how many bases are soft clipped on each end
    """
    c_left, c_right = 0, 0
    if r.is_unmapped or r.query_qualities is None:
        return 0, 0
    # 5' end
    if r.cigartuples[0][0] == BAM_SOFT_CLIP:
        Q = r.query_qualities[:r.cigartuples[0][1]]
        qc_fail_bases = 0
        for q in r.query_qualities:
            if q < min_qual:
                qc_fail_bases += 1
            else:
                break
        c_left = r.cigartuples[0][1] - qc_fail_bases
        c_left = max(c_left, 0)
    # 3' end
    if r.cigartuples[-1][0] == BAM_SOFT_CLIP:
        qc_fail_bases = 0
        for q in r.query_qualities[::-1]:
            if q < min_qual:
                qc_fail_bases += 1
            else:
                break
        c_right = r.cigartuples[-1][1] - qc_fail_bases
        c_right = max(c_right, 0)
    return c_left, c_right

def right_tail(r, c_right):
    """
    Right clipped tail, left breakpoint
    """
    L = len(r.query_sequence)
    p1 = L - r.cigartuples[-1][1]
    p2 = p1 + c_right
    clip_seq = r.query_sequence[p1:p2]
    if not clip_seq:
        clip_seq = "-"
        clip_qual = "-"
    else:
        Q = r.query_qualities[p1:p2]
        clip_qual = "".join(chr(q + Q_OFFSET) for q in Q)
    return clip_seq, clip_qual

def left_tail(r, c_left):
    """
    Left clipped tail, right breakpoint
    """
    p2 = r.cigartuples[0][1]
    p1 = p2 - c_left
    clip_seq = reverse_complement(r.query_sequence[p1:p2])
    if not clip_seq:
        clip_seq = "-"
        clip_qual = "-"
    else:
        Q = r.query_qualities[p1:p2]
        clip_qual = "".join(chr(q + Q_OFFSET) for q in Q)
        clip_qual = clip_qual[::-1]
    return clip_seq, clip_qual


# process the input sam file
samfile = pysam.AlignmentFile(sys.stdin if args.aln_file == "-"
    else args.aln_file)
for r in samfile:
    if r.is_unmapped:
        continue
    if r.mapping_quality < args.min_mapq:
        continue
    if r.reference_end - r.reference_start < args.min_map_length:
        continue
    c_left, c_right = soft_clip_length(r, args.min_end_qual)
    clip_seq = "-"
    clip_qual = "-"
    if args.as_clipped:
        if not c_left and not c_right:
            # clipped only for quality
            continue
        # note that seq and qual in SAM/BAM file are on the same strand
        # as the reference genome
        bkp = "R" if c_right < c_left else "L"
        if c_right == c_left:
            # tie -- assume the 5' end is the major clipped end
            bkp = "L" if r.is_reverse else "R"
        if bkp == "R":
            assert r.cigartuples[0][0] == BAM_SOFT_CLIP
            clip_seq, clip_qual = left_tail(r, c_left)
        else:
            assert r.cigartuples[-1][0] == BAM_SOFT_CLIP
            clip_seq, clip_qual = right_tail(r, c_right)
    elif r.is_reverse:
        bkp = "R"
        if c_left > 0 and r.cigartuples[0][0] == BAM_SOFT_CLIP:
            clip_seq, clip_qual = left_tail(r, c_left)
    else:
        bkp = "L"
        if c_right > 0 and r.cigartuples[-1][0] == BAM_SOFT_CLIP:
            clip_seq, clip_qual = right_tail(r, c_right)
    htsr = HTSRead.from_aligned_segment(r, c_left, c_right, bkp,
        args.sample_name, clip_seq, clip_qual)
    sys.stdout.write(f"{htsr}\n")
