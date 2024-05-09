"""
Alignment wrapper
"""

import sys
from math import log
from array import array

from Bio.Seq import reverse_complement

from TRHelper.Constants import (Q_OFFSET, MAX_Q, POLYA_5TRIM, POLYA_R, POLYA_Q,
    POLYA_FULL_Q)

from TRHelper import _TRHelper

PIDENT = 0.8
LEN_CUTOFF = 16
LEAD_GAPS = 2

def polya_score(seq, trim=POLYA_5TRIM, s_match=POLYA_R, s_mism=POLYA_Q):
    byte_seq = bytes(seq.upper(), "ASCII")
    return _TRHelper.polya_score(byte_seq, trim, s_match, s_mism)

def polya_score_qual(seq, qual, Smatch, Smism, trim=POLYA_5TRIM,
q_offset=Q_OFFSET):
    byte_seq, byte_qual = bytes(seq.upper(), "ASCII"), bytes(qual, "ASCII")
    return _TRHelper.polya_score_qual(byte_seq, byte_qual, Smatch,
        Smism, trim, q_offset)

def polya_score_full_length(seq, s_match=POLYA_R, s_mism=POLYA_FULL_Q):
    byte_seq = bytes(seq.upper(), "ASCII")
    return _TRHelper.polya_score_full_length(byte_seq, s_match, s_mism)

def polya_score_full_length_qual(seq, qual, Smatch, Smism, q_offset=Q_OFFSET):
    byte_seq, byte_qual = bytes(seq.upper(), "ASCII"), bytes(qual, "ASCII")
    return _TRHelper.polya_score_full_length_qual(byte_seq, byte_qual, Smatch,
        Smism, q_offset)


def compute_clipping_from_qual(bqual, len_L, len_R, min_qual=3):
    clip_L, clip_R = _TRHelper.compute_clipping_from_qual(bqual, min_qual,
        len_L, len_R)
    return clip_L, clip_R


def compute_polya_trim(rec):
    """
    compute the length of the poly(A) at the 3' end of the reference sequence
    to be trimmed
    """
    rid = rec.id.upper()
    transposon = rid.rsplit(sep = ":", maxsplit=1)[-1].upper()
    if transposon == "LTR":
        return 0
    elif transposon in ("SINE1", "SVA", "LINE1", "LINE1INV"):
        seq = str(rec.seq).upper()
        if transposon != "LINE1INV":
            seq = reverse_complement(seq)
        score = 0
        max_score = 0
        max_score_pos = 0
        for pos, ch in enumerate(seq, 1):
            if ch == "T":
                score += 1
                if score >= max_score:
                    max_score = score
                    max_score_pos = pos
            else:
                score -= 2
        return max_score_pos
    return -1


def cluster_clipped_tails(read_list, pident=PIDENT, len_cutoff=LEN_CUTOFF,
lead_gaps=LEAD_GAPS, max_cluster_input_len=None):
    """
    CD-HIT like clustering
    """
    # TODO: resolve equal lengths using base qualities
    read_list_sorted = sorted(read_list, key=lambda r:len(r.cseq)) # ascending
    if (max_cluster_input_len is not None and len(read_list_sorted)
    > max_cluster_input_len):
        read_list_sorted = read_list_sorted[-max_cluster_input_len:]
    repr_list = []
    read_cluster_list = []
    while read_list_sorted:
        r = read_list_sorted.pop()  # pop the last element
        byte_cseq = bytes(r.cseq, "ASCII")
        for repr_seq, cluster in zip(repr_list, read_cluster_list):
            if _TRHelper.check_identity(repr_seq, byte_cseq, pident, len_cutoff, 
            lead_gaps):
                cluster.append(r)
                break
        else:
            repr_list.append(byte_cseq)
            read_cluster_list.append([r])
    return read_cluster_list


def check_identity(seq1, seq2, pident=PIDENT, len_cutoff=LEN_CUTOFF,
lead_gaps=LEAD_GAPS):
    byte_seq1 = bytes(seq1, "ASCII")
    byte_seq2 = bytes(seq2, "ASCII")
    return _TRHelper.check_identity(byte_seq1, byte_seq2, pident, len_cutoff,
        lead_gaps)


def check_vs_normal(cseq, normal_tails, pident=PIDENT, len_cutoff=LEN_CUTOFF,
lead_gaps=LEAD_GAPS):
    byte_cseq = bytes(cseq, "ASCII")
    for normal_seq in normal_tails:
        byte_ncseq = bytes(normal_seq, "ASCII")
        if _TRHelper.check_identity(byte_cseq, byte_ncseq, pident, len_cutoff,  
            lead_gaps):
            return True
    return False


def calc_num_clusters(ref, read_list, pident=PIDENT, len_cutoff=LEN_CUTOFF,
lead_gaps=LEAD_GAPS, num_clusters_cutoff=None):
    """
    ref = clipped seqeunce at the breakpoint
    read_list = reads to cluster first against ref and then the other reads
    """
    read_list_sorted = sorted(read_list, key=lambda r:len(r), reverse=True)
    byte_ref_list = [bytes(ref, "ASCII")]
    for r in read_list_sorted:
        br = bytes(r, "ASCII")
        for bref in byte_ref_list:
            if _TRHelper.check_identity(bref, br, pident, len_cutoff, lead_gaps):
                break
        else:
            byte_ref_list.append(br)
            if (num_clusters_cutoff is not None and len(byte_ref_list)
            >= num_clusters_cutoff):
                break
    return len(byte_ref_list)


def _f(p, a, b):
    """
    auxiliary function to equate to 0
    """
    return a * log(4. / 3 * (1 - p)) + b * log(4 * p)

def _calc_p(a, b, tol=1e-4):
    """
    solves for p equating _f(p, a, b) to zero
    """
    p0 = b / (a + b)
    assert _f(p0, a, b) > 0
    while _f(p0, a, b) > 0:
        p1 = (1 + p0) / 2
        if _f(p1, a, b) < 0:
            break
        else:
            p0 = p1
    while p1 - p0 > tol:
        p = (p0 + p1) / 2
        fn = _f(p, a, b)
        if fn == 0:
            return p
        elif fn > 0:
            p0 = p
        else:
            p1 = p
    return (p0 + p1) / 2

def calc_p_lambda(a, b, bitscore_scale=False):
    p = _calc_p(a, b)
    L = log(4 * p) / a
    if bitscore_scale:
        # use bitscore scale
        return p, L / log(2)
    return p, L

def _calc_score_from_q(pident, L, Q):
    perr = 10**(-Q / 10)
    Rmatch = 4 * pident
    Rmism = 4 * (1 - pident) / 3
    Rmatch_Q = (1 - perr) * Rmatch + perr * Rmism
    Rmism_Q = perr * Rmatch + (1 + perr) * Rmism
    s_match_Q = log(Rmatch_Q) / L
    s_mism_Q = - log(Rmism_Q) / L
    return s_match_Q, s_mism_Q

def compute_qual_scores(s_match, s_mism, max_Q=MAX_Q):
    p, L = calc_p_lambda(s_match, s_mism, bitscore_scale=False)
    Smatch, Smism = zip(*(_calc_score_from_q(p, L, Q) for Q in range(max_Q + 1)))
    Smatch = [round(s) for s in Smatch]
    Smism = [round(s) for s in Smism]
    Smatch = bytes(array("h", Smatch))
    Smism = bytes(array("h", Smism))
    return Smatch, Smism
