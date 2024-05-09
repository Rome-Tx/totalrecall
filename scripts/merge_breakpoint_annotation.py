#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import Breakpoint, BPSAnnotation
from TRHelper.Align import polya_score, polya_score_full_length, calc_p_lambda
from TRHelper.Constants import (POLYA_R, POLYA_Q, POLYA_FULL_Q, POLYA_MINSCORE,
    POLYA_MINSCORE_FULL, POLYA_5TRIM)


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Auxiliary script to merge info regarding mapping to repeats '
		    'for clipped sequences at the breakpoints')
parser.add_argument('byid_file', help='breakpoints sorted by integer id')
parser.add_argument("annot_file", help="best mapping for each breakpoint "
    "sequence, sorted by integer id")

args = parser.parse_args()

#Smatch, Smism = compute_qual_scores(POLYA_R, POLYA_Q, MAX_Q)
#Smatch_FL, Smism_FL = compute_qual_scores(POLYA_R, POLYA_FULL_Q, MAX_Q)
p, L = calc_p_lambda(POLYA_R, POLYA_Q, bitscore_scale=True)

with open(args.byid_file) as fh_bp, open(args.annot_file) as fh_annot:
    bp = None
    for line in fh_annot:
        bts = BPSAnnotation.from_text_repr(line)
        while bp is None or bp.id < bts.numeric_id:
            try:
                line_bp = next(fh_bp)
            except StopIteration:
                raise ValueError("breakpoint/annotation files "
                    "not sorted numerically?")
            bp = Breakpoint.from_line(line_bp)
            score_polya = polya_score(bp.seq, POLYA_5TRIM, POLYA_R, POLYA_Q)
            bitscore_polya = L * score_polya
            if polya_score_full_length(bp.seq, POLYA_R, 
            POLYA_FULL_Q) >= POLYA_MINSCORE_FULL:
                # full length poly(A):
                bp.annotation = BPSAnnotation.dummy_polya_annotation(bp.id,
                    bitscore_polya)
            elif score_polya >= POLYA_MINSCORE:
                dpa = BPSAnnotation.dummy_polya_annotation(bp.id,
                    bitscore_polya)
                if bp.id == bts.numeric_id:
                    bp.annotation = (bts if bts.bitscore >= bitscore_polya
                        else dpa)
                else:
                    bp.annotation = dpa
            elif bp.id == bts.numeric_id:
                    bp.annotation = bts
            sys.stdout.write(f"{bp}\n")
    # now handle the remaining breakpoints which have no annotation
    while True:
        try:
            line_bp = next(fh_bp)
        except StopIteration:
            # end of the file
            break
        bp = Breakpoint.from_line(line_bp)
        sys.stdout.write(f"{bp}\n")
