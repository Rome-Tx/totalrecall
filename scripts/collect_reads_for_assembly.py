#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import (HTSRead, NamedBreakpoint, IntervalWrapper,
    intersect_multiple_iters, template_length_stats)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Annotate breakpoints using discordant read info')
parser.add_argument('breakpoint_file',
    help='TSV file with named breakpoint info')
parser.add_argument('breakpoint', help='breakpoint to work on (L/R)')
parser.add_argument('disc_tsv_file',
    help='TSV file with discordant read info')
parser.add_argument('clip_tsv_file',
    help='TSV file with clipped read info')
parser.add_argument('tlen_file', help='file with template length info')
parser.add_argument('known_file', help='file with list of annotated'
    'insertions')
parser.add_argument('-c', '--clip-distance', type=int, default=5,
    help='allowed distance between clipping coordinates')
parser.add_argument('-Q', '--min-qual', type=int, default=1,
    help='minmal MAPQ for a discordant read to be used as anchor')


args = parser.parse_args()

obs_tlen, obs_sd, ext_len = template_length_stats(args.tlen_file)

assert args.breakpoint in ("L", "R")
left_ext = 0 if args.breakpoint == "R" else ext_len
right_ext = 0 if args.breakpoint == "L" else - ext_len

with open(args.known_file) as fh:
    know_id_list = {line.split(maxsplit=1)[0] for line in fh}

with open(args.breakpoint_file) as fid, open(
args.disc_tsv_file) as pe_fid, open(args.clip_tsv_file) as cr_fid:
    bp_iter = (IntervalWrapper(bp.ch, bp.pos, bp.pos, bp) for bp in 
        (NamedBreakpoint(line) for line in fid) if bp.bkp==args.breakpoint 
        and bp.name.strip("@LR") in know_id_list)
    pe_iter = (IntervalWrapper(r.rname, r.ref_locus + left_ext,
        r.ref_locus + right_ext, r) for r in (HTSRead.from_line(line)
        for line in pe_fid) if r.bkp==args.breakpoint)
    cr_iter = (IntervalWrapper(r.rname, r.ref_locus - args.clip_distance,
        r.ref_locus + args.clip_distance, r)  for r in (HTSRead.from_line(line)
        for line in cr_fid) if r.bkp==args.breakpoint)
    for bp, pe_reads, c_reads in intersect_multiple_iters(bp_iter, pe_iter,
    cr_iter):
        fq_seq = set()
        for r in pe_reads:
            fq_seq.add((r.seq, r.qual))
            fq_seq.add((r.mate_seq, r.mate_qual))
        for r in c_reads:
            if r.check_against_breakpoint(bp,
            allowed_dist = args.clip_distance):
                fq_seq.add((r.seq, r.qual))
        name = bp.name.strip("@LR")
        fq_list = (f"{seq}|{qual}" for seq, qual in fq_seq)
        fq_repr = " ".join(fq_list)
        if not fq_repr:
            fq_repr = "."
        sys.stdout.write(f"{name}\t{fq_repr}\n")
