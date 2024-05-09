#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict, Counter
from TRHelper.Aux import (HTSRead, NamedBreakpoint, IntervalWrapper,
    MateRepeatAlnInfo, intersect_iter, template_length_stats)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Annotate breakpoints using discordant read info')
parser.add_argument('breakpoint_file',
    help='TSV file with named breakpoint info')
parser.add_argument('breakpoint', help='breakpoint to work on (L/R)')
parser.add_argument('disc_tsv_file',
    help='TSV file with discordant read info')
parser.add_argument('tlen_file', help='file with template length info')
parser.add_argument('-b', '--bed-l1-fli',
    help='BED file with coordinates of FLI LINE1 elements')
parser.add_argument('-c', '--clip-distance', type=int, default=5,
    help='allowed distance between clipping coordinates')
parser.add_argument('-f', '--flanking-length', type=int, default=2000,
    help='maximal flanking length past the missed poly(A) signal '
    'for LINE readthrough transcripts')
parser.add_argument('-Q', '--min-qual', type=int, default=1,
    help='minmal MAPQ for a discordant read to be used as anchor')


args = parser.parse_args()

obs_tlen, obs_sd, ext_len = template_length_stats(args.tlen_file)

assert args.breakpoint in ("L", "R")
left_ext = 0 if args.breakpoint == "R" else - ext_len
right_ext = 0 if args.breakpoint == "L" else ext_len

# read intact LINE locations
# FLI = defaultdict(list)
# if args.bed_l1_fli:
    # with open(args.bed_l1_fli) as fid:
        # for line in fid:
            # S = line.split()
            # ch = S[0]
            # name = S[3]
            # strand = S[5]
            # if strand == "+":
                # start = int(S[2])
                # end = start + args.flanking_length
            # else:
                # end = int(S[1])
                # start = end - args.flanking_length
            # FLI[ch].append((start, end, name, strand))


# def check_downstream_line(r, FLI):
    # for start, end, name, strand in FLI[r.mate_rname]:
        # if r.mate_is_mapped:
            # if r.mate_is_rc == (strand == "-") and (start 
            # <= r.mate_ref_start <= end):
                # return name
    # return ""


with open(args.breakpoint_file) as fid, open(
args.disc_tsv_file) as pe_fid:
    bp_iter = (IntervalWrapper(bp.ch, bp.pos + left_ext, bp.pos + 
        right_ext, bp) for bp in (NamedBreakpoint(line) for line in 
        fid) if bp.bkp==args.breakpoint)
    pe_iter = (IntervalWrapper(r.rname, r.ref_locus, r.ref_locus, r) 
        for r in (HTSRead.from_line(line) for line in pe_fid) 
        if r.bkp==args.breakpoint)
    for bp, pe_reads in intersect_iter(bp_iter, pe_iter, one_by_one=False):
        hc = list()
        lc = list()
        n_unknown = 0
        for r in pe_reads:
            # check_downstream_line(r, FLI)
            if r.mate_aln_info == ".":
                n_unknown += 1
            else:
                for x in r.mate_aln_info.split(","):
                    mri = MateRepeatAlnInfo.from_short_repr(x)
                    mri.distance_to_breakpoint = (bp.pos - r.ref_start + 1
                        ) if args.breakpoint == "L" else (r.ref_end - bp.pos + 1)
                    mri.mate_seq = r.mate_seq
                    mri.mate_qual = r.mate_qual
                    if r.check_against_breakpoint(bp, args.clip_distance):
                        hc.append(mri)
                    else:
                        lc.append(mri)
        name = bp.name.strip("@LR")
        info_hc = "~".join(str(mri) for mri in hc) if hc else "."
        info_lc = "~".join(str(mri) for mri in lc) if lc else "."
        #downstream_line = Counter(filter(None, downstream_line))
        sys.stdout.write(f"{name}\t{n_unknown}\t{info_hc}\t{info_lc}\n")
