#!/usr/bin/env python3

import sys
import argparse
from TRHelper.Aux import DOC, TSDOC, TargetSite, IntervalWrapper, intersect_iter

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='TODO')
parser.add_argument("doc_file", help="DOC file")
parser.add_argument("-e", "--ext-length", type=int, default=100, 
                    help="extension length for intervals near insertions")
parser.add_argument("-q", "--lower-quantile", type=int, default=20, 
                    help="number of bases to compute lower coverage at")
parser.add_argument("-Q", "--upper-quantile", type=int, default=20, 
                    help="number of bases to compute upper coverage at")

args = parser.parse_args()

EL = args.ext_length

with open("case.genome") as fid:
    contig_len = {rid: int(rlen) for rid, rlen in (line.split() for line in fid)}

with open("case.ins.cand.tsv") as fid, open(args.doc_file) as doc_fid:
    ts_iter = (IntervalWrapper(ts.ch, ts.estart(EL), ts.eend(EL, contig_len),
        ts) for ts in (TargetSite.from_line(line) for line in fid))
    doc_iter = (IntervalWrapper(doc.ch, doc.pos, doc.pos, 
        doc) for doc in (DOC(line) for line in doc_fid))
    for ts, doc_list in intersect_iter(ts_iter, doc_iter,
    one_by_one=False, yield_all=True):
        # list of DOC values
        doc_d = {doc.pos:doc.depth for doc in doc_list}
        doc_values = [doc_d.get(pos, 0) for pos in range(ts.estart(EL), 
            ts.eend(EL, contig_len))]
        doc_values.sort()
        # quantiles
        L = ts.eend(EL, contig_len) - ts.estart(EL)
        median_q = L // 2
        low_q = min(median_q, args.lower_quantile)
        upper_q = min(L, args.upper_quantile)
        low_d = doc_values[low_q]
        upper_d = doc_values[-upper_q]
        median_d = doc_values[median_q]
        try:
            doc_ratio = (upper_d+low_d) / (median_d + low_d)
        except ZeroDivisionError:
            doc_ratio = -1
        # compute DOC at breakpoints
        d_lbp = doc_d.get(ts.lbp, 0)
        d_lts = doc_d.get(ts.lbp + 1, 0)
        d_rbp = doc_d.get(ts.rbp, 0)
        d_rts = doc_d.get(ts.rbp - 1, 0)
        tsdoc = TSDOC.from_values(ts.name, d_lbp, d_lts, d_rbp, d_rts,
            doc_ratio, low_d, median_d, upper_d)
        sys.stdout.write(f"{tsdoc}\n")
