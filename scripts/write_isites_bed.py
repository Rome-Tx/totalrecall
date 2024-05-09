#!/usr/bin/env python3

import sys
import gzip
import argparse

from TRHelper.Aux import TargetSite, VCF, template_length_stats

ADDITIONAL_BASES = 5

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Write BED file from list of insertions')
parser.add_argument('genome', help='.genome file with contig lengths')
parser.add_argument('infile', help='input TSV file')
parser.add_argument("-e", "--ext_length", type=int, default=0, 
                    help="extension length for BED intervals")
parser.add_argument("-t", "--template-length",
                    help="file with template length stats")
parser.add_argument("-v", "--vcf-gz", action='store_true', 
                    help="input is gzipped VCF as opposed to tab delimited")


args = parser.parse_args()

with open(args.genome) as fid:
    contig_len = {rid: int(rlen) for rid, rlen in (line.split() 
       for line in fid)}

if args.template_length:
    obs_tlen, obs_sd, ext = template_length_stats(args.template_length)
else:
    ext = args.ext_length
ext += ADDITIONAL_BASES

fh = gzip.open(args.infile, "rt") if args.vcf_gz else open(args.infile)
for line in fh:
    if args.vcf_gz:
        if line.startswith("#"):
            continue
        vc = VCF.from_line(line)
        ch = vc.ch
        start = max(1, vc.pos - ext)
        end = min(contig_len[ch], vc.info.END + ext) + 1
        nm = "{}.{}.{}".format(vc.id, vc.info.TE,
            vc.filters_repr.replace(",", "|"))

    else:
        ts = TargetSite.from_line(line)
        ch = ts.ch
        start = ts.estart(ext)
        end = ts.eend(ext, contig_len)
        nm = ts.name
    sys.stdout.write(f"{ch}\t{start}\t{end-1}\t{nm}\n")
