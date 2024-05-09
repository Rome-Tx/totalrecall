#!/usr/bin/env python3

import sys, os
import gzip
from operator import attrgetter
from random import sample
import argparse

from TRHelper.Aux import VCF, read_genome_sequence_lengths
from TRHelper.Constants import MAX_SCREENSHOTS

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Write IGV batch file for taking screenshots of insertions')

parser.add_argument('vcf_gz', help='input VCF file')
parser.add_argument('genome', help='.genome file with contig lengths')
parser.add_argument("reference_fa_gz", help="reference genome")
parser.add_argument("-e", "--ext_length", type=int, default=100,
    help="extension length for intervals near insertions")

args = parser.parse_args()

with open(args.genome) as fid:
    contig_len = {rid: int(rlen) for rid, rlen in (line.split() for line in fid)}

glh = read_genome_sequence_lengths(args.reference_fa_gz)

ext = args.ext_length

try:
    os.mkdir("screenshots")
except FileExistsError:
    pass

sys.stdout.write("new\n")
sys.stdout.write(f"genome {args.reference_fa_gz}\n")
sys.stdout.write("preference SAM.FILTER_SECONDARY_ALIGNMENTS true\n")
sys.stdout.write("preference SAM.FILTER_SUPPLEMENTARY_ALIGNMENTS false\n")
sys.stdout.write("preference SAM.SHOW_CENTER_LINE false\n")
sys.stdout.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
sys.stdout.write("maxPanelHeight 10000\n")
sys.stdout.write("load case.bam name=CASE\n")
sys.stdout.write("load control.bam name=CONTROL\n")
sys.stdout.write("load is.annot.bed name=BREAKPOINTS\n")
sys.stdout.write("collapse\n")
sys.stdout.write("snapshotDirectory screenshots\n")

with gzip.open(args.vcf_gz, "rt") as fid:
    results = [VCF.from_line(line) for line in fid if not line.startswith("#")]

results = [vc for vc in results if vc.ch in glh] # remove variants on missing contigs

if len(results) > MAX_SCREENSHOTS:
    results = [vc for vc in results if vc.info.TE == "LINE1"]
    if len(results) > MAX_SCREENSHOTS:
        results = sample(results, MAX_SCREENSHOTS)

results.sort(key=attrgetter("pos"))
results.sort(key=attrgetter("ch"))

for vc in results:
    nm = "{}.{}".format(vc.id, vc.info.TE)
    ch = vc.ch
    start, end = sorted((vc.info.BP5, vc.info.BP3))
    start = max(0, start - ext)
    end = min(contig_len[ch], end + ext)
    sys.stdout.write("goto {} {} {}\n".format(ch, start, end))
    sys.stdout.write("maxPanelHeight 10000\n")
    sys.stdout.write("collapse\n")
    sys.stdout.write("snapshot {}.png\n".format(nm))
    sys.stdout.write("squish\n")
    sys.stdout.write("snapshot {}.svg\n".format(nm))

sys.stdout.write("exit\n")
