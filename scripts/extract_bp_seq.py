#!/usr/bin/env python3

import sys
import argparse
from Bio.Seq import reverse_complement
from TRHelper.Aux import NamedBreakpoint, read_gzipped_genome

EXT_LEN = 15

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='extract sequence near the breakpoint')
parser.add_argument('genome_gz_file', help='genome fa.gz file')
parser.add_argument('-l', '--ext-length', type=int, default=EXT_LEN,
    help='length of extension on each side of the breakpoint')
parser.add_argument('-n', '--nbp-file', default="case.ins.nbp.tsv",
    help='length of the sequence to extract')

args = parser.parse_args()

gh = read_gzipped_genome(args.genome_gz_file)

with open(args.nbp_file) as fid:
    for line in fid:
        nbp = NamedBreakpoint(line)
        try:
            gseq = gh[nbp.ch][nbp.pos-args.ext_length:nbp.pos+args.ext_length]
        except KeyError:
            continue
        if len(gseq) < 2 * args.ext_length:
            # breakpoint near the end of the contig, skip it
            continue
        if nbp.bkp == "R":
            gseq = reverse_complement(gseq)
            # tseq = reverse_complement(genome[nbp.ch][nbp.pos:nbp.pos+EXT_LEN]
            #                           + nbp.seq)
        # else:
        #     tseq = genome[nbp.ch][nbp.pos-EXT_LEN:nbp.pos] + nbp.seq
        sys.stdout.write(">{}\n{}\n".format(nbp.name, gseq.upper()))
