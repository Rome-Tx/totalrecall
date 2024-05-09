#!/usr/bin/python3

import sys
from operator import attrgetter
import argparse
from Bio import SeqIO
from TRHelper.Aux import BlastTab, BPSAnnotation
from TRHelper.Align import compute_polya_trim

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Auxiliary script to assemble/convert annotation format')
parser.add_argument('fasta_file', help='reference fasta file')

args = parser.parse_args()


MIN_PIDENT = 80

rl = list(SeqIO.parse(args.fasta_file, "fasta"))
ref_length_h = {r.id : len(r.seq) for r in rl}
ref_polya_h = {r.id : compute_polya_trim(r) for r in rl}

def alignment_iter():
    aln_list = []
    qid = None
    for line in sys.stdin:
        aln = BlastTab(line)
        if aln.qname == qid:
            aln_list.append(aln)
            continue
        if aln_list:
            yield aln_list
        qid = aln.qname
        aln_list = [aln]
    if aln_list:
        yield aln_list

for aln_list in alignment_iter():
    aln_list = list(filter(lambda aln: aln.pident >= MIN_PIDENT, aln_list))
    if not aln_list:
        continue
    aln_list.sort(key=attrgetter("rstart"))  # prefer leftmost alignment
    aln_list.sort(key=attrgetter("bitscore"))
    aln = aln_list.pop()
    bts = BPSAnnotation.from_blast_tab(aln, ref_length=ref_length_h[aln.rname],
        ref_polya=ref_polya_h[aln.rname])
    sys.stdout.write(f"{bts}\n")
