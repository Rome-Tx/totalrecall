#!/usr/bin/env python3

import sys
import argparse
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Assemble left/right/global contings')
parser.add_argument("known_tsv", help='File with known transductions')
parser.add_argument("left_reads", help='Fasta file with reads '
    'for the left breakpoint ')
parser.add_argument("right_reads", help='Fasta file with reads '
    'for the right breakpoint ')

args = parser.parse_args()

def parse_line_with_reads(line):
    S = line.split()
    name = S[0]
    fq_set = set()
    for r in S[1:]:
        if r != ".":
            seq, qual = r.split("|")
            fq_set.add((seq, qual))
    return name, fq_set

def assemble_contigs(read_set):
    if not read_set:
        return "."
    fq_string = "".join(f"@{n}\n{seq}\n+\n{qual}\n" for n, (seq, qual)
        in enumerate(read_set))
    child = subprocess.Popen(["fml-asm", "-t", "1", "-"],
        stdin=subprocess.PIPE, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    child.stdin.write(fq_string)
    child.stdin.close()
    ctg_list = list(SeqIO.parse(child.stdout, "fastq"))
    # TODO: report qualities as well
    contigs = "|".join(str(rec.seq) for rec in ctg_list) if ctg_list else "."
    return contigs

with open(args.known_tsv) as fknown, open(args.left_reads) as fl, open(
args.right_reads) as fr:
    for line_known, line_l, line_r in zip(fknown, fl, fr):
        name_known, info_known = line_known.split(maxsplit=1)
        info_known, filter_known = info_known.rsplit(maxsplit=1)
        name_l, r_l = parse_line_with_reads(line_l)
        name_r, r_r = parse_line_with_reads(line_r)
        assert name_known == name_l == name_r
        ctg_l = assemble_contigs(r_l)
        ctg_r = assemble_contigs(r_r)
        ctg = assemble_contigs(r_l | r_r)
        sys.stdout.write(f"{name_known}\t{info_known}\t{ctg_l}\t{ctg_r}\t"
            f"{ctg}\t{filter_known}\n")
