#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
from TRHelper.Aux import HTSRead, MateRepeatAlnInfo

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Auxiliary script to merge info regarding mapping to repeats for mates '
		    'in discordant read pairs')
parser.add_argument('annot_file', help='annotated reads mapping to repeats, '
                    'sorted by read name')
parser.add_argument("tsv_file", help="TSV file with discordant reads, "
                    "sorted by read name")
parser.add_argument("fq_file", help="FASTQ files with discordant reads, "
                    "sorted by read name")

args = parser.parse_args()

def parse_interleaved_fastq(fname):
    """
    Amended FASTQ parser, trims "/1" and "/2" from the read id
    and adds the "pe" attribute indicating whether it is R1 or R2
    """
    read_list = list()
    for r in SeqIO.parse(fname, "fastq"):
        rid, pe = r.id.rsplit("/", maxsplit=1)
        r.id = rid
        r.pe = int(pe)
        if not read_list or read_list[-1].id == r.id:
            read_list.append(r)
        else:
            read_list.sort(key = lambda x:x.pe)
            yield read_list[0].id, read_list[0], read_list[1]
            read_list = [r]
    if read_list:
        yield read_list[0].id, read_list[0], read_list[1]


class RepeatAnnotation:
    """
    Annotation for mates of discordant reads
    """
    def __init__(self, name):
        """
        R1, R2 - annotations for mates of R1 and R2
        """
        self.name = name
        self.A1 = list()
        self.A2 = list()
        self.R1 = ""
        self.Q1 = ""
        self.R2 = ""
        self.Q2 = ""

    def update_info(self, r):
        if r.pe == 1:
            self.A1.append(r)
        else:
            self.A2.append(r)
    
    def add_raw_sequence(self, r):
        if r.id != self.name:
            raise ValueError(f"mismatched id's: {self.name} vs {r.id}")
        seq = str(r.seq)
        qual = "".join(chr(q+33) for q in r.letter_annotations["phred_quality"])
        if r.pe == 1:
            self.R1, self.Q1 = seq, qual
        else:
            self.R2, self.Q2 = seq, qual

    def get_info(self, r):
        info = list()
        if r.pe == 1:
            mate_aln_info = [str(r) for r in self.A2]
            mate_seq = self.R2
            mate_qual = self.Q2
        else:
            mate_aln_info = [str(r) for r in self.A1]
            mate_seq = self.R1
            mate_qual = self.Q1
        mate_aln_info = ",".join(mate_aln_info) if mate_aln_info else "."
        return mate_aln_info, mate_seq, mate_qual


def repeat_annot_generator(fname):
    """
    Group reads with the same name and get the repeat annotation for R1/R2
    Return:
        read name, repeat for R1, repeat for R2
    """
    fq_iter = parse_interleaved_fastq(args.fq_file)
    try:
        fq_id, R1, R2 = next(fq_iter)
    except StopIteration:
        # empty FASTQ file
        return
    ra = None  # repeat annotation
    fh = open(fname)
    for line in fh:
        r = MateRepeatAlnInfo.from_LAST(line)
        if ra is None:
            ra = RepeatAnnotation(r.qname)
        if ra.name == r.qname:
            ra.update_info(r)
        else:
            while fq_id < ra.name:
                fq_id, R1, R2 = next(fq_iter)
            ra.add_raw_sequence(R1)
            ra.add_raw_sequence(R2)
            yield ra
            ra = RepeatAnnotation(r.qname)
            ra.update_info(r)
    if ra is not None:
        while fq_id < ra.name:
            fq_id, R1, R2 = next(fq_iter)
        ra.add_raw_sequence(R1)
        ra.add_raw_sequence(R2)
        yield ra


# both files need to be pre-sorted by the read name
# this is not checked here
fh1 = sys.stdin if args.tsv_file == "-" else open(args.tsv_file)
rai = repeat_annot_generator(args.annot_file)
try:
    ra = next(rai)
except StopIteration:
    # empty annotation file
    sys.exit(0)
for line in fh1:
    r = HTSRead.from_line(line)
    if ra.name > r.name:
        continue
    while ra.name < r.name:
        try:
            ra = next(rai)
        except StopIteration:
            break
    else:
        if ra.name==r.name:
            mate_info = ra.get_info(r)
            r.set_repeat_info(*mate_info)
            sys.stdout.write(f"{r}\n")
        continue
    # we have reached the end of the annotation file
    # and the while loop was ended by the break statement
    break
