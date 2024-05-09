#!/usr/bin/env python3

import sys
import os
import gzip
import numpy as np
from sklearn.naive_bayes import GaussianNB
from TRHelper import Aux

LINE3END = 6010  # rightmost allowed alignment start TODO: embed in aln info

class RepeatAlignInfo(Aux.MateRepeatAlnInfo):
    def __init__(self):
        pass

    @classmethod
    def from_line(self, line):
        if line == ".":
            return None
        rai = super().from_long_repr(line)
        return rai
    
    def check_three(self, max_fragment_length=10000, coord=-1, 
    aligned_fraction = 0.2, inversion=False):
        if not inversion and self.strand == "-":
            return False
        if coord == -1:
            coord = self.rlen
        if self.transposon == "LINE1" and self.rstart > LINE3END:
            return False
        if (coord - self.rstart + self.distance_to_breakpoint
        ) > max_fragment_length:
            return False
        return (self.rend - self.rstart) >= aligned_fraction * self.qlen
    
    def check_five(self, max_fragment_length=10000, coord=-1,
    aligned_fraction = 0.7, inversion=False):
        if not inversion and self.strand == "+":
            return False
        if coord == -1:
            coord = 0
        if self.rstart < coord:
            return False
        if (self.rend - coord + self.distance_to_breakpoint
        ) > max_fragment_length:
            return False
        return (self.rend - self.rstart) >= aligned_fraction * self.qlen

class RepeatSummary:
    def __init__(self, X):
        self.n_unknown = int(X[0])
        self.hc = [RepeatAlignInfo.from_line(x) for x in X[1].split("~")
            if X[1] != "."]
        self.lc = [RepeatAlignInfo.from_line(x) for x in X[2].split("~")
            if X[2] != "."]
        self.N = self.n_unknown + len(self.hc) + len(self.lc)

    def count_transposon_reads(self, transposon):
        tl = self.hc + self.lc
        num_match = sum(r.transposon.upper() == transposon.upper() for r in tl)
        return num_match

def read_vcf(fname, info_L, info_R):
    vcf_h = dict()
    with open(fname) as fh:
        for line in fh:
            S = line.split(maxsplit=2)
            cmp = S[0]
            old_id = S[1]
            v = Aux.VCF.from_line(S[2])
            v.lc = 1 if v.filters[0] != "PASS" else 0
            v.cmp = cmp
            v.peL = info_L.get((v.cmp, v.id), None)
            v.peR = info_R.get((v.cmp, v.id), None)
            vcf_h[(cmp, old_id)] = v
    return vcf_h

def pe_repr(v, info, five=True):
    """
    numeric representation of PE annotation
    """
    if info is None:
        return [0] * 5
    nhc = 0
    nlc = 0
    for pea in info.hc:
        nhc += (pea.transposon == "LINE1" and (
            pea.check_five() if five else pea.check_three()))
    for pea in info.lc:
        nlc += (pea.transposon == "LINE1" and (
            pea.check_five() if five else pea.check_three()))
    return [info.N, nhc, nlc, nhc/v.info.CASE_DOC_MEDIAN,
        nlc/v.info.CASE_DOC_MEDIAN]
    
def vcf2num(v):
    vi = v.info
    va5 = vi.BP5ANNOT
    va3 = vi.BP3ANNOT
    CC5 = int(vi.CASE_CC5.split(":", maxsplit=1)[0])
    CC3 = int(vi.CASE_CC3.split(":", maxsplit=1)[0])
    X = [int(vi.INVERSION), vi.SVLEN, vi.MQ5, vi.MQ3, vi.BP5NR, vi.BP3NR,
        vi.DEPTH_5UP, vi.DEPTH_3UP, vi.BP5NR/vi.DEPTH_5UP,
        vi.BP3NR/vi.DEPTH_3UP, vi.CN_RATIO, vi.PVALUE_CNV, CC5, CC3,
        vi.CASE_DOC_MEDIAN, vi.CASE_DOC_RATIO, va5.qstart,
        len(vi.BP5SEQ) - va5.qend, len(vi.BP5SEQ), va3.qstart, 
        len(vi.BP3SEQ) - va3.qend, len(vi.BP3SEQ), va5.pident, va5.bitscore,
        va3.pident, va3.bitscore, v.lc]
    X = [int(vi.INVERSION), vi.SVLEN, vi.MQ5, vi.MQ3, vi.BP5NR, # 0-4
        vi.BP3NR, vi.DEPTH_5UP, vi.DEPTH_3UP, vi.BP5NR/vi.DEPTH_5UP, # 5-8
        vi.BP3NR/vi.DEPTH_3UP, vi.CN_RATIO, vi.PVALUE_CNV, CC5, CC3, # 9-13
        vi.CASE_DOC_MEDIAN, vi.CASE_DOC_RATIO, va5.qstart, # 14-16
        len(vi.BP5SEQ) - va5.qend, len(vi.BP5SEQ), va3.qstart, # 17-19
        len(vi.BP3SEQ) - va3.qend, len(vi.BP3SEQ), va5.pident, # 20-22
        va5.bitscore, va3.pident, va3.bitscore, v.lc] # 23-26
    I5 = pe_repr(v, v.peL if v.info.STRAND=="+" else v.peR, five=True) # 27-31
    I3 = pe_repr(v, v.peR if v.info.STRAND=="+" else v.peL, five=False) # 32-36
    return X + I5 + I3

cand_vars = [23, 25, 16, 17, 4, 26, 30, 31, 35, 36]

training_dir = sys.argv[1]
sl_path = os.path.join(training_dir, "training_classification.txt")
x_path = os.path.join(training_dir, "training_data.txt")
sl = [line.strip() for line in open(sl_path)]
x = np.loadtxt(x_path)

# bayesian classifier
gnb = GaussianNB()
pred = gnb.fit(x, sl).predict(x)

# results
ofh = open("line.vcf", "w")
vl = []
for line in gzip.open("somatic.vcf.gz", "rt"):
    if line.startswith("#"):
        ofh.write(line)
        continue
    v = Aux.VCF.from_line(line)
    if v.alt == "<INS:ME:LINE1>" and v.info.CASE_DOC_MEDIAN:
        vl.append(v)
id_set = {v.id for v in vl}

def PEInfoReader(fname, id_set):
    rs = dict()
    for line in open(fname):
        S = line.split()
        bp = S[0]
        if bp in id_set:
            rs[bp] = RepeatSummary(S[1:])
    return rs

info_L = PEInfoReader("aux/case.pe.L.tsv", id_set)
info_R = PEInfoReader("aux/case.pe.R.tsv", id_set)
for v in vl:
    v.peL = info_L.get(v.id, None)
    v.peR = info_R.get(v.id, None)
    v.lc = 1 if v.filters[0] != "PASS" else 0

X2 = np.array([vcf2num(v) for v in vl])
x2 = X2[:, cand_vars]
pred = gnb.predict(x2)

for v, flag in zip(vl, pred):
    if flag == "OK":
        ofh.write(f"{v}\n")
