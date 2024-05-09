#!/usr/bin/env python3

import sys
import os
import argparse
from operator import attrgetter
from collections import defaultdict
import subprocess

import pysam
from Bio import SeqIO

from TRHelper.Aux import (template_length_stats, TEI, VCF, MateRepeatAlnInfo,
    BPSAnnotation, Evidence, read_gzipped_genome)
from TRHelper.Constants import NCC_CUTOFF, LTR_TRANSPOSON_LEN

HC_CUTOFF = 0.5
DISC_READ_FRAC = 0.5

LINE3END = 6010  # rightmost allowed alignment start TODO: embed in aln info

FULL_LENGTH_OFFSET = 50

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Annotate and filter insertion sites')
parser.add_argument('genome_gz_file', help='genome fa.gz file')
parser.add_argument('-c', '--min-coverage-fraction', type=float, default=0.05,
    help='minmal number of clipped reads at the breakpoint as a fraction '
    'of the depth of coverage')


args = parser.parse_args()

gh = read_gzipped_genome(args.genome_gz_file)

sys.stdout.write(VCF.header("case.header.sam"))

contig_order = defaultdict(int)
with pysam.AlignmentFile("case.header.sam", "r") as af:
    for n, r in enumerate(af.references, 1):
        contig_order[r] = n


class RepeatAlignInfo(MateRepeatAlnInfo):
    def __init__(self):
        pass

    @classmethod
    def from_line(self, line):
        if line == ".":
            return None
        rai = super().from_long_repr(line)
        return rai
    
    def check_three(self, max_fragment_length, coord=-1, 
    aligned_fraction = 0.7, inversion=False):
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
    
    def check_five(self, max_fragment_length, coord=-1,
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
    def __init__(self, X, ext_len):
        self.max_fragment_length = ext_len
        self.n_unknown = int(X[0])
        self.hc = [RepeatAlignInfo.from_line(x) for x in X[1].split("~")
            if X[1] != "."]
        self.lc = [RepeatAlignInfo.from_line(x) for x in X[2].split("~")
            if X[2] != "."]
        self.N = self.n_unknown + len(self.hc) + len(self.lc)

    def filter_repeat_info(self, prime_end, coord=-1, aligned_fraction=0.7,
       inversion=False):
        if prime_end == 3:
            self.hc = [rei for rei in self.hc if rei.check_three(
                self.max_fragment_length, coord, aligned_fraction, inversion)]
            self.lc = [rei for rei in self.lc if rei.check_three(
                self.max_fragment_length, coord, aligned_fraction, inversion)]
        elif prime_end == 5:
            self.hc = [rei for rei in self.hc if rei.check_five(
                self.max_fragment_length, coord, aligned_fraction, inversion)]
            self.lc = [rei for rei in self.lc if rei.check_five(
                self.max_fragment_length, coord, aligned_fraction, inversion)]
        else:
            raise ValueError("Only 3'/5' is valid for the transposon end")

    def check_transposon_evidence(self, transposon):
        def count_from_list(ra_list):
            num_match = sum(r.transposon.upper() == transposon.upper() 
                for r in ra_list)
            return num_match, len(ra_list)
        # end aux funciton
        num_match, num_all = count_from_list(self.hc)
        if num_all > 0 and num_match > DISC_READ_FRAC * num_all:
            return True
        num_match, num_all = count_from_list(self.hc + self.lc)
        return num_match > DISC_READ_FRAC * num_all

    def count_transposon_reads(self, transposon):
        tl = self.hc + self.lc
        num_match = sum(r.transposon.upper() == transposon.upper() for r in tl)
        return num_match

    def assemble_contigs(self):
        return "."
        mate_list = self.hc + self.lc
        if not mate_list:
            return "."
        read_set = set((rai.mate_seq, rai.mate_qual) for rai in mate_list)
        fq_string = "".join(f"@{n}\n{seq}\n+\n{qual}\n" for n, (seq, qual)
            in enumerate(read_set))
        child = subprocess.Popen(["fml-asm", "-t", "1", "-"],
            stdin=subprocess.PIPE, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        child.stdin.write(fq_string)
        child.stdin.close()
        ctg_list = list(SeqIO.parse(child.stdout, "fastq"))
        contigs = "|".join(str(rec.seq) for rec in ctg_list) if ctg_list else "."
        return contigs


class PEInfoReader:
    def __init__(self, fname, tlen_fname):
        """
        fname - annotation file
        tlen_fname - file with template length stats
        """
        self.obs_tlen, self.obs_sd, self.ext_len = template_length_stats(
            tlen_fname)
        self.fh = open(fname)
        self.bp = None
        try:
            line = next(self.fh)
        except StopIteration:
            pass
        else:
            S = line.split()
            self.bp = S[0]
            rs = RepeatSummary(S[1:], self.ext_len)
            #rs.filter_repeat_info()
            self.rs = rs
    
    def __call__(self, bp_name):
        """
        bp - breakpoint name
        Subsequent calls should be made with ordered bp 
        arguments to allow the self.fh iterator to work correctly
        """
        while self.bp is not None and self.bp < bp_name:
            try:
                line = next(self.fh)
            except StopIteration:
                self.bp = None
            else:
                S = line.split()
                self.bp = S[0]
                rs = RepeatSummary(S[1:], self.ext_len)
                #rs.filter_repeat_info()
                self.rs = rs
        if self.bp is None or self.bp > bp_name:
            return RepeatSummary(["0", ".", "."], self.ext_len)
        else:
            return self.rs


with open("case.nbases") as fid:
    n_bases_case = int(fid.readline())
with open("control.nbases") as fid:
    try:
        n_bases_ctrl = int(fid.readline())
    except ValueError:
        n_bases_ctrl = 0
ctrl_present = bool(n_bases_ctrl)

with open("aux/control.clipped.read") as fid:
    filter_read = set(line.strip() for line in fid)

with open("aux/case.low_complexity") as fid:
    filter_lc = set(line.strip() for line in fid)

pe_left = PEInfoReader("aux/case.pe.L.tsv", "case.tlen")
pe_right = PEInfoReader("aux/case.pe.R.tsv", "case.tlen")

vc_list = []
with open("case.ins.namesorted.tsv") as fid, open("aux/doc_case.tsv"
    ) as fid_case, open("aux/doc_control.tsv") as fid_ctrl, open(
    "aux/case.read.cluster") as fid_cluster, open("aux/control.read.cluster"
    ) as fid_case_cluster:
    for line in fid:
        tei = TEI.from_line(line, fid_case.readline(), fid_ctrl.readline()
            if ctrl_present else None, fid_cluster.readline(),
            fid_cluster.readline(), fid_case_cluster.readline(),
            fid_case_cluster.readline())
        tei.compute_p_cnv(n_bases_case, n_bases_ctrl)
        case_info = "{}\t{}".format(tei.case_doc.d_median, tei.case_doc.ratio)
        ctrl_info = ("{}\t{}".format(tei.ctrl_doc.d_median, tei.ctrl_doc.ratio)
            if tei.ctrl_doc else "-1\t-1")
        ca_left = BPSAnnotation.from_text_repr(tei.lannot)
        ca_right = BPSAnnotation.from_text_repr(tei.rannot)
        ctail_left = ca_left.evidence if ca_left is not None else None
        ctail_right = ca_right.evidence if ca_right is not None else None
        # read info related to clipped tail clustering
        nc_left, nc_right = tei.case_cc_L, tei.case_cc_R
        cluster_info = "{}|{}".format(nc_left[2], nc_right[2])
        # filters
        if tei.name in filter_read:
            tei.add_filter("CONTROL.CLIPPED.READ")
        if tei.name in filter_lc:
            tei.add_filter("LOW_COMPLEXITY")
        tei.compute_min_nr_filter(args.min_coverage_fraction)
        tei.compute_cc_filter(NCC_CUTOFF)
        tei.compute_ctrl_cc_filter(NCC_CUTOFF)
        # retrieve the PE annotation
        annot_left = pe_left(tei.name)
        annot_right = pe_right(tei.name)
        annot_clipped = None  # this value is never used
        te = None  # transposon type
        confidence = "LOW"
        int_length = -1
        strand = "."
        inversion = False
        full_length = False
        # TODO: refactor: code for non-LTR transposons repeats itself
        if {Evidence.LTR3, Evidence.LTR5} == {ctail_right, ctail_left}:
            # 3' end is specific to HERVK, 5' end may be seen in SVAs
            # ==> both end signatures are required for HERVK
            te = "HERVK" # TODO: general retroviruses, read name from reference
            annot3, annot5, strand = (annot_right,
                annot_left, "-") if \
                ctail_right == Evidence.LTR3 else (annot_left, annot_right, 
                "+")
            annot3.filter_repeat_info(3)
            annot5.filter_repeat_info(5)
            check3end = annot3.check_transposon_evidence("HERVK")
            check5end = annot5.check_transposon_evidence("HERVK")
            if check3end and check5end:
                confidence = "HIGH"
                int_length = LTR_TRANSPOSON_LEN
            elif check3end or check5end:
                confidence = "MEDIUM"
        elif {Evidence.SINE, Evidence.polyA} == {ctail_right, ctail_left}:
            te = "SINE1"
            confidence = "MEDIUM"
            annot3, annot5, annot_clipped, strand = (
                annot_left, annot_right,
                ca_right, "-") if ctail_right == Evidence.SINE else (annot_right,
                annot_left, ca_left, "+")
            int_length = max(annot_clipped.rlen + 1 - annot_clipped.rpolya 
                - annot_clipped.rstart, 0)
            annot3.filter_repeat_info(3)
            annot5.filter_repeat_info(5, annot_clipped.rstart)
            if annot5.check_transposon_evidence("SINE1"
            ) or annot3.check_transposon_evidence("SINE1"):
                confidence = "HIGH"
        elif {Evidence.SVA, Evidence.polyA} == {ctail_right, ctail_left}:
            te = "SVA"
            annot3, annot5, annot_clipped, strand = (
                annot_left, annot_right, 
                ca_right, "-") if ctail_right == Evidence.SVA else (annot_right,
                annot_left, ca_left, "+")
            int_length = max(annot_clipped.rlen + 1 - annot_clipped.rpolya 
                - annot_clipped.rstart, 0)
            annot3.filter_repeat_info(3)
            annot5.filter_repeat_info(5, annot_clipped.rstart)
            if annot5.check_transposon_evidence("SVA"
            ) or annot3.check_transposon_evidence("SVA"):
                confidence = "HIGH"
        elif {Evidence.LINE_INV, Evidence.polyA} == {ctail_right, ctail_left}:
            inversion = True
            te = "LINE1"
            confidence = "MEDIUM"
            annot3, annot5, annot_clipped, strand = (
                annot_left, annot_right, 
                ca_right, "-") if ctail_right == Evidence.LINE_INV else (annot_right,
                annot_left, ca_left, "+")
            # int_length here is the length from inversion to the 3' end
            int_length = max(annot_clipped.rstart - annot_clipped.rpolya, 0)
            annot3.filter_repeat_info(3, inversion=True)
            #annot5.filter_repeat_info(annot_clipped.rstart)
            if annot5.check_transposon_evidence("LINE1"
            ) or annot3.check_transposon_evidence("LINE1"):
                confidence = "HIGH"
            # if annot3.readthrough is not None:
            #     te += ".RT"
        elif {Evidence.LINE, Evidence.polyA} == {ctail_right, ctail_left}:
            te = "LINE1"
            confidence = "MEDIUM"
            annot3, annot5, annot_clipped, strand = (
                annot_left, annot_right,
                ca_right, "-") if ctail_right == Evidence.LINE else (annot_right,
                annot_left, ca_left, "+")
            int_length = max(annot_clipped.rlen + 1 - annot_clipped.rpolya 
                - annot_clipped.rstart, 0)
            if annot_clipped.rstart < FULL_LENGTH_OFFSET:
                full_length = True
            annot3.filter_repeat_info(3)
            annot5.filter_repeat_info(5, annot_clipped.rstart)
            if annot5.check_transposon_evidence("LINE"
            ) or annot3.check_transposon_evidence("LINE"):
                confidence = "HIGH"
            # if annot3.readthrough is not None:
            #     te += ".RT"
        elif ctail_right == Evidence.polyA != ctail_left:
            strand = "+"
            te = "ORPHAN"
        elif ctail_left == Evidence.polyA != ctail_right:
            strand = "-"
            te = "ORPHAN"
        # skip candidates missing the required signatures
        if te is None:
            continue
        # assemble contigs if needed
        contigs_L = annot_left.assemble_contigs()
        contigs_R = annot_right.assemble_contigs()
        # summarize filters and output the result
        filters = tei.filter_status
        annot_left = str(ca_left) if ca_left is not None else "."
        annot_right = str(ca_right) if ca_right is not None else "."
        # VCF entry
        ch = tei.ch
        pos, end = sorted((tei.lbp, tei.rbp - 1)) # works for all cases
        try:
            ref_allele = gh[ch][max(0, pos - 1)].upper()
        except:
            ref_allele = "N"
        filters = tei.filter_status
        extra_vcf_args = {"CTRL_DOC_MEDIAN": tei.ctrl_doc.d_median,
            "CTRL_DOC_RATIO": tei.ctrl_doc.ratio} if tei.ctrl_doc else {}
        if te != "ORPHAN":
            # TODO: compute the actual values
            extra_vcf_args["DRP5"] = 0
            extra_vcf_args["DRP5S"] = 0
            extra_vcf_args["DRP5C"] = 0
            extra_vcf_args["DRP3"] = 0
            extra_vcf_args["DRP3S"] = 0
            extra_vcf_args["DRP3C"] = 0
        info_repr = VCF.construct_info_repr(
            END = end,
            STRAND = strand,
            INVERSION = inversion,
            SOMATIC = not tei.has_filter("CONTROL.CLIPPED.READ"),
            FULL_LENGTH = full_length,
            TE = te,
            TS_SEQUENCE = tei.ts_seq,
            TS_TYPE=tei.itype,
            BP5=tei.lbp if strand == "+" else tei.rbp,
            BP3=tei.rbp if strand == "+" else tei.lbp,
            BP5NR = tei.lnr if strand == "+" else tei.rnr,
            BP3NR = tei.rnr if strand == "+" else tei.lnr, 
            MQ5 = tei.lmq if strand == "+" else tei.rmq,
            MQ3 = tei.rmq if strand == "+" else tei.lmq,
            BP5SEQ = tei.lseq if strand == "+" else tei.rseq,
            BP3SEQ = tei.rseq if strand == "+" else tei.lseq,
            BP5ANNOT = tei.lannot if strand == "+" else tei.rannot,
            BP3ANNOT = tei.rannot if strand == "+" else tei.lannot, 
            CN_RATIO = tei.cn_ratio,
            PVALUE_CNV = tei.p_cnv,
            SVLEN = int_length,

            DEPTH_5UP = tei.case_doc.d_lbp if strand == "+"
                else tei.case_doc.d_rbp,
            DEPTH_5DN = tei.case_doc.d_lts if strand == "+"
                else tei.case_doc.d_rts,
            DEPTH_3UP = tei.case_doc.d_rbp if strand == "+"
                else tei.case_doc.d_lbp,
            DEPTH_3DN = tei.case_doc.d_rts if strand == "+"
                else tei.case_doc.d_lts,
            CASE_CC5 = tei.case_cc_L if strand == "+" else tei.case_cc_R,
            CASE_CC3 = tei.case_cc_R if strand == "+" else tei.case_cc_L,
            CTRL_CC5 = tei.ctrl_cc_L if strand == "+" else tei.ctrl_cc_R,
            CTRL_CC3 = tei.ctrl_cc_R if strand == "+" else tei.ctrl_cc_L,
            CASE_DOC_MEDIAN = tei.case_doc.d_median,
            CASE_DOC_RATIO = tei.case_doc.ratio,
            **extra_vcf_args)
        vc = VCF.from_values(ch, pos, tei.name, ref_allele, te, filters,
            info_repr)
        vc_list.append(vc)

vc_list.sort(key = attrgetter("pos"))
vc_list.sort(key = lambda vc: contig_order[vc.ch])

for vc in vc_list:
    sys.stdout.write(f"{vc}\n")
