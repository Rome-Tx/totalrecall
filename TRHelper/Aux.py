from math import floor, ceil, sqrt
from statistics import stdev, median
from collections import deque, OrderedDict
from operator import itemgetter
from types import SimpleNamespace
from enum import Enum
import datetime
import gzip

from scipy.stats import chi2_contingency
import pysam
from Bio import SeqIO

from TRHelper.Align import check_identity
from TRHelper.Constants import FLAG_MUNMAP, FLAG_MREVERSE, Q_OFFSET

from TRHelper import _pysam_ext

Evidence = Enum("Evidence", ["LINE", "LINE_INV", "SINE", "SVA", "polyA",
    "LTR3", "LTR5"])


class BlastTab:
    """
    Alignment in Blast tabular format:
    query name, reference name, percent identity, alignment length, mismatches,
      gap opens, query start, query end, reference start, reference end,
      E-value, bit score
    Produced by lastal with the -f BlastTab option
    Coordinates are 1-based
    """
    def __init__(self, line):
        S = line.split()
        self.qname = S[0]
        self.rname = S[1]
        self.pident = float(S[2])
        self.alen = int(S[3])
        self.qstart = int(S[6])
        self.qend = int(S[7])
        self.rstart = int(S[8])
        self.rend = int(S[9])
        self.evalue = float(S[10])
        self.bitscore = float(S[11])
        self.strand = "+" if self.qstart < self.qend else "-"


class BPSAnnotation:
    def __init__(self):
        pass
        
    LINE_INV_POLYA_OFFSET = 20
    
    @classmethod
    def from_blast_tab(cls, bt, ref_length=-1, ref_polya=-1):
        bts = cls()
        bts.qname = bt.qname
        bts.rname, bts.subtype, bts._transposon = bt.rname.rsplit(sep=":",
            maxsplit=2)
        bts.qstart = bt.qstart
        bts.qend = bt.qend
        bts.rstart = bt.rstart
        bts.rend = bt.rend
        bts.rlen = ref_length
        bts.rpolya = ref_polya
        bts.pident = bt.pident
        bts.bitscore = bt.bitscore
        return bts
        
    @classmethod
    def from_text_repr(cls, line):
        if line == ".":
            return None
        bts = cls()
        S = line.strip().split("|")
        bts.qname = S[0]
        bts.rname, bts.subtype, bts._transposon = S[1].rsplit(sep=":",
            maxsplit=2)
        bts.qstart = int(S[2])
        bts.qend = int(S[3])
        bts.rstart = int(S[4])
        bts.rend = int(S[5])
        bts.rlen = int(S[6])
        bts.rpolya = int(S[7])
        bts.pident = float(S[8])
        bts.bitscore = float(S[9])
        return bts
    
    @classmethod
    def dummy_polya_annotation(cls, qname, bitscore):
        bts = cls()
        bts.qname = qname
        bts.rname = "polyA"
        bts.subtype = "polyA"
        bts._transposon = "polyA"
        bts.qstart = -1
        bts.qend = -1
        bts.rstart = -1
        bts.rend = -1
        bts.rlen = -1
        bts.rpolya = -1
        bts.pident = -1
        bts.bitscore = bitscore
        return bts

    @property
    def transposon(self):
        return self._transposon.upper()

    @property
    def numeric_id(self):
        nid = -1
        try:
            nid = int(self.qname)
        except ValueError:
            pass
        return nid
    
    @property
    def evidence(self):
        # TODO: use enum instead of strings
        if hasattr(self, "_evidence"):
            return self._evidence
        if self.transposon == "POLYA":
            self._evidence = Evidence.polyA
        elif self.transposon.startswith("LINE1"):
            if self.rname.endswith("INV"):
                if self.rstart < max(self.LINE_INV_POLYA_OFFSET, self.rpolya):
                    self._evidence = Evidence.polyA
                else:
                    self._evidence = Evidence.LINE_INV
            else:
                self._evidence = Evidence.LINE
        elif self.transposon == "SINE1":
            self._evidence = Evidence.SINE
        elif self.transposon == "LTR":
            if self.rname.lower().endswith("3end"):
                self._evidence = Evidence.LTR3
            elif self.rname.lower().endswith("5end"):
                self._evidence = Evidence.LTR5
            else:
                raise ValueError("cannot determine the end of the LTR: "
                    + self.rname.lower())
        elif self.transposon == "SVA":
            self._evidence = Evidence.SVA
        else:
            print(self.transposon)
            raise ValueError(f"unknown transposon: {self.transposon}")
        return self._evidence


    def __str__(self):
        return(f"{self.qname}|{self.rname}:{self.subtype}:{self._transposon}|"
            f"{self.qstart}|{self.qend}|{self.rstart}|{self.rend}|"
            f"{self.rlen}|{self.rpolya}|{self.pident}|{self.bitscore}")


class HTSRead:
    """
    Read as represented by discordant read TSV file
    """
    def __init__(self):
        pass

    @classmethod
    def from_line(cls, line):
        htsr = cls()
        htsr.line = line.strip()
        S = line.split()
        htsr.name = S[0]
        htsr.pe = int(S[1])
        htsr.is_rc = int(S[2])
        htsr.rname = S[3]
        htsr.ref_start = int(S[4])
        htsr.ref_end = int(S[5])
        htsr.mapq = int(S[6])
        htsr.flag = int(S[7])
        htsr.mate_rname = S[8]
        htsr.mate_ref_start = int(S[9])
        htsr.clip_left = S[10]
        htsr.clip_right = S[11]
        htsr.ref_locus = int(S[12])
        htsr.bkp = S[13]
        htsr.is_supplementary = int(S[14])
        htsr.sample = S[15]
        htsr.seq = S[16]
        htsr.qual = S[17]
        htsr.cseq = S[18]
        htsr.cqual = S[19]
        htsr.mate_aln_info = S[20]
        htsr.mate_seq = S[21]
        htsr.mate_qual = S[22]
        return htsr

    @classmethod
    def from_aligned_segment(cls, r, c_left, c_right, bkp, sample_name,
    clip_seq, clip_qual, store_full_seq=False):
        htsr = cls()
        if r.mate_is_unmapped:
            mate_rname = "."
            mate_start = -1
        else:
            mate_rname = r.next_reference_name
            mate_start = r.next_reference_start
        ref_locus = r.reference_start + 1 if bkp=="R" else r.reference_end
        seq = "."
        qual = "."
        if store_full_seq:
            seq = r.query_sequence
            qual = "".join(chr(q + Q_OFFSET) for q in r.query_qualities)
        htsr.line = (f"{r.query_name}\t{2 - r.is_read1}\t{int(r.is_reverse)}\t" 
            f"{r.reference_name}\t{r.reference_start + 1}\t{r.reference_end}\t"
            f"{r.mapping_quality}\t{r.flag}\t{mate_rname}\t{mate_start}\t"
            f"{c_left}\t{c_right}\t{ref_locus}\t{bkp}\t"
            f"{int(r.is_supplementary)}\t{sample_name}\t{seq}\t"
            f"{qual}\t{clip_seq.upper()}\t{clip_qual}\t.\t.\t.")
        return htsr

    def __str__(self):
        return self.line

    @property
    def mate_is_mapped(self):
        if self.flag & FLAG_MUNMAP:
            return False
        return True

    @property
    def mate_is_rc(self):
        if self.flag & FLAG_MREVERSE:
            return True
        return False

    def set_repeat_info(self, mate_aln_info, mate_seq, mate_qual):
        line_split = self.line.rsplit(maxsplit=3)[0]
        self.line = line_split + f"\t{mate_aln_info}\t{mate_seq}\t{mate_qual}"
        self.mate_aln_info = mate_aln_info
        self.mate_seq = mate_seq
        self.mate_qual = mate_qual

    def check_against_breakpoint(self, bp, allowed_dist=0):
        """
        Check that the read is clipped at the breakpoint
        and that the clipped sequences match
        """
        if self.bkp != bp.bkp or self.rname != bp.ch or abs(self.ref_locus -
        bp.pos) > allowed_dist:
            return False
        return check_identity(self.cseq, bp.seq)


class Breakpoint:
    """
    Breakpoint as represented in the text files
    """
    def __init__(self):
        pass
    
    @classmethod
    def from_line(cls, line):
        bp = cls()
        bp.line = line.strip()
        S = line.split()
        bp.id = int(S[0])
        bp.ch = S[1]
        bp.pos = int(S[2])
        bp.bkp = S[3]
        bp.nr = int(S[4])
        bp.seq = S[5]
        bp.mq = int(S[6])
        bp.annot = S[9]
        bp.used = False
        return bp

    @classmethod
    def from_cluster(cls, rl, bp_id, seq=None):
        bp = cls()
        #rl = sorted(rl, key=lambda r:len(r.cseq), reverse=True)  # not needed
        bp.id = bp_id
        ch = rl[0].rname
        bkp = rl[0].bkp
        if seq is None:
            seq = rl[0].cseq
        bp.seq = seq
        if bkp == "L":
            pos = floor(median([r.ref_locus for r in rl]))
        else:
            pos = ceil(median([r.ref_locus for r in rl]))
        pos_std = -1 if len(rl) == 1 else ceil(stdev(r.ref_locus for r in rl))
        Q = max(r.mapq for r in rl)
        seq_list = "|".join(r.cseq for r in rl)
        qual_list = "|".join(r.cqual for r in rl)
        bp.line = (f"{bp_id}\t{ch}\t{pos}\t{bkp}\t{len(rl)}\t{seq}\t{Q}\t"
            f"{seq_list}\t{qual_list}\t.")
        return bp

    @property
    def annotation(self):
        if hasattr(self, "_annotation"):
            return self._annotation
        elif hasattr(self, "annot"):
            self._annotation = (BPSAnnotation.from_text_repr(self.annot)
                if self.annot != "." else None)
            return self._annotation
        return None

    @annotation.setter
    def annotation(self, annot):
        if not isinstance(annot, BPSAnnotation):
            raise ValueError("BPSAnnotation instance required as annotation")
        self._annotation = annot
        ls = self.line.rsplit(maxsplit=1)[0]
        self.line = f"{ls}\t{annot}"

    def clipped_fasta(self):
        return f">{self.id}\n{self.seq}\n"

    def __str__(self):
        return self.line


class NamedBreakpoint:
    """
    Breakpoint of integration site
    """
    def __init__(self, line):
        self.line = line.strip()
        S = line.split()
        self.ch = S[0]
        self.pos = int(S[1])
        self.bkp = S[2]
        self.name = S[3]
        self.seq = S[4]

    def __str__(self):
        return self.line


class DOC:
    """
    DOC as computed by samtools depth
    """
    def __init__(self, line):
        S = line.split()
        self.ch = S[0]
        self.pos = int(S[1])
        self.depth = int(S[2])


class TargetSite:
    """
    Target site (insertion)
    """
    def __init__(self):
        pass
    
    @classmethod
    def from_line(cls, line):
        ts = cls()
        S = line.split()
        ts._split_line = S
        ts.name = S[0]
        ts.ch = S[1]
        ts.lbp = int(S[2])
        ts.rbp = int(S[3])
        if S[5] != ".":
            ts._ts_seq = S[5]
        ts.lnr = int(S[6])
        ts.lmq = int(S[7])
        ts.lseq = S[8]
        ts.lannot = S[9]
        ts.rnr = int(S[10])
        ts.rmq = int(S[11])
        ts.rseq = S[12]
        ts.rannot = S[13]
        return ts

    @classmethod
    def from_breakpoints(cls, name, bp1, bp2, gh=None):
        ts = cls()
        bpl, bpr = (bp1, bp2) if bp1.bkp == "L" else (bp2, bp1)
        ts.name = name
        ts.ch =bpl.ch
        ts.lbp = bpl.pos
        ts.rbp = bpr.pos
        ts.lnr = bpl.nr
        ts.lmq = bpl.mq
        ts.lseq = bpl.seq
        ts.lannot = bpl.annot
        ts.rnr = bpr.nr
        ts.rmq = bpr.mq
        ts.rseq = bpr.seq
        ts.rannot = bpr.annot
        if gh is not None:
            ts.compute_ts_seq(gh)
        return ts

    @property
    def itype(self):
        if self.rbp <= self.lbp:
            return "tsdup"
        elif self.rbp - self.lbp == 1:
            return "tsexact"
        else:
            return "tsdel"

    def compute_ts_seq(self, gh):
        TSS_FLANKING_LEN = 10
        try:
            seq = gh[self.ch]
        except KeyError:
            return
        # compute zero-based coordinates
        if self.rbp <= self.lbp:
            # TS duplication
            tss = self.rbp -1
            tse = self.lbp
        else:
            # TS deletion or an exact cut
            tss = self.lbp
            tse = self.rbp - 1
        start = max(0, tss - TSS_FLANKING_LEN)
        end = min(len(seq), tse + TSS_FLANKING_LEN)
        # extract the sequences
        fs_left = seq[start:tss].lower()
        fs_right = seq[tse:end].lower()
        tss = seq[tss:tse].upper()
        # set the attribute
        self._ts_seq = f"{fs_left}|{tss}|{fs_right}"

    @property
    def ts_seq(self):
        ts_seq = "."
        try:
            ts_seq = self._ts_seq
        except AttributeError:
            pass
        return ts_seq

    def estart(self, ext_length=0):
        """
        start = first basee left of the target site, 1-based
        """
        start = min(self.lbp, self.rbp - 1)
        return max(1, start - ext_length)

    def eend(self, ext_length=0, contig_length=None):
        """
        end  = first base right of the target site, 1-based
        """
        end = max(self.lbp + 1, self.rbp)
        if contig_length is not None:
            return min(contig_length[self.ch], end + ext_length)
        return end + ext_length

    def __str__(self):
        return (f"{self.name}\t{self.ch}\t{self.lbp}\t{self.rbp}\t"
            f"{self.itype}\t{self.ts_seq}\t{self.lnr}\t{self.lmq}\t"
            f"{self.lseq}\t{self.lannot}\t{self.rnr}\t{self.rmq}\t"
            f"{self.rseq}\t{self.rannot}")


class TSDOC:
    """
    Various DOC metrics at the target site
    """
    def __init__(self):
        pass

    @classmethod
    def from_line(cls, line):
        tsdoc = cls()
        S = line.split()
        tsdoc.ts_id = S[0]
        tsdoc.d_lbp = int(S[1])
        tsdoc.d_lts = int(S[2])
        tsdoc.d_rbp = int(S[3])
        tsdoc.d_rts = int(S[4])
        tsdoc.ratio = float(S[5])
        tsdoc.d_min = int(S[6])
        tsdoc.d_median = int(S[7])
        tsdoc.d_max = int(S[8])
        return tsdoc

    @classmethod
    def from_values(cls, ts_id, d_lbp, d_lts, d_rbp, d_rts, ratio, d_min,
    d_median, d_max):
        tsdoc = cls()
        tsdoc.ts_id = ts_id
        tsdoc.d_lbp = d_lbp
        tsdoc.d_lts = d_lts
        tsdoc.d_rbp = d_rbp
        tsdoc.d_rts = d_rts
        tsdoc.ratio = ratio
        tsdoc.d_min = d_min
        tsdoc.d_median = d_median
        tsdoc.d_max = d_max
        return tsdoc
    
    def __str__(self):
        return (f"{self.ts_id}\t{self.d_lbp}\t{self.d_lts}\t{self.d_rbp}\t"
            f"{self.d_rts}\t{self.ratio}\t{self.d_min}\t{self.d_median}\t"
            f"{self.d_max}")


class TEI(TargetSite):
    """
    Transposable element insertion
    """
    def __init__(self):
        pass

    @classmethod
    def from_line(cls, line, case_doc_line, ctrl_doc_line, cluster_line1,
        cluster_line2, ctrl_cluster_line1, ctrl_cluster_line2):
        tei = super().from_line(line)
        tei._src_line = line.strip()
        tei._filter = set()
        tsdoc = TSDOC.from_line(case_doc_line)
        if tsdoc.ts_id != tei.name:
            raise ValueError(f"case DOC id ({tsdoc.ts_id}) does not match "
            f"that of the insertion ({tei.name})")
        tei.case_doc = tsdoc
        if ctrl_doc_line:
            tsdoc = TSDOC.from_line(ctrl_doc_line)
            if tsdoc.ts_id != tei.name:
                raise ValueError(f"control DOC id ({tsdoc.ts_id}) does not match "
                f"that of the insertion ({tei.name})")
            tei.ctrl_doc = tsdoc
        else:
            tei.ctrl_doc = None
        # control clipped read clusters
        for ln in (cluster_line1, cluster_line2):
            S = ln.split()
            if S[0].strip("@LR") != tei.name:
                raise ValueError(f"Id in the cluster file ({S[0]}) "
                    f"does not match that of the insertion ({tei.name})")
            if S[0].endswith("R"):
                tei.nc_right = [int(x) for x in S[1:]]
            elif S[0].endswith("L"):
                tei.nc_left = [int(x) for x in S[1:]]
        # case clipped read clusters
        for ln in (ctrl_cluster_line1, ctrl_cluster_line2):
            S = ln.split()
            if S[0].strip("@LR") != tei.name:
                raise ValueError(f"Id in the cluster file ({S[0]}) "
                f"does not match that of the insertion ({tei.name})")
            if S[0].endswith("R"):
                tei.ctrl_nc_right = [int(x) for x in S[1:]]
            elif S[0].endswith("L"):
                tei.ctrl_nc_left = [int(x) for x in S[1:]]
        return tei

    def compute_p_cnv(self, n_bases_case, n_bases_ctrl):
        if n_bases_ctrl:
            try:
                stat, self.p_cnv, dof, expected = chi2_contingency([
                [self.ctrl_doc.d_median, self.case_doc.d_median],
                [n_bases_ctrl, n_bases_case]])
            except ValueError:
                self.p_cnv = -1
            try:
                self.cn_ratio = (self.case_doc.d_median * n_bases_case 
                    ) / (self.ctrl_doc.d_median * n_bases_ctrl)
            except ZeroDivisionError:
                self.cn_ratio = -1
        else:
            self.p_cnv = -1
            self.cn_ratio = -1

    def compute_min_nr_filter(self, min_clipped_read_fraction):
        if (self.lnr < min_clipped_read_fraction * self.case_doc.d_lbp and
        self.rnr < min_clipped_read_fraction * self.case_doc.d_rbp):
            self._filter.add("LOW_FRAC_OF_DEPTH")

    def compute_cc_filter(self, ncc_cutoff):
        for ncL, ncR, cutoff in zip(self.nc_left, self.nc_right, ncc_cutoff):
            if ncL > cutoff or ncR > cutoff:
                self._filter.add("CLIPPED.CLUSTERS")
                break
    def compute_ctrl_cc_filter(self, ncc_cutoff):
        for ncL, ncR, cutoff in zip(self.ctrl_nc_left, self.ctrl_nc_right,
        ncc_cutoff):
            if ncL > cutoff or ncR > cutoff:
                self._filter.add("CLIPPED.CLUSTERS")
                break

    def add_filter(self, filter_name):
        self._filter.add(filter_name)

    @property
    def filter_status(self):
        if self._filter:
            return ",".join(sorted(self._filter))
        return "PASS"
    
    def has_filter(self, filter_name):
        return filter_name in self._filter

    @property
    def case_cc_L(self):
        return ":".join(str(x) for x in self.nc_left)

    @property
    def case_cc_R(self):
        return ":".join(str(x) for x in self.nc_right)

    @property
    def ctrl_cc_L(self):
        return ":".join(str(x) for x in self.ctrl_nc_left)

    @property
    def ctrl_cc_R(self):
        return ":".join(str(x) for x in self.ctrl_nc_right)

    def __str__(self):
        return (f"{self.name}\t{self.ch}\t{self.lbp}\t{self.rbp}\t{self.itype}\t"
            f"{self.ts_seq}\t{self.lnr}\t{self.case_doc.d_lbp}\t"
            f"{self.case_doc.d_lts}\t{self.lseq}\t{self.lannot}\t"
            f"{self.rnr}\t{self.case_doc.d_rbp}\t{self.case_doc.d_rts}\t"
            f"{self.rseq}\t{self.rannot}")
        #return self._src_line


class VCF:
    """
    VCF representation
    """
    class info_field:
        def __init__(self, description, reader=lambda x:x, formatter=str):
            self.formatter = formatter
            self.reader = reader
            self.description = description
        
        @property
        def type(self):
            if self.reader is int:
                return "Integer"
            elif self.reader is float:
                return "Float"
            elif self.reader is bool:
                return "Flag"
            return "String"

        @property
        def number(self):
            if self.reader is bool:
                return 0
            return 1

    info_attr_list = OrderedDict([
        ("SOMATIC", info_field("Variant is somatic", bool)),
        ("INVERSION", info_field("LINE element contains inversion as a result "
            "of twin priming", bool)),
        ("FULL_LENGTH", info_field("Non-LTR Transposon is full length", bool)),
        ("STRAND", info_field("Strand")),
        ("END", info_field("Another breakpoint", int)),
        ("TE", info_field("Transposon")),
        ("TS_SEQUENCE", info_field("Target site sequence")),
        ("TS_TYPE",
            info_field("Target site type (duplication/deletion/exact cut)")),
        ("SVLEN", info_field("Inferred transposon length excluding the "
            "poly(A); for inverted LINEs distance from the 3' end to the twin "
            "priming locus", int)),
        ("BP5", info_field("Breakpoint for the 5' end of the transposon", int)),
        ("BP3", info_field("Breakpoint for the 3' end of the transposon", int)),
        ("MQ5", info_field("Mapping quality of the breakpoint for the 5' end "
            "of the transposon", int)),
        ("MQ3", info_field("Mapping quality of the breakpoint for the 3' end "
            "of the transposon", int)),
        ("BP5NR", info_field("Number of clipped reads for the 5' end "
            "of the transposon", int)),
        ("BP3NR", info_field("Number of clipped reads for the 3' end "
            "of the transposon", int)),
        ("BP5SEQ", info_field("Clipped sequence for the breakpoint for the 5' "
            "end of the transposon")),
        ("BP3SEQ", info_field("Clipped sequence for the breakpoint for the 3' "
            "end of the transposon")),
        ("BP5ANNOT", info_field("Alignment of the clipped sequence for the 5' "
            "end", BPSAnnotation.from_text_repr)),
        ("BP3ANNOT", info_field("Alignment of the clipped sequence for the 3' "
            "end", BPSAnnotation.from_text_repr)),
        ("DRP5", info_field("Number of DRPs at the 5' end", int)),
        ("DRP5S", info_field("Number of DRPs at the 5' end supporting the call",
            int)),
        ("DRP5C", info_field("Number of DRPs at the 5' end supporting the "
            "call which are consistently clipped at the 5' breakpoint", int)),
        ("DRP3", info_field("Number of DRPs at the 3' end", int)),
        ("DRP3S", info_field("Number of DRPs at the 3' end supporting the call",
            int)),
        ("DRP3C", info_field("Number of DRPs at the 3' end supporting the "
            "call which are consistently clipped at the 3' breakpoint", int)),
        ("CN_RATIO", info_field("Copy number ratio, case:control", float)),
        ("PVALUE_CNV", info_field("P-value for the CNV test", float)),
        ("DEPTH_5UP", info_field("Depth of coverage upstream of the 5' end",
            int)),
        ("DEPTH_5DN", info_field("Depth of coverage downstream of the 5' end",
            int)),
        ("DEPTH_3UP", info_field("Depth of coverage upstream of the 3' end",
            int)),
        ("DEPTH_3DN", info_field("Depth of coverage downstream of the 3' end",
            int)),
        ("CASE_CC5",  info_field("Number of clipped read clusters at various "
            "distances form the 5' end breakpoint, case")),
        ("CASE_CC3", info_field("Number of clipped read clusters at various "
            "distances form the 3' end breakpoint, case")),
        ("CTRL_CC5",  info_field("Number of clipped read clusters at various "
            "distances form the 5' end breakpoint, control")),
        ("CTRL_CC3", info_field("Number of clipped read clusters at various "
            "distances form the 3' end breakpoint, control")),
        ("CASE_DOC_MEDIAN", info_field("Median depth of coverage, case", int)),
        ("CASE_DOC_RATIO", info_field("High to low depth of coverage ratio, "
            "case", float)),
        ("CTRL_DOC_MEDIAN", info_field("Median depth of coverage, control",
            int)),
        ("CTRL_DOC_RATIO", info_field("High to low depth of coverage ratio, "
            "control", float))
    ])

    @classmethod
    def header(cls, sam_header_file=None):
        d = datetime.date.today()
        p1 = ("##fileformat=VCFv4.3\n"
        "##fileDate={:4d}{:02d}{:02d}\n"
        "##source=TotalReCall\n").format(d.year, d.month, d.day)
        p2 = ""
        if sam_header_file:
            with pysam.AlignmentFile(sam_header_file, "r") as af:
                p2 = "".join("""##contig=<ID={},length={}>\n""".format(
                    r, af.get_reference_length(r)) for r in af.references)
        p3 = "".join(f"""##INFO=<ID={k},Number={v.number},"""
            f"""Type={v.type},Description="{v.description}">\n"""
            for k, v in cls.info_attr_list.items())
        p4 =("""##ALT=<ID=INS:ME:LINE1,Description="Insertion of L1 element">\n"""
        """##ALT=<ID=INS:ME:SINE1,Description="Insertion of ALU element">\n"""
        """##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">\n"""
        """##ALT=<ID=INS:ME:RV,Description="Insertion of a retroviral element">\n""")
        p5 = ("""##FILTER=<ID=CONTROL.CLIPPED.READ,Description="Matching clipped read in the control">\n"""
        """##FILTER=<ID=CLIPPED.CLUSTERS,Description="Too many clusters of clipped reads in the vicinity">\n"""
        """##FILTER=<ID=LOW_COMPLEXITY,Description="TODO">\n"""
        """##FILTER=<ID=LOW_FRAC_OF_DEPTH,Description="Number of clipped reads is a small fraction of total depth">\n"""
        """#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n""")
        return p1 + p2 + p3 + p4 + p5

    def __init__(self):
        pass

    @classmethod
    def construct_info_repr(cls, **kwargs):
        flag_list = [k for k, v in cls.info_attr_list.items()
            if k in kwargs and v.reader is bool and kwargs[k]]
        kv_list = [f"{k}={v.formatter(kwargs[k])}"
            for k, v in cls.info_attr_list.items() if k in kwargs
            and v.reader is not bool]
        return ";".join(flag_list + kv_list)

    @classmethod
    def from_values(cls, ch, pos, v_id, ref, me_name, filters_repr, info_repr):
        vc = cls()
        vc.ch = ch
        vc.pos = pos
        vc.id = v_id
        vc.ref = ref
        vc.alt = f"<INS:ME:{me_name}>"
        vc.filters_repr = filters_repr
        vc.info_repr = info_repr
        return vc

    @classmethod
    def from_line(cls, line):
        if line.startswith("#"):
            raise ValueError("Comment/header line cannot be parsed as a "
                f"variant: {line.strip()}")
        vc = cls()
        S = line.split()
        vc.ch, vc.pos, vc.id, vc.ref = S[0], int(S[1]), S[2], S[3]
        vc.alt, vc.filters_repr, vc.info_repr = S[4], S[6], S[7] # S[5] is "."
        return vc

    @property
    def info(self):
        if not hasattr(self, "_info"):
            self._info = SimpleNamespace()
            for k, v in self.info_attr_list.items():
                if v.reader is bool:
                    setattr(self._info, k, False)
            for attr in self.info_repr.split(";"):
                if "=" in attr:
                    k, v = attr.split("=")
                    setattr(self._info, k, self.info_attr_list[k].reader(v))
                elif (attr in self.info_attr_list
                and self.info_attr_list[attr].reader is bool):
                    setattr(self._info, attr, True)
                else:
                    raise ValueError(f"malformed info field: {self.info_repr}")
        return self._info

    @property
    def filters(self):
        if not hasattr(self, "_filters"):
            self._filters = self.filters_repr.split(",")
        return self._filters

    def __str__(self):
        return (f"{self.ch}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t.\t"
            f"{self.filters_repr}\t{self.info_repr}")


class MateRepeatAlnInfo:
    def __init__(self):
        pass

    @classmethod
    def from_LAST(cls, processed_last_output):
        """
        Select output from LAST:
        read name, read 1/2, subject, pident, alignment length,
        sstart, send, query length, subject length, raw score, strand
        """
        mri = cls()
        S = processed_last_output.split()
        mri.qname = S[0]
        mri.pe = int(S[1])
        mri.rname, mri.subtype, mri._transposon = S[2].rsplit(sep=":",
            maxsplit=2)
        mri.pident = float(S[3])
        mri.alen = int(S[4])
        mri.rstart = int(S[5])
        mri.rend = int(S[6])
        mri.qlen = int(S[7])
        mri.rlen = int(S[8])
        mri.score = int(S[9])
        mri.strand = S[10]
        return mri
    
    @classmethod
    def from_short_repr(cls, line):
        mri = cls()
        S = line.strip().split(sep="|")
        mri._split_line = S
        mri.rname, mri.subtype, mri._transposon = S[0].rsplit(sep=":",
            maxsplit=2)
        mri.score = int(S[1])
        mri.pident = float(S[2])
        mri.rstart = int(S[3])
        mri.rend = int(S[4])
        mri.rlen = int(S[5])
        mri.qlen = int(S[6])
        mri.strand = S[7]
        return mri

    @classmethod
    def from_long_repr(cls, line):
        mri = cls.from_short_repr(line)
        mri.distance_to_breakpoint = int(mri._split_line[8])
        mri.mate_seq = mri._split_line[9]
        mri.mate_qual = mri._split_line[10]
        return mri

    @property
    def transposon(self):
        return self._transposon.upper()

    def _str_short(self):
        return "{}:{}:{}|{}|{}|{}|{}|{}|{}|{}".format(self.rname, 
            self.subtype, self._transposon, self.score,
            self.pident, self.rstart, self.rend, self.rlen, self.qlen,
            self.strand)

    def __str__(self):
        if hasattr(self, "mate_seq"):
            return (self._str_short() + f"|{self.distance_to_breakpoint}|"
                f"{self.mate_seq}|{self.mate_qual}")
        else:
            return f"{self._str_short()}"


class IntervalWrapper:
    """
    Wrapper for intersection function
    """
    def __init__(self, ch, start, end, obj):
        """
        ch, start, end - interval
        obj - wrapped object
        """
        self.ch, self.start, self.end = ch, start, end
        self.obj = obj

    def __lt__(self, other):
        return (self.ch, self.end) < (other.ch, other.start)

    def __gt__(self, other):
        return (self.ch, self.start) > (other.ch, other.end)

    def __and__(self, other):
        """
        Check for intersection
        """
        return not (self < other or self > other)

def intersect_iter(iter1, iter2, one_by_one=True, yield_all=False):
    """
    Uses sorted iter1 and iter2 and iterates through them 
    to find the intersection
    Iterators should yield IntervalWrapper instances which have
    ch, start, end attributes
    Iterators should be sorted by ch, then by start, then (optionally) by end
    one_by_one: yield pairs from i1, i2 (True) or i1, list_of_matches_from_i2
    (False)
    yield_all: when one_by_one=False, do not skip i1 even if there are 
    no intersecting elements from i2 (intersecting_i2 is empty)
    """
    d = deque()
    for i1 in iter1:
        while not d or not d[-1] > i1:
            try:
                i2 = next(iter2)
            except StopIteration:
                break
            if not i2 < i1:
                d.append(i2)
        while d and d[0] < i1:
            d.popleft()
        intersecting_i2 = tuple(i2.obj for i2 in d if i1 & i2)
        if one_by_one:
            for i2obj in intersecting_i2:
                yield i1.obj, i2obj
        elif yield_all or intersecting_i2:
            yield i1.obj, intersecting_i2

def intersect_multiple_iters(iter1, *iter_list):
    """
    Uses sorted iter1 and iterators form iter_list and iterates through them 
    to find the intersection
    Iterators should yield IntervalWrapper instances which have
    ch, start, end attributes
    Iterators should be sorted by ch, then by start, then (optionally) by end
    yield tuples (i1, matches_from_first_iter, matches_from_second_iter, ...)
    """
    deque_list = list(deque() for i in iter_list)
    for i1 in iter1:
        intersecting_elements = deque()
        for (d, iter2) in zip(deque_list, iter_list):
            while not d or not d[-1] > i1:
                try:
                    i2 = next(iter2)
                except StopIteration:
                    break
                if not i2 < i1:
                    d.append(i2)
            while d and d[0] < i1:
                d.popleft()
            intersecting_i2 = tuple(i2.obj for i2 in d if i1 & i2)
            intersecting_elements.append(intersecting_i2)
        yield (i1.obj,) + tuple(intersecting_elements)


def MaskedInterval(line):
    """
    Simple interval: contig, start, end
    """
    S = line.split()
    ch = S[0]
    start = int(S[1])
    end = int(S[2])
    return IntervalWrapper(ch, start, end, None)


def read_iter(fname, masked_fname, offset=10):
    """
    Iterator for reads excluding masked regions
    """
    if masked_fname:
        with open(fname) as fid, open(masked_fname) as masked_fid:
            masked_iter = (MaskedInterval(line) for line in masked_fid)
            htsr_iter = (IntervalWrapper(htsr.rname, htsr.ref_locus - offset,
                htsr.ref_locus + offset, htsr) for htsr in 
                (HTSRead.from_line(line) for line in fid))
            for htsr, masked in intersect_iter(htsr_iter, masked_iter,
            one_by_one=False, yield_all=True):
                if not masked:  # masked is empty tuple
                    yield htsr
    else:
        with open(fname) as fid:
            yield from (HTSRead.from_line(line) for line in fid)


def template_length_stats(fname, drop_frac=0.001):
    """
    fname = file with fragment length distribution
    drop_frac = percentile to drop when computing the cutoff length
    """
    tlen = 0
    tlen_sq = 0
    N = 0
    dist = list()
    with open(fname) as fid:
        for line in fid:
            S = line.split()
            L, n = int(S[0]), int(S[1])
            dist.append((L, n))
            N += n
            tlen += n * L
            tlen_sq += n * L**2
    N = sum(n for L, n in dist)
    tlen = sum(n * L for L, n in dist)
    tlen_sq = sum(n * L**2 for L, n in dist)
    obs_tlen = tlen / N
    obs_sd = sqrt(tlen_sq/N - obs_tlen**2)
    dist.sort(key=itemgetter(0), reverse=True)  # sort by length
    N_drop = N * drop_frac
    n_current = 0
    for L, n in dist:
        if n_current > N_drop:
            ext_len = L
            break
        n_current += n
    else:
        raise ValueError("cannot compute fragment length")
    return obs_tlen, obs_sd, ext_len


def read_gzipped_genome(genome_gz_file_name):
    """
    Read genome into memory from a gzipped fasta file
    """
    with gzip.open(genome_gz_file_name, "rt") as fh:
        gh = {rec.id:str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
    return gh

def read_genome_sequence_lengths(genome_gz_file_name):
    """
    Read genome into memory from a gzipped fasta file
    """
    with gzip.open(genome_gz_file_name, "rt") as fh:
        glh = {rec.id:len(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
    return glh


def split_bam(infile, outfile_clip, outfile_disc, args):
    return _pysam_ext.split_bam(infile, outfile_clip, outfile_disc,
        dist_conc = args.dist_conc, dist_disc = args.dist_disc,
        min_clip_length = args.min_clip_length, min_end_qual =
        args.min_end_qual, write_drp = args.write_drp)
