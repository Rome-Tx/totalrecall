# cython: language_level=3

from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcalignmentfile cimport (BAM_FPROPER_PAIR, BAM_FPAIRED,
    BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP, BAM_FUNMAP, BAM_FMUNMAP,
    BAM_FSUPPLEMENTARY)
from pysam.libcalignedsegment cimport (pysam_get_n_cigar, pysam_bam_get_cigar,
    pysam_bam_get_qual)
from pysam.libchtslib cimport (BAM_CIGAR_MASK, BAM_CIGAR_SHIFT, BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP, BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF)
from collections import defaultdict

#TODO: use hts_pos_t instead of uint64_t - how to import it?

cdef inline int check_concordant(AlignedSegment r, int64_t dist_conc):
    """
    check that both reads in the pair is mapped
    and check that they are next to each other
    """
    cdef uint16_t f = r._delegate.core.flag
    if f & (BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP):
        return 0
    if not (f & BAM_FPROPER_PAIR):
        return 0
    if r._delegate.core.tid == r._delegate.core.mtid and abs(
    r._delegate.core.pos - r._delegate.core.mpos) <= dist_conc:
        return 1
    return 0


cdef inline int check_discordant(AlignedSegment r, int64_t dist_disc):
    """
    check that at least one of the reads in the pair is mapped
    and check that they are discordant
    """
    cdef uint16_t f = r._delegate.core.flag
    if (f & BAM_FUNMAP) and (f & BAM_FMUNMAP):
        # both read and its mate are unmapped
        return 0
    cdef uint16_t unmap = 1 if (f & BAM_FUNMAP) else 0
    cdef uint16_t munmap = 1 if (f & BAM_FMUNMAP) else 0
    if unmap != munmap:
        # only one of the read and its mate is unmapped
        return 1
    # do other checks only after the check for unmapped since some other
    # flags are meaningless for unmapped reads
    if f & (BAM_FPROPER_PAIR | BAM_FSUPPLEMENTARY):
        # read is mapped as a proper pair or read is supplementary
        return 0
    if r._delegate.core.tid == r._delegate.core.mtid and abs(
    r._delegate.core.pos - r._delegate.core.mpos) <= dist_disc:
        # reads are not far enough
        return 0
    return 1


cdef inline int check_soft_clip(AlignedSegment r, uint32_t min_clip_length,
uint8_t min_end_qual):
    """
    check that the read is soft clipped and that it happens on one side only
    """
    cdef uint16_t f = r._delegate.core.flag
    if f & BAM_FUNMAP:
        return 0
    if f & BAM_FSUPPLEMENTARY:
        # all supplementary alignments are clipped
        return 1
    # see code for property cigartuples and property query_qualities
    # in pysam.libcalignedsegment.pyx for how to work 
    # with the cigar string and base qualities
    cdef uint32_t n_cigar = pysam_get_n_cigar(r._delegate)
    if n_cigar == 0:
        return 0
    cdef uint32_t * cigar_p  # integer representation of op, l in cigar
    cdef uint8_t* q  # qualities
    cigar_p = pysam_bam_get_cigar(r._delegate)
    cdef uint32_t op, l  # cigar entry: operation, length
    cdef int32_t L  # query length
    cdef uint32_t qc_fail_bases, pos
    # TODO: account for both hard and soft clips at the end
    # 5' end
    cdef uint8_t is_clipped_left = 0
    op = cigar_p[0] & BAM_CIGAR_MASK
    if op == BAM_CSOFT_CLIP:
        l = cigar_p[0] >> BAM_CIGAR_SHIFT
        qc_fail_bases = 0
        q = pysam_bam_get_qual(r._delegate)
        for pos in range(l):
            if q[pos] < min_end_qual:
                qc_fail_bases += 1
            else:
                break
        if l >= min_clip_length + qc_fail_bases:
            is_clipped_left = 1
    # 3' end
    cdef uint8_t is_clipped_right = 0
    op = cigar_p[n_cigar - 1] & BAM_CIGAR_MASK
    if op == BAM_CSOFT_CLIP:
        l = cigar_p[n_cigar - 1] >> BAM_CIGAR_SHIFT
        L = r._delegate.core.l_qseq
        qc_fail_bases = 0
        q = pysam_bam_get_qual(r._delegate)
        for pos in range(L - 1, L - 1 - l, -1):
            if q[pos] < min_end_qual:
                qc_fail_bases += 1
            else:
                break
        if l >= min_clip_length + qc_fail_bases:
            is_clipped_right = 1
    # check that exactly one end is clipped
    # TODO: specify why
    return 1 if is_clipped_left != is_clipped_right else 0


cdef inline int compute_mapped_bases(AlignedSegment r):
    cdef uint32_t n_cigar = pysam_get_n_cigar(r._delegate)
    if n_cigar == 0:
        return 0
    cdef uint32_t* cigar_p = pysam_bam_get_cigar(r._delegate)
    cdef uint32_t k, mapped_bases = 0
    cdef int op
    for k in range(n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        # TODO: count insertion as mapped bases?
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            mapped_bases += cigar_p[k] >> BAM_CIGAR_SHIFT
    return mapped_bases


def split_bam(AlignmentFile infile, AlignmentFile outfile_clip,
AlignmentFile outfile_disc, dist_conc, dist_disc, min_clip_length,
min_end_qual, write_drp):
    cdef int64_t dist_conc_c = dist_conc
    cdef int64_t dist_disc_c = dist_disc
    cdef uint32_t min_clip_length_c = min_clip_length
    cdef uint8_t min_end_qual_c = min_end_qual
    nBases = 0
    tlen_dist = defaultdict(int)
    rlen_dist = defaultdict(int)

    cdef AlignedSegment read
    cdef uint16_t f, skip_flag

    skip_flag = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    
    for read in infile:
        f = read._delegate.core.flag
        if f & skip_flag:
            # skip secondary, duplicate and qcfail alignments
            continue
        # check for paired in sequencing
        if not (f & BAM_FPAIRED):
            continue
        # count mapped bases/read length
        mapped_bases = compute_mapped_bases(read)
        if mapped_bases > 0:
            nBases += mapped_bases
            # read length distribution
            if not (f & BAM_FSUPPLEMENTARY):
                rlen_dist[read.infer_read_length()] += 1
        # check for clipping
        if check_soft_clip(read, min_clip_length_c, min_end_qual_c):
            outfile_clip.write(read)
        # process properly paired/discordant reads
        if check_concordant(read, dist_conc_c):
            # template length
            if read.template_length > 0:
                tlen_dist[read.template_length] += 1
        elif write_drp and check_discordant(read, dist_disc_c):
            outfile_disc.write(read)
    return nBases, tlen_dist, rlen_dist
