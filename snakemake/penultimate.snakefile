include: "variables.snakefile"

rule annotated_breakpoints:
    input:
        byid="{sample}.breakpoints.byid.tsv",
        annot="{sample}.clip_annot.tsv"
    output:
        "{sample}.breakpoints.tsv"
    shell:
        "merge_breakpoint_annotation.py {input.byid} {input.annot} | "
        "sort -T . -k 2,2 -k 3,3n -k 4,4 > {output}"

rule candidate_integration_sites:
    input:
        "{sample}.breakpoints.tsv"
    output:
        temp("{sample}.ins.cand.tsv")
    shell:
        "compute_candidate_sites.py -d {MAX_SITE_LEN} -r {MIN_NR_2P} "
        "{input} {IXDIR}/genome.fa.gz > {output}"

rule pe_info_annotation:
    input:
        bp="case.ins.nbp.tsv",
        pe="{sample}.disc.tsv",
        tlen="{sample}.tlen"
    output:
        "aux/{sample}.pe.{breakpoint}.tsv"
    params:
        lr="{breakpoint}"
    shell:
        "annotate_with_pe_info.py -Q 1 -b {IXDIR}/L1_FLI.bed {input.bp} {params.lr} "
        "{input.pe} {input.tlen} | sort -S 500M -T . -k 1,1 > {output}"

rule filter_clipped_reads_cluster:
    input:
        bp="case.ins.nbp.tsv",
        cclip="{sample}.clip.tsv",
        masked = "case.mask.tsv"
    output:
        "aux/{sample}.read.cluster"
    shell:
        """filter_clipped_reads_cluster.py -l """
        """{FILTER_CLUSTER_CLIP_LEN} -m {input.masked} """
        """{input.bp} {input.cclip} """
        """| sort -S 500M -T . -t "@" -k 1,1 > {output}"""

rule filter_control_reads:
    input:
        bp="case.ins.nbp.tsv",
        cclip="control.clip.tsv",
    output:
        "aux/control.clipped.read"
    shell:
        "filter_control_reads.py -d {FILTER_READ_EXT} -l "
        "{MIN_CONTROL_READ_LEN} {input.bp} {input.cclip} "
        "| sort -S 100M -T . -u > {output}"

rule isites_to_named_breakpoints:
    input:
        "{sample}.ins.cand.tsv"
    output:
        temp("{sample}.ins.nbp.tsv")
    shell:
        """cat {input} | awk '{{print ${TS_FIELD_CH} "\t" ${TS_FIELD_LBP} """
        """ "\tL\t" ${TS_FIELD_NAME} "@L\t" ${TS_FIELD_LSEQ}; """
        """print ${TS_FIELD_CH} "\t" ${TS_FIELD_RBP} "\tR\t" """
        """${TS_FIELD_NAME} "@R\t" ${TS_FIELD_RSEQ}}}' """
        """| sort -S 500M -T . -k 1,1 -k 2,2n -k 3,3 > {output}"""

rule isites_bed:
    input:
        is_tsv="case.ins.cand.tsv",
        genome="case.genome",
        tlen="case.tlen"
    output:
        is_bed=temp("is.bed"),
        drp_bed=temp("drp.bed")
    shell:
        "write_isites_bed.py -e {BED_EXT_LEN} {input.genome} {input.is_tsv} "
        "| sort -S 500M -T . -k 1,1 -k 2,2n > {output.is_bed}; "
        "write_isites_bed.py -t {input.tlen} {input.genome} {input.is_tsv} "
        "| sort -S 500M -T . -k 1,1 -k 2,2n > {output.drp_bed}"

## auxiliary file with insertion sites sorted by the name
rule insertions_name_sorted:
    input:
        "{sample}.ins.cand.tsv"
    output:
        temp("{sample}.ins.namesorted.tsv")
    shell:
        "sort -S 500M -T . -k 1,1 {input} > {output}"

# ==== filters based on low complexity ==== #

rule low_complexity_filter:
    input:
        nbp="{sample}.ins.nbp.tsv",
        dust="{sample}.bps.dust"
    output:
        "aux/{sample}.low_complexity"
    params:
        bn="{sample}"
    shell:
        "compute_lc_filter.py {params.bn} | sort -T . -u > {output}"

rule mask_breakpoint_sequence:
    input:
        "{sample}.bps.fa"
    output:
        "{sample}.bps.dust"
    shell:
        "dustmasker -in {input} -outfmt fasta > {output}"

rule genomic_breakpoint_sequence:
    input:
        "{sample}.ins.nbp.tsv"
    output:
        "{sample}.bps.fa"
    shell:
        "extract_bp_seq.py -n {input} {IXDIR}/genome.fa.gz > {output}"

# ==== repeat annotation for discordant reads ==== #

rule merge_discordant_repeat_info:
    input:
        bam="{sample}.drp.subset.bam",
        repeatinfo="{sample}.repeatinfo.tsv",
        fq="{sample}.disc.fq"
    output:
        tsv=temp("{sample}.disc.tsv"),
        tmp_tsv=temp("{sample}.disc.tmp.tsv")
    params:
        basename="{sample}"
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "bam2tsv.py -Q 0 -n {params.basename} -t {threads} {input.bam} | "
        "sort --parallel={threads} -S 500M -T . -k 1,1 > {output.tmp_tsv}; "
        "merge_repeat_info.py {input.repeatinfo} {output.tmp_tsv} {input.fq} | "
        "sort --parallel={threads} -S 500M -T . "
        "-k {HTSR_FIELD_CONTIG},{HTSR_FIELD_CONTIG} "
        "-k {HTSR_FIELD_LOCUS},{HTSR_FIELD_LOCUS}g "
        "-k {HTSR_FIELD_BREAKPOINT},{HTSR_FIELD_BREAKPOINT} "
        "> {output.tsv}"
