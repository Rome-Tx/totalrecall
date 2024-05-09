include: "variables.snakefile"

rule breakpoints_1st_pass_fastq:
    input:
        "{sample}.clip.tsv"
    output:
        fq=temp("{sample}.1p.fq"),
        polya=temp("{sample}.1p.polya.tsv"),
        mask_pre=temp("{sample}.mask.pre.tsv"),
        mask=temp("{sample}.mask.tsv")
    shell:
        "compute_breakpoint_fq.py -d {MAX_PRECLUSTER_DIST} -l 5 "
        "-n {MIN_NR_1P} -o {output.polya} -O {output.mask_pre} -q 1 {input} "
        "> {output.fq};"
        "sort -k 1,1 -k 2,2n -k 3,3n {output.mask_pre} > {output.mask}"

rule clipped_reads_tsv:
    input:
        "{sample}.clip.namesorted.bam"
    output:
        temp("{sample}.clip.tsv")
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "bam2tsv_hc.py -n {wildcards.sample} -t {threads} {input} | "
        "sort --parallel={threads} -S 1G -T . "
        "-k {HTSR_FIELD_CONTIG},{HTSR_FIELD_CONTIG} "
        "-k {HTSR_FIELD_LOCUS},{HTSR_FIELD_LOCUS}g "
        "-k {HTSR_FIELD_BREAKPOINT},{HTSR_FIELD_BREAKPOINT} "
        "> {output}"

