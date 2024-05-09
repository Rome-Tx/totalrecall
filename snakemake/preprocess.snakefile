include: "variables.snakefile"

# ==== init and split the input BAMs ==== #

## store the header of the BAM files for reference (genomic contigs used
## in the reference, alignment parameters, deduplication etc.
rule sam_header:
    input:
        "{sample}.bam"
    output:
        "{sample}.header.sam"
    shell:
        "samtools view -H {input} > {output}"

rule sort_clipped_bam:
    input:
        "{sample}.clip.bam"
    output:
        "{sample}.clip.namesorted.bam"
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "samtools sort -n -@ {threads} -m {BAM_SORT_MB_PERCORE}M -o {output} {input}"

rule split_case_bam:
    input:
        {CASE_ALN}
    output:
        temp("case.disc.bam"),
        temp("case.clip.bam"),
        "case.tlen",
        "case.nbases",
        "case.rlen"
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "split_bam.py --write-drp -r {IXDIR}/genome.fa.gz -t {threads} "
        "-o case {input}"

rule split_ctrl_bam:
    input:
        {CONTROL_ALN}
    output:
        "control.clip.bam",
        "control.tlen",
        "control.nbases",
        "control.rlen"
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "split_bam.py -r {IXDIR}/genome.fa.gz -t {threads} -o control {input}"

rule genome_file:
    input:
        "{sample}.clip.bam"
    output:
        "{sample}.genome"
    shell:
        "mkgenomefile.py {input} > {output}"

rule case_bam:
    input:
        {CASE_ALN}
    output:
        "case.bam"
    shell:
        "ln -s {CASE_ALN} case.bam"

rule control_bam:
    input:
        {CONTROL_ALN}
    output:
        "control.bam"
    shell:
        "ln -s {CONTROL_ALN} control.bam"

rule case_index:
    input:
        {CASE_INDEX}
    output:
        "case.bam.bai"
    shell:
        "ln -s {CASE_INDEX} case.bam.bai"

rule control_index:
    input:
        {CONTROL_INDEX}
    output:
        "control.bam.bai"
    shell:
        "ln -s {CONTROL_INDEX} control.bam.bai"
