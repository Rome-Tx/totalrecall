include: "variables.snakefile"

# ==== compute coverage ==== #

rule doc_distribution:
    input:
        doc="{sample}.doc",
        genome="case.genome",
        ins="case.ins.cand.tsv"
    output:
        "aux/doc_{sample}.tsv"
    shell:
        "compute_coverage_info.py {input.doc} | sort -S 500M -T . -k 1,1"
        " > {output}"

rule samtools_compute_doc:
    input:
        bam="{sample}.bam",
        bai="{sample}.bam.bai",
        bed="is.bed"
    output:
        temp("{sample}.doc")
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "samtools depth -b {input.bed} {input.bam} | "
        "sort -S 1G -T . -k 1,1 -k 2,2n --parallel={threads} > {output}"
