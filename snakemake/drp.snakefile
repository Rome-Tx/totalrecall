include: "variables.snakefile"

# ==== repeat annotation for discordant reads ==== #

rule drp_names:
    input:
        bam="{sample}.disc.bam",
        bed="drp.bed"
    output:
        temp("{sample}.drp_names")
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "samtools view -@ {threads} -L {input.bed} {input.bam} | awk '{{print $1}}' | "
        "sort -T . -u --parallel={threads} -S 500M > {output}"

rule drp_subset_bam:
    input:
        bam="{sample}.disc.bam",
        names="{sample}.drp_names"
    output:
        temp("{sample}.drp.subset.bam")
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "samtools view -@ {threads} --qname-file {input.names} "
        "-o {output} {input.bam}"

rule discordant_reads_fastq:
    input:
        "{sample}.drp.subset.bam"
    output:
        temp("{sample}.disc.fq")
    params:
        basename="{sample}"
    threads: max(1, min(workflow.cores / 2, 4))
    shell:
        "samtools sort -@ {threads} -u -m {BAM_SORT_MB_PERCORE}M -T "
        "./tmp.{params.basename}_ -n {input} | "
        "samtools fastq -o /dev/stdout -N -0 /dev/null -s /dev/null - | "
        "paste - - - - | sort -T . -t / -k 1,1 --parallel={threads} -S 500M | "
        """tr "\t" "\n" > {output}"""

rule discordant_reads_repeat_mapping:
    input:
        "{sample}.disc.fq"
    output:
        pre=temp("{sample}.repeatinfo.pre.tsv"),
        final=temp("{sample}.repeatinfo.tsv")
    params:
        basename="{sample}"
    threads: 50
    shell:
    # TODO: explain what BlastTab+ is, what field numbers are...
        "lastal -Q sanger -P {threads} -f BlastTab+  {IXDIR}/last/transposon "
        "-r {LAST_R} -q {LAST_Q} -a {LAST_A} -b {LAST_B} {input} | "
        # remove comments
        "awk '!/^#/' | "
        # add strand (+/-)
        """awk '{{if($7>$8){{s="-"}}else{{s="+"}}; print $0 "\t" s}}' | """
        # sort by read name, strand, score (decreasing)
        "sort -k 1,1 -k 16,16 -k 15,15gr --parallel={threads} -S 500M -T . | "
        # print the best alignment for each read/strand; cut the fields:
        # 1-4: query id, subject id, pident, alignment length
        # 9-10: subject start, subject end
        # 13-16: query length, subject length, raw score, strand
        "awk 'a !~ $1 $NF; {{a=$1 $NF}}' | cut -f 1-4,9-10,13-16 | "
        # replace /1 and /2 in read names with a separate field
        """awk 'BEGIN{{OFS="\t"}};
        {{$1=gensub("/([12])$", "\t\\\\1", "g", $1); print}}' > {output.pre};"""
        # use {output.pre} to avoid running to sort commands in the same dir
        # sort and write the final output
        "sort -k 1,1 --parallel={threads} -S 500M -T . {output.pre} "
        "> {output.final}"
