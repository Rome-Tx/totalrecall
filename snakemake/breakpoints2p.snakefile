include: "variables.snakefile"

rule breakpoints_2nd_pass:
    input:
        clip="{sample}.clip_subset.tsv",
        mask="{sample}.mask.tsv"
    output:
        bp=temp("{sample}.breakpoints.byid.tsv"),
        fa=temp("{sample}.clipped.fa")
    shell:
        "compute_breakpoints.py -d {MAX_PRECLUSTER_DIST} -l 1 "
        "-n 1 -q 0 -m {input.mask} {input.clip} {output.bp} {output.fa}"

rule subset_clipped_reads:
    input:
        clip="{sample}.clip.tsv",
        loci="{sample}.1p.loci.tsv"
    output:
        temp("{sample}.clip_subset.tsv")
    shell:
        "subset_clipped_reads.py {input.clip} {input.loci} -e {SEARCH_EXT_1P} "
        "> {output}"

rule loci_1st_pass:
    input:
        "{sample}.1p.other.tsv",
        "{sample}.1p.polya.tsv"
    output:
        temp("{sample}.1p.loci.tsv")
    shell:
        "cat {input} | sort -k 1,1 -k 2,2n -u -S 500M "
        "-T . > {output}"
