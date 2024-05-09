include: "variables.snakefile"

rule align_1st_pass_fastq:
    input:
        "{sample}.1p.fq"
    output:
        temp("{sample}.1p.other.tsv")
    threads: 50
    shell:
    # TODO: explain what is done here
        "lastal -P {threads} -Q fastx -r {LAST_R} -q {LAST_Q} -a {LAST_A} "
        "-b {LAST_B} -d {LAST_D_1p} -e {LAST_E_1p} -l 5 -m 1000 -s 1 -K 0 "
        "-f BlastTab {LASTDB_1P} {input} | awk '! /^#/' | cut -f1 | "
        "cut -d : -f 2- | rev | sed 's/:/\t/' | rev > {output}"

rule clipped_tail_te_ref_alignment:
    input:
        "{sample}.clipped.fa"
    output:
        temp("{sample}.clip_annot.tsv")
    threads: 50
    shell:
        "lastal -P {threads} -r {LAST_R} -q {LAST_Q} -a {LAST_A} "
        "-b {LAST_B} -d {LAST_D} -e {LAST_E} -l 5 -m 1000 -s 1 -K 0 "
        """-f BlastTab {IXDIR}/last/signatures {input} | grep -v "^#" | """
        "sort -S 500M --parallel={threads} -T . -k 1,1n | "
        "blastTab2annot.py {IXDIR}/signatures.fa > {output}"
