include: "variables.snakefile"

# ==== assemble the results ==== #

## create archive with the results
rule results_tgz:
    input:
        "case.header.sam",
        "control.header.sam",
        "case.tlen",
        "case.nbases",
        "case.rlen",
        "control.tlen",
        "control.nbases",
        "control.rlen",
        "aux/case.pe.L.tsv",
        "aux/case.pe.R.tsv",
        "screenshots.tgz",
        "all.vcf.gz",
        "somatic.vcf.gz",
        "line_hc.vcf",
        "line.vcf"
    output:
        "results.tgz"
    shell:
        "tar -czf {output} {input}"

# ==== annotate/filter integration sites and take screenshots ==== #

## run IGV in batch mode to take screenshots using xvfb X server emulator
## possible IGV failure errors is silenced to allow the workflow to continue
rule screenshots_tgz:
    input:
        batch="igv.batch",
        bed="is.annot.bed",
        case="case.bam",
        case_bai="case.bam.bai",
        control="control.bam",
        control_bai="control.bam.bai"
    output:
        "screenshots.tgz"
    shell:
         """timeout 5h xvfb-run --auto-servernum """
         """--server-args="-screen 0 1920x1080x24" {IGVDIR}/igv.sh """
         """--igvDirectory $(pwd) -b {input.batch} || mkdir -p screenshots; """
         """tar -czf {output} screenshots"""

## write IGV batch files for going to each identified insertion
## and taking a screenshot
rule igv_batch_file:
    input:
        tsv="somatic.vcf.gz",
        genome="case.genome"
    output:
        "igv.batch"
    shell:
        "mkigvbatch.py {input.tsv} {input.genome} {IXDIR}/genome.fa.gz "
        "-e {IGV_EXT_LEN} > {output}"

## create and index a BED file for identified insertion to show
## on the screenshots; two ends of the BED interval correspond to
## the two breakpoints
rule igv_bed_file:
    input:
        vcf_gz="somatic.vcf.gz",
        genome="case.genome"
    output:
        bed="{infile}.bed",
        idx="{infile}.bed.idx"
    shell:
        "write_isites_bed.py -v {input.genome} {input.vcf_gz} | "
        "sort -T . -k 1,1 -k 2,2n > {output.bed}; "
        "{IGVDIR}/igvtools index {output.bed}"

## classifier 
rule classifier:
    input:
        "somatic.vcf.gz",
        "aux/case.pe.L.tsv",
        "aux/case.pe.R.tsv"
    output:
        "line.vcf"
    shell:
        """classifier_pe.py {IXDIR}; sed -i '/CLIPPED.CLUSTERS/d' line.vcf"""

## high confidence LINEs
rule line_hc_vcf:
    input:
        "somatic.vcf.gz",
    output:
        "line_hc.vcf"
    shell:
        "zcat {input} | filter_vcf.py -l {MIN_LEN_FILTER} -L "
        "-n {MIN_NR_FILTER} > {output}"""

## remove germline variants
rule somatic_vcf:
    input:
        "all.vcf.gz"
    output:
        "somatic.vcf.gz"
    shell:
        """zcat {input} | awk '$7 !~ "CONTROL.CLIPPED.READ"' | """
        """gzip -c > {output}"""

## classify insertions based on the clipped reads at their ends
## and DRPs near them, apply the relevant filters
rule annotate_filter_insertions:
    input:
        "case.ins.namesorted.tsv",
        "case.header.sam",
        "case.nbases",
        "control.nbases",
        "case.tlen",
        "control.tlen",
        "aux/control.clipped.read",
        "aux/case.read.cluster",
        "aux/control.read.cluster",
        "aux/doc_case.tsv",
        "aux/doc_control.tsv",
        "aux/case.low_complexity",
        "aux/case.pe.L.tsv",
        "aux/case.pe.R.tsv"
    output:
        "all.vcf.gz"
    shell:
        "annotate_insertions.py -c {MIN_COVERAGE_FRACTION} "
        "{IXDIR}/genome.fa.gz | gzip > {output}"
