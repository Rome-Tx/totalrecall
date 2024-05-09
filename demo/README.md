##Example data
Copy `case.bam`, `case.bam.bai`,  `control.bam`,  `control.bam.bai`
to some directory

## Expected results
`line.vcf` - expected results

## Running
Run using Docker from the directory input files were copied to:
```
LC_ALL=C snakemake -j 1 -s /opt/totalrecall/Snakefile -p results.tgz
```
