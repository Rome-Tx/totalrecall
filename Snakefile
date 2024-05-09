## Snakemake workflow for transposon identification from WGS
## 

include: "snakemake/preprocess.snakefile"
include: "snakemake/drp.snakefile"
include: "snakemake/lastal.snakefile"
include: "snakemake/breakpoints.snakefile"
include: "snakemake/breakpoints2p.snakefile"
include: "snakemake/penultimate.snakefile"
include: "snakemake/depth.snakefile"
include: "snakemake/final.snakefile"
