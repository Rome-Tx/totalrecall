# HTSR class
HTSR_FIELD_CONTIG = 4
HTSR_FIELD_LOCUS = 13
HTSR_FIELD_BREAKPOINT = 14


# === workflow parameters === #
MIN_NR_1P = 2
MIN_NR_2P = 2
MAX_SITE_LEN = 50
SEARCH_EXT_1P = MAX_SITE_LEN + 20
MAX_PRECLUSTER_DIST = 5
FILTER_BP_EXT = 1
FILTER_READ_EXT = 5
MIN_CONTROL_READ_LEN = 3
FILTER_CLUSTER_CLIP_LEN = 10
MIN_COVERAGE_FRACTION = 0.05
BED_EXT_LEN = 250
IGV_EXT_LEN = 100

MIN_NR_FILTER = config.get("min_nr_final", 3)
MIN_LEN_FILTER = config.get("min_len_final", 50)

IGVDIR = config.get("igvdir", "/opt/igv")

# === data paths and input files === #
IXDIR = config.get("ixdir", "/opt/totalrecall")
LASTDB_1P = os.path.join(IXDIR, "last/first_pass")

CASE_ALN = config.get("case_alignment", "case.bam")
CASE_INDEX = config.get("case_index", "case.bam.bai")
CONTROL_ALN = config.get("control_alignment", "control.bam")
CONTROL_INDEX = config.get("control_index", "control.bam.bai")

# TargetSite class
TS_FIELD_NAME = 1
TS_FIELD_CH = 2
TS_FIELD_LBP = 3
TS_FIELD_RBP = 4
TS_FIELD_LSEQ = 9
TS_FIELD_RSEQ = 13

# === alignment parameters === #
LAST_A=30
LAST_B=20
LAST_R=10
LAST_Q=20
LAST_D_1p=60
LAST_E_1p=150
LAST_D=70
LAST_E=120

BAM_SORT_MB_PERCORE = config.get('bam_sort_mb', 4096)
