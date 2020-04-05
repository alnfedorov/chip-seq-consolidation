# Global options for fetching and processing encode data

NARROW_HISTONE_MARKS = [
    "H2AFZ",
    "H3ac",
    "H3K27ac",
    "H3K4me2",
    "H3K4me3",
    "H3K9ac"
]

BROAD_HISTONE_MARKS = [
    "H3F3A",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K79me2",
    "H3K79me3",
    "H3K9me1",
    "H3K9me2",
    "H3K9me3",
    "H4K20me1"
]

ORIGINAL_DIR = "original"

SOFT_FDR_CUTOFF = 0.1
IDR_CUTOFF = 0.05
CONSENSUS_FILENAME = "consensus.bed"

ASSEMBLY = 'hg19'
GENOME_FASTA = "/data/hg19/hg19.fasta"
BWA_GENOME_INDEX = "/data/hg19/bwa-index/hg19.fasta"

# See https://github.com/Boyle-Lab/Blacklist/blob/v2.0/lists/hg19-blacklist.v2.bed.gz
BLACKLISTED_REGIONS = "/data/hg19/hg19-blacklist.v2.bed.gz"

LOGS_DIR = "logs"
PEAK_CALLING_DIR = "peak-calling"
BAM_QC_DIR = "qc"

UNFILTERED_READS_PREFIX = 'unfiltered'
DUPLICATED_READS_PREFIX = 'duplicates'
UNIQUE_READS_PREFIX = 'unique'

UNFILTERED_READS_DIR = f'{UNFILTERED_READS_PREFIX}-bam'
DUPLICATED_READS_DIR = f'{DUPLICATED_READS_PREFIX}-bam'
UNIQUE_READS_DIR = f'{UNIQUE_READS_PREFIX}-bam'

SIMULATED_DIR = "simulated"
CHIPS_DIR = "chips"
SUBSAMPLE_DIR = "subsample"
CHIPS_MODELS_DIR = "models"

SIMULATED_READS_DIRS = ["q0.25", "q0.5", "q0.75"]
READS_TO_SIMULATE = {
    hm: {"q0.25": 14705741, "q0.5": 22154364, "q0.75": 31509056} for hm in NARROW_HISTONE_MARKS
}
READS_TO_SIMULATE["H3K27me3"] = {"q0.25": 14548131, "q0.5": 23287822, "q0.75": 41352521}
assert all(
    len(dispatch) == len(SIMULATED_READS_DIRS) and all(k in SIMULATED_READS_DIRS for k in dispatch.keys())
    for dispatch in READS_TO_SIMULATE.values()
)
