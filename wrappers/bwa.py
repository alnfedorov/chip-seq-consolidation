import os
import tempfile
import logging
from .utils import run

logger = logging.getLogger(__name__)


async def aln(fastq: str, index: str, saveto: str = None, quality_threshold: int = 5, seedlen: int = 32,
              max_seeddiff: int = 2, threads: int = 1) -> str:
    assert os.path.isfile(fastq) and seedlen > 0 and quality_threshold > 0 and \
           seedlen > 0 and max_seeddiff > 0 and threads >= 1
    saveto = saveto if saveto is None else tempfile.mkstemp()[1]
    cmd = ["bwa", "aln", "-q", quality_threshold, "-l", seedlen, "-k", max_seeddiff, "-t", threads, index, fastq]
    cmd = [str(x) for x in cmd]
    run(cmd, logger=logger, logbefore=f"Running bwa aln with cmd {' '.join(cmd)}", logafter="Alignment finished")
    return saveto


# ========================================
# Map reads to create raw SAM file
# ========================================
# SAI_FILE_1="${OFPREFIX}.sai"
# RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
# RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" # To be stored
# RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File
#
# module add bwa/0.7.13
# module add samtools/1.7
#
# bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} > ${SAI_FILE_1}
#
# bwa samse ${BWA_INDEX_NAME} ${SAI_FILE_1} ${FASTQ_FILE_1} | samtools view -Su - | samtools sort -o ${RAW_BAM_FILE} -
#
# rm ${SAI_FILE_1}
#
# # Use bwa-mem for reads >= 70 bp
# # bwa mem -M -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} | samtools sort -o ${RAW_BAM_FILE} -
#
# samtools sort -n --threads 10 ${RAW_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${RAW_BAM_FILE_MAPSTATS}
