# Single end algorithm
# # ========================================
# # Map reads to create raw SAM file
# # ========================================
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
# ==============================================================
# ==============================================================
# ==============================================================
# ==============================================================
# ==============================================================
# Paired end algorithm
# SAI_FILE_1="${OFPREFIX}_1.sai"
# SAI_FILE_2="${OFPREFIX}_2.sai"
# RAW_SAM_FILE="${OFPREFIX}.raw.sam.gz"
#
# bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} > ${SAI_FILE_1}
#
# bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_2} > ${SAI_FILE_2}
#
# bwa sampe ${BWA_INDEX_NAME} ${SAI_FILE_1} ${SAI_FILE_2} ${FASTQ_FILE_1} ${FASTQ_FILE_2} | gzip -nc > ${RAW_SAM_FILE}
#
# rm ${SAI_FILE_1} ${SAI_FILE_2}
#
# # Use bwa-mem for reads >= 70 bp
# # bwa mem -M -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} ${FASTQ_FILE_2} | gzip -nc > ${RAW_SAM_FILE}
#
# # ==============================================================
# # Remove read pairs with bad CIGAR strings and sort by position
# # ==============================================================
# RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
# RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" # To be stored
# BADCIGAR_FILE="${TMP}/badReads${RANDOM}.tmp"
# RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File
#
# # Find bad CIGAR read names
# zcat ${RAW_SAM_FILE} | awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && $6!="*" { cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) print $1”\t” ; }' | sort | uniq > ${BADCIGAR_FILE}
#
# # Remove bad CIGAR read pairs
# if [[ $(cat ${BADCIGAR_FILE} | wc -l) -gt 0 ]]
# then
#     zcat ${RAW_SAM_FILE} | grep -v -F -f ${BADCIGAR_FILE} | samtools view -Su - | samtools sort - ${RAW_BAM_PREFIX}
# else
#     samtools view -Su ${RAW_SAM_FILE} | samtools sort - ${RAW_BAM_PREFIX}
# fi
#
# rm ${BADCIGAR_FILE} ${RAW_SAM_FILE}
#
# samtools sort -n --threads 10 ${RAW_BAM_FILE} -O SAM | SAMstats --sorted_sam_file -  --outf ${RAW_BAM_FILE_MAPSTATS}
# samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}
import os
import tempfile
from wrappers import bwa, samtools, sambamba


class align:
    """Note: current implementation is not optimal, because bash pipes(or they alternatives) are not used."""
    @staticmethod
    async def send(bwaindex: str, fastq: str, readlen: int, threads: int = 1, saveto: str = None):
        """Align single end fastq files"""
        assert readlen > 26 and os.path.exists(fastq) and threads >= 1
        saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
        if readlen < 70:    # bwa aln/samse algorithm
            sai = await bwa.aln(fastq, bwaindex, threads=threads)
            sam = await bwa.samse(bwaindex, sai, fastq)
            await sambamba.fromsam(sam, threads, saveto)
            os.remove(sai)
            os.remove(sam)
        else:   # bwa mem algorithm
            raise NotImplementedError()

    @staticmethod
    async def pend():
        """Align paired-end fastq files"""
        raise NotImplementedError()
