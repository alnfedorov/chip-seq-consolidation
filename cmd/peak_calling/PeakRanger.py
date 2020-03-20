from ..run import run


def bcp():
    pass

# duplicates are ignored by default
def ranger(data: [str], control: [str], saveto: str,
           threads: int = 1, pcutoff: float = 0.0001, fdrcutoff: float = 0.01, format: str = "bam"):
    """           peakranger ranger --data=/data/encode/H3K4me3/ENCSR661MUS/data/duplicates-bam/duplicates-breplica-ENCFF755PXB.bam \
                  --control=/data/encode/H3K4me3/ENCSR661MUS/data/duplicates-bam/duplicates-control-ENCFF004LNU.bam \
                  --output=ranger.bed \
                  --format=bam \
                  --pval=0.001 \
                  --FDR=0.05 \
                  --verbose \
                  --thread=12"""
    pass


def bcp():
    """peakranger bcp --data=/data/encode/H3K4me3/ENCSR661MUS/data/unique-bam/unique-breplica-ENCFF755PXB.bam \
                  --control=/data/encode/H3K4me3/ENCSR661MUS/data/unique-bam/unique-control-ENCFF004LNU.bam \
                  --format=bam \
                  --output=bcp.bed \
                  --pval=0.001
"""
    pass