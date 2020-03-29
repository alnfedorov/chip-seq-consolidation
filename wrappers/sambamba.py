import os
import tempfile
import logging
from .run import run, _move


logger = logging.getLogger(__name__)
FILTER_KEEP_DUPS = "mapping_quality >= 1 " \
                   "and not (unmapped or secondary_alignment) " \
                   "and not ([SA] != null)"
FILTER_KEEP_UNIQUE = "mapping_quality >= 1 " \
                     "and not (unmapped or secondary_alignment or duplicate) " \
                     "and not ([SA] != null)"


def fromsam(path: str, threads: int = 1, saveto: str = None):
    assert os.path.isfile(path)
    assert threads >= 1

    saveto = saveto if saveto else tempfile.mktemp(dir=os.path.dirname(path))
    run(
        [
            "sambamba", "view", "--sam-input", "--with-header", "--show-progress", "--compression-level=9",
            f"--nthreads={threads}", "--format=bam", f"--output-filename={saveto}", path
        ], logger, logbefore=f"Start sambamba SAM->BAM transform for {path}", logafter="Transform finished"
    )
    return saveto


def sort(path: str, threads: int = 1, saveto: str = None, inplace: bool = True):
    assert os.path.isfile(path)
    assert threads > 0

    # sambamba sort -t 12 -l 9 --show-progress -o OUT.bam IN.bam
    saveto = saveto if saveto and not inplace else tempfile.mktemp(dir=os.path.dirname(path))
    run(
        [
            "sambamba", "sort", f"--nthreads={threads}", "--show-progress", "--compression-level=9",
            f"--out={saveto}", path
        ], logger, logbefore=f"Start sambamba sort for {path}", logafter="Sort finished"
    )
    if inplace:
        saveto = _move(saveto, path)
    return saveto


def markdup(path: str, threads: int = 1, saveto: str = None, inplace: bool = True):
    assert os.path.isfile(path)
    assert threads > 0

    # sambamba markdup -t 12 -l 9 --show-progress IN.bam OUT.bam
    saveto = saveto if saveto and not inplace else tempfile.mktemp(dir=os.path.dirname(path))
    run(
        [
            "sambamba", "markdup", f"--nthreads={threads}", "--show-progress", "--compression-level=9",
            path, saveto
        ], logger, logbefore=f"Start sambamba markdup for {path}", logafter="markdup finished"
    )
    if inplace:
        saveto = _move(saveto, path)
    return saveto


def filter(path: str, rule: str, threads: int = 1, saveto: str = None, inplace: bool = False):
    assert os.path.isfile(path)
    assert threads > 0

    # sambamba view -t 12 -h -f bam -F "" -o ENCFF496QVI-no-duplicates.bam ENCFF496QVI-with-duplicates.bam
    saveto = saveto if saveto and not inplace else tempfile.mktemp(dir=os.path.dirname(path))
    run(
        [
            "sambamba", "view", "--with-header", "--show-progress", "--compression-level=9",
            f"--nthreads={threads}", f"--filter={rule}", "--format=bam", f"--output-filename={saveto}", path
        ], logger, logbefore=f"Start sambamba view for {path} with rule {rule}", logafter="view finished"
    )

    if inplace:
        saveto = _move(saveto, path)
    return saveto



# Filter multimappers and other dummies
# sambamba view -t 12 -h --show-progress -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([SA] != null)" -o ENCFF074TRC.sorted.markdup.filtered.bams ENCFF074TRC.sorted.markdup.bam

# Pipeline:
# 1. Align reads
# (single end)
# "bowtie2 --mm --threads 12 -x /data/hg19/bowtie2-index/hg19 -U ENCFF900SUO.fastq.gz -S ENCFF900SUO.sam"

# 2. SAM -> BAM
# "sambamba view --sam-input --with-header --show-progress --compression-level=9 --nthreads=12 --format=bam --output-filename=OUT.bam IN.sam"

# 3. Sort and markdup bam files
# sambamba sort -t 12 -l 9 --show-progress -o OUT.bam IN.bam
# sambamba markdup -t 12 -l 9 --show-progress IN.bam OUT.bam

# 4. filter out dummies and [keep duplicates, filter duplicates]
# sambamba view --nthreads=12 --header --show-progress --format=bam -F "[XS] == null and not unmapped  and not duplicate" --output-filename=ENCFF074TRC.sorted.markdup.filtered.bams ENCFF074TRC.sorted.markdup.bam
#
# sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([SA] != null)" -o ENCFF496QVI-no-duplicates.bam ENCFF496QVI-with-duplicates.bam

