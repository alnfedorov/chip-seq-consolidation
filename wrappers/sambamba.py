import os
import shutil
import tempfile
import logging
from .utils import run, replace_bam
from .piping import pipe
from typing import Awaitable


logger = logging.getLogger(__name__)
FILTER_KEEP_DUPS = "mapping_quality >= 1 " \
                   "and not (supplementary or unmapped or secondary_alignment or failed_quality_control or chimeric) " \
                   "and not ([SA] != null)"
FILTER_KEEP_UNIQUE = FILTER_KEEP_DUPS + " and not duplicate"

PAIRED_FILTER_KEEP_DUPS = FILTER_KEEP_DUPS + " and proper_pair"
PAIRED_FILTER_KEEP_UNIQUE = FILTER_KEEP_UNIQUE + " and proper_pair"


@pipe(writearg=("saveto", "p"), bam="p")
async def subsample(bam: str, fraction: float = None, reads: int = None, saveto: str = None,
                    threads: int = 1) -> str:
    assert os.path.exists(bam) and threads >= 1
    if fraction is not None:
        assert reads is None and 0 < fraction < 1
    else:
        assert fraction is None
        fraction = reads / (await numreads(bam, threads))
        assert 0 < fraction < 1
    saveto = saveto if saveto else tempfile.mkstemp()[1]
    cmd = ["sambamba", "view", "-h", f"--nthreads={threads}", f"--subsample={fraction}",
           "--format=bam", "--subsampling-seed=123", f"--output-filename={saveto}", bam]
    await run(cmd, logger, f"subsample reads with cmd: {' '.join(cmd)}", "subsampling finished")
    return saveto


async def numreads(bam: str, threads: int = 1) -> int:
    assert os.path.isfile(bam) and threads >= 1
    cmd = ["sambamba", "view", "-c", f"--nthreads={threads}", bam]
    result = await run(cmd, logger, f"count reads with cmd: {' '.join(cmd)}", "numreads finished")
    return int(result.stdout.decode())


@pipe(writearg=("saveto", "p"), bam="p")
async def index(bam: str, threads: int = 1, saveto: str = None) -> str:
    assert os.path.exists(bam)
    assert threads >= 1
    cmd = ["sambamba", "index", "-c", "-t", str(threads), bam]
    await run(cmd, logger, logbefore=f"Indexing file {bam}", logafter=f"finished indexing")
    if saveto is not None:
        shutil.move(bam + ".bai", saveto)
    else:
        saveto = bam + ".bai"
    return saveto


@pipe(writearg=("saveto", "p"), sam="p")
async def fromsam(sam: str, threads: int = 1, saveto: str = None) -> str:
    assert os.path.exists(sam)
    assert threads >= 1

    saveto = saveto if saveto else tempfile.mkstemp()[1]
    await run(
        [
            "sambamba", "view", "--sam-input", "--with-header", "--show-progress", "--compression-level=9",
            f"--nthreads={threads}", "--format=bam", f"--output-filename={saveto}", sam
        ], logger, logbefore=f"Start sambamba SAM->BAM transform for {sam}", logafter="Transform finished"
    )
    return saveto


@pipe(writearg=("saveto", "p"), path="p")
async def sort(path: str, saveto: str = None, threads: int = 1,
               byname: bool = False, match_mates: bool = None, filter: str = None) -> str:
    assert os.path.exists(path)
    assert threads > 0

    saveto = saveto if saveto else tempfile.mkstemp()[1]
    # saveto = saveto if saveto and not inplace else tempfile.mkstemp()[1]
    cmd = ["sambamba", "sort", f"--nthreads={threads}", "--show-progress", "--compression-level=9", f"--out={saveto}"]
    if byname:
        cmd.append("--sort-by-name")
    if match_mates:
        cmd.append("--match-mates")
    if filter:
        cmd.append(f'--filter={filter}')
    cmd.append(path)
    await run(cmd, logger, logbefore=f"Start sambamba sort for {path}", logafter="Sort finished")
    # if inplace:
    #     saveto = replace_bam(saveto, path)
    return saveto


@pipe(writearg=("saveto", "p"), path="p")
async def markdup(path: str, threads: int = 1, saveto: str = None) -> str:
    assert os.path.exists(path)
    assert threads > 0
    saveto = saveto if saveto else tempfile.mkstemp()[1]
    await run(
        [
            "sambamba", "markdup", f"--nthreads={threads}", "--show-progress", "--compression-level=9",
            path, saveto
        ], logger, logbefore=f"Start sambamba markdup for {path}", logafter="markdup finished"
    )
    # if inplace:
    #     saveto = replace_bam(saveto, path)
    return saveto


@pipe(writearg=("saveto", "p"), path="p")
async def filter(path: str, rule: str, threads: int = 1, saveto: str = None) -> str:
    assert os.path.exists(path)
    assert threads > 0

    # saveto = saveto if saveto and not inplace else tempfile.mkstemp()[1]
    saveto = saveto if saveto else tempfile.mkstemp()[1]
    await run(
        [
            "sambamba", "view", "--with-header", "--show-progress", "--compression-level=9",
            f"--nthreads={threads}", f'--filter={rule}', "--format=bam", f"--output-filename={saveto}", path
        ], logger, logbefore=f"Start sambamba view for {path} with rule {rule}", logafter="view finished"
    )

    # if inplace:
    #     saveto = replace_bam(saveto, path)
    return saveto
