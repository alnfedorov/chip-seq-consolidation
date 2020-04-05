import os
import shutil
import tempfile
import logging
import asyncio
from .utils import run

logger = logging.getLogger(__name__)


async def faidx(file: str, saveto: str = None):
    assert os.path.exists(file)
    await run(["samtools", "faidx", file],
              logbefore=f"samtools faidx {file}", logafter="fasta index building is finished")
    if saveto:
        shutil.move(saveto, file + ".fai")
    else:
        saveto = file + ".fai"
    return saveto


async def merge(files: [str], threads: int = 1, saveto: str = None):
    assert threads > 0 and all(os.path.exists(f) for f in files)
    if len(files) == 1:
        return files[0]
    saveto = saveto if saveto else tempfile.mkstemp()[1]

    await run(["samtools", "merge", "-f", f"--threads={threads}", saveto, *files], logger,
              logbefore=f"start samtools merge for {files}, saved in {saveto}", logafter="samtools merge finished")
    assert os.path.exists(saveto)
    return saveto


# samtools flagstat in.bam -> simple stats mapped/unmapped/passed etc
# samtools stats -> statistics about alignment etc
async def stats(path: str, saveto: str = None):
    assert os.path.isfile(path) and ".bam" in path

    saveto = saveto if saveto else tempfile.mkstemp()[1]
    stats = (await run(
        ["samtools", "stats", path],
        logger, logbefore=f"Start samtools stats for {path}", logafter="samtools stats finished"
    )).stdout.decode()

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto


async def flagstat(path: str, saveto: str = None):
    assert os.path.isfile(path) and ".bam" in path

    saveto = saveto if saveto else tempfile.mkstemp()[1]

    stats = (await run(
        ["samtools", "flagstat", path],
        logger, logbefore=f"Start samtools flagstat for {path}", logafter="samtools flagstat finished"
    )).stdout.decode()

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto
