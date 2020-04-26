import os
import shutil
import tempfile
import logging
from .utils import run
from .piping import pipe

logger = logging.getLogger(__name__)


@pipe(writearg=("saveto", "p"), file="p")
async def faidx(file: str, saveto: str = None) -> str:
    assert os.path.exists(file)
    await run(["samtools", "faidx", file],
              logbefore=f"samtools faidx {file}", logafter="fasta index building is finished")
    if saveto:
        shutil.move(saveto, file + ".fai")
    else:
        saveto = file + ".fai"
    return saveto


@pipe(writearg=("saveto", "p"), files="p")
async def merge(files: [str], threads: int = 1, saveto: str = None) -> str:
    assert threads > 0 and all(os.path.exists(f) for f in files)
    if len(files) == 1 and saveto is None:
        return files[0]
    saveto = saveto if saveto else tempfile.mkstemp()[1]

    await run(["samtools", "merge", "-f", f"--threads={threads}", saveto, *files], logger,
              logbefore=f"start samtools merge for {files}, saved in {saveto}", logafter="samtools merge finished")
    assert os.path.exists(saveto)
    return saveto


# samtools stats -> statistics about alignment etc
@pipe(writearg=("saveto", "p"), file='p')
async def stats(file: str, saveto: str = None) -> str:
    assert os.path.exists(file)

    saveto = saveto if saveto else tempfile.mkstemp()[1]
    stats = (await run(
        ["samtools", "stats", file],
        logger, logbefore=f"Start samtools stats for {file}", logafter="samtools stats finished"
    )).stdout.decode()

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto


# samtools flagstat in.bam -> simple stats mapped/unmapped/passed etc
@pipe(writearg=("saveto", "p"), file='p')
async def flagstat(file: str, saveto: str = None) -> str:
    assert os.path.exists(file)
    saveto = saveto if saveto else tempfile.mkstemp()[1]

    stats = (await run(
        ["samtools", "flagstat", file],
        logger, logbefore=f"Start samtools flagstat for {file}", logafter="samtools flagstat finished"
    )).stdout.decode()

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto
