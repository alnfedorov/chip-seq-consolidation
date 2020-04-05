import os
import tempfile
import logging
from .utils import run

logger = logging.getLogger(__name__)


async def aln(fastq: str, index: str, saveto: str = None, quality_threshold: int = 5, seedlen: int = 32,
              max_seeddiff: int = 2, threads: int = 1) -> str:
    assert os.path.isfile(fastq) and seedlen > 0 and quality_threshold > 0 and \
           seedlen > 0 and max_seeddiff > 0 and threads >= 1
    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    cmd = [
        "bwa", "aln", "-q", quality_threshold, "-l", seedlen, "-k", max_seeddiff, "-t", threads, index, fastq
    ]
    cmd = [str(x) for x in cmd]
    with open(saveto, 'w') as file:
        await run(
            cmd, logger, f"Running bwa aln with cmd {' '.join(cmd)}", "Alignment finished",
            logstdout=False, stdout=file
        )
    return saveto


async def samse(index: str, sai: str, fastq: str, saveto: str = None):
    """Single end aligner, should be used on top of the aln command"""
    assert os.path.isfile(sai) and os.path.isfile(fastq)
    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    cmd = ["bwa", "samse", index, sai, fastq]
    with open(saveto, 'w') as file:
        await run(
            cmd, logger, f"Running bwa samse with cmd {' '.join(cmd)}", "Alignment finished",
            logstdout=False, stdout=file
        )
    return saveto


async def sampe(index: str, sai1: str, sai2: str, fastq1: str, fastq2: str, saveto: str = None):
    """Paired end aligner, should be used on top of the aln command"""
    assert os.path.isfile(sai1) and os.path.isfile(sai2) and os.path.isfile(fastq1) and os.path.isfile(fastq2)
    saveto = saveto if saveto is None else tempfile.mkstemp()[1]
    raise NotImplementedError()
    cmd = ["bwa", "sampe", index, sai1, sai2, fastq1, fastq2, saveto]
    await run(cmd, logger=logger, logbefore=f"Running bwa samse with cmd {' '.join(cmd)}",
              logafter="Alignment finished")
    return saveto


async def mem(index, fastq1: str, fastq2: str = None, threads: int = 1):
    assert os.path.exists(fastq1) and fastq2 is None or os.path.exists(fastq2)
    assert threads >= 1
    raise NotImplementedError()

__all__ = [aln, sampe, samse, mem]
