import os
import asyncio
from wrappers import sambamba, picard, samtools, fastqc


class filter:
    @staticmethod
    async def toduplicates(bam: str, saveto: str, maxthreads: int) -> str:
        duplicates = saveto
        assert await sambamba.filter(bam, sambamba.FILTER_KEEP_DUPS, maxthreads, saveto=duplicates) == duplicates
        assert await sambamba.sort(duplicates, threads=maxthreads, inplace=True) == duplicates
        assert await picard.MarkDuplicates(duplicates, inplace=True) == duplicates
        return duplicates

    @staticmethod
    async def tounique(bam: str, saveto: str, maxthreads: int) -> str:
        unique = saveto
        assert await sambamba.filter(bam, sambamba.FILTER_KEEP_UNIQUE, maxthreads, saveto=unique) == unique
        assert await sambamba.sort(unique, threads=maxthreads, inplace=True) == unique
        return saveto


async def qc(bam: str, saveto: str, doflagstats: bool = True, dostats: bool = True, dofastqc: bool = False):
    """
    Run QC pipeline for the given BAM file
    :param bam: path to the BAM file
    :param saveto: where to save qc results
    :param doflagstats: run samtools.flagstats or not
    :param dostats: run samtools.stats or not
    :param dofastqc: run fastqc or not
    :return: None
    """
    assert os.path.isdir(saveto) and bam.endswith("bam")
    filename = os.path.split(bam)[-1]
    coro = []
    if doflagstats:
        coro.append(samtools.flagstat(bam, os.path.join(saveto, filename.replace(".bam", ".samtools-flagstat"))))
    if dostats:
        coro.append(samtools.stats(bam, os.path.join(saveto, filename.replace(".bam", ".samtools-stats"))))
    if dofastqc:
        folder = os.path.join(saveto, filename.replace(".bam", "-fastqc"))
        os.mkdir(folder)
        coro.append(fastqc.fastqc(bam, folder))
    await asyncio.gather(*coro)


__all__ = [filter, qc]
