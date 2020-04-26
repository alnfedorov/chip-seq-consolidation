import os
import asyncio
from typing import Awaitable
from multiprocessing import cpu_count
from wrappers.peak_calling import macs2, epic2, PePr, PeakRanger, idr
from pybedtools import BedTool, Interval
from utils.bed import compute_conservative_regions
from .meta import ExperimentMeta
from .config import SOFT_FDR_CUTOFF, IDR_CUTOFF, CONSENSUS_FILENAME, BLACKLISTED_REGIONS, BROAD_HISTONE_MARKS


def is_broad(hm: str):
    return hm in BROAD_HISTONE_MARKS


async def call_peaks(pcalling: str, experiments: [ExperimentMeta], maxthreads: int, force: bool):
    assert len(experiments) == 2 and experiments[0].target == experiments[1].target
    consensus = os.path.join(pcalling, CONSENSUS_FILENAME)
    maxthreads = max(1, maxthreads // 3)
    coro = [asyncio.create_task(_run_pepr(experiments, pcalling, IDR_CUTOFF, maxthreads, force))]

    # macs2
    coro.append(
        _run_idr([
            _run_macs2(experiments[0], pcalling, SOFT_FDR_CUTOFF, force),
            _run_macs2(experiments[1], pcalling, SOFT_FDR_CUTOFF, force)
        ], pcalling, rankby="p.value", idrcutoff=IDR_CUTOFF, force=force)
    )

    # epic2
    coro.append(
        _run_idr([
            _run_epic2(experiments[0], pcalling, SOFT_FDR_CUTOFF, maxthreads, force),
            _run_epic2(experiments[1], pcalling, SOFT_FDR_CUTOFF, maxthreads, force)
        ], pcalling, rankby="p.value", idrcutoff=IDR_CUTOFF, force=force)
    )

    # peakranger
    # There is no p.value in the bcp algorithm, q.value only
    coro.append(
        _run_idr([
            _run_peakranger(experiments[0], pcalling, SOFT_FDR_CUTOFF, maxthreads, force),
            _run_peakranger(experiments[1], pcalling, SOFT_FDR_CUTOFF, maxthreads, force)
        ], pcalling, rankby="q.value" if is_broad(experiments[0].target) else "p.value",
            idrcutoff=IDR_CUTOFF, force=force)
    )

    results = await asyncio.gather(*coro)
    # also subtract blacklisted regions and presort
    blacklisted = BedTool(BLACKLISTED_REGIONS).sort()
    bed = []
    for r in results:
        b = BedTool(r).sort().subtract(blacklisted).sort()
        if len(b) == 0:
            import warnings
            warnings.warn(f"ERROR: {r} has no regions after idr")
            continue
        bed.append(b)

    if not os.path.exists(consensus) or force:
        # create consensus peak set
        threshold = len(results) // 2

        def merge(parents, consensus: Interval):
            score = f"{sum(float(i.fields[4]) for i in parents) / len(parents):.2f}"
            return Interval(consensus.chrom, consensus.start, consensus.end, score=score)

        bed = compute_conservative_regions(bed, threshold, merge)
        bed.sort().saveas(consensus)
    return consensus


async def _run_idr(files: [Awaitable[str]], saveto: str, rankby: str, idrcutoff: float, force: bool):
    files = await asyncio.gather(*files)
    fnames = [os.path.split(f)[-1] for f in files]
    assert fnames[0].split("-")[0] == fnames[1].split("-")[0] and fnames[0].split(".")[-1] == fnames[1].split(".")[-1]
    prefix, format = fnames[0].split("-")[0], fnames[0].split(".")[-1]
    saveto = os.path.join(saveto, f"idr+{prefix}.{format}")
    if not os.path.isfile(saveto) or force:
        await idr.idr(files, rankby, idrcutoff, format=format, saveto=saveto, plotto=saveto+".png")
    return saveto


async def _run_macs2(experiment: ExperimentMeta, saveto: str, fdrcutoff: float, force: bool):
    # mode = 'duplicated'     # macs2 automatically handles duplicates in a smart way
    mode = 'unique'
    isbroad = is_broad(experiment.target)
    postfix = "broadPeak" if isbroad else "narrowPeak"
    file = os.path.join(saveto, f"macs2-{experiment.name}.{postfix}")
    if not os.path.exists(file) or force:
        await macs2.callpeak(
            experiment.alltreatment(mode), experiment.allcontrol(mode), fdrcutoff=fdrcutoff,
            isbroad=isbroad, saveto=file, format='BAMPE' if experiment.paired_data else 'AUTO'
        )
    assert os.path.isfile(file)
    return file


async def _run_epic2(experiment: ExperimentMeta, saveto: str, fdrcutoff: float, threads: int, force: bool):
    mode = 'unique'
    file = os.path.join(saveto, f"epic2-{experiment.name}.broadPeak")
    if not os.path.isfile(file) or force:
        treatment, control = experiment.alltreatment(mode), experiment.allcontrol(mode)
        await epic2.epic2(treatment, control, experiment.paired_data, threads, saveto=file, fdrcutoff=fdrcutoff)
    assert os.path.isfile(file)
    return file


async def _run_peakranger(experiment: ExperimentMeta, saveto: str, fdrcutoff: float, threads: int, force: bool):
    mode = 'unique'
    treatment, control = experiment.alltreatment(mode), experiment.allcontrol(mode)

    if is_broad(experiment.target):
        file = os.path.join(saveto, f"bcp-{experiment.name}.broadPeak")
        if not os.path.exists(file) or force:
            await PeakRanger.bcp(treatment, control, saveto=file, maxthreads=threads, fdrcutoff=fdrcutoff)
    else:
        file = os.path.join(saveto, f"ranger-{experiment.name}.narrowPeak")
        if not os.path.exists(file) or force:
            await PeakRanger.ranger(treatment, control, saveto=file, maxthreads=threads, fdrcutoff=fdrcutoff)
    assert os.path.isfile(file)
    return file


async def _run_pepr(experiments: [ExperimentMeta], saveto: str, fdrcutoff: float, threads: int, force: bool):
    assert len(experiments) >= 2, "PePr works only with replicated ChIP-seq experiments"
    assert all(experiments[0].paired_data == e.paired_data and experiments[0].target == e.target for e in experiments)
    is_paired = experiments[0].paired_data
    threads = cpu_count() if threads < 0 else max(1, threads)
    peaktype = "broad" if is_broad(experiments[0].target) else "sharp"
    file = os.path.join(saveto, f"pepr-{'+'.join(e.name for e in experiments)}.broadPeak")
    if not os.path.exists(file) or force:
        mode = 'unique'
        treatment, control = sum([e.alltreatment(mode) for e in experiments], []), \
                             sum([e.allcontrol(mode) for e in experiments], [])
        treatment, control = list(set(treatment)), list(set(control))

        format = "bam" if not is_paired else "bampe"
        await PePr.PePr(treatment, control, peaktype, saveto=file, maxthreads=threads, format=format, fdrcutoff=fdrcutoff)
    assert os.path.exists(file)
    return file


__all__ = [call_peaks]
