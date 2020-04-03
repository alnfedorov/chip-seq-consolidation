import os
import asyncio
from typing import Awaitable
from multiprocessing import cpu_count
from wrappers.peak_calling import macs2, epic2, PePr, PeakRanger, idr
from pybedtools import BedTool, Interval
from utils.bed import compute_conservative_regions
from .utils import ExperimentMeta, BamMeta, is_broad
from .config import SOFT_FDR_CUTOFF, IDR_CUTOFF, CONSENSUS_FILENAME, BLACKLISTED_REGIONS


async def call_peaks(pcalling: str, experiments: [ExperimentMeta], maxthreads: int, force: bool):
    assert len(experiments) == 2
    consensus = os.path.join(pcalling, CONSENSUS_FILENAME)

    if os.path.exists(consensus) and not force:
        return consensus

    coro = [asyncio.create_task(_run_pepr(experiments, pcalling, IDR_CUTOFF, maxthreads))]

    # macs2
    coro.append(
        _run_idr([
            _run_macs2(experiments[0], pcalling, SOFT_FDR_CUTOFF),
            _run_macs2(experiments[1], pcalling, SOFT_FDR_CUTOFF)
        ], pcalling, rankby="p.value", idrcutoff=IDR_CUTOFF)
    )

    # epic2
    coro.append(
        _run_idr([
            _run_epic2(experiments[0], pcalling, SOFT_FDR_CUTOFF),
            _run_epic2(experiments[1], pcalling, SOFT_FDR_CUTOFF)
        ], pcalling, rankby="p.value", idrcutoff=IDR_CUTOFF)
    )

    # peakranger
    coro.append(
        _run_idr([
            _run_peakranger(experiments[0], pcalling, SOFT_FDR_CUTOFF, maxthreads),
            _run_peakranger(experiments[1], pcalling, SOFT_FDR_CUTOFF, maxthreads)
        ], pcalling, rankby="p.value", idrcutoff=IDR_CUTOFF)
    )

    await asyncio.gather(*coro)
    results = [BedTool(os.path.join(pcalling, file)) for file in os.listdir(pcalling) if file.endswith("Peak")]

    # subtract blacklisted regions and presort
    blacklisted = BedTool(BLACKLISTED_REGIONS).sort()
    results = [bed.sort().subtract(blacklisted).sort() for bed in results]

    # create consensus peak set
    threshold = len(results) // 2

    def merge(parents, consensus: Interval):
        score = f"{sum(float(i.fields[4]) for i in parents) / len(parents):.2f}"
        return Interval(consensus.chrom, consensus.start, consensus.end, score=score)

    bed = compute_conservative_regions(results, threshold, merge)
    bed.sort().saveas(consensus)
    return consensus


async def _run_idr(files: [Awaitable[str]], saveto: str, rankby: str, idrcutoff: float):
    files = await asyncio.gather(*files)
    fnames = [os.path.split(f)[-1] for f in files]
    assert fnames[0].split("-")[0] == fnames[1].split("-")[0] and fnames[0].split(".")[-1] == fnames[1].split(".")[-1]
    prefix, format = fnames[0].split("-")[0], fnames[0].split(".")[-1]
    saveto = os.path.join(saveto, f"idr+{prefix}.{format}")
    await idr.idr(files, rankby, idrcutoff, format=format, saveto=saveto, plotto=saveto+".png")
    return saveto


async def _run_macs2(experiment: ExperimentMeta, saveto: str, fdrcutoff: float):
    mode = 'duplicated'
    isbroad = is_broad(experiment.target)
    postfix = "broadPeak" if isbroad else "narrowPeak"
    file = os.path.join(saveto, f"macs2-{experiment.name}.{postfix}")
    await macs2.callpeak(
        experiment.alltreatment(mode), experiment.allcontrol(mode), fdrcutoff=fdrcutoff,
        isbroad=isbroad, saveto=file, format='AUTO' if not experiment.paired_data else 'BAMPE'
    )
    assert os.path.isfile(file)
    return file


async def _run_epic2(experiment: ExperimentMeta, saveto: str, fdrcutoff: float):
    assert not experiment.paired_data, "run_epic2 is not implemented to BAM->BEDPE on-the-fly conversion"
    mode = 'unique'
    file = os.path.join(saveto, f"epic2-{experiment.name}.broadPeak")
    await epic2.epic2(experiment.alltreatment(mode), experiment.allcontrol(mode), saveto=file, fdrcutoff=fdrcutoff)
    assert os.path.isfile(file)
    return file


async def _run_peakranger(experiment: ExperimentMeta, saveto: str, fdrcutoff: float, threads: int = -1):
    mode = 'unique'
    treatment, control = experiment.alltreatment(mode), experiment.allcontrol(mode)

    if is_broad(experiment.target):
        file = os.path.join(saveto, f"bcp-{experiment.name}.broadPeak")
        await PeakRanger.bcp(treatment, control, saveto=saveto, fdrcutoff=fdrcutoff)
    else:
        file = os.path.join(saveto, f"ranger-{experiment.name}.narrowPeak")
        threads = cpu_count() if threads < 0 else max(1, threads)
        await PeakRanger.ranger(
            treatment, control, saveto=file, threads=threads, fdrcutoff=fdrcutoff
        )

    assert os.path.isfile(file)
    return file


async def _run_pepr(experiments: [ExperimentMeta], saveto: str, fdrcutoff: float, threads: int = -1):
    assert len(experiments) >= 2, "PePr works only with replicated ChIP-seq experiments"
    assert all(experiments[0].paired_data == e.paired_data and experiments[0].hmod == e.hmod for e in experiments)
    threads = cpu_count() if threads < 0 else max(1, threads)
    peaktype = "broad" if is_broad(experiments[0].hmod) else "sharp"
    file = os.path.join(saveto, f"pepr-{'+'.join(e.name for e in experiments)}.broadPeak")

    mode = 'unique'
    treatment, control = sum([e.alltreatment(mode) for e in experiments], []), \
                         sum([e.allcontrol(mode) for e in experiments], [])
    treatment, control = list(set(treatment)), list(set(control))

    is_paired = experiments[0].paired_data
    format = "bam" if not is_paired else "bampe"
    await PePr.PePr(treatment, control, peaktype, saveto=file, threads=threads, format=format, fdrcutoff=fdrcutoff)
    assert os.path.exists(file)
    return file


__all__ = [call_peaks]
