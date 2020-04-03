import os
import asyncio
from .mktree import mktree
from . import encode, bam, simulations
from .peak_calling import call_peaks
from .utils import ExperimentMeta, BamMeta, make_filename, config_logging
from .config import ORIGINAL_DIR, UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR, ASSEMBLY, \
                    PEAK_CALLING_DIR, CHIPS_DIR, BAM_QC_DIR, GENOME_FASTA, LOGS_DIR


async def _process_bam(expfolder: str, experiment: str, meta: BamMeta, maxthreads: int, force: bool):
    unfiltered, duplicated, unique, qc = (
        os.path.join(expfolder, f) for f in [UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR, BAM_QC_DIR]
    )
    coro = []

    # Fetch aligned but not filtered reads, bam format
    meta.unfiltered = make_filename(unfiltered, mode="unfiltered", name=meta.name, accession=meta.accession)
    if not os.path.exists(meta.unfiltered) or force:
        await encode.fetch(meta.accession, experiment, meta.target == "control", ASSEMBLY, saveto=meta.unfiltered)
        coro.append(asyncio.create_task(
            bam.qc(meta.unfiltered, saveto=qc, doflagstats=True, dostats=True, dofastqc=True)
        ))

    # Mark duplicates, sort and filter. (Duplicates are kept)
    meta.duplicated = make_filename(duplicated, mode="duplicated", name=meta.name, accession=meta.accession)
    if not os.path.exists(meta.duplicated) or force:
        await bam.filter.toduplicates(meta.unfiltered, saveto=meta.duplicated, maxthreads=maxthreads, postsort=True)
        coro.append(asyncio.create_task(
            bam.qc(meta.duplicated, saveto=qc, doflagstats=True, dostats=True, dofastqc=False)
        ))

    # Remove duplicates and keep only unique reads
    meta.unique = make_filename(unique, mode="unique", name=meta.name, accession=meta.accession)
    if not os.path.exists(meta.unique) or force:
        await bam.filter.tounique(meta.duplicated, saveto=meta.unique, maxthreads=maxthreads, postsort=False)
        coro.append(asyncio.create_task(
            bam.qc(meta.unique, saveto=qc, doflagstats=True, dostats=True, dofastqc=False)
        ))

    await asyncio.gather(*coro)


async def _process_experiments(expfolder: str, experiments: [ExperimentMeta], maxthreads: int, force: bool):
    # Filter and markdup files
    coro = [(exp.accession, meta) for exp in experiments for meta in exp.treatment + exp.control]
    coro = set(coro)    # filter possible duplicates in BamMeta objects
    threads = max(1, int(round(len(coro) / maxthreads)))
    coro = [_process_bam(expfolder, accession, meta, maxthreads=threads, force=force) for accession, meta in coro]
    await asyncio.gather(*coro)

    # Call peaks
    pcalling = os.path.join(expfolder, PEAK_CALLING_DIR)
    consensus = await call_peaks(pcalling, experiments, maxthreads, force)
    return consensus


async def run(root: str, experiments: [ExperimentMeta],
              genome: str = GENOME_FASTA, maxthreads: int = -1, force: bool = False):
    if maxthreads < 0:
        maxthreads = os.cpu_count()
    else:
        maxthreads = max(1, maxthreads)

    # 1. create experiment folders tree
    mktree(root)
    logs_dir = os.path.join(root, ORIGINAL_DIR, LOGS_DIR)
    config_logging(logs_dir)

    # 2. Generate ground-true consensus labels
    expfolder = os.path.join(root, ORIGINAL_DIR)
    consensus = await _process_experiments(expfolder, experiments, maxthreads, force)

    # 3. Simulate chip-seq pipeline and generate a consensus for synthetic data.
    # chipsfolder = os.path.join(root, CHIPS_DIR)
    # async for (expfolder, experiments) in simulations.simulate(chipsfolder, genome,
    #                                                            consensus, experiments, maxthreads=maxthreads):
    # logs_dir = os.path.join(expfolder, LOGS_DIR)
    # config_logging(LOGS_DIR)
    #     await _process_experiments(expfolder, experiments, maxthreads)
