import os
import asyncio
from .path import mktree
from . import encode, bam, simulations
from wrappers import cookbook
from .peak_calling import call_peaks
from .utils import ExperimentMeta, BamMeta, make_filename, config_logging
from .config import ORIGINAL_DIR, UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR, ASSEMBLY, CHROMINFO, \
                    PEAK_CALLING_DIR, SIMULATED_DIR, BIGWIG_DIR, BAM_QC_DIR, GENOME_FASTA, LOGS_DIR, SUBSAMPLE_DIR


async def _process_bam(expfolder: str, experiment: str, meta: BamMeta, maxthreads: int, force: bool):
    unfiltered, duplicated, unique, qc = (
        os.path.join(expfolder, f) for f in [UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR, BAM_QC_DIR]
    )
    coro = []

    # if needed, fetch aligned but not filtered reads, bam format
    if meta.unfiltered is None:
        meta.unfiltered = make_filename(unfiltered, mode="unfiltered", name=meta.name, accession=meta.accession,
                                        reads=meta.reads)
    if not meta.simulated:
        existed = os.path.exists(meta.unfiltered)
        await encode.fetch(meta.accession, experiment, meta.target == "control", ASSEMBLY, saveto=meta.unfiltered)
        if not existed or force:
            coro.append(asyncio.create_task(
                bam.qc(meta.unfiltered, saveto=qc, doflagstats=True, dostats=True, dofastqc=True)
            ))

    # Mark duplicates, sort and filter.
    if meta.duplicated is None:
        meta.duplicated = make_filename(duplicated, mode="duplicated", name=meta.name, accession=meta.accession,
                                        reads=meta.reads)
    if not os.path.exists(meta.duplicated) or force:
        await bam.filter.toduplicates(meta.unfiltered, meta.paired, saveto=meta.duplicated, maxthreads=maxthreads)
        coro.append(asyncio.create_task(
            bam.qc(meta.duplicated, saveto=qc, doflagstats=True, dostats=True, dofastqc=False)
        ))

    # Remove duplicates and keep only unique reads
    if meta.unique is None:
        meta.unique = make_filename(unique, mode="unique", name=meta.name, accession=meta.accession, reads=meta.reads)
    if not os.path.exists(meta.unique) or force:
        await bam.filter.tounique(meta.duplicated, meta.paired, saveto=meta.unique, maxthreads=maxthreads)
        coro.append(asyncio.create_task(
            bam.qc(meta.unique, saveto=qc, doflagstats=True, dostats=True, dofastqc=False)
        ))

    # Make bigwig file from the unique bam
    bigwig = make_filename(expfolder, BIGWIG_DIR, mode="unique", name=meta.name, accession=meta.accession, format=".bw")
    if not os.path.exists(bigwig) or force:
        await cookbook.bam_to_bigwig(meta.unique, CHROMINFO, bigwig)

    await asyncio.gather(*coro)


async def _process_experiments(expfolder: str, experiments: [ExperimentMeta], maxthreads: int, force: bool):
    # Filter and markdup files
    coro = [(exp.accession, meta) for exp in experiments for meta in exp.treatment + exp.control]
    coro = set(coro)    # filter possible duplicates in BamMeta objects
    threads = int(round(maxthreads / len(coro))) + 2
    # threads = maxthreads
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

    subsfolder = os.path.join(root, SIMULATED_DIR, SUBSAMPLE_DIR)
    async for (expfolder, experiments) in simulations.subsample_simulations(subsfolder, experiments,
                                                                            maxthreads=maxthreads, force=force):
        logs_dir = os.path.join(expfolder, LOGS_DIR)
        config_logging(logs_dir)
        await _process_experiments(expfolder, experiments, maxthreads, force)

    # 3. Simulate chip-seq pipeline and generate a consensus for synthetic data.
    # chipsfolder = os.path.join(root, SIMULATED_DIR, CHIPS_DIR)
    # async for (expfolder, experiments) in simulations.chips_simulations(chipsfolder, genome, consensus, experiments,
    #                                                                     maxthreads=maxthreads, force=force):
    #     logs_dir = os.path.join(expfolder, LOGS_DIR)
    #     config_logging(logs_dir)
    #     await _process_experiments(expfolder, experiments, maxthreads, force)
