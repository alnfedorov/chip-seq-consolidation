import os
from inspect import isawaitable
import asyncio
from wrappers import chips, sambamba
from copy import deepcopy
from typing import Awaitable
from .fastq import align
from .utils import BamMeta, ExperimentMeta, make_filename, batched_gather
from .config import CHIPS_MODELS_DIR, READS_TO_SIMULATE, UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, BWA_GENOME_INDEX

# Question, should I use unique reads or something else?
# There is no reason to use unfiltered,
# because I would end up with the same reads after processing but with less number of them
async def subsample_simulations(subsfolder: str, experiments: [ExperimentMeta],
                                maxthreads: int = 1, force: bool = False) -> [ExperimentMeta]:
    """
    Simulate ChIP-seq experiment from the existing reads by naive subsampling of the unique reads.
    """
    assert len(experiments) == 2 and all(exp.target == experiments[0].target for exp in experiments)
    hmod = experiments[0].target
    for folder, reads in READS_TO_SIMULATE[hmod].items():
        folder = os.path.join(subsfolder, folder)
        subsampled = deepcopy(experiments)
        allbam = set(sum([exp.treatment + exp.control for exp in subsampled], []))
        for meta in allbam:
            meta.reads = reads
            meta.unfiltered = None    # keep dummy pointer to avoid reprocessing
            # meta.unique = None
            meta.simulated = True
            duplicated = make_filename(folder, DUPLICATED_READS_DIR, mode="duplicated",
                                       name=meta.name, accession=meta.accession, reads=meta.reads)
            if not os.path.isfile(duplicated) or force:
                await sambamba.subsample(
                    meta.duplicated, reads=reads, saveto=duplicated, threads=maxthreads
                )
                await sambamba.sort(duplicated, inplace=True, threads=maxthreads)
            meta.duplicated = duplicated
        yield folder, subsampled


# async def chips_simulations(chipsfolder: str, genome: str, consensus: str,
#                             experiments: [ExperimentMeta], maxthreads: int = 1, force: bool = False) -> [ExperimentMeta]:
#     assert len(experiments) == 2 and all(exp.target == experiments[0].target for exp in experiments)
#     # 1. learn chips model for treatment files.
#     score_column = 5  # default for the narrowPeak and boradPeak bed files
#     models_dir = os.path.join(chipsfolder, CHIPS_MODELS_DIR)
#     treatment = set(sum([exp.treatment for exp in experiments], []))
#     models = {}
#     for meta in treatment:
#         saveto = os.path.join(models_dir, f"{meta.accession}.json")
#         if os.path.exists(saveto) and not force:
#             models[meta.accession] = saveto
#         else:
#             models[meta.accession] = asyncio.create_task(chips.learn(
#                 meta.duplicated, meta.paired, consensus, score_column, saveto
#             ))
#     models = {key: await model if isawaitable(model) else model for key, model in models.items()}
#
#     # 2. Run simulations for all treatment and control files, align results.
#     #    asyncio is not used, because BWA and chips perform well with multiple threads
#     hmod = experiments[0].target
#     alltreatment = set(sum([exp.treatment for exp in experiments], []))
#     allcontrol = set(sum([exp.control for exp in experiments], []))
#     allbam = alltreatment.union(allcontrol)
#     bam_meta: {str, BamMeta} = {meta.accession: meta for meta in allbam}
#
#     for folder, reads in READS_TO_SIMULATE[hmod].items():
#         expfolder = os.path.join(chipsfolder, folder)
#         fastqc = {}  # accession -> raw fastqc reads
#         bam = {   # accession -> unfiltered bam file
#             meta.accession: make_filename(expfolder, UNFILTERED_READS_DIR, mode='unfiltered', name=meta.name,
#                                           accession=meta.accession, reads=reads) for meta in allbam
#         }
#         # simulate reads if bam files doesn't exist
#         adjusted_reads = int(1.2 * reads)   # adjust number of reads to address possible duplicates and unmapped
#         for meta in alltreatment:
#             if not os.path.exists(bam[meta.accession]) or force:
#                 # Don't use BAM scoring, it sucks
#                 fastqc[meta.accession] = chips.simreads(
#                     genome, models[meta.accession], peaks=consensus, numreads=adjusted_reads, readlen=meta.readlen,
#                     paired=meta.paired, threads=maxthreads, score_column=score_column, scale_outliers=True
#                 )
#         for meta in allcontrol:
#             if not os.path.exists(bam[meta.accession]) or force:
#                 fastqc[meta.accession] = chips.wce(
#                     genome, numreads=adjusted_reads, readlen=meta.readlen, paired=meta.paired,
#                     threads=maxthreads, bam=meta.duplicated
#                 )
#
#         # align results if needed
#         async def job(meta: BamMeta, fastqc: Awaitable[str]):
#             fastqc = await fastqc
#             if meta.paired:
#                 await align.pend()
#             else:
#                 await align.send(BWA_GENOME_INDEX, fastqc, meta.readlen, maxthreads, saveto=bam[meta.accession])
#             os.remove(fastqc)
#
#         await batched_gather([job(bam_meta[accession], file) for accession, file in fastqc.items()], batch_size=2)
#
#         # construct experiments meta
#         simulated_experiments = deepcopy(experiments)
#         simulated_bam = set(sum([exp.treatment + exp.control for exp in simulated_experiments], []))
#         for meta in simulated_bam:
#             meta.unfiltered = bam[meta.accession]
#             meta.duplicated = None
#             # meta.unique = None
#             meta.reads = reads
#             meta.simulated = True
#
#         yield expfolder, simulated_experiments
