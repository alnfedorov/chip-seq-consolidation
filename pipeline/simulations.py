import os
import asyncio
from wrappers import chips
from copy import deepcopy
from .utils import BamMeta, ExperimentMeta, make_filename
from .config import CHIPS_MODELS_DIR, READS_TO_SIMULATE, UNFILTERED_READS_DIR, UNFILTERED_READS_PREFIX

# How to properly separate the workload?
# The best option is to learn/and simulate in the same time for all experiments.
# and return new ExperimentMeta's tuple


async def simulate(chipsfolder: str, genome: str, consensus: str,
                   experiments: [ExperimentMeta], maxthreads: int = 1) -> [ExperimentMeta]:
    assert len(experiments) == 2 and all(exp.hmod == experiments[0].hmod for exp in experiments)
    # TODO: ensure duplicates are working fine
    # 1. learn chips model for treatment files
    score_column = 5
    treatment = sum([exp.treatment for exp in experiments], [])
    models_dir = os.path.join(chipsfolder, CHIPS_MODELS_DIR)
    models = [
        asyncio.create_task(chips.learn(meta.duplicated, meta.paired, consensus,
                                        score_column, os.path.join(models_dir, f"{meta.accession}.json")))
        for meta in treatment
    ]
    models = {meta.accession: await model for meta, model in zip(treatment, models)}

    # 2. Run simulations for all treatment and control files, align results.
    #    asyncio is not used, because BWA and chips perform well with multiple threads
    hmod = experiments[0].hmod
    files = {}  # accession -> (raw fastqc reads -> unfiltered bam)
    alltreatment = sum([exp.treatment for exp in experiments], [])
    allcontrol = sum([exp.control for exp in experiments], [])
    for folder, reads in READS_TO_SIMULATE[hmod].items():
        expfolder = os.path.join(chipsfolder, folder)
        # simulate reads
        for meta in alltreatment:
            files[meta.accession] = await chips.simreads(
                consensus, genome, models[meta.accession], numreads=reads, readlen=meta.readlen,
                paired=meta.paired, threads=maxthreads, bam=meta.duplicated, score_column=score_column
            )
        for meta in allcontrol:
            files[meta.accession] = await chips.wce(
                genome, numreads=reads, readlen=meta.readlen, paired=meta.paired,
                threads=maxthreads, bam=meta.duplicated
            )

        # align results
        unfiltered_dir = os.path.join(expfolder, UNFILTERED_READS_DIR)
        pass

        # construct experiments meta
        simulated_experiments = deepcopy(experiments)
        bam_meta = sum([exp.treatment + exp.control for exp in simulated_experiments], [])
        for meta in bam_meta:
            meta.unfiltered = files[meta.accession]
            meta.duplicated = None
            meta.unique = None

        yield expfolder, simulated_experiments
