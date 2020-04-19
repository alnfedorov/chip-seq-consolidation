import os
from wrappers import sambamba
from copy import deepcopy
from .utils import ExperimentMeta, make_filename
from .config import READS_TO_SIMULATE, DUPLICATED_READS_DIR


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
