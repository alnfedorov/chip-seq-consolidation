import os
from wrappers import sambamba
from copy import deepcopy
from .meta import ExperimentMeta
from .path import make_filename
from .config import READS_TO_SIMULATE, UNIQUE_READS_DIR


# There is no reason to use unfiltered reads during subsampling.
# However, sampling must be done based on the file with duplicates because it improves the possibility to include
# duplicated reads into the resulted file. How to determine the exact number of reads after that? Great question. I don't know.
# Good point, however, we never can be sure is it a PCR duplicate or library problem. Hence, for the lack of the simplicity we shall use unique reads.
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
            meta.simulated = True
            unique = make_filename(folder, UNIQUE_READS_DIR, mode="unique",
                                   name=meta.name, accession=meta.accession, reads=meta.reads)
            if not os.path.isfile(unique) or force:
                file = sambamba.subsample(
                    meta.unique, reads=reads, threads=maxthreads
                )
                await sambamba.sort(file, saveto=unique, threads=maxthreads)
            meta.unfiltered = unique    # keep dummy pointer to avoid reprocessing
            meta.unique = unique
            meta.duplicated = unique
        yield folder, subsampled
