import numpy as np
import utils
from pybedtools import BedTool
from data.pod import IntervalMeta
from data.dataset import ChIPseq


def by_ambiguity(regions: BedTool, ambiguous: BedTool, ambiguity_thr: float) -> BedTool:
    # filter ambiguous regions
    # windows = bed.window_maker(b=bed, w=site_size)
    ambiguity = utils.bed.coverage(regions, ambiguous, presorted=True)
    assert len(ambiguity) == len(regions)
    result = BedTool([reg for ambg, reg in zip(ambiguity, regions) if ambg <= ambiguity_thr]).sort()

    # check for sanity. Note that composition of good windows not necessarily imply good regions:
    # all inter_i / len_i < threshold != sum(inter) / sum(len) < threshold
    intersection = result.intersect(ambiguous).sort()
    ambiguity = intersection.total_coverage() / result.total_coverage()
    assert ambiguity <= ambiguity_thr + 0.05, \
        f"Filtering ambiguous regions failed: {ambiguity}, threshold {ambiguity_thr}"
    return result


def by_meta(dataset: ChIPseq, annotation: {str: BedTool}, verbose: bool = True) -> [IntervalMeta]:
    """Drop regions that doesn't pass QC in-place"""
    metas = list(dataset.itermeta(annotation))

    # drop first 25 percent with the lowest std
    std = np.asarray([m.ratio_std for m in metas])
    threshold = np.quantile(std, 0.25)
    if verbose:
        print(f"Thresholding std by 0.25 quantile: {threshold}")

    todrop = [ind for ind, m in enumerate(metas) if m.ratio_std <= threshold]
    dataset.remove(todrop)
    metas = np.delete(np.asarray(metas, dtype=np.object), todrop)
    return metas.tolist()
