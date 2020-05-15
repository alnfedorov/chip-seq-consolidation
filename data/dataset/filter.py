from typing import List

import numpy as np
import utils
from pybedtools import BedTool
from data.dataset import ChIPseqReplicaTrainDataset
from data.pod import IntervalMeta
from tqdm import tqdm


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


def by_meta(dataset: ChIPseqReplicaTrainDataset, annotation: {str: str}, verbose: bool = True) -> [IntervalMeta]:
    """Drop regions that are not 'difficult' in-place"""
    # 1. compute per region intersection with annotations
    intervals = BedTool(dataset.intervals).sort()
    annocov = {key: utils.bed.coverage(intervals, BedTool(fn), presorted=True) for key, fn in annotation.items()}
    # 2. compute intersection between all peaks and intervals
    peakscov = utils.bed.coverage(intervals, BedTool(dataset.reference.peaks_index.bedfn), presorted=True)
    # 3. iterate and build interval meta
    metas: List[IntervalMeta, ...] = []
    assert len(intervals) == len(dataset) == len(peakscov) and all(len(x) == len(intervals) for x in annocov.values())
    for idx, interval in tqdm(enumerate(intervals), total=len(intervals)):
        enrichment = dataset.subsampled.enrichment(interval=interval)
        assert enrichment.ndim == 1
        anno = {key: value[idx] for key, value in annocov.items()}
        metas.append(
            IntervalMeta(peakscov[idx], np.nanstd(enrichment, ddof=1), dataset.target, dataset.sampling, anno)
        )
    # drop first 25 percent with the lowest std
    std = np.asarray([m.ratio_std for m in metas if not np.isnan(m.ratio_std)])
    threshold = np.quantile(std, 0.25)
    if verbose:
        print(f"Thresholding std by 0.25 quantile: {threshold}")
    todrop = [ind for ind, m in enumerate(metas) if m.ratio_std <= threshold or np.isnan(m.ratio_std)]
    return metas, todrop
