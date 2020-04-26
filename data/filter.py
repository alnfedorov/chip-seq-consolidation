import numpy as np
from pybedtools import BedTool, Interval
from .misc import IntervalMeta
from .dataset import ChIPseq


def by_ambiguity(regions: BedTool, ambiguous: BedTool, ambiguity_thr: float) -> BedTool:
    # filter ambiguous regions
    # windows = bed.window_maker(b=bed, w=site_size)
    regions = regions.intersect(ambiguous, wao=True)  # type: BedTool
    filtered = []
    for w in regions:
        # assert len(w.fields) == 7, f"{w.fields}"
        if int(w.fields[-1]) / w.length < ambiguity_thr:
            filtered.append(Interval(w.chrom, w.start, w.end))
    result = BedTool(filtered).sort()

    # check for safety. Note that composition of good windows not necessarily imply good regions:
    # all inter_i / len_i < threshold != sum(inter) / sum(len) < threshold
    intersection = result.intersect(ambiguous).sort()
    ambiguity = intersection.total_coverage() / result.total_coverage()
    assert ambiguity <= ambiguity_thr + 0.05, \
        f"Filtering ambiguous regions failed: {ambiguity}, threshold {ambiguity_thr}"
    return result


def by_meta(dataset: ChIPseq, annotation: {str: BedTool}, ) -> [IntervalMeta]:
    """Drop regions that doesn't pass QC in-place"""
    metas = []
    todrop = []
    for ind, meta in enumerate(dataset.itermeta(annotation)):
        metas.append(meta)
        # We are not interested in regions where ratio std is equal to 0, there are either constant or no reads at all
        # More checks might be added later
        if meta.ratio_std < 1e-6:
            todrop.append(ind)
    dataset.remove(todrop)
    metas = np.delete(np.asarray(metas, dtype=np.object), todrop)
    return metas.tolist()
