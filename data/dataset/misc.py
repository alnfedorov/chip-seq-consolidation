from collections import defaultdict
from typing import Tuple, Dict, Optional

import numpy as np
import pyBigWig
from pybedtools import BedTool, Interval

from data.pod import BigWigMeta, ChIPseqReplicaMeta


class PeaksIndex:
    """Helper class to make fast in-memory intersections between peaks and arbitrary intervals"""
    def __init__(self, bed: BedTool, presorted_and_merged=True):
        if not presorted_and_merged:
            bed = bed.sort().merge().sort()
        index = defaultdict(list)
        for inter in bed:
            index[inter.chrom].append(inter)
        self.index = {
            chrom: (np.asarray([x.start for x in intervals]), np.asarray([x.end for x in intervals]))
            for chrom, intervals in index.items()
        }
        self.bedfn = bed.saveas()

    def intersect(self, interval: Interval):
        intersected = []
        starts, ends = self.index[interval.chrom]
        ind = np.searchsorted(starts, interval.start)
        if ind > 0:
            # check to avoid missing hits with the end of the previous interval
            interstrt, interend = max(starts[ind - 1], interval.start), min(ends[ind - 1], interval.end)
            if interstrt < interend:
                intersected.append(Interval(interval.chrom, interstrt, interend))
        while ind < len(starts):
            interstrt, interend = max(starts[ind], interval.start), min(ends[ind], interval.end)
            if interstrt >= interend:
                break
            intersected.append(Interval(interval.chrom, interstrt, interend))
            ind += 1
        return intersected


# @dataclass(frozen=True)
# class ItemMeta:
#     target: str
#     sampling: str
#     peaks: PeaksIndex
#     treatment: Tuple[BigWigMeta, ...]
#     control: Tuple[BigWigMeta, ...]
#     reftreatment: Tuple[BigWigMeta, ...]
#     refcontrol: Tuple[BigWigMeta, ...]
#
#     @staticmethod
#     def make(sampling, replica: ChIPseqReplicaMeta):
#         bymode = lambda mode, files: tuple(f[mode] for f in files)
#         peaks = BedTool(replica.peaks).sort().merge().sort()
#         kwargs = {
#             "target": replica.target, "sampling": sampling, "peaks": PeaksIndex(peaks),
#             "reftreatment": bymode("original", replica.treatment),
#             "refcontrol": bymode("original", replica.control),
#             "treatment": bymode(sampling, replica.treatment),
#             "control": bymode(sampling, replica.control),
#         }
#         return ItemMeta(**kwargs)


# dummy enrichment computation. It is NOT read length / fragment size aware
def enrichment(treatment: Tuple[BigWigMeta, ...], control: Tuple[BigWigMeta, ...],
               reads: Dict[BigWigMeta, np.ndarray]) -> np.ndarray:
    # cat and normalize reads
    catnormalize = lambda bwmeta: np.sum([reads[bw] for bw in bwmeta], axis=0).astype(np.float64) / \
                                  sum(bw.readlen * bw.numreads for bw in bwmeta)
    eps = np.finfo(np.float64).eps
    treatment, control = catnormalize(treatment), catnormalize(control)
    enrichment = treatment / (treatment + control + eps)
    return enrichment.astype(np.float32)


class ChIPseqReplica:
    def __init__(self, enrichment: BigWigMeta, peaks: Optional[str] = None, uid: str = ""):
        if peaks is not None:
            self.peaks_index: Optional[PeaksIndex] = PeaksIndex(
                BedTool(peaks), presorted_and_merged=False
            )
        else:
            self.peaks_index: Optional[PeaksIndex] = None
        self.uid = uid
        self._enrichmentfn = enrichment.path
        self._enrichment = None
        self._reread = False
        self.reread_bigwig_on_query()

    @staticmethod
    def frommeta(meta: ChIPseqReplicaMeta, sampling: str) -> 'ChIPseqReplica':
        assert sampling in meta.enrichment
        peaks = meta.peaks if sampling == "original" else None
        uid = sorted(meta.treatment)
        uid = uid[0] if len(uid) == 1 else ','.join(uid)
        return ChIPseqReplica(meta.enrichment[sampling], peaks, uid)

    def _reread_bigwig(self):
        self._enrichment = pyBigWig.open(self._enrichmentfn)

    def reread_bigwig_on_query(self):
        self._enrichment = None
        self._reread = True

    def enrichment(self, interval: Interval) -> np.ndarray:
        if self._reread:
            self._reread_bigwig()
            self._reread = False

        chrom, start, end = interval.chrom, interval.start, interval.end
        return self._enrichment.values(chrom, start, end, numpy=True)

    def peaks(self, interval: Interval) -> np.ndarray:
        if self.peaks_index is None:
            raise ValueError("Intersection with peaks is requested, but peaks file wasn't provided in the constructor")
        intersection = self.peaks_index.intersect(interval)
        peaks = np.zeros(interval.length, dtype=np.float32)
        if len(intersection) == 0:
            return peaks

        start = interval.start
        for inter in intersection:
            b, e = inter.start - start, inter.end - start
            peaks[b:e] = 1

        return peaks

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_enrichment"] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.reread_bigwig_on_query()
