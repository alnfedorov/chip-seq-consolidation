from collections import defaultdict
from itertools import chain
from typing import Dict, Optional

import numpy as np
from pybedtools import BedTool, Interval

from data.pod import SeparatedReadsMeta, IntervalReads
from utils import enrichment


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
        if interval.chrom not in self.index:
            return intersected
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


class ChIPseqReplica:
    def __init__(self, treatment: SeparatedReadsMeta, control: SeparatedReadsMeta,
                 chrominfo: Dict[str, int], effective_genome_size: int, peaks: Optional[str] = None):
        if peaks is not None:
            self.peaks_index: Optional[PeaksIndex] = PeaksIndex(
                BedTool(peaks), presorted_and_merged=False
            )
        else:
            self.peaks_index: Optional[PeaksIndex] = None
        self.treatment = treatment
        self.control = control
        self.chrominfo = chrominfo
        self.effective_genome_size = effective_genome_size
        self._ndarrays: Dict[str, np.ndarray] = {}

    def reread(self):
        for path in chain(self.treatment.forward.values(), self.treatment.reverse.values(),
                          self.control.forward.values(), self.control.reverse.values()):
            # self._ndarrays[path] = np.load(path)
            self._ndarrays[path] = np.memmap(path, dtype=np.int32, mode='r')

    def _slice(self, forward, reverse, interval: Interval, five_slope: int, three_slope: int):
        # TODO: coordinates correction for 1
        # That is rather strange, but pure python bisect is faster than np.searchsorted
        # Perhaps, searchsorted is heavy optimized for the batched search
        from bisect import bisect_left, bisect_right

        beginf = bisect_left(forward, interval.start - five_slope)
        endf = bisect_right(forward, interval.end + three_slope, lo=beginf)
        # forward = np.asarray(forward[beginf: endf])
        forward = forward[beginf: endf]

        beginr = bisect_left(reverse, interval.start - three_slope)
        endr = bisect_right(reverse, interval.end + five_slope, lo=beginr)
        reverse = reverse[beginr: endr]
        return np.asarray(forward), np.asarray(reverse)
        # return forward, reverse
        # beginf = np.searchsorted(forward, interval.start - five_slope, side="left")
        # endf = np.searchsorted(forward[beginf:], interval.end + three_slope, side="right") + beginf
        # # forward = np.asarray(forward[beginf: endf])
        # forward = forward[beginf: endf]
        #
        # beginr = np.searchsorted(reverse, interval.start - three_slope, side="left")
        # endr = np.searchsorted(reverse[beginr:], interval.end + five_slope, side="right") + beginr
        # reverse = reverse[beginr: endr]
        # beginf = np.searchsorted(forward, interval.start - five_slope, side="left")
        # endf = np.searchsorted(forward, interval.end + three_slope, side="right")
        # # forward = np.asarray(forward[beginf: endf])
        # forward = forward[beginf: endf]
        #
        # beginr = np.searchsorted(reverse, interval.start - three_slope, side="left")
        # endr = np.searchsorted(reverse, interval.end + five_slope, side="right")
        # reverse = reverse[beginr: endr]
        # reverse = np.asarray(reverse[beginr: endr])
        # return forward, reverse

    def reads(self, bamreads: SeparatedReadsMeta, interval: Interval, five_slope: int,
              three_slope: int) -> IntervalReads:
        chr = interval.chrom
        if len(self._ndarrays) == 0:
            self.reread()

        forward, reverse = self._ndarrays[bamreads.forward[chr]], self._ndarrays[bamreads.reverse[chr]]
        forward, reverse = self._slice(forward, reverse, interval, five_slope, three_slope)
        return IntervalReads(forward, reverse)

    def enrichment(self, interval: Interval, fragment_size: int):
        chr = interval.chrom
        chromlen = self.chrominfo[chr]

        treatment = self.reads(self.treatment, interval, fragment_size, 0)
        control = self.reads(self.control, interval, 5000, 5000)
        enrich_ends, enrich_values = enrichment.endtoend(chromlen, fragment_size, self.effective_genome_size,
                                                         self.treatment.numreads, self.control.numreads,
                                                         treatment.forward, treatment.reverse,
                                                         control.forward, control.reverse)
        result = enrichment.todense(enrich_ends, enrich_values, np.int32(interval.start), np.int32(interval.end))
        return treatment, control, result

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
        del state["_ndarrays"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._ndarrays = {}
