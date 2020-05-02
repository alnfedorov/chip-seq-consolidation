import numpy as np
import tempfile
from typing import Tuple
from pybedtools import BedTool, Interval
from data.pod import BigWigMeta, ReplicaMeta
from collections import defaultdict
from dataclasses import dataclass


class PeaksIndex:
    """Helper class to make fast in-memory intersections between peaks and arbitrary region"""
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
        self.bedfn = bed.saveas(tempfile.mkstemp()[-1]).fn

    def intersect(self, interval: Interval):
        intersected = []
        starts, ends = self.index[interval.chrom]
        ind = np.searchsorted(starts, interval.start)
        if ind > 0:
            # check to avoid missing hits with the end of the previous interval
            interstrt, interend = max(starts[ind-1], interval.start), min(ends[ind-1], interval.end)
            if interstrt < interend:
                intersected.append(Interval(interval.chrom, interstrt, interend))
        while ind < len(starts):
            interstrt, interend = max(starts[ind], interval.start), min(ends[ind], interval.end)
            if interstrt >= interend:
                break
            intersected.append(Interval(interval.chrom, interstrt, interend))
            ind += 1
        return intersected


@dataclass(frozen=True)
class ItemMeta:
    target: str
    sampling: str
    peaks: PeaksIndex
    treatment: Tuple[BigWigMeta, ...]
    control: Tuple[BigWigMeta, ...]
    fulltreatment: Tuple[BigWigMeta, ...]
    fullcontrol: Tuple[BigWigMeta, ...]

    @staticmethod
    def make(sampling, replica: ReplicaMeta):
        bymode = lambda mode, files: tuple(f[mode] for f in files)
        peaks = BedTool(replica.peaks).sort().merge().sort()
        kwargs = {
            "target": replica.target, "sampling": sampling, "peaks": PeaksIndex(peaks),
            "fulltreatment": bymode("original", replica.treatment),
            "fullcontrol": bymode("original", replica.control),
            "treatment": bymode(sampling, replica.treatment),
            "control": bymode(sampling, replica.control),
        }
        return ItemMeta(**kwargs)
