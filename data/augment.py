import numpy as np
from typing import Tuple
from .dataset.misc import BigWigMeta
from pybedtools import Interval
from .hg19.annotation import CHROMOSOME_SIZE


def shift(interval: Interval, limits: (float, float)) -> Interval:
    mins, maxs = limits[0] * interval.length, limits[1] * interval.length
    shift = np.random.randint(int(mins), int(maxs))
    # dont get off the chromosome
    chrom = interval.chrom
    chromsize = CHROMOSOME_SIZE[chrom]
    shift = min(shift, chromsize - interval.end)
    shift = max(shift, -interval.start)
    result = Interval(chrom, interval.start + shift, interval.end + shift)
    assert result.start >= 0 and result.end <= chromsize
    return result


def randreads(treatment: Tuple[BigWigMeta, ...], control: Tuple[BigWigMeta, ...],
              fulltreatment: Tuple[BigWigMeta, ...], fullcontrol: Tuple[BigWigMeta, ...],
              reads: {BigWigMeta: np.ndarray}, random_reads_fraction: float) -> {BigWigMeta: np.ndarray}:
    for meta in treatment + control:
        oreads = reads[meta]
        readscount = int(oreads.sum() / meta.readlen)
        togenerate = int(round(readscount * random_reads_fraction))
        if togenerate == 0:
            continue
        for strt in np.random.randint(0, len(oreads), size=togenerate):
            oreads[strt: strt + meta.readlen] += 1
    return reads
