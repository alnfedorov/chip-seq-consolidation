import numpy as np
from typing import Tuple
from pybedtools import Interval
from .hg19.annotation import CHROMOSOME_SIZE


# def shift(interval: Interval, limits: (float, float)) -> Interval:
#     mins, maxs = limits[0] * interval.length, limits[1] * interval.length
#     shift = np.random.randint(int(mins), int(maxs))
#     # dont get off the chromosome
#     chrom = interval.chrom
#     chromsize = CHROMOSOME_SIZE[chrom]
#     shift = min(shift, chromsize - interval.end)
#     shift = max(shift, -interval.start)
#     result = Interval(chrom, interval.start + shift, interval.end + shift)
#     assert result.start >= 0 and result.end <= chromsize
#     return result
# TODO: horizontal flipping
