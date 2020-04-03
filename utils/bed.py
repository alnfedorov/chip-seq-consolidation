from pybedtools import BedTool, Interval
from typing import Iterable, Callable


def compute_conservative_regions(bed: Iterable[BedTool], threshold: int,
                                 merge: Callable[[Iterable[Interval], Interval], Interval]) -> BedTool:
    """
    :param bed: Original bed files
    :param threshold: How many votes do we need to consider region a conservative one?
    :param merge: Function to transfer meta information from voted intervals -> new interval
    :return:
    """
    # 1. Create intervals that have the same number of full hits with respect to the original bed files.
    # 2. Threshold regions by the number of hits
    # |▓▓▓| |▓▓▓▓▓▓▓▓▓|
    # |▓▓|    |▓▓▓▓▓▓▓▓▓|
    # |▓|         |▓▓▓▓▓▓▓▓▓|
    # -----------------------
    # |▓|▓| |▓▓|▓▓|▓▓▓|▓▓|▓▓|
    # |3|2| |1 |2 |3  |2 |1 |
    result = []

    # In parallel loop over bed intervals
    # At each step select subinterval and count hits

    iterators = [iter(b.sort()) for b in bed]
    intervals = [(ind, next(b_iter)) for ind, b_iter in enumerate(iterators)]

    while intervals:
        boundaries = sorted(set(
            [(inter.chrom, inter.start) for _, inter in intervals] + [(inter.chrom, inter.end) for _, inter in
                                                                      intervals]
        ))  # [(chr, boundary]), ...]
        schrom, start = boundaries[0]
        echrom, end = boundaries[1]
        assert schrom == echrom

        hits = []
        for _, inter in intervals:
            if inter.chrom == schrom and inter.start <= start and end <= inter.end:
                hits.append(inter)
        assert len(hits) > 0
        if len(hits) >= threshold:
            consolidated = Interval(schrom, start, end)
            consolidated = merge(hits, consolidated)
            result.append(consolidated)

        # Push interval forward
        new_intervals = []
        for ind, inter in intervals:
            if inter.chrom != schrom:
                new_intervals.append((ind, inter))
                continue
            assert inter.start >= start

            if start == inter.start:
                assert end <= inter.end
                inter.start = end

            # request next interval
            if inter.start == inter.end:
                try:
                    inter = next(iterators[ind])
                    new_intervals.append((ind, inter))
                except StopIteration:
                    continue
            else:
                new_intervals.append((ind, inter))
        intervals = new_intervals
    return BedTool(result).sort()
