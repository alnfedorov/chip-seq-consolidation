from pybedtools import BedTool, Interval
from typing import Iterable, Callable
from math import log10


def _default_merge(intervals: Iterable[Interval], new_interval: Interval) -> Interval:
    return new_interval


def coverage(source: BedTool, to_intersect: BedTool, presorted: bool = True) -> [float]:
    """
    For each interval in the source compute coverage by to_intersect features. Source intervals must be non overlapping
    """
    if not presorted:
        source, to_intersect = source.sort(), to_intersect.sort()
    intersection = source.intersect(to_intersect, wao=True)
    intersection, source = list(intersection), list(source)

    ind = 0
    curinter = intersection[ind]
    curcov = float(curinter.fields[-1])
    coverage = []
    assert source[ind] == curinter, f"Fail: expected {source[ind]}, got {curinter}"
    for inter in intersection:
        if inter == curinter:
            curcov += float(inter.fields[-1])
        else:
            coverage.append(curcov / curinter.length)
            assert coverage[-1] <= 1
            curinter = inter
            curcov = float(curinter.fields[-1])
            ind += 1
            assert curinter == source[ind], f"Fail: expected {source[ind]}, got {curinter}"
    coverage.append(curcov / curinter.length)
    assert len(coverage) == len(source)
    return coverage


def consensus(bed: Iterable[BedTool], weights: [int], threshold: float,
              merge: Callable[[Iterable[Interval], Interval], Interval] = _default_merge) -> BedTool:
    """
    :param bed: Original bed files
    :param weights: vote's weights, pass ones if you are not sure
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
    assert len(weights) == len(iterators), "len(regions) != len(weights)"
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
        for ind, inter in intervals:
            if inter.chrom == schrom and inter.start <= start and end <= inter.end:
                hits.append(weights[ind])
        assert len(hits) > 0
        if sum(hits) >= threshold:
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
    return BedTool(result).sort().merge()


def threshold_qvalue(bed: BedTool, qvalue: float = 0.05, column=8, verbose: bool = True):
    """
    Threshold bed file entries by qvalue.
    NOTE: bed files must be in the standard narrow or broad peak format and include -log10(qvalue) column
    """
    nlog10_qvalue = -log10(qvalue)

    lenbefore = len(bed)
    filtered = []
    for interval in bed:
        # assert len(interval.fields) in (20, 18, 10, 9), f"{bed.fn}-{len(interval.fields)}"
        value = float(interval.fields[column])
        assert value >= 0
        if value >= nlog10_qvalue:
            filtered.append(interval)
    filtered = BedTool(filtered).sort()
    lenafter = len(filtered)
    if verbose:
        print(f"Filtered {lenbefore - lenafter} intervals")
    return filtered
