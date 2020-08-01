import numpy as np

from data.pod import IntervalReads


def randreads(*, minreads=1e6, maxreads=30e6):
    return np.random.randint(minreads, maxreads + 1)


def subsample(reads: IntervalReads, p: float, replace=False):
    assert p <= 1 or replace, f"p={p}, replace={replace}"
    forward = np.random.choice(reads.forward, int(round(reads.forward.size * p)), replace=replace)
    reverse = np.random.choice(reads.reverse, int(round(reads.reverse.size * p)), replace=replace)
    return IntervalReads(forward, reverse)


def lowenrichment(treatment: IntervalReads, control: IntervalReads, enrichfrac: float):
    # treatment / (treatment + control) = x
    # treatment = x * treatment + x * control
    # treatment * (1-x) = x * control
    # control = treatment * (1 - x) / x
    if treatment.numreads == 0 or control.numreads == 0:
        return treatment

    tosample = treatment.numreads * (1 - enrichfrac) / enrichfrac
    control = subsample(control, p=tosample / control.numreads, replace=control.numreads < tosample)

    lowenrich = IntervalReads(
        forward=np.concatenate([treatment.forward, control.forward]),
        reverse=np.concatenate([treatment.reverse, control.reverse])
    )
    if treatment.numreads > 100 and control.numreads > 100:
        assert abs(treatment.numreads / lowenrich.numreads - enrichfrac) < 0.01

    # subsample to the original depth
    lowenrich = subsample(lowenrich, p=treatment.numreads / lowenrich.numreads)
    return lowenrich


def shift(reads: IntervalReads, limits):
    reads.forward += np.random.randint(*limits, size=reads.forward.size, dtype=np.int32)
    reads.reverse += np.random.randint(*limits, size=reads.reverse.size, dtype=np.int32)
    return reads
