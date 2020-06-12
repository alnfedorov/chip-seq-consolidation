import numpy as np
from data.pod import IntervalReads


def randreads(*, minreads=1e6, maxreads=None):
    return np.random.randint(minreads, maxreads + 1)


def subsample(reads: IntervalReads, before: int, after: int, replace=False):
    p = after / before
    assert p <= 1 or replace
    forward = np.random.choice(reads.forward, int(round(reads.forward.size * p)), replace=replace)
    reverse = np.random.choice(reads.reverse, int(round(reads.reverse.size * p)), replace=replace)
    return IntervalReads(forward, reverse)


def lowenrichment(treatment: IntervalReads, treatment_numreads: int,
                  control: IntervalReads, control_numreads: int, noise: float):
    noise_numreads = noise * treatment_numreads / (1 - noise)
    noise_numreads = int(round(noise_numreads))
    maxreads_by_noise = noise_numreads + treatment_numreads

    total_numreads = randreads(maxreads=maxreads_by_noise)
    signal_numreads = int(round(total_numreads * (1 - noise)))
    noise_numreads = int(round(total_numreads * noise))

    signal = subsample(treatment, before=treatment_numreads, after=signal_numreads)
    if noise_numreads > control_numreads:
        noise = subsample(control, before=control_numreads, after=noise_numreads, replace=True)
    else:
        noise = subsample(control, before=control_numreads, after=noise_numreads, replace=False)

    lowenrich = IntervalReads(
        forward=np.concatenate([signal.forward, noise.forward]),
        reverse=np.concatenate([signal.reverse, noise.reverse])
    )
    return total_numreads, lowenrich


def shift(reads: IntervalReads, limits):
    reads.forward += np.random.randint(*limits, size=reads.forward.size, dtype=np.int32)
    reads.reverse += np.random.randint(*limits, size=reads.reverse.size, dtype=np.int32)
    return reads
