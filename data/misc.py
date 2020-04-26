from dataclasses import dataclass


@dataclass
class IntervalMeta:
    # intersection with consensus peaks in this experiment
    peaks: float
    # std of ratio between reads in the sampled treatment and control
    ratio_std: float
    # target histone modification
    target: str
    # sampling quartile
    quartile: str
    # intersection with genomic annotation(exons, UTR, introns - etc)
    annotation: {str: float}
