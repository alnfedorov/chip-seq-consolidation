import os
from typing import TypedDict, Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class BigWigMeta:
    path: str

    def __post_init__(self):
        assert os.path.isfile(self.path)


SimulatedBigWig = TypedDict(
    "SimulatedBigWig", {
        "q0.25": BigWigMeta,
        "q0.5": BigWigMeta,
        "q0.75": BigWigMeta,
        "original": BigWigMeta
    }
)


@dataclass(frozen=True)
class ChIPseqReplicaMeta:
    experiment_accession: str
    target: str
    enrichment: SimulatedBigWig
    peaks: str
    treatment: Tuple[str, ...]
    control: Tuple[str, ...]

    def __post_init__(self):
        assert os.path.exists(self.peaks)


@dataclass
class IntervalMeta:
    # intersection with consensus peaks in this experiment
    peaks: float
    # std of ratio between reads in the sampled treatment and control
    ratio_std: float
    # target histone modification
    target: str
    # sampling quartile
    sampling: str
    # intersection with genomic annotation(exons, UTR, introns - etc)
    annotation: {str: float}
