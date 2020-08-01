import os
import numpy as np
from typing import Tuple, Dict
from dataclasses import dataclass


@dataclass(frozen=True)
class BamMeta:
    path: str
    estimated_fragment_size: int
    numreads: int
    readlen: int

    def __post_init__(self):
        assert os.path.isfile(self.path)

    def __repr__(self):
        return f"BamMeta(numreads={self.numreads}, readlen={self.readlen}, path={os.path.split(self.path)[-1]})"


@dataclass(frozen=True)
class SeparatedReadsMeta:
    forward: Dict[str, str]
    reverse: Dict[str, str]
    readlen: int
    numreads: int

    def __post_init__(self):
        assert set(self.forward.keys()) == set(self.reverse.keys())
        assert all(os.path.isfile(x) for x in set(self.forward.values()).union(self.reverse.values()))

    def __repr__(self):
        return f"SeparatedReadsMeta(readlen={self.readlen}, numreads={self.numreads})"


@dataclass(frozen=False)
class IntervalReads:
    forward: np.ndarray
    reverse: np.ndarray

    @property
    def numreads(self):
        return self.forward.size + self.reverse.size


@dataclass(frozen=True)
class AccessionInfo:
    experiment: str
    treatment: Tuple[str, ...]
    control: Tuple[str, ...]


@dataclass(frozen=True)
class ChIPseqReplicaMeta:
    accesions: AccessionInfo
    target: str
    peaks: str
    # sampling -> bam meta
    samplemods: Tuple[str, ...]
    treatment: Dict[str, Tuple[BamMeta, ...]]
    control: Dict[str, Tuple[BamMeta, ...]]

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
