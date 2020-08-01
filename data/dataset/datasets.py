from itertools import groupby
from typing import Dict, List
from typing import Tuple

import numpy as np
from pybedtools import BedTool, Interval
from torch.utils.data import Dataset

import utils
from data.pod import SeparatedReadsMeta
from utils import augreads
from .misc import ChIPseqReplica


class ChIPseqReplicaTrainDataset(Dataset):
    def __init__(self, target: str, intervals: BedTool, peaks: str, treatment: SeparatedReadsMeta,
                 control: SeparatedReadsMeta, chromnfo: Dict[str, int], fragment_sizes: Tuple[int, ...],
                 binsize: int, effective_genome_size: int):
        assert binsize >= 1 and effective_genome_size > 1
        self.target = target
        self.intervals: List[Interval, ...] = list(intervals.sort())
        self._dataset = ChIPseqReplica(treatment, control, chromnfo, effective_genome_size, peaks)
        self.fragments_sizes = fragment_sizes
        self.binsize = binsize

    def __len__(self):
        return len(self.intervals)

    def remove(self, indices):
        lenbefore = len(self)
        for ind in sorted(indices, reverse=True):
            del self.intervals[ind]
        assert lenbefore - len(self) == len(indices)

    def __getitem__(self, idx):
        assert 0 <= idx < len(self)
        interval = self.intervals[idx]

        # 1. reference enrichment
        reffragment_size = sum(self.fragments_sizes) // len(self.fragments_sizes)
        reftreatment, refcontrol, refenrichment = self._dataset.enrichment(interval, reffragment_size)

        # 2. augment treatment reads
        shift_limits = (-75, 75)  # randomly shift by half of the nucleosome size
        subsampleto = augreads.randreads(maxreads=self._dataset.treatment.numreads)
        treatment = augreads.subsample(reftreatment, p=subsampleto / self._dataset.treatment.numreads)

        # mixin treatment and control to have treatment / (treatment + control) = enrichfrac
        enrichfrac = np.random.uniform(0.65, 1)
        treatment = augreads.lowenrichment(treatment, refcontrol, enrichfrac)
        treatment = augreads.shift(treatment, shift_limits)

        # 3. augment control reads
        if subsampleto < self._dataset.control.numreads:
            subsampleto = augreads.randreads(minreads=subsampleto, maxreads=self._dataset.control.numreads)
            control = augreads.subsample(refcontrol, p=subsampleto / self._dataset.control.numreads)
        else:
            control = augreads.subsample(refcontrol, p=subsampleto / self._dataset.control.numreads, replace=True)
        control = augreads.shift(control, shift_limits)

        # 4. calculate the enrichment
        # TODO: Use log-normal distribution
        # https://en.wikipedia.org/wiki/Log-normal_distribution
        # https://github.com/fnaumenko/isChIP
        # https://github.com/fnaumenko/bioStat/blob/master/pict/PEdistribs_medium.png
        minfs, maxfs = min(self.fragments_sizes), max(self.fragments_sizes)
        fragment_size = np.random.randint(min(minfs, 150), max(250 + 1, maxfs + 1))

        chromlen = self._dataset.chrominfo[interval.chrom]
        enrich_ends, enrich_values = utils.enrichment.endtoend(chromlen, fragment_size,
                                                               self._dataset.effective_genome_size,
                                                               self._dataset.treatment.numreads,
                                                               self._dataset.control.numreads,
                                                               treatment.forward, treatment.reverse,
                                                               control.forward, control.reverse)
        enrichment = utils.enrichment.todense(enrich_ends, enrich_values, interval.start, interval.end)

        assert enrichment.shape == refenrichment.shape and refenrichment.ndim == 1
        nbins = interval.length // self.binsize

        dobinning = lambda x: x.reshape(nbins, -1).mean(axis=-1)
        enrichment, refenrichment = dobinning(enrichment), dobinning(refenrichment)

        peaks = dobinning(self._dataset.peaks(interval))

        assert enrichment.shape == refenrichment.shape == peaks.shape and \
               enrichment.dtype == refenrichment.dtype == peaks.dtype == np.float32

        ohepeaks = np.empty((2, peaks.size), dtype=np.float32)
        ohepeaks[0, :] = 1 - peaks
        ohepeaks[1, :] = peaks
        return {"enrichment": enrichment, "refenrichment": refenrichment, "peaks": ohepeaks}


class ChIPseqReplicaValDataset(Dataset):
    def __init__(self, uid: str, experiment_accession: str, target: str, sampling: str, intervals: BedTool, peaks: str,
                 reftreatment: SeparatedReadsMeta, refcontrol: SeparatedReadsMeta,
                 reffragment_size: int,
                 sampledtreatment: SeparatedReadsMeta, sampledcontrol: SeparatedReadsMeta,
                 sampledfragment_size: int,
                 chrominfo: Dict[str, int],
                 binsize: int,
                 effective_genome_size: int):
        assert binsize >= 1 and effective_genome_size > 1
        self.uid = uid
        self.experiment_accession = experiment_accession
        self.target = target
        self.sampling = sampling
        self.intervals: List[Interval, ...] = list(intervals.sort())
        self._reference = ChIPseqReplica(reftreatment, refcontrol, chrominfo, effective_genome_size, peaks)
        self._sampled = ChIPseqReplica(sampledtreatment, sampledcontrol, chrominfo, effective_genome_size)
        self._reffragment_size = reffragment_size
        self._sampledfragment_size = sampledfragment_size
        self.binsize = binsize

    def __len__(self):
        return len(self.intervals)

    def remove(self, indices):
        lenbefore = len(self)
        for ind in sorted(indices, reverse=True):
            del self.intervals[ind]
        assert lenbefore - len(self) == len(indices)

    def __getitem__(self, idx):
        assert 0 <= idx < len(self)
        interval = self.intervals[idx]

        *_, refenrichment = self._reference.enrichment(interval, self._reffragment_size)
        *_, enrichment = self._sampled.enrichment(interval, self._sampledfragment_size)

        assert enrichment.shape == refenrichment.shape and refenrichment.ndim == 1
        nbins = interval.length // self.binsize

        dobinning = lambda x: x.reshape(nbins, -1).mean(axis=-1)
        enrichment, refenrichment = dobinning(enrichment), dobinning(refenrichment)

        peaks = dobinning(self._reference.peaks(interval))

        assert enrichment.shape == refenrichment.shape == peaks.shape and \
               enrichment.dtype == refenrichment.dtype == peaks.dtype == np.float32

        ohepeaks = np.empty((2, peaks.size), dtype=np.float32)
        ohepeaks[0, :] = 1 - peaks
        ohepeaks[1, :] = peaks
        return {"enrichment": enrichment, "refenrichment": refenrichment, "peaks": ohepeaks}


class ConcatChIPseqDataset(Dataset):
    def __init__(self, datasets):
        super().__init__()
        assert len(datasets) > 0 and all(isinstance(dst, ChIPseqReplicaTrainDataset) for dst in datasets)
        self.datasets: List[ChIPseqReplicaTrainDataset, ...] = list(datasets)
        self.cumsizes = None
        self._reindex()

    def _reindex(self):
        lengths = [len(dst) for dst in self.datasets]
        # drop dataset if no instances are left
        todrop = [ind for ind, length in enumerate(lengths) if length == 0]
        if len(todrop) != 0:
            self.datasets = np.delete(self.datasets, todrop)
            lengths = np.delete(self.datasets, todrop)
        assert all(length > 0 for length in lengths)
        self.cumsizes = np.cumsum(lengths, dtype=np.uint32)

    def unroll_idx(self, idx):
        if idx < 0:
            if -idx > len(self):
                raise ValueError("absolute value of index should not exceed dataset length")
            idx = len(self) + idx

        dataset_idx = np.searchsorted(self.cumsizes, idx, side="right")
        if dataset_idx == 0:
            sample_idx = idx
        else:
            sample_idx = idx - self.cumsizes[dataset_idx - 1]
        return dataset_idx, sample_idx

    def remove(self, indices):
        unrolled = [self.unroll_idx(idx) for idx in indices]
        unrolled = sorted(unrolled, key=lambda x: x[0])
        visited_dst = set()
        for dataset_idx, dataset_unrolled in groupby(unrolled, key=lambda x: x[0]):
            assert dataset_idx not in visited_dst
            visited_dst.add(dataset_idx)
            samples_idx = [x[1] for x in indices]
            self.datasets[dataset_idx].remove(samples_idx)
        self._reindex()

    def __len__(self):
        return self.cumsizes[-1]

    def __getitem__(self, idx):
        dataset_idx, sample_idx = self.unroll_idx(idx)
        return self.datasets[dataset_idx][sample_idx]
