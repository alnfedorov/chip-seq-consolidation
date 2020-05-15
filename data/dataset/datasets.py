import numpy as np
import pyBigWig
import utils
from torch.utils.data import Dataset, ConcatDataset
from data.pod import IntervalMeta, ChIPseqReplicaMeta, BigWigMeta
from .misc import ChIPseqReplica
from tqdm import tqdm
from typing import Callable, Generator, Tuple, Optional
from pybedtools import BedTool, Interval
from itertools import groupby, chain
from typing import Dict, List


class ChIPseqReplicaTrainDataset(Dataset):
    def __init__(self, intervals: BedTool, meta: ChIPseqReplicaMeta, binsize: int, sampling: str):
        assert binsize >= 1
        self.intervals: [Interval] = list(intervals.sort())
        self.binsize = binsize
        self.target = meta.target
        self.experiment_accession = meta.experiment_accession
        self.sampling = sampling
        self.subsampled = ChIPseqReplica.frommeta(meta, sampling)
        self.reference = ChIPseqReplica.frommeta(meta, "original")

    def reread_bigwig_on_query(self):
        self.subsampled.reread_bigwig_on_query()
        self.reference.reread_bigwig_on_query()

    # def add_augmentation(self, intervals: Callable, bigwig: Callable):
    #     self.intaug = intervals
    #     self.bigwigaug = bigwig

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

        # if self.intaug:
        #     interval = self.intaug(interval)

        # reads = self.subsampled.reads(interval)
        # refreads = self.reference.reads(interval)

        # if self.bigwigaug:
        #     assert set(reads.keys()) != set(refreads.keys())
        #     totalreads = {k}

        enrichment = self.subsampled.enrichment(interval)
        refenrichment = self.reference.enrichment(interval)
        assert enrichment.shape == refenrichment.shape and refenrichment.ndim == 1
        nbins = enrichment.size // self.binsize

        dobinning = lambda x: x.reshape(nbins, -1).mean(axis=-1)
        enrichment, refenrichment = dobinning(enrichment), dobinning(refenrichment)

        peaks = dobinning(self.reference.peaks(interval))

        assert enrichment.shape == refenrichment.shape == peaks.shape and \
               enrichment.dtype == refenrichment.dtype == peaks.dtype == np.float32

        tozero = np.isnan(enrichment) | np.isnan(refenrichment)
        enrichment[tozero] = 0
        refenrichment[tozero] = 0
        peaks[tozero] = 0

        ohepeaks = np.empty((2, peaks.chromlen), dtype=np.float32)
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



# 1. For each experiment, what do we want to do?
# Список регионов -> (мета региона, метка -> (мета эксперимента, отсортированные регионы))
# Как мне сохранить такую же укладку данных, но выкинуть какие-либо индексы экспериментов?
# Пусть тупо массив булов - включен/выключен массив в этом эксперименте
# Как по нему индексировать? Хороший вопрос. Хочу бинарный поиск, но так просто это не зайдет.
# 1111100001111111100000111111100011111100111111111 - нужно взять n-yю единицу за минимум времени
# (0,5,5)(10,17,12) - да, можно переводить из такого представления и в него обратно.
# Но мы что? А мы ничего. Мы будем тупо брать массив интов и использовать его в качестве отображения(!!!!!)
# for exp in experiments:
#   for q in ["q0.25", "q0.5", "q0.75"]:
#       for

# Alright, there is some bool arrays. How to fast index it then?
# Find experiment should be easy - just binary search cumsum
# How to index inside experiment

# """
# Внешние интервалы передаются в датасет
# Датасет НЕ делает аугментацию, ее делают ребята в Pipeline. Или все-таки датасет это все вместе?
#
# Да плевать на этот факт!
# Что делаем? Просто датасет. Просто один класс. Который делает следующее:
# 1. Получает интервалы и BigWig файлы
# 2. Получает аугментацию отдельным методом
# 3. Имеет методы для удаления окон по регионам
#
# Кто хранит данные?
# * Датасет
#
# Кто фильтрует данные?
# * Датасет
#
# Кто сэмплирует данные?
# * Сэмплер)
#
# Интервалы -> аугментация -> BigWigInjector -> ReadsNoise ->
#
# Интервалы можно вообще просто передать в виде BedTool, а окна делать в makewindows каком-нибудь
# Где происходит аугментация? Хороший вопрос. В самом датасете, вероятно.
#
# Ок, тогда исходные окна просто передаются датасету. Или он знает о мете?
# Она мне нужна только чтобы подсчитать единожды(!) веса для каждого интервала -> мета не должна храниться в датасете.
#
# Окна создаются внешней функцией и просто передаются в датасет, у которого есть аугментер
# Аугментер, следовательно, можно будет выключать в какой-то момент.
#
# Датасет в свою очередь НЕ знает ни о какой мете. И просто применяет аугментор и посылает их дальше.
# НЕ знает он ни о какой мете
#
#
# Как тогда классы устроить?
# Я хочу простой класс - провайдер интервалов
# Простой класс с мета информацией
# Простой класс, который считает
# """
