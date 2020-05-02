import numpy as np
import pyBigWig
import utils
from torch.utils.data import Dataset
from data.pod import IntervalMeta, ReplicaMeta, BigWigMeta
from .misc import ItemMeta
from tqdm import tqdm
from typing import Callable, Generator, Tuple
from pybedtools import BedTool, Interval
from itertools import groupby, chain


class ChIPseq(Dataset):
    def __init__(self, intervals: BedTool, replicas: Tuple[ReplicaMeta, ...], binsize: int,
                 simmods: Tuple[str, ...] = ("q0.25", "q0.5", "q0.75")):
        assert binsize >= 1
        self.intaug = None
        self.bigwigaug = None
        self.intervals: [Interval] = list(intervals.sort())
        self.binsize = binsize

        # mapping = [[0, 1, 2, 3], [1,3,5], [1,2,5]] for each subsampled experiment -> indices of the original intervals
        # meta = [('q0.25', BigWigMeta), ('q0.5', BigWigMeta), ('q0.75', BigWigMeta)]
        self.mapping: [np.ndarray] = []
        self.meta: [ItemMeta] = []
        for r in replicas:
            for q in simmods:
                self.meta.append(ItemMeta.make(q, r))
                self.mapping.append(np.arange(len(intervals)))

        # bigwig = {"filename1": bigwig", "filename2": bigwig, ...}
        bigwig = set(chain(*[x.treatment + x.control + x.fulltreatment + x.fullcontrol for x in self.meta]))
        self.bigwig: {BigWigMeta: pyBigWig.pyBigWig} = {bw: pyBigWig.open(bw.path) for bw in bigwig}

        # for fast search of the flattened index, cumsum of lengths
        self.index = np.empty(len(self.mapping), dtype=np.uint32)
        self.reindex()
        self.reread_bigwig = False

    def reindex(self):
        newmeta, newmapping = [], []
        for meta, mapping in zip(self.meta, self.mapping):
            if len(mapping) > 0:
                newmeta.append(meta)
                newmapping.append(mapping)
        self.meta, self.mapping = newmeta, newmapping
        self.index[:] = [len(x) for x in self.mapping]
        self.index = np.cumsum(self.index, dtype=np.uint32)

    def add_augmentation(self, intervals: Callable[[Interval], Interval],
                                bigwig):
        self.intaug = intervals
        self.bigwigaug = bigwig

    def __len__(self):
        return self.index[-1]

    def reread_bigwig_on_fetch(self):
        self.reread_bigwig = True

    def remove(self, indices):
        # For a given experiment all indices MUST be removed at once
        indices = [self._get_interval(ind) for ind in indices]
        indices = sorted(indices, key=lambda x: x[0])
        visited = set()
        for exp_ind, indices in groupby(indices, key=lambda x: x[0]):
            assert exp_ind not in visited
            visited.add(exp_ind)
            interval_indices = [x[1] for x in indices]
            self.mapping[exp_ind] = np.delete(self.mapping[exp_ind], interval_indices)
        self.reindex()

    def _get_interval(self, global_index):
        """global index -> experiment index, experiment interval index"""
        ind = np.searchsorted(self.index, global_index, side='right')
        interval_ind = global_index - self.index[ind - 1] if ind != 0 else global_index
        return ind, interval_ind

    def _get_reads(self, filesmeta: Tuple[BigWigMeta, ...], interval: Interval) -> {BigWigMeta: np.ndarray}:
        chrom, start, end = interval.chrom, interval.start, interval.end
        return {meta: self.bigwig[meta].values(chrom, start, end, numpy=True) for meta in filesmeta}

    def _cat_normalize(self, reads: {BigWigMeta: np.ndarray}, bigwigs: [BigWigMeta]):
        # It is not clear how to do this properly, hence use simple approximation (x+y) / (X+Y)
        total = sum(bw.readlen * bw.numreads for bw in bigwigs)
        return np.sum([reads[bw] for bw in bigwigs], axis=0) / total

    def itermeta(self, annotation: {str: BedTool}) -> Generator[IntervalMeta, None, None]:
        # 1. compute per region intersection with annotations
        intervals = BedTool(self.intervals).saveas()  # already sorted
        annocov = {key: utils.bed.coverage(intervals, bed, presorted=True) for key, bed in annotation.items()}

        # 2. compute intersection between all peaks and intervals
        peakscov = set(x.peaks.bedfn for x in self.meta)
        peakscov = {fn: utils.bed.coverage(intervals, BedTool(fn), presorted=True) for fn in peakscov}

        # 3. iterate and build interval meta
        for item in tqdm(range(len(self)), total=len(self)):
            ind, interval_ind = self._get_interval(item)
            meta = self.meta[ind]
            mapped_raw_interval_ind = self.mapping[ind][interval_ind]
            interval = self.intervals[mapped_raw_interval_ind]
            cov = {key: cov[mapped_raw_interval_ind] for key, cov in annocov.items()}
            pcov = peakscov[meta.peaks.bedfn][mapped_raw_interval_ind]

            # fetch raw reads, cat them and
            allbw = meta.treatment + meta.control
            reads = self._get_reads(allbw, interval)
            treatment, control = self._cat_normalize(reads, meta.treatment), \
                                 self._cat_normalize(reads, meta.control)
            ratio = treatment / (treatment + control + np.finfo(np.float32).eps)
            std = ratio.std()
            yield IntervalMeta(pcov, std, meta.target, meta.sampling, cov)

    def __getitem__(self, item):
        if self.reread_bigwig:
            self.bigwig = {meta: pyBigWig.open(meta.path) for meta in self.bigwig.keys()}
            self.reread_bigwig = False

        assert 0 <= item < self.index[-1]
        ind, interval_ind = self._get_interval(item)
        meta = self.meta[ind]
        interval = self.intervals[self.mapping[ind][interval_ind]]

        if self.intaug:
            interval = self.intaug(interval)

        allbw = tuple(set(meta.treatment + meta.control + meta.fullcontrol + meta.fulltreatment))
        reads = self._get_reads(allbw, interval)

        if self.bigwigaug:
            reads = self.bigwigaug(meta.treatment, meta.control, meta.fulltreatment, meta.fullcontrol, reads)

        # 1. Cat and normalize each track by the number of total reads * readlen.
        treatment, fulltreatment, control, fullcontrol = self._cat_normalize(reads, meta.treatment), \
                                                         self._cat_normalize(reads, meta.fulltreatment), \
                                                         self._cat_normalize(reads, meta.control), \
                                                         self._cat_normalize(reads, meta.fullcontrol)

        # 3. Calculate the ratio to predict
        eps = np.finfo(np.float32).eps
        ratio = treatment / (treatment + control + eps)
        refratio = fulltreatment / (fulltreatment + fullcontrol + eps)

        # 4. Binning
        ratio, refratio = ratio.astype(np.float32), refratio.astype(np.float32)
        assert ratio.shape == refratio.shape and ratio.ndim == 1
        nbins = ratio.size // self.binsize
        ratio, refratio = ratio.reshape(nbins, -1).mean(axis=-1), refratio.reshape(nbins, -1).mean(axis=-1)

        # 5. reference peaks
        intersection = meta.peaks.intersect(interval)
        if len(intersection) == 0:
            peaks = np.zeros_like(ratio)
        else:
            peaks = np.zeros(interval.length, dtype=np.float32)
            start = interval.start
            for inter in intersection:
                b, e = inter.start - start, inter.end - start
                peaks[b:e] = 1
            peaks = peaks.reshape(nbins, -1).mean(axis=-1)

        assert ratio.shape == refratio.shape == peaks.shape and \
               ratio.dtype == refratio.dtype == peaks.dtype == np.float32
        return {"ratio": ratio, "refratio": refratio, "peaks": peaks}

    def __getstate__(self):
        state = self.__dict__.copy()
        state["bigwig"] = list(state["bigwig"].keys())
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.bigwig = {bw: pyBigWig.open(bw.path) for bw in self.bigwig}




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
