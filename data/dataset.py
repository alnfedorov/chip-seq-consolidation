import numpy as np
import pyBigWig
from torch.utils.data import Dataset
from .misc import IntervalMeta
from pybedtools import BedTool, Interval
from typing import TypedDict, Callable, Generator
from itertools import groupby

SubsampledBigWig = TypedDict("SubsampledBigWig",
                             {"q0.25": str, "q0.5": str, "q0.75": str, "original": str})


class BigWigMeta(TypedDict):
    name: str
    target: str
    treatment: SubsampledBigWig
    control: SubsampledBigWig
    peaks: BedTool


class ChIPseq(Dataset):
    def __init__(self, intervals: BedTool, experiments: [BigWigMeta], binsize: int, pseudocounts: int = 1):
        assert binsize >= 1
        self.intaug = None
        self.bigwigaug = None
        self.intervals: [Interval] = list(intervals.sort())
        self.binsize = binsize
        self.pseudocounts: int = pseudocounts

        # mapping = [[0, 1, 2, 3], [1,3,5], [1,2,5]] for each subsampled experiment -> indices of the original intervals
        # meta = [('q0.25', BigWigMeta), ('q0.5', BigWigMeta), ('q0.75', BigWigMeta)]
        self.mapping: [np.ndarray] = []
        self.meta: [(str, BigWigMeta)] = []
        for exp in experiments:
            for q in ("q0.25", "q0.5", "q0.75"):
                assert q in exp['treatment'] and exp['treatment'][q] is not None
                self.meta.append((q, exp))
                self.mapping.append(np.arange(len(intervals)))

        # bigwig = {"filename1": bigwig", "filename2": bigwig, ...}
        self.bigwig: {str: pyBigWig.pyBigWig} = {}
        for exp in experiments:
            for q in ("q0.25", "q0.5", "q0.75", "original"):
                for mode in ("treatment", "control"):
                    filename = exp[mode][q]
                    self.bigwig[filename] = pyBigWig.open(filename)

        # for fast search of the flattened index, cumsum of lengths
        self.index = np.empty(len(self.mapping), dtype=np.int32)
        self.reindex()

    def reindex(self):
        newmeta, newmapping = [], []
        for meta, mapping in zip(self.meta, self.mapping):
            if len(mapping) > 0:
                newmeta.append(meta)
                newmapping.append(mapping)
        self.meta, self.mapping = newmeta, newmapping
        self.index[:] = [len(x) for x in self.mapping]
        self.index = np.cumsum(self.index)

    def augmentation(self, intervals: Callable[[Interval], Interval],
                     bigwig: Callable[[np.ndarray, np.ndarray], (np.ndarray, np.ndarray)]):
        self.intaug = intervals
        self.bigwigaug = bigwig

    def __len__(self):
        return self.index[-1]

    def remove(self, indices):
        # For a given experiment all indices MUST be removed at once
        indices = [self._get_interval(ind) for ind in indices]
        indices = sorted(indices, key=lambda x: x[0])
        visited = set()
        for exp_ind, interval_indices in groupby(indices, key=lambda x: x[0]):
            assert exp_ind not in visited
            visited.add(exp_ind)
            interval_indices = [x[1] for x in interval_indices]
            self.mapping[exp_ind] = np.delete(self.mapping[exp_ind], interval_indices)
        self.reindex()

    def _get_interval(self, global_index):
        """global index -> experiment index, experiment interval index"""
        exp_ind = np.searchsorted(self.index, global_index, side='right')
        exp_interval_ind = global_index - self.index[exp_ind - 1] if exp_ind != 0 else global_index
        return exp_ind, exp_interval_ind

    def _get_reads(self, filenames: str, interval: Interval) -> [np.ndarray]:
        # BigWig uses half open intervals, we assumes closed intervals
        chrom, start, end = interval.chrom, interval.start, interval.end + 1
        nbins = (end - start) // self.binsize
        return [self.bigwig[filename].stats(chrom, start, end, type="mean", nBins=nbins) for filename in filenames]

    def itermeta(self, annotation: {str: BedTool}) -> Generator[IntervalMeta, None, None]:
        # 1. compute per region intersection with annotations
        intervals = BedTool(self.intervals)  # already sorted
        intersections = {key: list(intervals.intersect(bed, wao=True)) for key, bed in annotation.items()}
        coverage: [{str: float}] = []
        for ind, interval in enumerate(intervals):
            cov = {key: intervals[ind] for key in intersections}
            assert all(interval.chrom == x.chrom and
                       interval.start == x.start and
                       interval.end == x.end
                       for x in cov.values())
            coverage.append({key: int(w.fields[-1]) / w.length for key, w in cov.items()})

        # 2. iterate
        for ind in range(len(self)):
            exp_ind, exp_interval_ind = self._get_interval(ind)
            q, exp = self.meta[exp_ind]
            interval = self.intervals[self.mapping[exp_ind][exp_interval_ind]]
            treatment, control = [
                self._get_reads(filename, interval) for filename in (exp["treatment"][q], exp["control"][q])
            ]
            std = ((treatment + self.pseudocounts) / (control + self.pseudocounts)).std()
            meta = IntervalMeta(0., std, target=exp['target'], quartile=q,
                                annotation=coverage[self.mapping[exp_ind][exp_interval_ind]])
            yield meta

    def __getitem__(self, item):
        assert 0 < item < self.index[-1]
        exp_ind, exp_interval_ind = self._get_interval(item)
        q, exp = self.meta[exp_ind]
        interval = self.intervals[self.mapping[exp_ind][exp_interval_ind]]

        if self.intaug:
            interval = self.intaug(interval)
        treatment, control, fulltreatment, fullcontrol = [
            self._get_reads(filename, interval) for filename in
            (exp["treatment"][q], exp["control"][q], exp["treatment"]["orignal"], exp["control"]["orignal"])
        ]

        if self.bigwigaug:
            treatment, control = self.bigwigaug(treatment, control)

        refratio = (fulltreatment + self.pseudocounts) / (fullcontrol + self.pseudocounts)
        ratio = (treatment + self.pseudocounts) / (control + self.pseudocounts)

        # Peaks
        # peaks = exp['peaks'].intersect(BedTool(interval))
        print("PEAKS ARE NOT READY YET!!!!")
        return {"ratio": ratio, "refratio": refratio}


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
