import numpy as np
from random import choice
from data.pod import IntervalMeta
from collections import defaultdict
from torch.utils.data import BatchSampler


class BalancedBatchSampler(BatchSampler):
    def __init__(self, meta: [IntervalMeta],
                 batch_size: int = 128, batches_per_epoch: int = None, drop_last: bool = True):
        """
        Sampler that tries to address all the diversity of the DNA intervals in a balanced way.
        First, every interval is associated with the following labels:
            * top most annotation group
            * peaks coverage quartile
            * reads ratio std quartile
            * target (given)
            * sampling quartile (given)
        peaks coverage and reads ratio quartiles are computed based on all the available data
        After that all intervals with the same set of labels is grouped together.
        Sampling is a 2 step procedure: 1-sample group, 2-sample interval inside group.
        """
        batches_per_epoch = batches_per_epoch if batches_per_epoch else len(meta) // batch_size
        assert batch_size >= 1 and batches_per_epoch >= 1
        self.batch_size = batch_size
        self.batches_per_epoch = batches_per_epoch
        self.drop_last = drop_last

        index: {tuple: [int]} = self._make_balanced_index(meta)
        ulabels = tuple(index.keys())
        assert len(ulabels) == len(set(ulabels))
        self.weights = self._rebalance(ulabels)
        self.groups = [index[lbls] for lbls in ulabels]
        self.total_samples = len(meta)

        # Workaround to support dummy Engine in the ignite
        self.sampler = None

    def _make_balanced_index(self, meta: [IntervalMeta]):
        # features to its quantiles
        quantiles = (0.25, 0.5, 0.75, 1)
        def toquantiles(values) -> np.ndarray:
            q = np.quantile(values, quantiles)
            q = np.asarray(sorted((set(q))))    # exclude repeated quantiles
            qv = np.searchsorted(q, values, side="left")
            assert qv.max() < len(q)
            return qv

        peaksq = toquantiles(np.asarray([x.peaks for x in meta]))
        readsq = toquantiles(np.asarray([x.ratio_std for x in meta]))

        annkeys = list(meta[0].annotation.keys())
        annotation = {k: np.asarray([x.annotation[k] for x in meta]) for k in annkeys}
        annotation = {k: toquantiles(v) for k, v in annotation.items()}
        assert peaksq.size == readsq.size and all(annotation[k].size == peaksq.size for k in annkeys)

        index = defaultdict(list)
        for ind, m in enumerate(meta):
            key = (
                peaksq[ind], readsq[ind], *(annotation[k][ind] for k in annkeys), m.sampling, m.target
            )
            index[key].append(ind)
        return dict(index)

    def _rebalance(self, ulabels):
        # How to balance this stuff?
        # I want for each subcategory of the single class to be equally probable, no matter other categories.
        # P(x==a) == P(x==b) == P(x==c) ....
        # And each multi-class category looks like P(x==a)*P(y==b)*P(z==c)..etc
        # Hence, simply compute sampling P for each class and multiply probabilities in the end
        assert all(len(u) == len(ulabels[0]) for u in ulabels)
        counts = [defaultdict(int) for _ in range(len(ulabels[0]))]
        for lbls in ulabels:
            for cls, subclass in enumerate(lbls):
                counts[cls][subclass] += 1
        probs = [{k: v / sum(cls.values()) for k, v in cls.items()} for cls in counts]
        assert all(abs(sum(cls.values()) - 1) < 1e-8 for cls in probs)

        newweights = []
        for lbls in ulabels:
            p = 1
            for cls, subclass in enumerate(lbls):
                p *= probs[cls][subclass]
            newweights.append(p)
        # newweights are NOT expected to sum to 1 as not all combinations are covered
        # reweight for simplicity
        newweights = np.asarray(newweights)
        newweights = newweights / np.linalg.norm(newweights, ord=1)
        return newweights

    def __iter__(self):
        grps = np.random.choice(np.arange(len(self.groups)), self.total_samples,
                                replace=True, p=self.weights)
        assert grps.max() < len(self.groups)
        for ind in range(self.batches_per_epoch):
            batch_grps = grps[ind*self.batch_size: (ind+1)*self.batch_size]
            batch_indices = [choice(self.groups[grp]) for grp in batch_grps]
            yield batch_indices

        if len(self.groups) / self.batch_size != self.batches_per_epoch and not self.drop_last:
            batch_grps = grps[self.batches_per_epoch*self.batch_size:]
            assert len(batch_grps) < self.batch_size
            batch_indices = [choice(self.groups[grp]) for grp in batch_grps]
            yield batch_indices

    def __len__(self):
        if self.drop_last:
            return self.total_samples // self.batch_size
        else:
            return (self.total_samples + self.batch_size - 1) // self.batch_size
