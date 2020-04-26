import numpy as np
from random import choice
from .misc import IntervalMeta
from collections import defaultdict
from torch.utils.data import Dataset, Sampler


class BalancedBatchSampler(Sampler):
    def __init__(self, data_source: Dataset, meta: [IntervalMeta],
                 batch_size: int = 128, batches_per_epoch: int = 10**6):
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
        super().__init__(data_source)
        assert len(meta) == len(data_source)
        assert batch_size >= 1 and batches_per_epoch >= 1
        self.batch_size = batch_size
        self.batches_per_epoch = batches_per_epoch
        self.balanced_index: {tuple: [int]} = self._make_balanced_index(meta)
        self.unique_labels = tuple(set(self.balanced_index.keys()))

    def _make_balanced_index(self, meta: [IntervalMeta]):
        # compute quantiles
        quantiles = (0.25, 0.5, 0.75, 1)
        peakq = np.quantile(np.asarray([x.peaks for x in meta]), quantiles, overwrite_input=True)
        readstdq = np.quantile(np.asarray([x.ratio_std for x in meta]), quantiles, overwrite_input=True)

        def quantile(value, quantiles):
            for q in quantiles:
                if value <= q:
                    return q
            assert False, "np.quantile is broken"

        # meta -> set of labels
        index = defaultdict(list)
        for ind, m in meta:
            key = (
                quantile(m.peaks, peakq), quantile(m.ratio_std, peakq),
                max(m.annotation.items(), key=lambda x: x[1])[0], m.quartile, m.target
            )
            index[key].append(ind)
        return dict(index)

    def __iter__(self):
        for _ in range(self.batches_per_epoch):
            batch_indices = []
            for bb in range(self.batch_size):
                group = choice(self.unique_labels)
                batch_indices.append(choice(self.balanced_index[group]))
            yield batch_indices

    def __len__(self):
        return self.batches_per_epoch
