import os
import pysam
import pickle
from hashlib import md5
from itertools import chain
from multiprocessing import cpu_count
from typing import Dict, Tuple
from dataclasses import dataclass

import numpy as np
import torch
from joblib import Parallel, delayed
from pybedtools import BedTool, Interval
from torch.utils.data import DataLoader

import data
from data import hg19, dataset
from data.pod import ChIPseqReplicaMeta, IntervalMeta


def make_regions(region_size: int, ambiguity_thr: float = 0.5):
    path = os.path.join("/tmp/", f"cached.make_regions({region_size}, {ambiguity_thr}).bed")

    if not os.path.exists(path):
        regions = BedTool().window_maker(w=region_size, genome="hg19")  # type: BedTool

        # pybedtools is not very consistent when dealing with Interval objects.
        # For instance, str(interval) in some cases will return only
        # 3 fields (chr, start, end).Another time, when fields are specified explicitly, 4 fields and more are
        # printed.It is possible to invoke intersection with '3-field' Intervals and receive mixed Intervals.
        # Intersected are wrongly read as one with strand and name being interval-hit information, and
        # non-intersected are turned into intervals with default additional fields.Workaround is to recreate each
        # interval to include the same number of fields; strand must be included.
        regions = BedTool([Interval(x.chrom, x.start, x.end) for x in regions]).saveas()

        regions = regions.filter(
            lambda r: r.length == region_size and r.chrom in hg19.annotation.CHROMOSOME_SIZE
        ).saveas()
        regions = dataset.filter.by_ambiguity(regions, BedTool(hg19.annotation.AMBIGUOUS), ambiguity_thr)
        regions.saveas(path, compressed=False)
    return BedTool(path)


def parse(experiments: Tuple[str, ...]):
    assert all(os.path.isdir(s) for s in experiments)
    experiments = sorted(set(experiments))
    hsh = md5('-'.join(experiments).encode()).hexdigest()
    path = os.path.join("/tmp/", f"cached.parse(hash={hsh}).pkl")

    if not os.path.exists(path):
        replicas = Parallel(n_jobs=-1)(
            delayed(data.parse.fromjson)(root) for root in experiments
        )
        replicas = list(chain(*replicas))
        with open(path, 'wb') as file:
            pickle.dump(replicas, file)

    with open(path, 'rb') as file:
        replicas = pickle.load(file)
    return replicas


def make_dataset(intervals: str, replica: ChIPseqReplicaMeta, binsize: int,
                 sampling: str) -> Tuple[ChIPseqReplicaMeta, IntervalMeta]:
    hsh = f"{intervals}{replica}{binsize}{sampling}"
    path = os.path.join("/tmp/", f"cached.make_dataset(hash={md5(hsh.encode()).hexdigest()}).pkl")

    if not os.path.exists(path):
        dst = dataset.ChIPseqReplicaTrainDataset(BedTool(intervals), replica, binsize, sampling)

        annotation = ("3'utr", "5'utr", "exons", "introns", "upstream1k", "downstream1k", "intergenic")
        annotation = {k: hg19.annotation.REGIONS[k].fn for k in annotation}
        meta, todrop = dataset.filter.by_meta(dst, annotation, True)

        dst.remove(todrop)
        meta = np.delete(meta, todrop).tolist()

        with open(path, 'wb') as file:
            pickle.dump((dst, meta), file)

    with open(path, 'rb') as file:
        dst, meta = pickle.load(file)
        dst.reread_bigwig_on_query()
    return dst, meta


def make_trainloader(intervals: BedTool, replicas: Tuple[ChIPseqReplicaMeta, ...], binsize: int,
                     batch_size: int, allsampling: Tuple[str, ...], threads: int) -> DataLoader:
    allsampling = sorted(set(allsampling))
    hsh = intervals.fn + ''.join(sorted(str(r) for r in replicas)) + str(binsize) + str(allsampling)
    path = os.path.join("/tmp/", f"cached.make_trainloader(hash={md5(hsh.encode()).hexdigest()}).pkl")

    threads = threads if threads >= 0 else cpu_count() + threads + 1
    if not os.path.exists(path):
        result = Parallel(n_jobs=max(1, threads), verbose=1, batch_size=1)(
            delayed(make_dataset)(intervals.fn, r, binsize, s) for r in replicas for s in allsampling
        )
        dsts, metas = zip(*result)
        concatdst = dataset.ConcatChIPseqDataset(dsts)
        sampler = torch.utils.data.BatchSampler(
            torch.utils.data.RandomSampler(concatdst), batch_size=batch_size, drop_last=False
        )
        # sampler = data.sampling.BalancedBatchSampler(metas, drop_last=False)

        with open(path, 'wb') as file:
            pickle.dump((concatdst, sampler), file)

    with open(path, 'rb') as file:
        concatdst, sampler = pickle.load(file)
        for d in concatdst.datasets:  # type: dataset.ChIPseqReplicaTrainDataset
            d.reread_bigwig_on_query()

    # dst.add_augmentation(
    #     partial(data.augment.shift, limits=(-0.5, 0.5)),
    #     partial(data.augment.randreads, random_reads_fraction=0.15)
    # )
    return DataLoader(concatdst, batch_sampler=sampler, num_workers=threads, pin_memory=True)


def make_valloaders(intervals: BedTool, replicas: Tuple[ChIPseqReplicaMeta, ...], binsize: int,
                    batch_size: int, allsampling: Tuple[str, ...], threads: int) -> Dict[str, DataLoader]:
    allsampling = sorted(set(allsampling))
    hsh = f"{intervals.fn}{' '.join(sorted(str(r) for r in replicas))}{binsize}{batch_size}{allsampling}"
    path = os.path.join("/tmp/", f"cached.make_valloaders(hash={md5(hsh.encode()).hexdigest()}).pkl")

    threads = threads if threads >= 0 else cpu_count() + threads + 1
    if not os.path.exists(path):
        result = Parallel(n_jobs=max(1, threads), verbose=1, batch_size=1)(
            delayed(make_dataset)(intervals.fn, r, binsize, s) for r in replicas for s in allsampling
        )
        dsts, _ = zip(*result)
        with open(path, 'wb') as file:
            pickle.dump(dsts, file)

    with open(path, 'rb') as file:
        valdatasets = pickle.load(file)
        for dst in valdatasets:
            dst.reread_bigwig_on_query()

    valloaders = {
        dst: DataLoader(dst, sampler=torch.utils.data.SequentialSampler(dst), batch_size=batch_size,
                        drop_last=False, num_workers=threads, pin_memory=True) for dst in valdatasets
    }
    return valloaders


@dataclass
class _ReadsSeparationContext:
    start_positions: Dict[str, list]

    def add(self, chrom: str, position: int):
        # 0 based indices
        self.start_positions[chrom].append(position)

    def serialize(self, path):
        paths = {}
        for chrom, positions in self.start_positions.items():
            chrompath = os.path.join(path, chrom + ".npy")
            positions = np.asarray(positions, dtype=np.int32)
            np.save(chrompath, positions)
            paths[chrom] = chrompath
        return paths


def separate_reads(file):
    assert os.path.isfile(file)
    targetdir = os.path.join("/tmp/", f"cached.separate_reads({md5(file.encode()).hexdigest()})")

    forwarddir = os.path.join(targetdir, "forward")
    reversedir = os.path.join(targetdir, "reverse")
    chrominfo = os.path.join(targetdir, "chromosomes.pkl")

    if not os.path.exists(targetdir):
        samfile = pysam.AlignmentFile(file, "rb")

        chromosomes = {chr: samfile.get_reference_length(chr) for chr in samfile.references}

        forward, reverse = _ReadsSeparationContext({chr: [] for chr in chromosomes}), \
                           _ReadsSeparationContext({chr: [] for chr in chromosomes})

        for read in samfile.fetch():    # type: pysam.AlignedRead
            if read.is_reverse:
                reverse.add(read.reference_name, read.reference_end)
            else:
                forward.add(read.reference_name, read.reference_start)

        os.makedirs(forwarddir, exist_ok=False)
        forward.serialize(forwarddir)

        os.makedirs(reversedir, exist_ok=False)
        reverse.serialize(reversedir)

        with open(chrominfo, 'wb') as file:
            pickle.dump(chromosomes, file)

    assert os.path.isdir(forwarddir) and os.path.isdir(reversedir) and os.path.isfile(chrominfo)

    with open(chrominfo, 'rb') as file:
        chromosomes = pickle.load(file)

    result = {"forward": {}, "reverse": {}, "chromosomes": chromosomes}
    for key, folder in [("forward", forwarddir), ("reverse", reversedir)]:
        for file in os.listdir(folder):
            chrom = file.replace(".npy", "")
            result[key][chrom] = os.path.join(folder, file)
    assert set(result["forward"].keys()) == set(result["reverse"].keys()) == set(result["chromosomes"].keys())
    return result
