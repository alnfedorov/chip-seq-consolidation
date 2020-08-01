import os
from collections import defaultdict

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
from data.pod import ChIPseqReplicaMeta, IntervalMeta, SeparatedReadsMeta


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


# def make_dataset(intervals: str, replica: ChIPseqReplicaMeta, binsize: int,
#                  sampling: str) -> Tuple[ChIPseqReplicaMeta, IntervalMeta]:
#     hsh = f"{intervals}{replica}{binsize}{sampling}"
#     path = os.path.join("/tmp/", f"cached.make_dataset(hash={md5(hsh.encode()).hexdigest()}).pkl")
#
#     if not os.path.exists(path):
#         dst = dataset.ChIPseqReplicaTrainDataset(BedTool(intervals), replica, binsize, sampling)
#
#         annotation = ("3'utr", "5'utr", "exons", "introns", "upstream1k", "downstream1k", "intergenic")
#         annotation = {k: hg19.annotation.REGIONS[k].fn for k in annotation}
#         meta, todrop = dataset.filter.by_meta(dst, annotation, True)
#
#         dst.remove(todrop)
#         meta = np.delete(meta, todrop).tolist()
#
#         with open(path, 'wb') as file:
#             pickle.dump((dst, meta), file)
#
#     with open(path, 'rb') as file:
#         dst, meta = pickle.load(file)
#         dst.reread_bigwig_on_query()
#     return dst, meta


def make_trainloader(intervals: BedTool, replicas: Tuple[ChIPseqReplicaMeta, ...], binsize: int,
                     batch_size: int, threads: int) -> DataLoader:
    hsh = intervals.fn + ''.join(sorted(str(r) for r in replicas)) + str(binsize)
    path = os.path.join("/tmp/", f"cached.make_trainloader(hash={md5(hsh.encode()).hexdigest()}).pkl")

    threads = threads if threads >= 0 else cpu_count() + threads + 1
    if not os.path.exists(path):
        # 1. Turn bam files into separate reads coordinates
        allbam_groups = set(tuple(r.treatment["original"] for r in replicas) +
                            tuple(r.control["original"] for r in replicas))
        separated_reads = Parallel(n_jobs=6, verbose=1, batch_size=1)(
            delayed(separate_reads)([f.path for f in files]) for files in allbam_groups
        )
        bamreads = {bgroup: reads for bgroup, reads in zip(allbam_groups, separated_reads)}

        # 2. make dataset for each replica(original data, sampled version will be generated on the fly)
        datasets = []
        for replica in replicas:
            treatment, tchrominfo = bamreads[replica.treatment["original"]]
            control, cchrominfo = bamreads[replica.control["original"]]
            assert treatment.numreads == sum(t.numreads for t in replica.treatment["original"]) and \
                   control.numreads == sum(c.numreads for c in replica.control["original"])
            assert set(tchrominfo.keys()) == set(cchrominfo.keys()) and \
                   all(tchrominfo[k] == cchrominfo[k] for k in tchrominfo)

            fragment_size = tuple(x.estimated_fragment_size for x in replica.treatment["original"])
            datasets.append(dataset.ChIPseqReplicaTrainDataset(
                replica.target, intervals, replica.peaks, treatment, control, tchrominfo, fragment_size,
                binsize, hg19.annotation.EFFECTIVE_GENOME_SIZE
            ))

        # 3. Make single dataset for simplicity
        concatdst = dataset.ConcatChIPseqDataset(datasets)
        sampler = torch.utils.data.BatchSampler(
            torch.utils.data.RandomSampler(concatdst), batch_size=batch_size, drop_last=False
        )
        # sampler = data.sampling.BalancedBatchSampler(metas, drop_last=False)

        with open(path, 'wb') as file:
            pickle.dump((concatdst, sampler), file)

    with open(path, 'rb') as file:
        concatdst, sampler = pickle.load(file)
        # for d in concatdst.datasets:  # type: dataset.ChIPseqReplicaTrainDataset
        #     d.reread_bigwig_on_query()
    return DataLoader(concatdst, batch_sampler=sampler, num_workers=threads, pin_memory=True)


def make_valloaders(intervals: BedTool, replicas: Tuple[ChIPseqReplicaMeta, ...], binsize: int,
                    batch_size: int, allsampling: Tuple[str, ...], threads: int) -> Dict[str, DataLoader]:
    allsampling = sorted(set(allsampling))
    hsh = f"{intervals.fn}{' '.join(sorted(str(r) for r in replicas))}{binsize}{batch_size}{allsampling}"
    path = os.path.join("/tmp/", f"cached.make_valloaders(hash={md5(hsh.encode()).hexdigest()}).pkl")

    threads = threads if threads >= 0 else cpu_count() + threads + 1
    if not os.path.exists(path):
        # 1. Turn bam files into separate reads coordinates
        allbam_groups = set(chain(*[r.treatment.values() for r in replicas], *[r.control.values() for r in replicas]))
        separated_reads = Parallel(n_jobs=6, verbose=1, batch_size=1)(
            delayed(separate_reads)([f.path for f in files]) for files in allbam_groups
        )
        bamreads = {bgroup: reads for bgroup, reads in zip(allbam_groups, separated_reads)}

        # 2. make val dataset for each replica
        datasets = []
        for replica in replicas:
            reftreatment, reftchrominfo = bamreads[replica.treatment["original"]]
            refcontrol, refcchrominfo = bamreads[replica.control["original"]]
            assert reftreatment.numreads == sum(t.numreads for t in replica.treatment["original"]) and \
                   refcontrol.numreads == sum(c.numreads for c in replica.control["original"])
            assert set(reftchrominfo.keys()) == set(refcchrominfo.keys()) and \
                   all(reftchrominfo[k] == refcchrominfo[k] for k in reftchrominfo)

            reffragment_size = int(np.mean([x.estimated_fragment_size for x in replica.treatment["original"]]))
            for sampling in allsampling:
                sampledtreatment, sampledtchrominfo = bamreads[replica.treatment[sampling]]
                sampledcontrol, sampledcchrominfo = bamreads[replica.control[sampling]]
                assert sampledtreatment.numreads == sum(t.numreads for t in replica.treatment[sampling]) and \
                       sampledcontrol.numreads == sum(c.numreads for c in replica.control[sampling])
                assert set(sampledtchrominfo.keys()) == set(sampledcchrominfo.keys()) and \
                       all(sampledtchrominfo[k] == sampledcchrominfo[k] for k in reftchrominfo)

                assert set(sampledtchrominfo.keys()) == set(reftchrominfo.keys()) and \
                       all(sampledtchrominfo[k] == reftchrominfo[k] for k in reftchrominfo)

                sampledfragment_size = int(np.mean([x.estimated_fragment_size for x in replica.treatment[sampling]]))

                taccessions = replica.accesions.treatment
                uid = ",".join(taccessions) if len(taccessions) > 1 else taccessions[0]
                datasets.append(dataset.ChIPseqReplicaValDataset(
                    uid, replica.accesions.experiment,
                    replica.target, sampling, intervals, replica.peaks,
                    reftreatment, refcontrol, reffragment_size,
                    sampledtreatment, sampledcontrol, sampledfragment_size,
                    reftchrominfo, binsize, hg19.annotation.EFFECTIVE_GENOME_SIZE
                ))

        with open(path, 'wb') as file:
            pickle.dump(datasets, file)

    with open(path, 'rb') as file:
        valdatasets = pickle.load(file)

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
            positions = np.sort(positions)
            np.save(chrompath, positions)
            paths[chrom] = chrompath
        return paths


def separate_reads(files: Tuple[str, ...]) -> Tuple[SeparatedReadsMeta, Dict[str, int]]:
    assert all(os.path.isfile(f) for f in files)
    files = sorted(files)
    targetdir = os.path.join("/tmp/", f"cached.separate_reads({md5(''.join(files).encode()).hexdigest()})")

    forwarddir = os.path.join(targetdir, "forward")
    reversedir = os.path.join(targetdir, "reverse")
    chrominfo = os.path.join(targetdir, "chromosomes.pkl")
    meta = os.path.join(targetdir, "meta.pkl")

    if not os.path.exists(targetdir):
        chromosomes = {}
        numreads = 0
        readlen = None
        forward, reverse = _ReadsSeparationContext(defaultdict(list)), \
                           _ReadsSeparationContext(defaultdict(list))

        # concat reads from all files
        for file in files:
            samfile = pysam.AlignmentFile(file, "rb")
            for chr in samfile.references:
                length = samfile.get_reference_length(chr)
                if chr in chromosomes:
                    assert chromosomes[chr] == length
                else:
                    chromosomes[chr] = length

            for read in samfile.fetch():    # type: pysam.AlignedRead
                numreads += 1
                if readlen is None:
                    readlen = read.query_length
                assert readlen == read.query_length

                if read.is_reverse:
                    reverse.add(read.reference_name, read.reference_end)
                else:
                    forward.add(read.reference_name, read.reference_start)

        assert readlen is not None and numreads > 0

        os.makedirs(forwarddir, exist_ok=False)
        forward.serialize(forwarddir)

        os.makedirs(reversedir, exist_ok=False)
        reverse.serialize(reversedir)

        with open(chrominfo, 'wb') as file:
            pickle.dump(chromosomes, file)

        with open(meta, "wb") as file:
            pickle.dump((numreads, readlen), file)

    assert os.path.isdir(forwarddir) and os.path.isdir(reversedir) and \
           os.path.isfile(chrominfo) and os.path.isfile(meta)

    with open(chrominfo, 'rb') as file:
        chromosomes = pickle.load(file)

    with open(meta, 'rb') as file:
        numreads, readlen = pickle.load(file)

    result = {"forward": {}, "reverse": {}, "chromosomes": chromosomes}
    for key, folder in [("forward", forwarddir), ("reverse", reversedir)]:
        for file in os.listdir(folder):
            chrom = file.replace(".npy", "")
            result[key][chrom] = os.path.join(folder, file)

    assert set(result["forward"].keys()) == set(result["reverse"].keys()) == set(result["chromosomes"].keys())
    return SeparatedReadsMeta(result["forward"], result["reverse"], readlen, numreads), result["chromosomes"]
