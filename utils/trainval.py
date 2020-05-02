import os
import data
import pickle
import torch
from functools import partial
from hashlib import md5
from data import hg19, dataset
from pybedtools import BedTool, Interval
from ignite.engine import Engine
from typing import Tuple
from joblib import Parallel, delayed
from ignite.metrics import MeanSquaredError
from itertools import chain
from data.pod import ReplicaMeta
from torch.utils.data import DataLoader
from multiprocessing import cpu_count
from data.dataset import ChIPseq
from data.sampling import BalancedBatchSampler
from metrics import NonZeroMeanSquaredError, EpochMetricsContainer, auroc, auprc, average_precision


def run_iteration(engine: Engine, batch, *, lossratio, losspeaks, model, device, optimizer=None):
    fetch = lambda key: batch.pop(key).to(device).unsqueeze(1)
    ratio, refratio, peaks = fetch("ratio"), fetch("refratio"), fetch("peaks")
    predratio, predpeaks = model(ratio)
    assert ratio.shape == refratio.shape == peaks.shape == predratio.shape == predpeaks.shape

    lratio, lpeaks = lossratio(predratio, refratio).mean(), losspeaks(predpeaks, peaks).mean()
    if optimizer:
        loss = (lratio + lpeaks)
        loss.backward()
        optimizer.step()
    batchsize = ratio.shape[0]
    return {
        "loss-peaks": lpeaks.item(), "loss-ratio": lratio.item(),
        "refratio": refratio.reshape(batchsize, -1),
        "ratio": ratio.reshape(batchsize, -1),
        "predratio": predratio.reshape(batchsize, -1),
        "peaks": peaks.reshape(batchsize, -1),
        "predpeaks": predpeaks.reshape(batchsize, -1)
    }


    # partial(run_validation, val_engine=evaluator, val_dataloader=val_dataloader, model=MODEL),
    # trainer.add_event_handler(Events.EPOCH_COMPLETED, h)
# What do I want to do during validation:
# 1. For each dataloader: fetch new


# def run_validation(engine: Engine, val_engine: Engine, val_dataloader: DataLoader, model):
#     torch.cuda.empty_cache()
#     model.eval()
#
#     with torch.no_grad():
#         metrics = val_engine.run(val_dataloader, max_epochs=1).metrics
#
#     engine.state.metrics.update(metrics)
#     model.train()
#     torch.cuda.empty_cache()


# def log_metrics(engine: Engine, cls_to_lbl):
#     print("Metrics: ")
#     for key, value in engine.state.metrics.items():
#         if isinstance(value, (float, int)) or value.ndim == 0:
#             value = float(value)
#             wandb.log({key: value}, commit=False)
#         else:
#             assert value.ndim == 2 and value.shape[0] == value.shape[1] and "Confusion" in key
#             value = value.cpu().numpy()
#             n_classes = value.shape[0]
#             wandb_m = wandb.Table(columns=["cls\\cls"]+[cls_to_lbl[lbl] for lbl in range(n_classes)])
#             for cls, row in enumerate(value):
#                 row = [cls_to_lbl[cls]] + [str(round(v, 5)) for v in row]
#                 wandb_m.add_data(*row)
#             wandb.log({"confusion_matrix": wandb_m}, commit=False)
#         print(f"{key}: {value}", end=' \n')
#     wandb.log(commit=True)



# def focal_loss_multiclass(logits, targets, alpha=1, gamma=2):
#     ce = F.binary_cross_entropy_with_logits(logits, targets, reduce=False)
#     ce = ce.reshape(-1)
#     pt = torch.exp(-ce)
#     # reweight by probability
#     loss = alpha * torch.pow((1 - pt), gamma) * ce


# def focal_loss(logits, targets, alpha=1, gamma=2):
#     ce = F.cross_entropy(logits, targets, reduce=False)
#     ce = ce.reshape(-1)
#     pt = torch.exp(-ce)
#     # reweight by probability
#     loss = alpha * torch.pow((1 - pt), gamma) * ce
#     return loss


def attach_metrics(engine: Engine):
    # Peaks related metrics
    peaks_container = EpochMetricsContainer(output_transform=lambda x: (x["predpeaks"], x["peaks"]))
    for name, metric in (("AUROC", auroc), ("AUPRC", auprc), ("AveragePrecision", average_precision)):
        peaks_container.add_metric(name, metric)
    peaks_container.attach(engine)

    # Reads ratio related metrics
    transform = lambda x: (x["predratio"], x["refratio"])
    MeanSquaredError(transform).attach(engine, "MSE")
    NonZeroMeanSquaredError(transform).attach(engine, "NonZeroMSE")
    return engine


class cached:
    @staticmethod
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

    @staticmethod
    def parse(experiments: Tuple[str, ...]):
        assert all(os.path.isdir(s) for s in experiments)

        hsh = md5('-'.join(sorted(experiments)).encode()).hexdigest()
        path = os.path.join("/tmp/", f"cached.parse(hash={hsh}).pkl")

        if not os.path.exists(path):
            replicas = Parallel(n_jobs=-1)(
                delayed(data.parse.fromjson)(root, hg19.annotation.BLACKLISTED.fn) for root in experiments
            )
            replicas = list(chain(*replicas))
            with open(path, 'wb') as file:
                pickle.dump(replicas, file)

        with open(path, 'rb') as file:
            replicas = pickle.load(file)
        return replicas

    @staticmethod
    def make_dataset_and_sampler(regions: BedTool, replicas: Tuple[ReplicaMeta, ...], binsize: int,
                                 batch_size: int, simmods: Tuple[str, ...]) -> Tuple[ChIPseq, BalancedBatchSampler]:
        hsh = regions.fn + ''.join(sorted(str(r) for r in replicas)) + str(binsize) + str(simmods)
        hsh = md5(hsh.encode()).hexdigest()
        path = os.path.join("/tmp/", f"cached.make_dataset_sampler(hash={hsh}).pkl")

        if not os.path.exists(path):
            dst = dataset.ChIPseq(regions, replicas, binsize)
            annotation = ("3'utr", "5'utr", "exons", "introns", "upstream1k", "downstream1k", "intergenic")
            annotation = {k: hg19.annotation.REGIONS[k] for k in annotation}
            metas = dataset.filter.by_meta(dst, annotation)
            sampler = data.sampling.BalancedBatchSampler(metas, drop_last=False)
            with open(path, 'wb') as file:
                pickle.dump((dst, sampler), file)

        with open(path, 'rb') as file:
            dst, sampler = pickle.load(file)
            dst.reread_bigwig_on_fetch()
        return dst, sampler

    @staticmethod
    def make_trainloader(regions: BedTool, replicas: Tuple[ReplicaMeta, ...], binsize: int,
                         batch_size: int, simmods: Tuple[str, ...], threads: int) -> DataLoader:
        threads = threads if threads >= 0 else cpu_count() + threads + 1
        dst, sampler = cached.make_dataset_and_sampler(regions, replicas, binsize, batch_size, simmods)
        dst.add_augmentation(
            partial(data.augment.shift, limits=(-0.5, 0.5)),
            partial(data.augment.randreads, random_reads_fraction=0.15)
        )
        return DataLoader(dst, batch_sampler=sampler, num_workers=threads, pin_memory=True)

    @staticmethod
    def make_valloaders(regions: BedTool, replicas: Tuple[ReplicaMeta, ...],
                        batch_size: int, binsize: int, simmods: Tuple[str, ...], threads: int) -> DataLoader:
        pass
