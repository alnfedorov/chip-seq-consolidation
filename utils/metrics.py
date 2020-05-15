import torch
from ignite.engine import Engine
from ignite.metrics import MeanSquaredError, confusion_matrix, Loss, RunningAverage

from metrics import EpochMetricsContainer, average_precision, auprc, auroc
from metrics import NonZeroMeanSquaredError


def _attach_losses(engine: Engine, prefix: str = "", running_average: bool = False):
    loss_peaks = Loss(
        output_transform=lambda x: (x["loss-peaks"], x["loss-peaks"]), loss_fn=lambda *x: x[0].mean()
    )
    loss_peaks.attach(engine, prefix + 'loss-peaks')
    loss_enrichment = Loss(
        output_transform=lambda x: (x["loss-enrichment"], x["loss-enrichment"]), loss_fn=lambda *x: x[0].mean()
    )
    loss_enrichment.attach(engine, prefix + 'loss-enrichment')

    if running_average:
        RunningAverage(loss_peaks).attach(engine, prefix + 'ra-loss-peaks')
        RunningAverage(loss_enrichment).attach(engine, prefix + 'ra-loss-enrichment')


def _attach_cm_related(engine: Engine, prefix: str = ""):
    transform = lambda x: (x["predpeaks"], x["peaks"].argmax(dim=1))
    cmatrix = confusion_matrix.ConfusionMatrix(num_classes=2, output_transform=transform)
    confusion_matrix.IoU(cmatrix, ignore_index=0)[0].attach(engine, prefix + "iou")

    precision = confusion_matrix.cmPrecision(cmatrix, average=False)[1]
    precision.attach(engine, prefix + "prec")

    recall = confusion_matrix.cmRecall(cmatrix, average=False)[1]
    recall.attach(engine, prefix + "recall")

    f1 = precision * recall * 2 / (precision + recall + 1e-20)
    f1.attach(engine, prefix + "f1")


def _attach_peaks_related(engine: Engine, prefix: str = ""):
    transform = lambda x: (x["refenrichment"].flatten(), x["predenrichment"].flatten())
    MeanSquaredError(transform).attach(engine, prefix + "mse")
    NonZeroMeanSquaredError(transform).attach(engine, prefix + "non_zero_mse")


def _attach_aucs(engine: Engine, prefix: str = ""):
    # Note: tracking AUROC and AUPRC works but requires enormous amount of memory.
    def transform(output):
        predpeaks, peaks = output["predpeaks"], output["peaks"]
        peaks = peaks.argmax(dim=1, keepdim=True)
        predpeaks = predpeaks.gather(dim=1, index=peaks)
        assert peaks.shape == predpeaks.shape
        return predpeaks.flatten(), peaks.flatten()
    peaks_container = EpochMetricsContainer(output_transform=transform)
    for name, metric in (("auroc", auroc), ("auprc", auprc), ("average_prec", average_precision)):
        peaks_container.add_metric(name, metric)
    peaks_container.attach(engine)


def attach(engine: Engine, mode: str):
    assert mode in ("train", "val")
    if mode == "train":
        _attach_losses(engine, running_average=True)
        _attach_cm_related(engine)
        _attach_peaks_related(engine)
    else:
        _attach_losses(engine, running_average=False)
        _attach_cm_related(engine)
        _attach_peaks_related(engine)
        # _attach_aucs(engine)
