import warnings
import torch

from data import dataset, parse
from data import hg19
from utils import metrics, trainval, cached
from data.pod import *
from models import NaiveCNN
from ignite.engine import Engine, Events
from functools import partial
from ignite.contrib.handlers import ProgressBar
from ignite.metrics import RunningAverage
from utils.losses import dice_loss, focal_loss

warnings.simplefilter("ignore")

ALL_SIMMODS = ("q0.25", "q0.5", "q0.75")
BIN_SIZE = 25
BATCH_SIZE = 256
REGION_SIZE = 25 * 1024
THREADS = 11
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# 1. make
regions = cached.make_regions(REGION_SIZE)

experiments = (
    # H3K4me3
    # "/data/encode/H3K4me3/ENCSR849YFO/",
    "/data/encode/H3K4me3/ENCSR956CTX",
)

replicas = cached.parse(experiments)
trainloader = cached.make_trainloader(regions, replicas, BIN_SIZE, BATCH_SIZE, ALL_SIMMODS, THREADS)
valloaders = cached.make_valloaders(regions, replicas, BIN_SIZE, BATCH_SIZE * 2, ALL_SIMMODS, THREADS)
model = NaiveCNN().to(DEVICE)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)


# What do I want to do?
# I want per-dataset metric averaged
# over metric + sampling, over metric, over sampling, over all
# How to understand/manage this problem? The problem is unique identification for the replicas

# def losspeaks(y_pred, y):
#     assert y_pred.shape == y.shape
#     loss = -y * y_pred   # y_pred is a log(q) already
#
#     # in the usual case I would simply increase the weight of the corrseponding positions by some constant factor
#     # Is it worse it? To always increase by constant factor? I would say otherwise. It is better to sample with balance
#     # - each batch has a positive and negative samples and to reweight within a batch to balance contribution of the
#     # corresponding positions.
#
#     # https://arxiv.org/pdf/1901.05555.pdf - balancing by the inverse root of class frequency. perhaps. Or by the effective number of samples
#
#     with torch.no_grad():
#         assert loss.ndim == 3
#         perclass_items = y.sum(dim=[0, 2])
#         mask = perclass_items < 1e-3
#         balanced_perclass_items = loss.shape[0] * loss.shape[2] // (y.shape[1] - mask.sum())
#         weights = balanced_perclass_items / perclass_items
#         weights[mask] = 0
#     return loss * weights.reshape(1, 2, 1)

kldivloss = torch.nn.KLDivLoss(reduction='none')
# losspeaks = lambda y_pred, y: (focal_loss(y_pred, y).mean(dim=2) + dice_loss(y_         pred, y)).mean(dim=1)
losspeaks = lambda y_pred, y: dice_loss(y_pred, y) / 2 + (kldivloss(y_pred, y) * 50).mean(dim=[1, 2])
# losspeaks = lambda y_pred, y: focal_loss(y_pred, y)
loss_enrichment = torch.nn.MSELoss(reduction='none')
trainer = Engine(partial(
    trainval.doiteration, loss_enrichment=loss_enrichment, loss_peaks=losspeaks, model=model, device=DEVICE, optimizer=optimizer
))
trainval.attach_validation(
    trainer, partial(trainval.doiteration, loss_enrichment=loss_enrichment, loss_peaks=losspeaks, model=model, device=DEVICE),
    model, valloaders
)
metrics.attach(trainer, "train")
trainval.attach_logger(trainer)
ProgressBar(ncols=100).attach(trainer, metric_names='all')
trainer.add_event_handler(Events.EPOCH_COMPLETED, lambda engine: engine.state.metrics.clear())
trainer.run(trainloader, max_epochs=100)



# # There is a problem! I need BALANCED big-wig files.
# # What is a model input?
# # What is a model output?  called-peaks + enrichment over control
#
# # How to normalize enrichment over control? Number of reads is variable - for sure.
# # What is about computing probability density functions for bckg and fg reads?
# # bckg = bckg / (total_bckg * readlen), fg = fg / (total_fg * readlen). (Perhaps it must be a bit more complicated)
# # Ok, it normalizes read count. In a sense that it is a valid probability density function now and doesn't depend on the number of reads.
# # What do I wan't to predict? Called peaks, for sure. What is else? Something that allows me to easily predict reads distribution
# # 5'--------------------------3'
# #       5'------------------------------------3'
# #               5'----------------3'
# #           5'--------------3'
# # What are this probabilities about? It is a probability for the given position to get sequenced during ChIP-seq(control ChIP-seq)
