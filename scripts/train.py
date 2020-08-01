import warnings
from functools import partial

import torch
from ignite.contrib.handlers import ProgressBar
from ignite.engine import Engine, Events

from models import NaiveCNN
from utils import metrics, trainval, cached
from utils.losses import dice_loss

warnings.simplefilter("ignore")

ALL_SAMPLING = ("original", "1m", "5m", "10m", "25m")
BIN_SIZE = 25
BATCH_SIZE = 256
REGION_SIZE = 25 * 1024
THREADS = -2
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# 1. make
regions = cached.make_regions(REGION_SIZE)

train_experiments = (
    # H3K4me3
    "/data/encode/H3K4me3/ENCSR153NDQ",
    "/data/encode/H3K4me3/ENCSR206STN",
    "/data/encode/H3K4me3/ENCSR361FWQ",
)
trainloader = cached.make_trainloader(regions, cached.parse(train_experiments), BIN_SIZE, BATCH_SIZE, THREADS)

val_experiments = (
    # H3K4me3
    "/data/encode/H3K4me3/ENCSR652QNW",
    "/data/encode/H3K4me3/ENCSR930HLX",
    "/data/encode/H3K4me3/ENCSR956CTX",
)
valloaders = cached.make_valloaders(regions, cached.parse(val_experiments), BIN_SIZE, BATCH_SIZE * 2, ALL_SAMPLING,
                                    THREADS)

model = NaiveCNN().to(DEVICE)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=15, gamma=0.1)

def losspeaks2(y_pred, y):
    assert y_pred.shape == y.shape
    loss = -y * y_pred  # y_pred is a log(q) already

    # in the usual case I would simply increase the weight of the corrseponding positions by some constant factor
    # Is it worse it? To always increase by constant factor? I would say otherwise. It is better to sample with balance
    # - each batch has a positive and negative samples and to reweight within a batch to balance contribution of the
    # corresponding positions.

    # https://arxiv.org/pdf/1901.05555.pdf - balancing by the inverse root of class frequency. perhaps. Or by the effective number of samples

    with torch.no_grad():
        assert loss.ndim == 3
        perclass_items = y.sum(dim=[0, 2])
        mask = perclass_items < 1e-3
        balanced_perclass_items = loss.shape[0] * loss.shape[2] // (y.shape[1] - mask.sum())
        weights = balanced_perclass_items / perclass_items
        weights[mask] = 0
    return loss * weights.reshape(1, 2, 1)


kldivloss = torch.nn.KLDivLoss(reduction='none')
# losspeaks = lambda y_pred, y: (focal_loss(y_pred, y).mean(dim=2) + dice_loss(y_         pred, y)).mean(dim=1)
losspeaks = lambda y_pred, y: dice_loss(y_pred, y) / 2 + (kldivloss(y_pred, y) * 50).mean(dim=[1, 2])
# losspeaks = lambda y_pred, y: focal_loss(y_pred, y)
loss_enrichment = torch.nn.MSELoss(reduction='none')
trainer = Engine(partial(
    trainval.doiteration, loss_enrichment=loss_enrichment, loss_peaks=losspeaks, model=model, device=DEVICE,
    optimizer=optimizer
))
metrics.attach(trainer, "train")
ProgressBar(ncols=100).attach(trainer, metric_names='all')

trainval.attach_validation(
    trainer,
    partial(trainval.doiteration, loss_enrichment=loss_enrichment, loss_peaks=losspeaks, model=model, device=DEVICE),
    model, valloaders
)
trainval.attach_logger(trainer)

trainer.add_event_handler(Events.EPOCH_COMPLETED, lambda engine: engine.state.metrics.clear())
trainer.add_event_handler(Events.EPOCH_COMPLETED, lambda engine: scheduler.step())

trainer.run(trainloader, max_epochs=120)
torch.save(model, "/data/model-final.pth")
