import warnings
import torch

from data import dataset, parse
from data import hg19
from data.pod import *
from models import NaiveCNN
from data.sampling import BalancedBatchSampler
from utils import trainval
from ignite.engine import Engine
from functools import partial
from ignite.contrib.handlers import ProgressBar
from ignite.metrics import RunningAverage

warnings.simplefilter("ignore")

REGION_SIZE = 25000
ALL_SIMMODS = ("q0.25", "q0.5", "q0.75")
BIN_SIZE = 25
BATCH_SIZE = 128

# 1. make
regions = trainval.cached.make_regions(REGION_SIZE)

experiments = (
    # H3K4me3
    "/data/encode/H3K4me3/ENCSR849YFO/",
    # "/data/encode/H3K4me3/ENCSR956CTX"
)

replicas = trainval.cached.parse(experiments)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
loader = trainval.cached.make_trainloader(regions, replicas, BIN_SIZE, BATCH_SIZE, ALL_SIMMODS, 0)
model = NaiveCNN().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
losspeaks = torch.nn.BCELoss()
lossratio = torch.nn.MSELoss()

def run_iteration(engine: Engine, batch, *,
                  model, device, optimizer=None):
    ratio, refratio, peaks = batch["ratio"].to(device), batch["refratio"].to(device), batch["peaks"].to(device)
    ratio, refratio, peaks = ratio.unsqueeze(1), refratio.unsqueeze(1), peaks.unsqueeze(1)
    predratio, predpeaks = model(ratio)

    lratio, lpeaks = lossratio(predratio, refratio).mean(), losspeaks(predpeaks, peaks).mean()
    if optimizer:
        loss = (lratio + lpeaks)
        loss.backward()
        optimizer.step()
    return {"loss-peaks": lpeaks.item(), "loss-ratio": lratio.item()}


trainer = Engine(partial(run_iteration, model=model, device=device, optimizer=optimizer))
RunningAverage(output_transform=lambda x: x["loss-peaks"]).attach(trainer, 'avg-loss-peaks')
RunningAverage(output_transform=lambda x: x["loss-ratio"]).attach(trainer, 'avg-loss-ratio')
ProgressBar(position=0, persist=True).attach(trainer, 'all')
trainer.run(loader, max_epochs=100)
#
# # # There is a problem! I need BALANCED big-wig files.
# # # What is a model input?
# # # What is a model output?  called-peaks + enrichment over control
# #
# # # How to normalize enrichment over control? Number of reads is variable - for sure.
# # # What is about computing probability density functions for bckg and fg reads?
# # # bckg = bckg / (total_bckg * readlen), fg = fg / (total_fg * readlen). (Perhaps it must be a bit more complicated)
# # # Ok, it normalizes read count. In a sense that it is a valid probability density function now and doesn't depend on the number of reads.
# # # What do I wan't to predict? Called peaks, for sure. What is else? Something that allows me to easily predict reads distribution
# # # 5'--------------------------3'
# # #       5'------------------------------------3'
# # #               5'----------------3'
# # #           5'--------------3'
# # # What are this probabilities about? It is a probability for the given position to get sequenced during ChIP-seq(control ChIP-seq)
