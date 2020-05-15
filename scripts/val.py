# I have a set of files -> big wig files, treatment + control
# I wan't to turn them into the model input. How to do so?
# And, I want to predict the whole genome.
# How that would looks like?
# Working with a set of files -> for sure.
# I get experiment json files, split genome into overlapping windows and lets go

# 1. Parse json experiments
import warnings
import torch

from data import dataset, parse
from data import hg19
from data.pod import *
from models import NaiveCNN
from data.sampling import BalancedBatchSampler
from utils import trainval
from pybedtools import BedTool
from ignite.engine import Engine, Events
from functools import partial
from ignite.contrib.handlers import ProgressBar
from ignite.metrics import RunningAverage

warnings.simplefilter("ignore")

EXPERIMENT = "/data/encode/H3K4me3/ENCSR849YFO/"
REGION_SIZE = 25000


# 1. parse json experiment meta
replicas = parse.fromjson(EXPERIMENT)

# 2. split genome into regions
regions = BedTool().window_maker(w=REGION_SIZE, genome="hg19")

# 3. for each replica and for each mode I need to save original ratio, reference ratio and predict result\


# Right now each replica is a separate sample in the training process
# How to fix it? Is it worth it? No, it is not. I might want to combine replicas in a smart way. For sure. But not now.
# Hence, the only option is to do the following -> For each replica:
# 1. Predict consensus and enrichment
# 2. Save prediction and reference files





REGION_SIZE = 25000
ALL_SIMMODS = ("q0.25", "q0.5", "q0.75")
BIN_SIZE = 25
BATCH_SIZE = 128
THREADS = -1
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")





regions = trainval.cached.make_regions(REGION_SIZE)

experiments = (
    # H3K4me3
    "/data/encode/H3K4me3/ENCSR849YFO/",
    # "/data/encode/H3K4me3/ENCSR956CTX"
)

replicas = trainval.cached.parse(experiments)
trainloader = trainval.cached.make_trainloader(regions, replicas, BIN_SIZE, BATCH_SIZE, ALL_SIMMODS, THREADS)
valloaders = trainval.cached.make_valloaders(regions, replicas, BIN_SIZE, BATCH_SIZE * 2, ALL_SIMMODS, THREADS)
model = NaiveCNN().to(DEVICE)

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

