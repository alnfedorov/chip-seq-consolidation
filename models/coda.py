import torch
from torch import nn


# It is not clear how to reuse current data producers and reproduce coda model in a simple way
# coda is a simple CNN, Conv1d(filters=6, length=51) -> Relu -> FC()
# Hence, coda was trained in a bin-by-bin bases, not seq-by-seq.
