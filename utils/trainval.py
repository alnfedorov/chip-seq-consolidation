from itertools import chain
from typing import Dict

import numpy as np
import torch
from ignite.engine import Engine, Events
from torch.utils.data import DataLoader

from data.dataset import ChIPseqReplicaValDataset
from . import metrics


def doiteration(engine: Engine, batch, *, loss_enrichment, loss_peaks, device, model: torch.nn.Module, optimizer=None):
    fetch = lambda key: batch.pop(key).to(device).unsqueeze(1)
    enrichment, refenrichment, peaks = fetch("enrichment"), fetch("refenrichment"), batch.pop("peaks").to(device)

    if optimizer:
        optimizer.zero_grad()

    predenrichment, predpeaks = model(enrichment)
    assert enrichment.shape == refenrichment.shape == predenrichment.shape
    assert predpeaks.shape == (enrichment.shape[0], 2, enrichment.shape[-1])
    batchsize = enrichment.shape[0]
    enrichment, refenrichment, predenrichment = enrichment.reshape(batchsize, -1), \
                                                refenrichment.reshape(batchsize, -1), \
                                                predenrichment.reshape(batchsize, -1)

    lenrich, lpeaks = loss_enrichment(predenrichment, refenrichment), loss_peaks(predpeaks, peaks)
    if optimizer:
        (lenrich.mean() + lpeaks.mean()).backward()
        optimizer.step()

    return {
        "loss-peaks": lpeaks, "loss-enrichment": lenrich,
        "refenrichment": refenrichment, "enrichment": enrichment, "predenrichment": predenrichment,
        "peaks": peaks, "predpeaks": predpeaks
    }


def attach_validation(engine: Engine, doiter, model: torch.nn.Module, dataloaders: Dict[str, DataLoader]):
    def doval(engine: Engine, valengines: Dict[ChIPseqReplicaValDataset, Engine],
              valloaders: Dict[ChIPseqReplicaValDataset, DataLoader]):
        if engine.state.epoch % 10 != 0 or engine.state.epoch == 0:
            return

        torch.cuda.empty_cache()
        model.eval()

        with torch.no_grad():
            # 1. compute metrics for each dataset separately
            metrics_per_dst = {
                dst: valengines[dst].run(valloaders[dst], max_epochs=1, seed=123).metrics for dst in valengines
            }
            metricnames = set(chain(*[metrics.keys() for metrics in metrics_per_dst.values()]))

            # 2. report per replica metrics with key: val-target-sampling-experiment_accession(treatment accession)
            result = {}
            for dst, dstmetrics in metrics_per_dst.items():
                base_key = f"val-{dst.target}-{dst.sampling}-{dst.experiment_accession}({dst.uid})"
                assert all(f"{base_key}-{k}" not in result for k in dstmetrics)
                result.update({f"{base_key}-{k}": v for k, v in dstmetrics.items()})

            # 3. report metrics averaged over target
            utargets = set(x.target for x in metrics_per_dst.keys())
            for target in utargets:
                # fetch all relevant datasets
                dstmetrics = [dstmetrics for dst, dstmetrics in metrics_per_dst.items() if dst.target == target]
                for mname in metricnames:
                    values = [metrics[mname] for metrics in dstmetrics]
                    key = f"val-avg-{target}-{mname}"
                    assert key not in result
                    result[key] = np.mean(values)

            # 4. report metrics averaged over sampling method
            usamplings = set(x.sampling for x in metrics_per_dst.keys())
            for sampling in usamplings:
                # fetch all relevant datasets
                dstmetrics = [dstmetrics for dst, dstmetrics in metrics_per_dst.items() if dst.sampling == sampling]
                for mname in metricnames:
                    values = [metrics[mname] for metrics in dstmetrics]
                    key = f"val-avg-{sampling}-{mname}"
                    assert key not in result
                    result[key] = np.mean(values)

            # 5. report all combinations sampling/target
            usamtar = set((x.sampling, x.target) for x in metrics_per_dst.keys())
            for sampling, target in usamtar:
                # fetch all relevant datasets
                dstmetrics = [dstmetrics for dst, dstmetrics in metrics_per_dst.items() if
                              dst.sampling == sampling and dst.target == target]
                for mname in metricnames:
                    values = [metrics[mname] for metrics in dstmetrics]
                    key = f"val-avg-{sampling}/{target}-{mname}"
                    assert key not in result
                    result[key] = np.mean(values)

            assert all(k not in engine.state.metrics for k in result)
            engine.state.metrics.update(result)

        model.train()
        torch.cuda.empty_cache()

    valengines = {}
    for dst in dataloaders:
        valengine = Engine(doiter)
        metrics.attach(valengine, 'val')
        valengines[dst] = valengine

    engine.add_event_handler(Events.EPOCH_COMPLETED, doval, valengines, dataloaders)
    return engine


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


def attach_logger(engine: Engine):
    def dolog(engine: Engine):
        print("Metrics: ")
        for key, value in engine.state.metrics.items():
            if isinstance(value, (float, int)) or value.ndim == 0:
                value = float(value)
            print(f"{key}: {value}", end=' \n')

    engine.add_event_handler(Events.EPOCH_COMPLETED, dolog)
    return engine
