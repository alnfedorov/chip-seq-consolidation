import warnings
import torch
import gc

from typing import Callable
from ignite.engine import Events, Engine


class EpochMetricsContainer:
    """
    Container for the network results that invokes attached metrics in the very end, when epoch is completed.
    Use output transform to turn engine.state.output -> (y_pred, y_true)
    Use add_metrics to attach new metrics, each metric is a callable(y_pred, y_true) -> float
    """
    def __init__(self, output_transform=lambda x: x):
        self._output_transform = output_transform
        self._predictions: [torch.Tensor] = []
        self._targets: [torch.Tensor] = []
        self._metrics: {str: Callable[[torch.Tensor, torch.Tensor], float]} = {}

    def add_metric(self, name: str, metric: Callable[[torch.Tensor, torch.Tensor], float]) -> "EpochMetricsContainer":
        self._metrics[name] = metric
        return self

    def reset(self):
        self._predictions = []
        self._targets = []
        gc.collect()

    def started(self, engine):
        self.reset()

    def update(self, output):
        y_pred, y_true = output

        if y_pred.ndimension() not in (1, 2):
            raise ValueError("Predictions should be of shape (batch_size, n_classes) or (batch_size, ).")

        if y_true.ndimension() not in (1, 2):
            raise ValueError("Targets should be of shape (batch_size, n_classes) or (batch_size, ).")

        if y_true.ndimension() == 2:
            if not torch.equal(y_true ** 2, y_true):
                raise ValueError("Targets should be binary (0 or 1).")

        if y_pred.ndimension() == 2 and y_pred.shape[1] == 1:
            y_pred = y_pred.squeeze(dim=-1)

        if y_true.ndimension() == 2 and y_true.shape[1] == 1:
            y_true = y_true.squeeze(dim=-1)

        y_pred = y_pred.to(torch.float32)
        y_true = y_true.to(torch.long)

        self._targets.append(y_true.cpu())
        self._predictions.append(y_pred.cpu())

    @torch.no_grad()
    def iteration_completed(self, engine):
        output = self._output_transform(engine.state.output)
        self.update(output)

    def compute(self):
        self._predictions = [torch.cat(self._predictions, dim=0)]
        self._targets = [torch.cat(self._targets, dim=0)]
        y_pred, y_true = self._predictions[0], self._targets[0]
        metrics = {}
        for name, metric in self._metrics.items():
            metrics[name] = metric(y_pred, y_true)
        return metrics

    def completed(self, engine):
        metrics = self.compute()
        for k, v in metrics:
            if torch.is_tensor(v) and len(v.shape) == 0:
                v = v.item()
            if k in engine.state.metrics:
                warnings.warn(f"Overwriting engine.state.metrics key {k} with value {v}")
            engine.state.metrics[k] = v

    def attach(self, engine: Engine):
        engine.add_event_handler(Events.EPOCH_COMPLETED, self.completed)
        if not engine.has_event_handler(self.started, Events.EPOCH_STARTED):
            engine.add_event_handler(Events.EPOCH_STARTED, self.started)
        if not engine.has_event_handler(self.iteration_completed, Events.ITERATION_COMPLETED):
            engine.add_event_handler(Events.ITERATION_COMPLETED, self.iteration_completed)
