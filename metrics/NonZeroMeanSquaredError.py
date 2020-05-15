from ignite.metrics import MeanSquaredError
from ignite.metrics.metric import reinit__is_reduced


class NonZeroMeanSquaredError(MeanSquaredError):
    @reinit__is_reduced
    def update(self, output):
        y_pred, y = output
        mask = (y_pred != 0) | (y != 0)
        return super().update((y_pred[mask], y[mask]))
