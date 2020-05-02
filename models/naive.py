from torch import nn


class NaiveCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Conv1d(1, 64, kernel_size=7, padding=3), nn.BatchNorm1d(64), nn.ReLU(),
            nn.Conv1d(64, 128, kernel_size=7, padding=3), nn.BatchNorm1d(128), nn.ReLU(),

            nn.Conv1d(128, 256, kernel_size=7, padding=3), nn.BatchNorm1d(256), nn.ReLU(),
            nn.Conv1d(256, 128, kernel_size=7, padding=3), nn.BatchNorm1d(128), nn.ReLU(),

            nn.Conv1d(128, 64, kernel_size=7, padding=3), nn.BatchNorm1d(64), nn.ReLU(),
            nn.Conv1d(64, 2, kernel_size=7, padding=3), nn.Sigmoid()
        )

    def forward(self, x, **kwargs):
        assert x.shape[-2] == 1, f"Model works only with raw reads ratio input, got tensor with shape {x.shape}"
        predict = self.layers(x)
        predratio, predpeaks = predict[:, 0].unsqueeze(1), predict[:, 1].unsqueeze(1)
        assert x.shape == predratio.shape == predpeaks.shape
        return predratio, predpeaks
