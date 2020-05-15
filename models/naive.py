import torch
from torch import nn


class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, padding):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.seq = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=kernel_size, padding=padding),
            nn.BatchNorm1d(out_channels),
            nn.ReLU()
        )
        self.adjust = nn.Sequential(
                nn.Conv1d(in_channels, out_channels, kernel_size=1),
                nn.BatchNorm1d(out_channels),
                nn.ReLU()
        ) if in_channels != out_channels else lambda x: x
        self.bn = nn.BatchNorm1d(out_channels)

    def forward(self, x):
        x = self.seq(x) + self.adjust(x)
        x = self.bn(x).relu()
        return x


class NaiveCNN(nn.Module):
    def __init__(self):
        super().__init__()
        mul = 2
        self.encoder_block1 = nn.Sequential(
            ResidualBlock(1, int(8 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(8 * mul), int(16 * mul), kernel_size=3, padding=1)
        )
        self.pool1 = nn.MaxPool1d(kernel_size=4)
        self.encoder_block2 = nn.Sequential(
            ResidualBlock(int(16 * mul), int(16 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(16 * mul), int(32 * mul), kernel_size=3, padding=1)
        )
        self.pool2 = nn.MaxPool1d(kernel_size=4)
        self.encoder_block3 = nn.Sequential(
            ResidualBlock(int(32 * mul), int(32 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(32 * mul), int(64 * mul), kernel_size=3, padding=1)
        )
        self.pool3 = nn.MaxPool1d(kernel_size=4)
        self.encoder_block4 = nn.Sequential(
            ResidualBlock(int(64 * mul), int(64 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(64 * mul), int(128 * mul), kernel_size=3, padding=1)
        )

        self.upsample1 = nn.Upsample(scale_factor=4)
        self.decoder_block1 = nn.Sequential(
            nn.Conv1d(int(128 * mul) + int(64 * mul), int(128 * mul), kernel_size=1),
            nn.BatchNorm1d(int(128 * mul)),
            ResidualBlock(int(128 * mul), int(64 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(64 * mul), int(64 * mul), kernel_size=3, padding=1)
        )
        self.upsample2 = nn.Upsample(scale_factor=4)
        self.decoder_block2 = nn.Sequential(
            nn.Conv1d(int(64 * mul) + int(32 * mul), int(64 * mul), kernel_size=1),
            nn.BatchNorm1d(int(64 * mul)),
            ResidualBlock(int(64 * mul), int(32 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(32 * mul), int(32 * mul), kernel_size=3, padding=1)
        )
        self.upsample3 = nn.Upsample(scale_factor=4)
        self.decoder_block3 = nn.Sequential(
            nn.Conv1d(int(32 * mul) + int(16 * mul), int(32 * mul), kernel_size=1),
            nn.BatchNorm1d(int(32 * mul)),
            ResidualBlock(int(32 * mul), int(16 * mul), kernel_size=3, padding=1),
            ResidualBlock(int(16 * mul), int(8 * mul), kernel_size=3, padding=1),
            nn.Conv1d(int(8 * mul), 3, kernel_size=3, padding=1)
        )

    def forward(self, x, **kwargs):
        assert x.shape[-2] == 1, f"Model works only with raw reads ratio input, got tensor with shape {x.shape}"

        # down pass
        stage1 = self.encoder_block1(x)
        stage2 = self.encoder_block2(self.pool1(stage1))
        stage3 = self.encoder_block3(self.pool2(stage2))
        predict = self.encoder_block4(self.pool3(stage3))

        predict = self.upsample1(predict)
        predict = self.decoder_block1(torch.cat([predict, stage3], dim=1))
        predict = self.upsample2(predict)
        predict = self.decoder_block2(torch.cat([predict, stage2], dim=1))
        predict = self.upsample3(predict)
        predict = self.decoder_block3(torch.cat([predict, stage1], dim=1))

        predenrichment, predpeaks = predict[:, 0].unsqueeze(1), predict[:, 1:]
        predenrichment, predpeaks = predenrichment.relu(), predpeaks.log_softmax(dim=1)

        assert x.shape == predenrichment.shape
        return predenrichment, predpeaks
