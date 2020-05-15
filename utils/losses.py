import torch
from torch import nn


class DiceLoss(nn.Module):
    def __init__(self, pmode: str = "log-prob", smooth: float = 1e-3):
        super(DiceLoss, self).__init__()
        assert pmode in ("log-prob",)
        self.pmode = pmode
        self.smooth = smooth

    def forward(self, predict: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        '''
        :param predict: B,C,*dims: probabilities in format pmode
        :param target: B,C,*dims: one-hot encoded ground true labels
        :return: tensor: B - loss value
        '''
        assert predict.shape == target.shape
        B, C = predict.shape[:2]

        if self.pmode == "log-prob":
            predict = torch.exp(predict)
        else:
            raise NotImplementedError()

        predict, target = predict.view(B, C, -1), target.view(B, C, -1)
        intersection = (predict * target).sum(dim=-1)  # B,C

        predict_area, target_area = predict.sum(dim=-1), target.sum(dim=-1)  # B,C
        loss = 1 - (2.0 * intersection + self.smooth) / (predict_area + target_area + self.smooth)
        return loss.mean(dim=1)


class FocalLoss(nn.Module):
    def __init__(self, pmode: str = "prob", gamma: float = 0):
        super(FocalLoss, self).__init__()
        assert pmode in ("log-prob",)
        self.pmode = pmode
        self.gamma = gamma

    def forward(self, predict: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        '''
        :param predict: N,C,*dims: probabilities in format pmode
        :param target: N,C,*dims: one-hot encoded ground true labels (might be soft encoded)
        :return: tensor: N - loss value
        '''
        if self.pmode != "log-prob":
            raise NotImplementedError()

        B, C = predict.shape[:2]
        logprob = predict.view(B, C, -1)  # N,C,H*W
        logprob = logprob.transpose(1, 2)  # N,H*W,C
        logprob = logprob.contiguous().view(-1, C)  # N,H*W,C => N*H*W,C

        target = target.view(-1, 1).type(torch.int64)  # Indices must be long
        logprob = logprob.gather(dim=1, index=target) + 1e-6  # 1e-6 for numerical stability

        loss = -1 * torch.pow((1 - logprob), self.gamma) * torch.log(logprob)
        loss = loss.view(B, 1, *logprob.shape[-2:])
        return loss


def focal_loss(logits, targets, gamma=2):
    ce = -logits * targets
    pt = torch.exp(logits)
    # reweight by probability
    loss = torch.pow(torch.abs(targets - pt), gamma) * ce
    return loss


def dice_loss(logits: torch.Tensor, targets: torch.Tensor):
    assert logits.ndim == targets.ndim == 3 and logits.shape == targets.shape
    assert logits.shape[1] >= 2
    dim = [2]
    p = torch.exp(logits)
    loss = 1 - (2 * p * targets).sum(dim=dim) / (p.sum(dim=dim) + targets.sum(dim=dim) + 1e-3)
    return loss.mean(dim=1)
