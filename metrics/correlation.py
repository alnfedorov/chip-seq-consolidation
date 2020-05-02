from scipy import stats


def pearsonr(y_true, y_pred):
    y_true, y_pred = y_true.numpy(), y_pred.numpy()
    return stats.pearsonr(y_true, y_pred)[0]


def spearmanr(y_true, y_pred):
    y_true, y_pred = y_true.numpy(), y_pred.numpy()
    return stats.spearmanr(y_true, y_pred).correlation
