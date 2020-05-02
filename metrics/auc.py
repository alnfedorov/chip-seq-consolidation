from sklearn.metrics import roc_auc_score, average_precision_score, auc, precision_recall_curve, f1_score


def auprc(y_pred, y_true):
    precision, recall, thresholds = precision_recall_curve(y_true.numpy(), y_pred.numpy())
    return auc(recall, precision)


def auroc(y_pred, y_true):
    return roc_auc_score(y_true.numpy(), y_pred.numpy())


def average_precision(y_pred, y_true):
    return average_precision_score(y_true.numpy(), y_pred.numpy())
