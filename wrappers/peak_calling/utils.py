from pybedtools import BedTool
from math import log10


# Used to replace 0 pvalue or 0 fdr with max(500, max(-log10(values))
def _zero_significance_filler(lines, ind, nlog_infinity: int = 500):
    key = lambda x: -log10(float(x[ind])) if float(x[ind]) != 0 else -1
    filler = max(lines, key=key)
    return max(nlog_infinity, -log10(float(filler[ind])))
