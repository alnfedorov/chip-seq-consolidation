import os
from pybedtools import BedTool
from pathlib import Path


_DEFAULT_CACHE_FOLDER = "/tmp/preprocessed_bed/"


def _filterqvalue(bed: BedTool, qvalue_thr: float = 10):
    if "narrow" in bed.fn:
        fields, qvalue_ind = 10, 8
    elif "broad" in bed.fn:
        fields, qvalue_ind = 9, 8
    else:
        assert "gapped" in bed.fn
        fields, qvalue_ind = 15, 14

    filtered = []
    for interval in bed:
        assert len(interval.fields) == fields
        if float(interval.fields[qvalue_ind]) >= qvalue_thr:
            filtered.append(interval)
    bed = BedTool(filtered).sort().merge().sort()
    return bed


def _load_file_cached(fn: str, preprocessing=lambda x: x.sort().merge().sort(),
                      cache_folder=_DEFAULT_CACHE_FOLDER) -> BedTool:
    path = Path(cache_folder + fn)
    # print(path.parent.as_posix())
    if path.exists():
        return BedTool(path.as_posix())
    bed = BedTool(fn)
    bed = preprocessing(bed)  # type: BedTool
    os.makedirs(path.parent.as_posix(), exist_ok=True)
    bed.saveas(path.as_posix(), compressed=False)
    return bed


def _compute_file_cached(file: str, func, cache_folder=_DEFAULT_CACHE_FOLDER):
    path = Path(cache_folder + file)
    if path.exists():
        return BedTool(path.as_posix())
    fn = func()
    fn = BedTool([inter for inter in fn])
    os.makedirs(path.parent.as_posix(), exist_ok=True)
    fn.saveas(path.as_posix(), compressed=False)
    return fn
