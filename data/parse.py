import os
import json
import tempfile
import numpy as np
from pybedtools import BedTool
from utils.bed import compute_conservative_regions
from data.pod import SimulatedBigWig, BigWigMeta, ReplicaMeta


def fromjson(root: str, blacklisted: str):
    file = os.path.join(root, "meta.json")
    assert os.path.isfile(file) and os.path.isfile(blacklisted), f"Files not found {file}, {blacklisted}"
    with open(file, 'r') as file:
        meta = json.load(file)
    blacklisted = BedTool(blacklisted)

    # 1. parse big-wig files
    bigwigs = {}
    for accession, info in meta['files'].items():
        readlen = info["readlen"]
        files = {
            k.replace("subsampled-", ""): BigWigMeta(os.path.join(root, info[k]["bigwig"]), readlen, info[k]["numreads"])
            for k in ["original", "subsampled-q0.25", "subsampled-q0.5", "subsampled-q0.75"]
        }
        assert all(os.path.isfile(f.path) for f in files.values())
        bigwigs[accession] = SimulatedBigWig(**files)

    # 2. make consensus bed files, simple strategy -> weighted majority voting
    mode = "original"
    weighting = {"pepr": len(meta["replicated"]), "macs2": 1, "epic2": 1, "peakranger": 1}

    bedfiles = [meta["replicated-peak-calling"][mode]["pepr"]]
    bedweights = [weighting["pepr"]]

    for r in meta["replicated"]:
        peaks = r["peak-calling"][mode]
        for k, path in peaks.items():
            bedfiles.append(path)
            bedweights.append(weighting[k])

    bedfiles = [BedTool(os.path.join(root, f)) for f in bedfiles]
    # l1 normalization for simplicity
    bedweights = np.asarray(bedweights, dtype=np.float32)
    bedweights = bedweights / np.linalg.norm(bedweights, ord=1)

    # compute consensus and save as a file
    consensus = compute_conservative_regions(bedfiles, bedweights, threshold=0.5)
    fname = tempfile.mkstemp()[1]
    consensus = consensus.subtract(
        blacklisted.sort()
    ).sort().saveas(fname).fn

    replicas = []
    for r in meta["replicated"]:
        replicas.append(
            ReplicaMeta(
                target=meta["target"].lower(),
                treatment=tuple(bigwigs[acc] for acc in r["treatment"]),
                control=tuple(bigwigs[acc] for acc in r["control"]),
                peaks=consensus
            )
        )
    return replicas
