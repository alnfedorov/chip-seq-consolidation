import os
import json
import tempfile
import numpy as np
from typing import Callable, Tuple
from pybedtools import BedTool
from utils import bed
from data.pod import SimulatedBigWig, BigWigMeta, ChIPseqReplicaMeta
from data.hg19.annotation import BLACKLISTED


def fromjson(root: str, blacklisted: str = BLACKLISTED.fn,
             consensusfn=bed.consensus) -> Tuple[ChIPseqReplicaMeta, ...]:
    file = os.path.join(root, "meta.json")
    assert os.path.isfile(file) and os.path.isfile(blacklisted), f"Files not found {file}, {blacklisted}"
    with open(file, 'r') as file:
        meta = json.load(file)
    blacklisted = BedTool(blacklisted)

    # 1. make consensus bed files, simple strategy -> weighted majority voting
    sampling = "original"
    weighting = {"pepr": len(meta["replicates"]), "macs2": 1, "epic2": 1, "peakranger": 1}

    bedfiles = []
    bedweights = []

    pepr = os.path.join(root, meta["replicated-peak-calling"][sampling]["pepr"])
    pepr = bed.threshold_qvalue(BedTool(pepr)).sort()
    if len(pepr) < 20:
        print(f"Skipped pepr, not enough peaks")
    else:
        bedfiles.append(pepr)
        bedweights.append(weighting["pepr"])

    for r in meta["replicates"]:
        peaks = r["peak-calling"][sampling]
        for k, path in peaks.items():
            bedfile = bed.threshold_qvalue(BedTool(os.path.join(root, path))).sort()
            if len(bedfile) < 20:
                print(f"Skipped bed file {path}, not enough peaks")
                continue
            bedfiles.append(bedfile)
            bedweights.append(weighting[k])

    # l1 normalization for simplicity
    bedweights = np.asarray(bedweights, dtype=np.float32)
    bedweights = bedweights / np.linalg.norm(bedweights, ord=1)

    # compute consensus and save as a file
    consensus = consensusfn(bedfiles, bedweights, threshold=0.5)
    fname = tempfile.mkstemp()[1]
    consensus = consensus.subtract(
        blacklisted.sort()
    ).sort().saveas(fname).fn

    # 2. parse replicates meta information
    replicas = []
    for r in meta["replicates"]:
        # parse big-wig enrichment files
        enrichment = {k: BigWigMeta(os.path.join(root, r["enrichment"][k]))
                      for k in ["original", "subsampled-q0.25", "subsampled-q0.5", "subsampled-q0.75"]}
        enrichment = {k.replace("subsampled-", ""): v for k, v in enrichment.items()}
        enrichment = SimulatedBigWig(**enrichment)
        # enrichment = {}
        # for accession, info in meta['files'].items():
        #     readlen = info["readlen"]
        #     files = {
        #         k: BigWigMeta(os.path.join(root, info[k]["bigwig"]), readlen, info[k]["numreads"], accession)
        #         for k in ["original", "subsampled-q0.25", "subsampled-q0.5", "subsampled-q0.75"]
        #     }
        #     files = {k.replace("subsampled-", ""): v for k, v in files.items()}
        #     assert all(os.path.isfile(f.path) for f in files.values())
        #     bigwigs[accession] = SimulatedBigWig(**files)

        replicas.append(
            ChIPseqReplicaMeta(
                experiment_accession=meta["accession"],
                target=meta["target"].lower(),
                enrichment=enrichment,
                peaks=consensus,
                treatment=r["treatment"],
                control=r["control"]
            )
        )
    return tuple(replicas)
