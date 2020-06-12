import json
import os
import tempfile
from collections import defaultdict
from typing import Tuple

import numpy as np
from pybedtools import BedTool

from data.hg19.annotation import BLACKLISTED
from data.pod import BamMeta, ChIPseqReplicaMeta, AccessionInfo
from utils import bed


def _consensus(meta: dict, root: str, balcklisted: str, consensusfn, minpeaks: int = 20):
    # 1. make consensus bed files, simple strategy -> weighted majority voting
    sampling = "original"
    weighting = {"pepr": len(meta["replicates"]), "macs2": 1, "epic2": 1, "peakranger": 1}

    bedfiles = []
    bedweights = []

    pepr = os.path.join(root, meta["replicated-peak-calling"][sampling]["pepr"])
    pepr = BedTool(pepr).sort()
    if len(pepr) < minpeaks:
        print(f"Skipped pepr, not enough peaks")
    else:
        bedfiles.append(pepr)
        bedweights.append(weighting["pepr"])

    for r in meta["replicates"]:
        peaks = r["peak-calling"][sampling]
        for k, path in peaks.items():
            bedfile = BedTool(os.path.join(root, path)).sort()
            if len(bedfile) < minpeaks:
                print(f"Skipped bed file {path}, not enough peaks")
                continue
            bedfiles.append(bedfile)
            bedweights.append(weighting[k])

    # l1 normalization for simplicity
    bedweights = np.asarray(bedweights, dtype=np.float32)
    bedweights = bedweights / np.linalg.norm(bedweights, ord=1)

    # compute consensus and save as a file
    consensus = consensusfn(bedfiles, bedweights, threshold=0.75)
    fname = tempfile.mkstemp()[1]
    consensus = consensus.subtract(
        BedTool(balcklisted).sort()
    ).sort().saveas(fname).fn
    return consensus


def fromjson(root: str, blacklisted: str = BLACKLISTED.fn,
             consensusfn=bed.consensus) -> Tuple[ChIPseqReplicaMeta, ...]:
    file = os.path.join(root, "meta.json")
    assert os.path.isfile(file) and os.path.isfile(blacklisted), f"Files not found {file}, {blacklisted}"
    with open(file, 'r') as file:
        meta = json.load(file)

    # 1. make consensus peaks set
    consensus = _consensus(meta, root, blacklisted, consensusfn)

    # 2. parse bam files information
    bam = defaultdict(dict)  # sampling -> accession -> bam meta
    for accession, bammeta in meta["files"].items():
        readlen = bammeta["readlen"]
        samplings = [k for k in bammeta.keys() if k != "readlen"]
        for sampling in samplings:
            bam[sampling][accession] = BamMeta(
                path=os.path.join(root, bammeta[sampling]["filtered"]),
                estimated_fragment_size=bammeta[sampling]["fragment_size"],
                numreads=bammeta[sampling]["numreads"],
                readlen=readlen
            )

    # 3. parse replicas information
    replicas = []
    for replica in meta["replicates"]:
        accessions = AccessionInfo(meta["accession"], replica["treatment"], replica["control"])
        treatment, control = defaultdict(list), defaultdict(list)
        for sampling in bam.keys():
            treatment[sampling] = tuple(bam[sampling][t] for t in accessions.treatment)
            control[sampling] = tuple(bam[sampling][c] for c in accessions.control)
        samplemods = tuple(sorted(treatment.keys()))
        assert samplemods == tuple(sorted(control.keys()))

        replicas.append(
            ChIPseqReplicaMeta(
                accesions=accessions,
                target=meta["target"],
                peaks=consensus,
                treatment=dict(treatment),
                control=dict(control),
                samplemods=samplemods
            )
        )
    return tuple(replicas)
