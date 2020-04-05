import asyncio
from pipeline import endtoend
from pipeline.utils import BamMeta, ExperimentMeta

ROOT = '/data/encode/H3K4me3/ENCSR956CTX'
TARGERT = "H3K4me3"
ACCESSION = "ENCSR956CTX"

bam = [
    # biological replicates
    BamMeta("breplica1", TARGERT, "ENCFF578CVN", paired=False, readlen=36),
    BamMeta("breplica2", TARGERT, "ENCFF074TRC", paired=False, readlen=36),

    # control
    BamMeta("control1", "control", "ENCFF151EVR", paired=False, readlen=36),
    BamMeta("control2", "control", "ENCFF468JNY", paired=False, readlen=36)
]
bam = {meta.accession: meta for meta in bam}

experiments = [
    ExperimentMeta("breplica1", TARGERT, ACCESSION, [bam["ENCFF578CVN"]], [bam["ENCFF151EVR"], bam["ENCFF468JNY"]]),
    ExperimentMeta("breplica2", TARGERT, ACCESSION, [bam["ENCFF074TRC"]], [bam["ENCFF151EVR"], bam["ENCFF468JNY"]]),
]
asyncio.run(endtoend.run(ROOT, experiments))
