import asyncio
from pipeline import endtoend
from pipeline.meta import BamMeta, ExperimentMeta

ROOT = '/data/encode/H3K4me3/ENCSR849YFO'
TARGERT = "H3K4me3"
ACCESSION = "ENCSR849YFO"

bam = [
    # biological replicates
    BamMeta("breplica1", TARGERT, "ENCFF132AQO", paired=True, readlen=76),
    BamMeta("breplica2", TARGERT, "ENCFF754VTC", paired=True, readlen=76),

    # control
    BamMeta("control1", "control", "ENCFF458IRD", paired=True, readlen=76),
    BamMeta("control2", "control", "ENCFF626LZE", paired=True, readlen=76)
]
bam = {meta.accession: meta for meta in bam}

experiments = [
    ExperimentMeta("breplica1", TARGERT, ACCESSION, [bam["ENCFF132AQO"]], [bam["ENCFF458IRD"], bam["ENCFF626LZE"]]),
    ExperimentMeta("breplica2", TARGERT, ACCESSION, [bam["ENCFF754VTC"]], [bam["ENCFF458IRD"], bam["ENCFF626LZE"]]),
]
asyncio.run(endtoend.run(ROOT, experiments))
