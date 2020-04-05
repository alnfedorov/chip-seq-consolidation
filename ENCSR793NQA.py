import asyncio
from pipeline import endtoend
from pipeline.utils import BamMeta, ExperimentMeta

ROOT = '/data/encode/H3K27me3/ENCSR793NQA'
TARGERT = "H3K27me3"
ACCESSION = "ENCSR793NQA"

bam = [
    # biological replicates
    BamMeta("breplica1", TARGERT, "ENCFF933CPQ", paired=True, readlen=101),
    BamMeta("breplica2", TARGERT, "ENCFF825SKJ", paired=True, readlen=101),

    # control
    BamMeta("control1", "control", "ENCFF941KVW", paired=True, readlen=101),
    BamMeta("control2", "control", "ENCFF379ZDL", paired=True, readlen=101)
]
bam = {meta.accession: meta for meta in bam}

experiments = [
    ExperimentMeta("breplica1", TARGERT, ACCESSION, [bam["ENCFF933CPQ"]], [bam["ENCFF941KVW"], bam["ENCFF379ZDL"]]),
    ExperimentMeta("breplica2", TARGERT, ACCESSION, [bam["ENCFF825SKJ"]], [bam["ENCFF941KVW"], bam["ENCFF379ZDL"]]),
]
asyncio.run(endtoend.run(ROOT, experiments))
