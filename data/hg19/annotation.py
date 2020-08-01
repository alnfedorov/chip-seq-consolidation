import pickle
import numpy as np
import os
from pybedtools import BedTool
from data.caching import _load_file_cached, _compute_file_cached
from pathlib import Path
from collections import defaultdict

__all__ = [
    "CHROMOSOME_SIZE", "UTR3", "UTR5", "AMBIGUOUS",
    "EXONS", "DOWNSTREAM_1K", "INTRONS", "UPSTREAM_1K"
]

# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
CHROMOSOME_SIZE = {
    "chr1":	249250621,
    "chr2":	243199373,
    "chr3":	198022430,
    "chr4":	191154276,
    "chr5":	180915260,
    "chr6":	171115067,
    "chr7":	159138663,
    "chr8":	146364022,
    "chr9":	141213431,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chrY":	59373566,
    "chrX":	155270560,
}

EFFECTIVE_GENOME_SIZE = 2.7e9

root = Path("/data/hg19/")
annotation = root.joinpath("annotation")

UTR3 =          _load_file_cached(annotation.joinpath("hg19_3UTR_exons.bed.gz").as_posix())
UTR5 =          _load_file_cached(annotation.joinpath("hg19_5UTR_exons.bed.gz").as_posix())
EXONS =         _load_file_cached(annotation.joinpath("hg19_coding_exons.bed.gz").as_posix())
INTRONS =       _load_file_cached(annotation.joinpath("hg19_introns.bed.gz").as_posix())
UPSTREAM_1K =   _load_file_cached(annotation.joinpath("hg19_upstream_1kb.bed.gz").as_posix())
DOWNSTREAM_1K = _load_file_cached(annotation.joinpath("hg19_downstream_1kb.bed.gz").as_posix())
AMBIGUOUS =     _load_file_cached(root.joinpath("hg19_ambiguous_regions.bed").as_posix())
BLACKLISTED =   _load_file_cached(root.joinpath("hg19-blacklist.v2.bed.gz").as_posix())

INTERGENIC = _compute_file_cached(
    file=Path.joinpath(annotation, "hg19_intergenic.bed").as_posix(),
    func=lambda: BedTool([]).cat(UTR3, UTR5, EXONS, INTRONS, UPSTREAM_1K, DOWNSTREAM_1K)
                            .sort().complement(genome='hg19').sort()
)

REGIONS = {
    "3'utr": UTR3,
    "5'utr": UTR5,
    "exons": EXONS,
    "introns": INTRONS,
    "upstream1k": UPSTREAM_1K,
    "downstream1k": DOWNSTREAM_1K,
    "intergenic": INTERGENIC,
    "ambiguous": AMBIGUOUS
}
