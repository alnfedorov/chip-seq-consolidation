import os
from dataclasses import dataclass, field
from .config import UNFILTERED_READS_PREFIX, UNIQUE_READS_PREFIX, DUPLICATED_READS_PREFIX, BROAD_HISTONE_MARKS


def is_broad(hm: str):
    return hm in BROAD_HISTONE_MARKS


def make_filename(*path: [str], mode: str, name: str, accession: str, format: str = "bam", reads: int = None):
    prefix = {"unique": UNIQUE_READS_PREFIX,
              "duplicated": DUPLICATED_READS_PREFIX,
              "unfiltered": UNFILTERED_READS_PREFIX}[mode]
    filename = [prefix]
    if reads is not None:
        filename.append(f"{reads // 10**6}mln")
    filename += [name, accession]
    filename = "-".join(filename) + f".{format}"
    return os.path.join(*path, filename)


@dataclass(frozen=False, unsafe_hash=True)
class BamMeta:
    """
    Note, BamMeta with the same name-target-accession-paired will always be considered identical,
    no matter actual paths to the bam files
    """
    name: str = field(hash=True, compare=True)
    target: str = field(hash=True, compare=True)
    accession: str = field(hash=True, compare=True)
    paired: bool = field(hash=True, compare=True)
    readlen: int = field(hash=True, compare=True)
    simulated: bool = False
    reads: int = None
    # paths to the bam files
    unfiltered: str = None
    duplicated: str = None
    unique: str = None


@dataclass(frozen=False)
class ExperimentMeta:
    name: str
    target: str
    accession: str
    treatment: [BamMeta]
    control: [BamMeta]
    paired_data: bool = None

    def __post_init__(self):
        # all reads must either paired or not
        assert all(t.paired for t in self.treatment) or all(not t.paired for t in self.treatment)
        assert all(c.paired for c in self.control) or all(not c.paired for c in self.control)
        assert self.control[0].paired == self.treatment[0].paired
        self.paired_data = self.control[0].paired

    def _paths(self, bam_meta: [BamMeta], mode: str) -> [str]:
        if mode == "duplicated":
            return [m.duplicated for m in bam_meta]
        elif mode == "unique":
            return [m.unique for m in bam_meta]
        else:
            raise ValueError(f"mode expected to be one of the (duplicated, unique), got {mode}")

    def allcontrol(self, mode: str) -> [str]:
        return self._paths(self.control, mode)

    def alltreatment(self, mode: str) -> [str]:
        return self._paths(self.treatment, mode)


__all__ = [is_broad, BamMeta, ExperimentMeta]
