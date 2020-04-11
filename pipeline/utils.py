import os
import sys
import logging
import asyncio
import itertools
from typing import Awaitable
from dataclasses import dataclass, field
from .config import UNFILTERED_READS_PREFIX, UNIQUE_READS_PREFIX, DUPLICATED_READS_PREFIX, BROAD_HISTONE_MARKS
from logging import Handler, LogRecord, Formatter


def is_broad(hm: str):
    return hm in BROAD_HISTONE_MARKS


async def batched_gather(tasks: [Awaitable], batch_size: int, verbose=False):
    results = []
    while tasks:
        if verbose:
            print("tasks left ", len(tasks))
        batch = tasks[:batch_size]
        if len(batch) == 0:
            break
        batch = await asyncio.gather(*batch)
        tasks = tasks[batch_size:]
        results.extend(batch)
    return results


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


def config_logging(logfolder: str, level: str = "NOTSET"):
    assert os.path.isdir(logfolder)
    logging.basicConfig(level=level, format="%(levelname)s: %(name)s: %(message)s")
    logger = logging.getLogger()

    handler = PerFileHandler(logfolder)
    handler.setLevel(level)
    handler.setFormatter(Formatter("%(asctime)s: %(levelname)s: %(message)s"))
    logger.addHandler(handler)


class PerFileHandler(Handler):
    def __init__(self, logfolder, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logfolder = logfolder
        self.__files = {}

    def emit(self, record: LogRecord) -> None:
        file = os.path.join(self.logfolder, record.name)
        if file not in self.__files:
            self.__files[file] = open(file, 'a')
        print(self.format(record), file=self.__files[file], flush=True)

    def flush(self) -> None:
        for f in self.__files.values():
            f.close()
        self.__files.clear()


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
    # unique: str = None


@dataclass(frozen=False)
class ExperimentMeta:
    name: str
    target: str
    accession: str
    treatment: [BamMeta]
    control: [BamMeta]

    # paired_data: bool

    def __post_init__(self):
        # all reads must either paired or not
        assert all(t.paired for t in self.treatment) or all(not t.paired for t in self.treatment)
        assert all(c.paired for c in self.control) or all(not c.paired for c in self.control)
        assert self.control[0].paired == self.treatment[0].paired
        object.__setattr__(self, 'paired_data', self.control[0].paired)

    def _paths(self, bam_meta: [BamMeta], mode: str) -> [str]:
        if mode == "duplicated":
            return [m.duplicated for m in bam_meta]
        # elif mode == "unique":
        #     return [m.unique for m in bam_meta]
        else:
            raise ValueError(f"mode expected to be one of the (duplicated, unique), got {mode}")

    def allcontrol(self, mode: str) -> [str]:
        return self._paths(self.control, mode)

    def alltreatment(self, mode: str) -> [str]:
        return self._paths(self.treatment, mode)


__all__ = [is_broad, BamMeta, ExperimentMeta]
