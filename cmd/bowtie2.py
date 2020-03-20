import os
import tempfile
import logging as log
import subprocess
from subprocess import run

# Buggy, not tested

def build_index(genome: str, threads: int = 1, seed: int = 0, saveto: str = None):
    assert os.path.isfile(genome) and genome.endswith(".fasta"), "Genome is expected to be a single fasta file"
    assert threads > 0
    assert saveto is None or os.path.isdir(saveto)

    threads, seed = str(threads), str(seed)
    dir = tempfile.mkdtemp() if saveto is None else saveto
    dir = dir + '/'

    log.debug(f"Build bowtie2 genome index for {genome}")
    result = run(["bowtie2-build", "--verbose", "--threads", threads, "--seed", seed, genome, dir], check=True,
                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.debug(result.stdout)
    log.debug("Build finished")

    return dir


class align:
    @staticmethod
    def single(fastq: [str], index: str, threads: int = 1, saveto: str = None):
        assert all(os.path.isfile(f) for f in fastq)
        assert threads > 0

        file = saveto if saveto else tempfile.mktemp() + ".sam"
        fastq = ','.join(fastq)
        log.debug(f"Start bowtie2 aligning {fastq}")
        result = run(["bowtie2", "--mm", "--threads", threads, "-x", index, "-U", fastq, "-S", file],
                     check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.debug(result.stdout)
        log.debug("Build finished")
        return file

    @staticmethod
    def paired(mate1: str, mate2: str, index: str, threads: int = 1, saveto: str = None):
        assert all(os.path.isfile(f) for f in [mate1, mate2])
        assert threads > 0

        file = saveto if saveto else tempfile.mktemp() + ".sam"
        log.debug(f"Start bowtie2 aligning {mate1}, {mate2}")
        result = run(["bowtie2", "--mm", "--threads", threads, "-x", index, "-1", mate1, "-2", mate2, "-S", file],
                     check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.debug(result.stdout)
        log.debug("Build finished")
        return file


