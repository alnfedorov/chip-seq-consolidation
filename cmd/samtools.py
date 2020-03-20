import os
import sys
import tempfile
import logging as log
import subprocess
from subprocess import run


def qc_report(files: [str], saveto: str):
    """Run samtools stats and samtools flagstats for each file. saveto a is destination folder"""
    assert os.path.exists(saveto) and os.path.isdir(saveto)
    for file in files:
        assert ".bam" in file
        to = os.path.split(file)[-1]
        stats(file, os.path.join(saveto, to.replace(".bam", ".samtools-stats")))
        flagstat(file, os.path.join(saveto, to.replace(".bam", ".samtools-flagstat")))


# samtools flagstat in.bam -> simple stats mapped/unmapped/passed etc
# samtools stats -> statistics about alignment etc
def stats(path: str, saveto: str = None):
    assert os.path.isfile(path) and ".bam" in path

    saveto = saveto if saveto else tempfile.mktemp()
    log.debug(f"Start samtools stats for {path}")
    result = run(["samtools", "stats", path],
                 check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stats = result.stdout.decode()
    log.debug(stats)
    log.debug("samtools stats finished")

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto


def flagstat(path: str, saveto: str = None):
    assert os.path.isfile(path) and ".bam" in path

    saveto = saveto if saveto else tempfile.mktemp()
    log.debug(f"Start samtools flagstat for {path}")
    result = run(["samtools", "flagstat", path],
                 check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stats = result.stdout.decode()
    log.debug(stats)
    log.debug("samtools flagstat finished")

    with open(saveto, 'w') as file:
        file.write(stats)
    return saveto
