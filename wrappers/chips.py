import os
import shutil
import logging
import tempfile
from wrappers import sambamba, samtools
from .utils import run

logger = logging.getLogger(__name__)


async def learn(bam: str, is_paired: bool, peaks: str, score_column: int, saveto: str = None,
                format: str = "bed", scale_outliers: bool = False):
    assert os.path.isfile(bam) and os.path.isfile(peaks)
    assert score_column > 0 and format in ("bed", "homer")
    if not os.path.exists(bam + ".bai"):
        await sambamba.index(bam)

    output = tempfile.mkstemp()[1]
    cmd = ["chips", "learn", "-b", bam, "-p", peaks, "-t", format, "-c", str(score_column), "-o", output]
    if is_paired:
        cmd.append("--paired")
    if scale_outliers:
        cmd.append("--scale-outliers")
    await run(cmd, logger, logbefore=f"chips learn with cmd {' '.join(cmd)}",
              logafter="finished learning ChIP-seq experiment parameters")
    if saveto is not None:
        shutil.move(output + ".json", saveto)
    else:
        saveto = output + ".json"
    return saveto


async def simreads(fagenome: str, model: str, peaks: str = None, numcopies: int = 100, numreads: int = 20*10**6,
                   readlen: int = 36, paired: bool = False, saveto: str = None,  threads: int = 1,
                   bam: str = None, score_column: int = None, scale_outliers: bool = False, sequencer: str = None):
    assert peaks is None or os.path.isfile(peaks) and os.path.isfile(fagenome) and os.path.isfile(model)
    assert numcopies > 1 and numreads > 1 and readlen > 26 and threads > 1
    assert bam is None and (peaks is not None and score_column is not None) or \
           bam is not None and (peaks is not None and score_column is None)
    if not os.path.exists(fagenome + ".fai"):
        await samtools.faidx(fagenome)
    output = tempfile.mkstemp()[1]
    cmd = ["chips", "simreads", "-p", peaks, "-f", fagenome, "-o", output, "--thread", threads, "-t", "bed",
           "--numcopies", numcopies, "--numreads", numreads, "--readlen", readlen, "--model", model]
    if paired:
        cmd.append("--paired")
    if bam:
        if not os.path.exists(bam + ".bai"):
            await sambamba.index(bam)
        cmd += ["-b", bam]
    else:
        assert score_column is not None
        cmd += ["-c", score_column]
    if scale_outliers:
        cmd.append("--scale-outliers")
    if sequencer:
        cmd += ["--sequencer", sequencer]
    cmd = [str(e) for e in cmd]
    await run(cmd, logger, logbefore=f"chips simreads with cmd: {' '.join(cmd)}",
              logafter="finished learning ChIP-seq experiment parameters")
    if saveto is not None:
        shutil.move(output + ".fastq", saveto)
    else:
        saveto = output + ".fastq"
    return saveto


async def wce(fagenome: str, numcopies: int = 100, numreads: int = 20*10**6, readlen: int = 36, paired: bool = False,
              saveto: str = None,  threads: int = 1, bam: str = None, sequencer: str = None):
    assert os.path.isfile(fagenome)
    assert numcopies > 1 and numreads > 1 and readlen > 26 and threads > 1
    if not os.path.exists(fagenome + ".fai"):
        await samtools.faidx(fagenome)
    output = tempfile.mkstemp()[1]
    cmd = ["chips", "simreads", "-t", "wce", "-f", fagenome, "-o", output, "--thread", threads,
           "--numcopies", numcopies, "--numreads", numreads, "--readlen", readlen]
    if paired:
        cmd.append("--paired")
    if bam:
        if not os.path.exists(bam + ".bai"):
            await sambamba.index(bam)
        cmd += ["-b", bam]
    # if scale_outliers:
    #     cmd.append("--scale-outliers")
    if sequencer:
        cmd += ["--sequencer", sequencer]

    cmd = [str(e) for e in cmd]
    await run(cmd, logger, logbefore=f"chips simreads with cmd: {' '.join(cmd)}",
              logafter="finished learning ChIP-seq experiment parameters")
    if saveto is not None:
        shutil.move(output + ".fastq", saveto)
    else:
        saveto = output + ".fastq"
    return saveto
