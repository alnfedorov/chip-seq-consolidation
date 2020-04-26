import os
import tempfile
import logging
from .utils import run
from .piping import pipe

logger = logging.getLogger(__name__)


async def fastqc(path: str, saveto: str = None, threads: int = 1):
    assert os.path.isfile(path)
    odir = tempfile.mkdtemp()[1] if saveto is None else saveto
    await run(["fastqc", f"--outdir={odir}", f"--threads={threads}", path], logger,
              logbefore=f"fastqc for file {path}, saveto {odir}", logafter=f"fastqc finished")
    return saveto
