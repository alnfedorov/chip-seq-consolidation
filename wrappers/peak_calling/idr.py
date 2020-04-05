import os
import shutil
import logging
import tempfile
from wrappers.utils import run

logger = logging.getLogger(__name__)


# Works only with 2 files, see:
# https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L523
async def idr(files: [str], rankby: str, idrcutoff: float, format: str = None, saveto: str = None,
              plotto: str = None, logto: str = None):
    assert rankby in ("signal.value", "p.value", "q.value", "columnIndex")
    assert len(files) == 2, "idr is implemented only for 2 replicates"
    assert idrcutoff > 0
    if format is None:
        f1_format, f2_format = (f.split('.')[-1] for f in files)
        if f1_format == f2_format:
            format = f1_format
    assert format in ("narrowPeak", "broadPeak")
    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    cmd = ["idr", "--samples", *files, f"--rank={rankby}", f"--output-file={saveto}",
           f"--idr-threshold={idrcutoff}", f"--input-file-type={format}"]
    if plotto is not None:
        cmd.append("--plot")
    if logto is not None:
        cmd.append(f"--log-output-file={logto}")
    await run(cmd, logger, logbefore=f"Starting idr with cmd {' '.join(cmd)}", logafter="idr finished")

    if plotto:
        shutil.move(saveto+".png", plotto)
    return saveto
