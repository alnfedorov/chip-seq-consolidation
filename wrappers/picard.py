import os
import tempfile
import logging
from .utils import run, _move


logger = logging.getLogger(__name__)


async def MarkDuplicates(path: str, saveto: str = None, inplace: bool = False):
    assert os.path.isfile(path)

    saveto = saveto if saveto and not inplace else tempfile.mkstemp(dir=os.path.dirname(path))[1]
    await run([
            "picard", "MarkDuplicates", f"I={path}", f"O={saveto}", f"M=/dev/null", "VALIDATION_STRINGENCY=LENIENT"
        ], logger, logbefore=f"Start picard MarkDuplicates for {path}", logafter="MarkDuplicates finished"
    )
    if inplace:
        saveto = _move(saveto, path)
    return saveto
