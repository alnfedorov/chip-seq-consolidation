import os
import tempfile
import logging
from .utils import run, replace_bam
from .piping import pipe

logger = logging.getLogger(__name__)


@pipe(writearg=("saveto", "p"), path="p")
async def MarkDuplicates(path: str, saveto: str = None) -> str:
    assert os.path.exists(path)

    saveto = saveto if saveto else tempfile.mkstemp()[1]
    # saveto = saveto if saveto and not inplace else tempfile.mkstemp()[1]
    await run([
            "picard", "MarkDuplicates", f"I={path}", f"O={saveto}", f"M=/dev/null", "VALIDATION_STRINGENCY=LENIENT",
        ], logger, logbefore=f"Start picard MarkDuplicates for {path}", logafter="MarkDuplicates finished"
    )
    # if inplace:
    #     saveto = replace_bam(saveto, path)
    return saveto
