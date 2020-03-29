import os
import tempfile
import logging
from .run import run, _move


logger = logging.getLogger(__name__)


def MarkDuplicates(path: str, saveto: str = None, inplace: bool = True):
    assert os.path.isfile(path)

    # picard MarkDuplicates \
    #       I=input.bam \
    #       O=marked_duplicates.bam
    saveto = saveto if saveto and not inplace else tempfile.mktemp(dir=os.path.dirname(path))
    run(
        [
            "picard", "MarkDuplicates", f"I={path}", f"O={saveto}", f"M=/dev/null", "VALIDATION_STRINGENCY=LENIENT"
        ], logger, logbefore=f"Start picard MarkDuplicates for {path}", logafter="MarkDuplicates finished"
    )
    if inplace:
        saveto = _move(saveto, path)
    return saveto
