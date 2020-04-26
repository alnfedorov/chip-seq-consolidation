import tempfile
import logging
from .utils import run
from .piping import pipe
from typing import Awaitable

logger = logging.getLogger(__name__)


@pipe(writearg=("saveto", "f"))
async def wget(url: str, saveto: str = None) -> str:
    file = tempfile.mkstemp()[1] if saveto is None else saveto
    await run(["wget", "--continue", "--retry-connrefused", "--tries=0", "--timeout=5", "-O", file, url], logger,
              logbefore=f"wget url: {url}", logafter="wget finished")
    return file
