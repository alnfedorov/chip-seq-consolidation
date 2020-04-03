import tempfile
import logging
from .utils import run

logger = logging.getLogger(__name__)


async def wget(url: str, saveto: str = None):
    file = tempfile.mkstemp()[1] if saveto is None else saveto
    await run(["wget", "--continue", "--retry-connrefused", "--tries=0", "--timeout=5", "-O", file, url], logger,
              logbefore=f"wget url: {url}", logafter="wget finished")
    return file
