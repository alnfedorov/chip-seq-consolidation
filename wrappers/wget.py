import os
import tempfile
import logging
import requests
import subprocess
from subprocess import run

logger = logging.getLogger(__name__)


def wget(url: str, saveto: str = None):
    file = tempfile.mkstemp()[1] if saveto is None else saveto

    result = run(["wget", "--continue", "--retry-connrefused", "--tries=0", "--timeout=5", "-O", file, url],
                 check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logger.debug(result.stdout.decode())
    return file

# logging.basicConfig(level=logging.NOTSET)
#
# a = download("ENCFF534UUQ", "ENCFF534UUQ.bed.gz", "b3ecbb186d2712af0541a2605259457c", "/tmp/tmp.bed")
