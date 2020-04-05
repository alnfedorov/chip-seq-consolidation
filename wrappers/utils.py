import os
import sys
import shutil
import logging
import asyncio
import subprocess
from subprocess import CalledProcessError, STDOUT, PIPE


async def run(cmd: [str], logger: logging.Logger = None, logbefore: str = None, logafter: str = None,
              logstdout: bool = True, stdout=PIPE, stderr=PIPE) -> subprocess.CompletedProcess:
    """Run command cmd and return stdout. With possible logging before and after cmd"""
    if logger is None:
        logger = logging.getLogger(__name__)
    if logbefore:
        logger.debug(logbefore)
    try:
        result = await asyncio.get_event_loop().run_in_executor(
            None, lambda: subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr))
    except CalledProcessError as e:
        if e.stdout is not None:
            logger.error(e.stdout.decode())
        raise e

    if logstdout and result.stdout is not None:
        logger.debug(result.stdout.decode())

    if logafter:
        logger.debug(logafter)
    return result


def replace_bam(src, dst):
    os.remove(dst)
    if os.path.exists(dst+'.bai'):
        os.remove(dst + '.bai')  # clear index, it is not valid anymore
    shutil.move(src, dst)
    return dst
