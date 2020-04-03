import os
import logging
import asyncio
import subprocess
from subprocess import CalledProcessError, STDOUT, PIPE


async def run(cmd: [str], logger: logging.Logger = None, logbefore: str = None, logafter: str = None) -> str:
    """Run command cmd and return stdout. With possible logging before and after cmd"""
    if logger is None:
        logger = logging.getLogger(__name__)

    if logbefore:
        logger.debug(logbefore)
    try:
        result = await asyncio.get_event_loop().run_in_executor(
            None, lambda: subprocess.run(cmd, check=True, stdout=PIPE, stderr=STDOUT))
    except CalledProcessError as e:
        logger.error(e.stdout.decode())
        raise e

    stdout = result.stdout.decode()
    logger.debug(stdout)

    if logafter:
        logger.debug(logafter)
    return stdout


def _move(src, dst):
    os.remove(dst)
    if os.path.exists(dst+'.bai'):
        os.remove(dst + '.bai')  # clear index, it is not valid anymore
    os.renames(src, dst)
    return dst
