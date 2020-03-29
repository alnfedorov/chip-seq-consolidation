import os
import logging
import subprocess


def run(cmd: [str], logger: logging.Logger = None, logbefore: str = None, logafter: str = None) -> str:
    """Run command cmd and return stdout. With possible logging before and after cmd"""
    if logger is None:
        logger = logging.getLogger(__name__)

    if logbefore:
        logger.debug(logbefore)
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = result.stdout.decode()
    logger.debug(result)

    if logafter:
        logger.debug(logafter)
    return result


def _move(src, dst):
    os.remove(dst)
    if os.path.exists(dst+'.bai'):
        os.remove(dst + '.bai')  # clear index, it is not valid anymore
    os.renames(src, dst)
    return dst
