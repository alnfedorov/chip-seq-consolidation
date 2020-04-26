import os
import sys
import asyncio
import tempfile
from .piping import pipe
from asyncio.subprocess import create_subprocess_shell


# It should be piped really but doesn't work for some reason...
@pipe(writearg=("saveto", "f"), bam="p")
async def bam_to_bedpe(bam: str, saveto: str = None, maxthreads: int = 1) -> str:
    assert maxthreads >= 1
    saveto = saveto if saveto else tempfile.mkstemp()[1]
    cmd = f"sambamba sort -t {maxthreads} -p -n -M {bam} -o /dev/stdout | " \
          f"samtools fixmate --threads {maxthreads} -r -O bam /dev/stdin /dev/stdout | " \
          f"bedtools bamtobed -i /dev/stdin -bedpe > {saveto}"
    process = await create_subprocess_shell(cmd, stderr=sys.stdout)
    assert await process.wait() == 0
    assert os.path.exists(saveto)
    return saveto


@pipe(writearg=("saveto", "p"), bam="p")
async def bam_to_bigwig(bam: str, chrominfo: str, saveto: str) -> str:
    tmp = tempfile.mkstemp()[1]
    cmd = f'LC_COLLATE=C bedtools genomecov -bga -ibam "{bam}" | ' \
          f'bedtools slop -i /dev/stdin -g "{chrominfo}" -b 0 | ' \
          f'bedClip /dev/stdin "{chrominfo}" /dev/stdout | ' \
          f'sort -k1,1 -k2,2n /dev/stdin > "{tmp}" && ' \
          f'bedGraphToBigWig "{tmp}" "{chrominfo}" "{saveto}"'
    process = await asyncio.subprocess.create_subprocess_shell(cmd)
    await process.wait()
    os.remove(tmp)

    assert process.returncode == 0
    assert os.path.exists(saveto)
    return saveto
