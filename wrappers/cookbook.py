import os
import sys
import asyncio
import tempfile
from asyncio.subprocess import create_subprocess_shell


async def bam_to_bedpe(bam: str, saveto: str = None, maxthreads: int = 1):
    assert maxthreads >= 1
    saveto = saveto if saveto else tempfile.mkstemp()[1]
    cmd = f"sambamba sort -t {maxthreads} -p -n -M {bam} -o /dev/stdout | " \
          f"samtools fixmate --threads {maxthreads} -r -O bam /dev/stdin /dev/stdout | " \
          f"bedtools bamtobed -i /dev/stdin -bedpe > {saveto}"
    process = await create_subprocess_shell(cmd, stderr=sys.stdout)
    assert await process.wait() == 0
    assert os.path.exists(saveto)
    return saveto


# LC_COLLATE=C bedtools genomecov -bga -ibam tmps.sorted.bam | \
# bedtools slop -i /dev/stdin -g /data/hg19/chromInfo.txt -b 0 | \
# bedClip /dev/stdin /data/hg19/chromInfo.txt /dev/stdout | \
# sort -k1,1 -k2,2n /dev/stdin > tmp.txt && \
# bedGraphToBigWig tmp.txt /data/hg19/chromInfo.txt tmps.bw

# bedtools genomecov -bga -ibam tmps.sorted.bam > tmp1.bdg && \
# bedtools slop -i tmp1.bdg -g /data/hg19/chromInfo.txt -b 0 > tmp2.bdg && \
# bedClip tmp2.bdg /data/hg19/chromInfo.txt tmp3.bdg && \
# LC_COLLATE=C sort -k1,1 -k2,2n tmp3.bdg > tmp4.bdg && \
# bedGraphToBigWig tmp4.bdg /data/hg19/chromInfo.txt tmps2.bw
def bam_to_bigwig(bam: str, chrominfo: str, saveto: str):
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
