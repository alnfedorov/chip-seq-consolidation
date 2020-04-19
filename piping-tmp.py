import os
from asyncio.subprocess import create_subprocess_shell, PIPE, create_subprocess_exec
ROOT = "/data/encode/H3K27me3/ENCSR793NQA/original/duplicates-bam/"


read, write = os.pipe()
process_1 = await create_subprocess_exec("bedtools", *(f"genomecov -bga -ibam {ROOT}tmps.sorted.bam".split(' ')), stdout=write)
os.close(write)

read2, write2 = os.pipe()
process_2 = await create_subprocess_exec('bedtools', *("slop -i /dev/stdin -g /data/hg19/chromInfo.txt -b 0").split(' '),
                                         stdin=read, stdout=write2)
os.close(read)
os.close(write2)

read3, write3 = os.pipe()
process_3 = await create_subprocess_exec('bedClip', *("/dev/stdin /data/hg19/chromInfo.txt /dev/stdout").split(' '),
                                         stdin=read2, stdout=write3)
os.close(read2)
os.close(write3)

file = open(f"{ROOT}tmp.python", 'w')
p4 = await create_subprocess_exec("sort", "-k1,1", "-k2,2", "/dev/stdin", stdin=read3, stdout=file)
os.close(read3)
await p4.communicate()
file.close()

# read, write = os.pipe()
# process_1 = await create_subprocess_exec('ls', stdout=write)
# os.close(write)
# process_2 = await create_subprocess_shell('sleep 3 && wc', stdin=read, stdout=PIPE)
# os.close(read)
# await process_2.stdout.read()

# class stderr by defaults is piped as well
# Надо ли это на самом деле? Вроде бы да. Так быстрее как минимум и не создаются промежуточные файлы
# class PipeableCommand():
#   stdin
# Пофиг. Пока без этого обойдемся. Но на будующее сделать будет полезно.


read, write = os.pipe()
p1 = await create_subprocess_shell(f"bedtools genomecov -bga -ibam {ROOT}tmps.sorted.bam", stdout=PIPE)
os.close(write)

read2, write2 = os.pipe()
p2 = await create_subprocess_shell(f"bedtools slop -i /dev/stdin -g /data/hg19/chromInfo.txt -b 0",
                                   stdin=read, stdout=write2)
os.close(read)
os.close(write2)

p3 = await create_subprocess_shell("bedClip /dev/stdin /data/hg19/chromInfo.txt /dev/stdout",
                                   stdin=read2, stdout=PIPE)
os.close(read2)

p4 = await create_subprocess_shell(
    f"sort -k1,1 -k2,2n /dev/stdin > {ROOT}tmp.python",
    stdin=p3.stdout)
await p4.wait()

# LC_COLLATE=C bedtools genomecov -bga -ibam tmps.sorted.bam | \
# bedtools slop -i /dev/stdin -g /data/hg19/chromInfo.txt -b 0 | \
# bedClip /dev/stdin /data/hg19/chromInfo.txt /dev/stdout | \
# sort -k1,1 -k2,2n /dev/stdin > tmp.txt && \
# bedGraphToBigWig tmp.txt /data/hg19/chromInfo.txt tmps.bw


