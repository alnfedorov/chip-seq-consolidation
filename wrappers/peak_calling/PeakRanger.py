import os
import tempfile
import logging
import asyncio
from math import log10
from .utils import _zero_significance_filler
from ..utils import run
from ..samtools import merge

logger = logging.getLogger(__name__)


async def ranger(treatment: [str], control: [str], saveto: str = None, threads: int = 1, pcutoff: float = None,
                 fdrcutoff: float = 0.05, format: str = "bam"):
    assert all(os.path.isfile(f) for f in treatment + control)
    assert threads > 0

    # merge treatment and control as a single file - peakranger requirement
    mtreatment, mcontrol = await asyncio.gather(merge(treatment, threads), merge(control, threads))
    tmpfile = os.path.join(tempfile.mkdtemp(), "ranger.bed")

    await _ranger(mtreatment, mcontrol, tmpfile, threads=threads, pcutoff=pcutoff, fdrcutoff=fdrcutoff, format=format)

    # free intermediates
    for f in [mtreatment, mcontrol]:
        if f not in treatment and f not in control:
            os.remove(f)

    logger.debug("Building narrow bed report file...")

    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]

    # We want to transform peakranger format to the
    # ENCODE narrow peak format https://genome.ucsc.edu/FAQ/FAQformat.html#format12
    # region_chr region_start region_end nearby_genes(6kbp) region_ID region_summits region_pvalue region_FDR region_strand region_treads region_creads
    # ------------->
    # chrom chromStart chromEnd name score strand signalValue pValue qValue source
    def converter(data):
        pvalue_filler = str(_zero_significance_filler(data, 6))
        qvalue_filler = str(_zero_significance_filler(data, 7))
        data = [[
            line[0],  # chrom
            line[1],  # chromStart
            line[2],  # chromEnd
            ".",  # name
            '0',  # score
            '.',  # strand
            '0',  # signalValue
            str(-log10(float(line[6]))) if float(line[6]) != 0 else pvalue_filler,  # pValue
            str(-log10(float(line[7]))) if float(line[7]) != 0 else qvalue_filler,  # qValue
            str(int(line[5]) - int(line[1]))  # source
        ] for line in data]
        return data
    _parse_ranger(converter, tmpfile, saveto, cleanup=True)

    logger.debug("Finished narrow bed report file...")
    return saveto


async def bcp(treatment: [str], control: [str], saveto: str = None, threads: int = 1, pcutoff: float = None,
              fdrcutoff: float = 0.05, format: str = "bam"):
    assert all(os.path.isfile(f) for f in treatment + control)
    assert threads > 0

    # merge treatment and control as a single file - peakranger requirement
    mtreatment, mcontrol = await asyncio.gather(merge(treatment, threads), merge(control, threads))
    tmpfile = os.path.join(tempfile.mkdtemp(), "bcp.bed")

    await _bcp(mtreatment, mcontrol, saveto=tmpfile, pcutoff=pcutoff, format=format)

    # free intermediates
    for f in [mtreatment, mcontrol]:
        if f not in treatment and f not in control:
            os.remove(f)

    logger.debug("Building narrow bed report file...")

    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]

    # We want to transform bcp format into the
    # ENCODE broad peak format https://genome.ucsc.edu/FAQ/FAQformat.html#format13
    # region_chr region_start region_end nearby_genes(6kbp) region_ID region_summits region_fdr region_strand region_treads region_creads
    # ------------->
    # chrom chromStart chromEnd name score strand signalValue pValue qValue
    def converter(data):
        qvalue_filler = str(_zero_significance_filler(data, 6))
        data = []
        for line in data:
            qvalue = float(line[6])
            # fdr control is not implemented in the bcp, do it here
            if qvalue >= fdrcutoff:
                continue
            data.append([
                line[0],  # chrom
                line[1],  # chromStart
                line[2],  # chromEnd
                ".",      # name
                '0',      # score
                '.',      # strand
                '0',      # signalValue
                "-1",     # pValue
                str(-log10(qvalue)) if qvalue != 0 else qvalue_filler,  # qValue
            ])
        return data
    _parse_ranger(converter, tmpfile, saveto, cleanup=True)

    logger.debug("Finished narrow bed report file...")
    return saveto


# PeakRanger -> ENCODE bed format
def _parse_ranger(converter, path: str, saveto: str, cleanup: bool = True):
    with open(path + "_details", 'r') as file:
        data = file.readlines()

    # preprocess and skip meta and header
    data = [line.strip() for line in data]
    begin = [ind for ind, line in enumerate(data) if line == ""]
    data = data[begin[0]+2:]
    data = [line.split('\t') for line in data if line != ""]

    data = converter(data)

    info = ['\t'.join(line) + "\n" for line in data]

    with open(saveto, 'w') as file:
        file.writelines(info)

    if cleanup:
        for postfix in ("_details", "_region.bed",  "_summit.bed"):
            if os.path.exists(path + postfix):
                os.remove(path + postfix)
        os.rmdir(os.path.dirname(path))


# duplicates are filtered by default
# paired end reads are automatically detected(not sure)
async def _ranger(treatment: str, control: str, saveto: str, threads: int,
                  pcutoff: float, fdrcutoff: float, format: str):
    cmd = [
        "peakranger", "ranger", f"--data={treatment}", f"--control={control}", f"--output={saveto}",
        f"--format={format}", f"--FDR={fdrcutoff}", "--verbose", f"--thread={threads}"
    ]
    if pcutoff is not None:
        cmd.append(f"--pval={pcutoff}")
    await run(cmd, logger, logbefore=f"Start peakranger ranger for {treatment} with control {control}",
              logafter="peak calling finished")
    return saveto


# duplicates are filtered by default
# paired end reads are automatically detected(not sure)
async def _bcp(treatment: str, control: str, saveto: str, pcutoff: float, format: str = "bam"):
    cmd = [
        "peakranger", "bcp", f"--data={treatment}", f"--control={control}", f"--output={saveto}",
        f"--format={format}", "--verbose"
    ]
    if pcutoff is not None:
        cmd.append(f"--pval={pcutoff}")
    await run(cmd, logger, logbefore=f"Start peakranger bcp for {treatment} with control {control}",
              logafter="peak calling finished")
    return saveto
