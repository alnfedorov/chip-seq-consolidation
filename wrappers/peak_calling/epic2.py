import os
import logging
import tempfile
import asyncio
from math import log10
from wrappers.peak_calling.macs2 import predictd
from wrappers.utils import run
from wrappers.peak_calling.utils import _zero_significance_filler

logger = logging.getLogger(__name__)


async def epic2(treatment: [str], control: [str], ispaired: bool, binsize: int = 200, gapsallowed: int = 3,
                fragment_size: int = None, saveto: str = None, fdrcutoff: float = 0.05, genome: str = 'hg19'):
    if ispaired:
        assert fragment_size is None
    elif fragment_size is None:
        # estimate as mean fragment size across libraries
        size = [predictd(file) for file in treatment]
        size = await asyncio.gather(*size)
        logger.debug(f"SICER(epic2) estimated fragment sizes for chip-seq files {size}")
        fragment_size = int(round(sum(size) / len(size)))
        assert fragment_size > 20

    output = tempfile.mkstemp()[1]
    await epic2_(treatment, control, binsize, gapsallowed, fragment_size, output, fdrcutoff, genome)

    # Convert to the broad bed format
    # Chromosome    Start	End	PValue	Score	Strand	ChIPCount	InputCount	FDR	log2FoldChange
    # ------------->
    # chrom chromStart chromEnd name score strand signalValue pValue qValue
    logger.debug("Started format conversion epic2 -> broad bed")
    with open(output, 'r') as file:
        data = file.readlines()
    os.remove(output)

    assert data[0].startswith("#")
    data = [line.strip().split("\t") for line in data[1:]]
    pvalue_filler = _zero_significance_filler(data, 3)
    qvalue_filler = _zero_significance_filler(data, 8)

    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    with open(saveto, 'w') as file:
        for line in data:
            pvalue, qvalue = float(line[3]), float(line[8])
            pvalue, qvalue = -log10(pvalue) if pvalue != 0 else pvalue_filler, \
                             -log10(qvalue) if qvalue != 0 else qvalue_filler
            string = "\t".join([
                line[0], line[1], line[2],  # chrom start end
                ".",    # name
                str(min(1000, max(0, float(line[4])))),      # score (0, 1000)
                ".",    # strand
                str(float(line[-1])),  # signal value filler
                str(pvalue),            # pValue
                str(qvalue)             # qValue
            ]) + "\n"
            file.write(string)
    logger.debug("Finished format conversion")
    return saveto


# duplicates -> auto filter to 1
# paired-end -> pass bedpe, no other way.
async def epic2_(treatment: [str], control: [str], binsize: int, gapsallowed: int, fragment_size: int, saveto: str,
                 fdrcutoff: float, genome: str):
    cmd = [
        "epic2", "-t", *treatment, "-c", *control, f"--genome={genome}", f"--false-discovery-rate-cutoff={fdrcutoff}",
        f"--output={saveto}",
    ]
    if binsize:
        cmd.append(f"--bin-size={binsize}")
    if gapsallowed:
        cmd.append(f"--gaps-allowed={gapsallowed}")
    if fragment_size:
        cmd.append(f"--fragment-size={fragment_size}")
    await run(cmd, logger, logbefore=f"Starting SICER(epic2) with cmd {' '.join(cmd)}", logafter="Finished SICER")
    return saveto
