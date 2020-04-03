import os
import tempfile
import logging
import asyncio
from wrappers.utils import run
from wrappers import sambamba
from wrappers.peak_calling.utils import _zero_significance_filler
from math import log10

logger = logging.getLogger(__name__)


# convert PePr peaks to Broad Bed format
def _pepr_to_broad_bed(folder: str, saveto: str, name="NA", fdrcutoff: float = 0.05, cleanup: bool = True):
    # Nothing to do, just don't forget to recompute p-value and fdr(qvalue) as -log10!
    # chrom chromStart chromEnd name score strand signalValue pValue qValue
    path = os.path.join(folder, f"{name}__PePr_peaks.bed")
    with open(path, 'r') as file:
        data = file.readlines()
    data = [line.strip().split('\t') for line in data]
    pind, qind = 7, 8
    pvalue_filler = _zero_significance_filler(data, pind)
    qvalue_filler = _zero_significance_filler(data, qind)

    for line in data:
        pvalue, qvalue = float(line[pind]), float(line[qind])
        if qvalue >= fdrcutoff:
            continue
        pvalue = -log10(pvalue) if pvalue != 0 else pvalue_filler
        qvalue = -log10(qvalue) if qvalue != 0 else qvalue_filler
        line[pind], line[qind] = str(pvalue), str(qvalue)
    data = ["\t".join(line) + "\n" for line in data]

    with open(saveto, 'w') as file:
        file.writelines(data)

    if cleanup:
        for file in os.listdir(folder):
            os.remove(os.path.join(folder, file))
        os.rmdir(folder)


async def PePr(treatment: [str], control: [str], peaktype: str, saveto: str = None, pcutoff: float = None,
               fdrcutoff: float = 0.05, threads: int = 1, format: str = 'bam'):
    assert len(treatment) >= 2 and (len(treatment) == len(control) or len(control) == 1), \
        "All experiments MUST be matched."
    assert peaktype in ("broad", "sharp")
    assert threads > 0
    outdir = tempfile.mkdtemp()

    if format.lower() == "bampe":
        # PePr makes some assumptions about reads order for paired-end reads - sorting by name is required.
        treatment = await asyncio.gather(*[sambamba.sort(t, threads=threads, byname=True) for t in treatment])
        control = await asyncio.gather(*[sambamba.sort(c, threads=threads, byname=True) for c in control])

    await PePr_(treatment, control, peaktype, outdir, pcutoff, threads, format)

    # PePr postprocessing - skip for now
    _pepr_to_broad_bed(outdir, saveto, cleanup=True, fdrcutoff=fdrcutoff)

    # cleanup if needed
    if format.lower() == 'bampe':
        for f in treatment + control:
            os.remove(f)

    assert not os.path.exists(outdir) and os.path.exists(saveto)
    return saveto


# duplicates: ok, keep everything by default
# paired-end: bampe AND MUST BE SORTED BY NAME
async def PePr_(treatment: [str], control: [str], peaktype: str, saveto: str, pcutoff: float,
                threads: int, format: str):
    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    cmd = [
        "PePr", "-c", ",".join(treatment), "-i", ",".join(control), "-f", format,
        f"--peaktype={peaktype}", f"--num-processors={threads}", f"--output-directory={saveto}"
    ]
    if pcutoff is not None:
        cmd.append(f"--threshold={pcutoff}")
    await run(cmd, logger, logbefore=f"Start PePr for data {treatment} and control {control}, peaks {peaktype}",
              logafter="PePr finished")
    return saveto
