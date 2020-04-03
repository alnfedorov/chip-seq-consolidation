import os
import logging
import tempfile
from wrappers.utils import run

logger = logging.getLogger(__name__)


async def predictd(file: str, rfile: str = "/dev/null", gsize: str = "hs"):
    cmd = ["macs2", "predictd", f"--gsize={gsize}", f"--rfile={rfile}", "-i", file]
    results = await run(
        cmd, logger, logbefore=f"Starting macs2 predictd for {file}", logafter=f"Finished macs2 predictd"
    )

    # Search in the output for the # predicted fragment length is 215 bps
    results = results.split('\n')
    fragment = [line.strip() for line in results if "# predicted fragment length is" in line]
    assert len(fragment) == 1, "macs2 predictd changed its output format"
    fragment = fragment[0].split(" ")
    assert fragment[-1] == "bps"

    fragment = int(fragment[-2])
    assert fragment > 0
    return fragment


# duplicates -> OK
# paired-end -> 5` end or specify format BAMPE
async def callpeak(data: [str], control: [str], pcutoff: float = None, fdrcutoff: float = None, isbroad: bool = False,
                   saveto: str = None, gsize: str = "hs", keepdup: str = 'auto', format: str = "AUTO"):
    assert pcutoff is not None or fdrcutoff is not None
    saveto = saveto if saveto is not None else tempfile.mkstemp()[1]
    output = tempfile.mkdtemp(dir=os.path.abspath(os.path.dirname(saveto)))
    cmd = [
        "macs2", "callpeak", "-t", *data, "-c", *control, "-g", gsize, f"--outdir={output}",
        "--seed=123", f"--keep-dup={keepdup}", f"--format={format}"
    ]
    if isbroad:
        cmd.append("--broad")
    if pcutoff is not None:
        cmd += ["-p", str(pcutoff)]
    else:
        cmd += ["-q", str(fdrcutoff)]
    await run(cmd, logger, logbefore=f"Starting macs2 callpeack with cmd {' '.join(cmd)}", logafter="Finished macs2")

    if isbroad:
        os.renames(os.path.join(output, "NA_peaks.broadPeak"), saveto)
    else:
        os.renames(os.path.join(output, "NA_peaks.narrowPeak"), saveto)

    # cleanup
    for file in os.listdir(output):
        os.remove(os.path.join(output, file))
    os.removedirs(output)

    return saveto
