import os
import tempfile
import logging
import requests
import subprocess
from subprocess import run

logger = logging.getLogger(__name__)
ENCODE_URL = "https://www.encodeproject.org/"


def get_unfiltered_bam(accession: str, saveto: str = None, assembly: str = 'hg19'):
    url = f"{ENCODE_URL}/files/{accession}/?format=json"
    response = requests.get(url, headers={'accept': 'application/json'}).json()
    assert response['accession'] == accession
    assert response['file_type'] == 'bam'
    assert "chip-seq" in response['assay_term_name'].lower()
    assert response['status'] == 'released'
    assert response['assembly'] == assembly
    assert response['output_type'] == 'unfiltered alignments'
    path = wget(accession, f'{accession}.bam', saveto=saveto)
    assert os.path.exists(path) and (saveto is None or saveto == path)
    result = run(["md5sum", path], check=True,
                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = result.stdout.decode()

    md5sum_download = result.split('  ')[0]
    assert response['md5sum'] == md5sum_download, "md5 sum wasn't matched"
    return path


def wget(accession: str, file: str, saveto: str = None):
    url = f"{ENCODE_URL}/{accession}/@@download/{file}"
    file = tempfile.mkstemp()[1] if saveto is None else saveto

    result = run(["wget", "--continue", "--retry-connrefused", "--tries=0", "--timeout=5", "-O", file, url],
                 check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logger.debug(result.stdout.decode())
    return file

# logging.basicConfig(level=logging.NOTSET)
#
# a = download("ENCFF534UUQ", "ENCFF534UUQ.bed.gz", "b3ecbb186d2712af0541a2605259457c", "/tmp/tmp.bed")
