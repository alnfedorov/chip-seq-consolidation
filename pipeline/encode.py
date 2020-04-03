import os
import asyncio
import requests
from wrappers import wget
from wrappers.utils import run

ENCODE_URL = "https://www.encodeproject.org/"


async def rget(*args, **kwargs):
    """async requests-get"""
    return await asyncio.get_event_loop().run_in_executor(None, lambda: requests.get(*args, **kwargs))


async def fetch(accession: str, experiment: str, is_control: bool, assembly: str, saveto: str):
    """
    Fetch unfiltered bam from the ENCODE data base
    :param accession: ENCODE accession number for the file
    :param experiment: ENCODE experiment accession
    :param is_control: control or not
    :param assembly: human genome assembly version
    :param saveto: file to save donwloaded bam file
    :return: saveto - path to the file
    """
    url = f"{ENCODE_URL}/files/{accession}/?frame=object&format=json&limit=all"
    response = (await rget(url, headers={'accept': 'application/json'})).json()
    assert response['accession'] == accession
    assert response['file_type'] == 'bam'
    assert "chip-seq" in response['assay_term_name'].lower()
    assert response['status'] == 'released'
    assert response['assembly'] == assembly
    assert response['output_type'] == 'unfiltered alignments'
    fileid = response['@id']
    md5sum = response['md5sum']

    # check accession is among experiment's files, or it's control files
    url = f"{ENCODE_URL}/experiments/{experiment}/?frame=object&format=json&limit=all"
    response = (await rget(url, headers={'accept': 'application/json'})).json()
    if not is_control:
        is_valid = fileid in response['files']
    else:
        is_valid = False
        for control in response['possible_controls']:
            control = control.split('/')[-2]
            url = f"{ENCODE_URL}/experiments/{control}/?frame=object&format=json&limit=all"
            response = (await rget(url, headers={'accept': 'application/json'})).json()
            if fileid in response['files']:
                is_valid = True
                break
    if not is_valid:
        raise ValueError(f"File {accession} doesn't belong to the experiment {experiment} files.")

    url = f"{ENCODE_URL}/{accession}/@@download/{accession}.bam"
    path = await wget.wget(url, saveto=saveto)
    assert os.path.exists(path) and (saveto is None or saveto == path)

    md5sum_downloaded = await run(["md5sum", path])
    md5sum_downloaded = md5sum_downloaded.split('  ')[0]
    assert md5sum == md5sum_downloaded, "md5 sum wasn't matched"
    return path

__all__ = [fetch]
