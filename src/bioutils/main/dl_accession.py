import json
import sys
from typing import List

import requests

from commonutils.libdl import download_file
from commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)


def dl_accession(acc: str):
    lh.info(f"Getting EBI REST API for {acc}")
    json_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=fastq_md5,fastq_ftp&format=json&download=true&limit=0"
    request_handler = requests.get(json_url, stream=False)
    try:
        json_content = request_handler.json()
    except json.decoder.JSONDecodeError as e:
        raise ValueError(f"Illegal accession {acc}") from e
    lh.info(f"Getting EBI REST API for {acc} FIN")
    for item in json_content:
        for accession_url, accession_md5 in zip(item['fastq_ftp'].split(';'), item['fastq_md5'].split(";")):
            download_file(accession_url, md5=accession_md5)



def main(args:List[str]):
    for name in args:
        dl_accession(name)

