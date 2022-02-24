import json
import logging
import sys

import requests

from commonutils.libdl import download_file

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()


def dl_accession(acc: str):
    logger_handler.info(f"Getting EBI REST API for {acc}")
    json_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=fastq_md5,fastq_ftp&format=json&download=true&limit=0"
    request_handler = requests.get(json_url, stream=False)
    try:
        json_content = request_handler.json()
    except json.decoder.JSONDecodeError as e:
        raise ValueError(f"Illegal accession {acc}") from e
    logger_handler.info(f"Getting EBI REST API for {acc} FIN")
    for item in json_content:
        accession_url = "ftp://" + item['fastq_ftp']
        accession_md5 = item['fastq_md5']
        download_file(accession_url, md5=accession_md5)


def main():
    dl_accession(sys.argv[1])


if __name__ == "__main__":
    main()
