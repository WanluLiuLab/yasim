import ftplib
import hashlib
import logging
import math
import os
import re
import warnings
from typing import Optional

import requests
import tqdm

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()

warnings.warn("This module is deprecated, so not updated.", DeprecationWarning, stacklevel=2)


def ftp_get_download(
        url: str,
        dest_filename: str,
        username: str = "anonymous",
        password: str = ""
):
    rm = re.match(r"^ftp://([^/]+)/(.*)$", url)
    if not rm:
        raise ValueError(f"Un-parsable FTP address: {url}")

    hostname = rm.group(1)
    filename = rm.group(2)
    logger_handler.info(f"Retrieving {url} -> {dest_filename}: HEADER")
    ftpclient = ftplib.FTP(hostname, username, password)
    try:
        file_size = ftpclient.size(filename)
    except ftplib.error_perm:
        file_size = math.inf
    logger_handler.info(f"Retrieving {url} -> {dest_filename}: HEADER FIN")
    with open(dest_filename, 'wb') as writer:
        with tqdm.tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024) as pbar:
            def write_binary(data: bytes):
                pbar.update(len(data))
                writer.write(data)

            ftpclient.retrbinary(f'RETR {filename}', write_binary)
    ftpclient.close()
    logger_handler.info(f"Retrieving {url} -> {dest_filename} FIN")


def http_get_download(url: str, dest_filename: str):
    logger_handler.info(f"Retrieving {url} -> {dest_filename}: HEADER")
    request_handler = requests.get(url, stream=True)
    file_size = int(request_handler.headers['content-length'])
    chunk_size = 1024
    logger_handler.info(f"Retrieving {url} -> {dest_filename}: HEADER FIN")
    with open(dest_filename, 'wb') as output_file:
        for chunk in tqdm.tqdm(iterable=request_handler.iter_content(chunk_size=chunk_size), total=file_size, unit='B',
                               unit_scale=True, unit_divisor=1024):
            if chunk:
                output_file.write(chunk)
    request_handler.close()
    logger_handler.info(f"Retrieving {url} -> {dest_filename} FIN")


def compute_md5(filename: str):
    logger_handler.info(f"MD5 {filename}")
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    logger_handler.info(f"MD5 {filename} FIN")
    return hash_md5.hexdigest()


def download_file(url: str, dest_filename: Optional[str] = None, md5: Optional[str] = None):
    """
    Download some file from url to dest_filename
    This code uses codes in https://www.simplifiedpython.net/python-download-file/
    """
    global logger_handler
    if dest_filename is None:
        dest_filename = os.path.basename(url)
    logger_handler.info(f"Retrieving {url} -> {dest_filename} START")
    if os.path.exists(dest_filename):
        if md5 is not None:
            actual_md5 = compute_md5(dest_filename)
            if actual_md5 == md5:
                logger_handler.info(f"Retrieving {url} -> {dest_filename} FIN")
                return
        logger_handler.info(f"Existing {dest_filename} will be removed")
        os.remove(dest_filename)
    if url.startswith("http://") or url.startswith("https://"):
        http_get_download(url, dest_filename)
    else:
        ftp_get_download(url, dest_filename)
    if md5 is not None:
        actual_md5 = compute_md5(dest_filename)
        if not actual_md5 == md5:
            logger_handler.warning(f"MD5 mismatch! Actual={actual_md5}, Get={md5}")
    logger_handler.info(f"Retrieving {url} -> {dest_filename} FIN with correct md5")
