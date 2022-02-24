import gzip
import shutil

import commonutils.io.file_system
import commonutils.shutil


def gz_compress(in_file: str, out_file: str, keep_in_file: bool = False, compresslevel: int = 9):
    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb', compresslevel=compresslevel) as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        commonutils.shutil.rm_rf(in_file)


def gz_decompress(in_file: str, out_file: str, keep_in_file: bool = False):
    with gzip.open(in_file, 'rb') as f_in, open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        commonutils.shutil.rm_rf(in_file)
