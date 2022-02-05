# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: ioctl.py -- A rewritten of important os.path modules and various sysadmin commands.
#
#  VERSION HISTORY:
#  2021-08-10 0.1  : Purposed and added by YU Zhejian, support get_abspath, readlink_f, mkdir_p, touch and rm_rf.
#  2021-08-11 0.1  : is_windows added.
#  2021-08-14 0.1  : Added get_realname (check_input_existence) and open_gz_txt.
#  2021-08-15 0.1  : Added check_output_exists.
#  2021-08-16 0.1  : OS detection upgraded. Add file_exists and is_user_admin.
#  2021-08-19 0.1  : Added is_64bit and is_little_endian.
#  2021-08-23 0.1  : Text merged back. Some system control functions moved to sysctl.py.
#
#  VERSION HISTORY (Text.py):
#  2021-08-16 0.1  : Migrated from ioctl.py, support txt, bgz and gz.
#  2021-08-16 0.1  : Added support for bz2 and xz.
#  2021-08-17 0.1  : Added suffix detection.
#  2021-08-23 0.1  : Added writer function.
# ==============================================================================

"""
ioctl.py -- A rewritten of important os.path modules.

Due to the issues caused by default os.path modules, ww rewrote some of them to make it behave as is used in GNU Core
Utils. Reasons are provided.

This module is intended to work on Microsoft Windows.
"""

import bz2
import gzip
import lzma
import os
import shutil
import stat
from typing import IO, Any

from commonutils import logger
from commonutils.logger import chronolog

lh = logger.get_logger(__name__)

# The following line is the name of the operating systems we will support.
# The following line is the name of operating system feature we will support.

_HAS_BIOPYTHON = True
try:
    import Bio
except ImportError:
    _HAS_BIOPYTHON = False

# Common file suffixes.
common_suffixes = {
    "GTF": ['.gtf.gz', '.gtf', '.gtf.bgz'],
    "GFF": ['.gff.gz', '.gff', '.gff.bgz'],
    "GFF3": ['.gff3.gz', '.gff3', '.gff3.bgz'],
    "BED": ['.bed.gz', '.bed', '.bed.bgz'],
    "RMSK": ['rmsk.txt', 'rmsk.txt.gz'],
    "FASTA": ['.fasta.gz', '.fasta', '.fasta.bgz', '.fa.gz', '.fa', '.fa.bgz'],
}


def get_file_type_from_suffix(filename: str) -> str:
    filename = filename.lower()
    for standard_suffix, real_suffixes in common_suffixes.items():
        for real_suffix in real_suffixes:
            if filename.endswith(real_suffix):
                return standard_suffix
    return 'UNKNOWN'


def _get_opener(filename: str, is_output: bool = False, **kwargs):
    if is_output:
        realname = ensure_output_existence(filename)
        mod = 'wt'
    else:
        realname = ensure_input_existence(filename)
        mod = 'rt'
    if 'encoding' not in kwargs:
        kwargs['encoding'] = 'utf-8'
    if 'newline' not in kwargs:
        kwargs['newline'] = '\n'
    reader = open(filename, mod, **kwargs)  # Buffer=256MiB
    if realname.endswith('.gz'):
        reader = gzip.open(filename, mod, **kwargs)
        # Test biopython
        # if _HAS_BIOPYTHON:
        #     try:
        #         Bio.bgzf.open(realname).close()
        #         reader = Bio.bgzf.open
        #     except Exception:
        #         pass
    elif realname.endswith('.xz') or realname.endswith('.lzma'):
        reader = lzma.open(filename, mod, **kwargs)
    elif realname.endswith('.bz2'):
        reader = bz2.open(filename, mod, **kwargs)
    return reader


def get_reader(filename: str, **kwargs):
    """
    Get a textIO reader for multiple format.
    """
    return _get_opener(filename, is_output=False, **kwargs)


def get_writer(filename: str, **kwargs):
    """
    Get a textIO writer for multiple format.
    """
    return _get_opener(filename, is_output=True, **kwargs)


@chronolog(display_time=True)
def get_abspath(path: str) -> str:
    """
    Get the absolute path of a given relative path. Can be a non-existent file or directory.

    The default one provided by os.path.abspath() or os.path.realpath()
    cannot perform tide expand (That is, expand ~ to your HOME directory).

    :param path: The relative path
    :return: The absolute path
    """
    if path == '':
        return path
    path = os.path.abspath(os.path.expanduser(path))
    return path


@chronolog(display_time=True)
def ensure_input_existence(input_path: str):
    """
    Ensure the existence of an input file.

    :param input_path: Path of the file. Can support special paths like ``/dev/stdin``.
    :return: The absolute path of the file.
    """
    input_path = get_abspath(input_path)
    open(input_path, 'r').close()
    return input_path


@chronolog(display_time=True)
def ensure_output_existence(output_path: str) -> str:
    """
    Ensure the existence of an output file. If it do not exist, try to create one.

    :param output_path: Path of the file. Can support special paths like ``/dev/stdin``.
    :return: The absolute path of the file.
    """
    output_path = get_abspath(output_path)
    touch(output_path)
    open(output_path, 'w').close()
    return output_path


@chronolog(display_time=True)
def readlink_f(path: str) -> str:
    """
    Remove soft links out of the path and return its abstract form, just like what is done by GNU CoreUtils readlink -f.

    Can be used to trace symlink of symlink.

    :param path: Input relative path
    :return: What you get from readlink -f

    FIXME: Errors in this function!
    FIXME: Sometimes you need to use os.readlink()
    """
    path = path.rstrip(os.sep)
    if path == '':
        return path
    is_link = stat.S_ISLNK(os.stat(path, follow_symlinks=False)[0])
    while is_link:
        path = get_abspath(os.path.realpath(path))
        is_link = stat.S_ISLNK(os.stat(path, follow_symlinks=False)[0])
    return path


@chronolog(display_time=True)
def wc_l(path: str, opener:Any=None) -> int:
    """
    Count lines in a file.

    :param path: File path.
    :return: Line number.
    """
    if opener is None:
        fd = get_reader(path)
    else:
        fd=opener(path)
    return wc_l_io(fd)


@chronolog(display_time=True)
def wc_c(path: str, opener:Any=None) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param path: File to open.
    :return: File length.
    """
    if opener is None:
        fd = get_reader(path)
    else:
        fd=opener(path)
    return wc_c_io(fd)

@chronolog(display_time=True)
def wc_l_io(underlying_fd:IO) -> int:
    """
    Count lines in a file.

    :param path: File path.
    :return: Line number.
    """
    if underlying_fd.seekable():
        curr_pos = underlying_fd.tell()
    else:
        return -1
    reti = 0
    underlying_fd.seek(0)
    while underlying_fd.readline():
        reti += 1
    underlying_fd.seek(curr_pos)
    return reti



@chronolog(display_time=True)
def wc_c_io(underlying_fd:IO) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param path: File to open.
    :return: File length.
    """
    if underlying_fd.seekable():
        curr_pos = underlying_fd.tell()
        underlying_fd.seek(0, 2)
        reti = underlying_fd.tell()
        underlying_fd.seek(curr_pos)
        return reti
    else:
        return -1


@chronolog(display_time=True)
def mkdir_p(path: str) -> str:
    """
    Protected mkdir, as is done by mkdir in GNU CoreUtils. It will:

    1) create the directory tree if not exist, and
    2) not complain about existing target directory.

    :param path: The directory you wish to create.
    :return: Absolute path of directory created.
    """
    abspath = get_abspath(path)
    os.makedirs(path, exist_ok=True)
    return abspath


def file_exists(path: str, allow_special_paths: bool = True) -> bool:
    """
    Check the existence of a file. Can check regular and special files.

    :param path: The path you wish to examine.
    :param allow_special_paths: Whether this is a special file, may be block device, etc.
    """
    if allow_special_paths:
        return os.path.exists(path) and not os.path.isdir(path)
    else:
        return os.path.isfile(path)


@chronolog(display_time=True)
def touch(filename: str):
    """
    touch: ensure the existence of a file, just like GNU CoreUtils touch.
    Please note that this version of touch CANNOT change file attributes like times.

    :param filename: The filename you wish to touch
    """
    if file_exists(filename):
        lh.debug(f"File '{filename}' exists")
    elif os.path.isdir(filename):
        raise IsADirectoryError(f"File '{filename}' is a directory")
    else:
        mkdir_p(os.path.dirname(filename))
        open(filename, mode="a").close()


@chronolog(display_time=True)
def rm_rf(path: str):
    """
    Remove path recursively from the filesystem, just like rm -rf
    Should not complain on non-existing files.

    :param path: The path you wish to remove
    """
    dbg_head = "rm(path='" + path + "')"
    if os.path.isdir(path) and not os.path.islink(path):
        lh.debug(f"{dbg_head} is a directory")
        shutil.rmtree(path)
    elif os.path.exists(path):
        lh.debug(f"{dbg_head} is a file")
        os.remove(path)
    else:
        lh.debug(f"{dbg_head} not exist")
    return 0


def write_line(f, *args, delimiter='\t', newline='\n'):
    f.write(delimiter.join(str(arg) for arg in args) + newline)
