import gzip
import os
import shutil
from typing import IO, Callable

from commonutils.io import get_reader
from commonutils.io.file_system import get_abspath, file_exists, is_soft_link
from commonutils.logger import chronolog, get_logger

lh = get_logger(__name__)


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
    while is_soft_link(path):
        path = get_abspath(os.path.realpath(path))
    return path


@chronolog(display_time=True)
def wc_l(filename: str, opener: Callable[[str], IO] = None) -> int:
    """
    Count lines in a file.

    :param filename: Input filename
    :param opener: Function to open this file. I.e., return an IO object.
    :return: Line number.
    """
    if opener is None:
        fd = get_reader(filename)
    else:
        fd = opener(filename)
    return wc_l_io(fd)


@chronolog(display_time=True)
def wc_c(filename: str, opener: Callable[[str], IO] = None) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param filename: Input filename
    :param opener: Function to open this file. I.e., return an IO object.
    :return: File length.
    """
    if opener is None:
        fd = get_reader(filename)
    else:
        fd = opener(filename)
    return wc_c_io(fd)


@chronolog(display_time=True)
def wc_l_io(fd: IO) -> int:
    """
    Count lines in a file.

    :param fd: An IO object.
    :return: Line number.
    """
    if fd.seekable():
        curr_pos = fd.tell()
    else:
        return -1
    reti = 1
    fd.seek(0)
    while fd.readline():
        reti += 1
    fd.seek(curr_pos)
    return reti


@chronolog(display_time=True)
def wc_c_io(fd: IO) -> int:
    """
    Count the number of chars inside a file, i.e. File length.

    :param fd: An IO object.
    :return: File length.
    """
    if fd.seekable():
        curr_pos = fd.tell()
        fd.seek(0, 2)
        reti = fd.tell()
        fd.seek(curr_pos)
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


def gz_compress(in_file: str, out_file: str, keep_in_file: bool = False, compresslevel: int = 9):
    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb', compresslevel=compresslevel) as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        rm_rf(in_file)


def gz_decompress(in_file: str, out_file: str, keep_in_file: bool = False):
    with gzip.open(in_file, 'rb') as f_in, open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_in_file:
        rm_rf(in_file)
