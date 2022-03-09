"""
commonutils.io -- Enhanced Python IO Functions/Classes
======================================================

This module includes some enhanced IO functions,
like IO that automatically creates missing intermediate directories or IO with tqdm progress bar.
It also supports Python standard archive types like ``bz2``, ``gzip`` or ``lzma``.

For those who wish to get a ``gzip`` command-like utilities,
please visit :py:mod:`commonutils.shell_utils` for more information
"""

import bz2
import gzip
import io
import lzma
import os
from typing import IO, AnyStr, Iterator, Iterable, Callable, Optional, Type, List, Union

from commonutils.stdlib_helper.docstring_helper import copy_doc

__all__ = (
    "determine_line_endings",
    "determine_file_line_endings",
    "ArchiveBaseIO",
    "SequentialReader",
    "get_reader",
    "get_writer",
    "get_appender"
)

from commonutils.typing import PathOrFDType


def determine_line_endings(fd: IO):
    """
    Determine line endings. If failed, will return OS default.

    """
    find_cr = False
    find_lf = False
    while True:
        c = fd.read(1)
        if c is None:
            break
        elif c == '\r':
            find_cr = True
            if find_lf:
                return '\n\r'
        elif c == '\n':
            find_lf = True
            if find_cr:
                return '\r\n'
        else:
            if find_cr:
                return '\r'
            if find_lf:
                return '\n'
            find_cr = False
            find_lf = False
    return os.linesep


def determine_file_line_endings(filename: str):
    return determine_line_endings(ArchiveBaseIO(filename))


class ArchiveBaseIO(IO):
    """
    An IO implementation which may read/write archives or texts.

    This IO wrapper will open archives of ``.gz``, ``.xz``, ``.lzma`` or ``.bz2`` suffix,
    with standard IO functions like ``read``,
    with iteration and context management (aka., ``__enter__`` or ``with`` statement) support.
    """

    _fd: IO
    """
    The underlying file descriptor.
    """

    def __init__(self,
                 path_or_fd:PathOrFDType,
                 *args,
                 **kwargs):
        """
        Open a file.

        :param path_or_fd: The filename to be opened.
        :param args: Positional arguments passed to underlying opener.
        :param kwargs: Keyword arguments passed to underlying opener.
        TODO: Example using zstd library.
        """
        if isinstance(path_or_fd, str):
            if path_or_fd.endswith('.gz'):
                self._fd = gzip.open(path_or_fd, *args, **kwargs)
            elif path_or_fd.endswith('.xz') or path_or_fd.endswith('.lzma'):
                self._fd = lzma.open(path_or_fd, *args, **kwargs)
            elif path_or_fd.endswith('.bz2'):
                self._fd = bz2.open(path_or_fd, *args, **kwargs)
            else:
                self._fd = open(path_or_fd, *args, **kwargs)
        elif isinstance(path_or_fd, IO) or isinstance(path_or_fd, io.IOBase):
            self._fd = path_or_fd
        else:
            raise TypeError(f"Type {type(path_or_fd)} not supported!")



    @property
    def mode(self) -> str:
        return self._fd.mode

    @property
    def name(self) -> str:
        return self._fd.name

    @property
    def closed(self) -> bool:
        return self._fd.closed

    @copy_doc(io.RawIOBase.close)
    def close(self) -> None:
        return self._fd.close()

    @copy_doc(io.RawIOBase.fileno)
    def fileno(self) -> int:
        return self._fd.fileno()

    @copy_doc(io.RawIOBase.flush)
    def flush(self) -> None:
        self._fd.flush()

    @copy_doc(io.RawIOBase.isatty)
    def isatty(self) -> bool:
        return self._fd.isatty()

    @copy_doc(io.RawIOBase.read)
    def read(self, size: int = -1) -> AnyStr:
        return self._fd.read(size)

    @copy_doc(io.RawIOBase.readable)
    def readable(self) -> bool:
        return self._fd.readable()

    @copy_doc(io.RawIOBase.readline)
    def readline(self, limit: int = -1) -> AnyStr:
        return self._fd.readline(limit)

    @copy_doc(io.RawIOBase.readlines)
    def readlines(self, hint: int = -1) -> List[AnyStr]:
        return self._fd.readlines(hint)

    @copy_doc(io.RawIOBase.seek)
    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        return self._fd.seek(offset, whence)

    @copy_doc(io.RawIOBase.seekable)
    def seekable(self) -> bool:
        return self._fd.seekable()

    @copy_doc(io.RawIOBase.tell)
    def tell(self) -> int:
        return self._fd.tell()

    @copy_doc(io.RawIOBase.truncate)
    def truncate(self, size: Optional[int]) -> int:
        return self._fd.truncate(size)

    @copy_doc(io.RawIOBase.writable)
    def writable(self) -> bool:
        return self._fd.writable()

    @copy_doc(io.RawIOBase.write)
    def write(self, s: AnyStr) -> int:
        return self._fd.write(s)

    @copy_doc(io.RawIOBase.writelines)
    def writelines(self, lines: Iterable[AnyStr]) -> None:
        self._fd.writelines(lines)

    @copy_doc(io.RawIOBase.__next__)
    def __next__(self) -> AnyStr:
        return self._fd.__next__()

    @copy_doc(io.RawIOBase.__iter__)
    def __iter__(self) -> Iterator[AnyStr]:
        return self._fd.__iter__()

    @copy_doc(io.RawIOBase.__enter__)
    def __enter__(self) -> IO[AnyStr]:
        try:
            return self._fd.__enter__()
        except AttributeError:
            pass

    @copy_doc(io.RawIOBase.__exit__)
    def __exit__(self,
                 exc_type: Optional[Type[BaseException]],
                 exc_val: Optional[BaseException],
                 exc_tb: Optional[BaseException]):
        try:
            self._fd.__exit__(exc_type, exc_val, exc_tb)
        except AttributeError:
            return


class SequentialReader(ArchiveBaseIO):
    """
    A sequential reader that is based on :py:class:`ArchiveBaseIO`
    which does not support random-access functions like :py:func:`seek`
    or writing functions like :py:func:`write`
    """

    def seek(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def seekable(self) -> bool:
        """False"""
        return False

    def truncate(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def writable(self) -> bool:
        """False"""
        return False

    def write(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def writelines(self, *args, **kwargs) -> None:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")


def get_reader(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a reader for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.

    :param filename: Filename to be opened.
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.

    .. warning::
        Do NOT specify ``mode`` keyword arguments!
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)


def get_writer(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a writer for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.
    
    :param filename: Filename to be opened.
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.
    
    :: warning..
        Do NOT specify ``mode`` keyword arguments!
    """
    if is_binary:
        mode = "wb"
    else:
        mode = "wt"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)


def get_appender(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get an appender for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.
    
    :param filename: Filename to be opened.
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.
    
    :: warning..
        Do NOT specify ``mode`` keyword arguments!
    """
    if is_binary:
        mode = "ab"
    else:
        mode = "at"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)
