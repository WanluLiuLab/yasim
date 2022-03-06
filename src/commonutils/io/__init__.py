"""
commonutils.io -- Enhanced Python IO Functions/Classes
======================================================

TODO
"""

import bz2
import gzip
import io
import lzma
import os
from typing import IO, AnyStr, Iterator, Iterable, Callable, Optional, Type, List

from commonutils.stdlib_helper.docstring_helper import copy_doc


def determine_line_endings(fd: IO):
    """
    Determine line endings. If failed, will return OS default.

    """
    FIND_CR = False
    FIND_LF = False
    while True:
        c = fd.read(1)
        if c is None:
            break
        elif c == '\r':
            FIND_CR = True
            if FIND_LF:
                return '\n\r'
        elif c == '\n':
            FIND_LF = True
            if FIND_CR:
                return '\r\n'
        else:
            if FIND_CR:
                return '\r'
            if FIND_LF:
                return '\n'
            FIND_CR = False
            FIND_LF = False
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
                 filename: str,
                 opener: Optional[Callable[..., IO]] = None,
                 *args,
                 **kwargs):
        """
        Open a file.

        :param filename: The filename to be opened.
        :param opener: The underlying opener you wish to use. Should return some IO.
        :param args: Positional arguments passed to underlying opener.
        :param kwargs: Keyword arguments passed to underlying opener.
        TODO: Example using zstd library.
        """
        if opener is not None:
            self._fd = opener(filename, *args, **kwargs)
        if filename.endswith('.gz'):
            self._fd = gzip.open(filename, *args, **kwargs)
        elif filename.endswith('.xz') or filename.endswith('.lzma'):
            self._fd = lzma.open(filename, *args, **kwargs)
        elif filename.endswith('.bz2'):
            self._fd = bz2.open(filename, *args, **kwargs)
        else:
            self._fd = open(filename, *args, **kwargs)

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
    or writting functions like :py:func:`write`
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
