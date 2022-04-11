import io
from typing import IO, Union, AnyStr, List, Optional, Iterable, Iterator, Type

from commonutils.io.rules import FileRuleRing
from commonutils.io.std_rules import use_python_std_rules
from commonutils.io.typing import IOProxyType
from commonutils.stdlib_helper.docstring_helper import copy_doc
from commonutils.typing import PathOrFDType, FDType

use_python_std_rules()


class IOProxy(IOProxyType):
    _fd: Union[IO, io.IOBase]
    """
    The underlying file descriptor.
    """

    def __init__(self, fd: FDType, *args, **kwargs):
        self._fd = fd

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
    def __enter__(self) -> IOProxyType:
        try:
            self._fd.__enter__()
        except AttributeError:
            pass
        return self

    @copy_doc(io.RawIOBase.__exit__)
    def __exit__(self,
                 exc_type: Optional[Type[BaseException]],
                 exc_val: Optional[BaseException],
                 exc_tb: Optional[BaseException]):
        try:
            self._fd.__exit__(exc_type, exc_val, exc_tb)
        except AttributeError:
            return


class RuleBasedIOProxy(IOProxy):
    """
    An IO implementation which may read/write archives or texts.

    This IO wrapper will open archives of ``.gz``, ``.xz``, ``.lzma`` or ``.bz2`` suffix,
    with standard IO functions like ``read``,
    with iteration and context management (aka., ``__enter__`` or ``with`` statement) support.
    """

    def __init__(self,
                 path_or_fd: PathOrFDType,
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
            super().__init__(FileRuleRing.open(path_or_fd, *args, **kwargs))
        elif isinstance(path_or_fd, IO) or isinstance(path_or_fd, io.IOBase):
            super().__init__(path_or_fd)
        else:
            raise TypeError(f"Type {type(path_or_fd)} not supported!")


class SequentialReader(RuleBasedIOProxy):
    """
    A sequential reader that is based on :py:class:`RuleBasedIOProxy`
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
