import bz2
import gzip
import lzma
from typing import IO, AnyStr, Iterator, Iterable, Callable, Optional


class ArchiveBaseIO(IO):
    fd: IO

    def __init__(self,
                 filename: str,
                 opener: Optional[Callable[..., IO]] = None,
                 *args,
                 **kwargs):
        if opener is not None:
            self.fd = opener(filename, *args, **kwargs)
        if filename.endswith('.gz'):
            self.fd = gzip.open(filename, *args, **kwargs)
        elif filename.endswith('.xz') or filename.endswith('.lzma'):
            self.fd = lzma.open(filename, *args, **kwargs)
        elif filename.endswith('.bz2'):
            self.fd = bz2.open(filename, *args, **kwargs)
        else:
            self.fd = open(filename, *args, **kwargs)

    @property
    def mode(self) -> str:
        return self.fd.mode

    @property
    def name(self) -> str:
        return self.fd.name

    @property
    def closed(self) -> bool:
        return self.fd.closed

    def close(self) -> None:
        return self.fd.close()

    def fileno(self) -> int:
        return self.fd.fileno()

    def flush(self) -> None:
        self.fd.flush()

    def isatty(self) -> bool:
        return self.fd.isatty()

    def read(self, *args, **kwargs) -> AnyStr:
        return self.fd.read(*args, **kwargs)

    def readable(self) -> bool:
        return self.fd.readable()

    def readline(self, *args, **kwargs) -> AnyStr:
        return self.fd.readline(*args, **kwargs)

    def readlines(self, *args, **kwargs) -> list[AnyStr]:
        return self.fd.readlines(*args, **kwargs)

    def seek(self, *args, **kwargs) -> int:
        return self.fd.seek(*args, **kwargs)

    def seekable(self) -> bool:
        return self.fd.seekable()

    def tell(self) -> int:
        return self.fd.tell()

    def truncate(self, *args, **kwargs) -> int:
        return self.fd.truncate(*args, **kwargs)

    def writable(self) -> bool:
        return self.fd.writable()

    def write(self, s: AnyStr) -> int:
        return self.fd.write(s)

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        self.fd.writelines(lines)

    def __next__(self) -> AnyStr:
        return self.fd.__next__()

    def __iter__(self) -> Iterator[AnyStr]:
        return self.fd.__iter__()

    def __enter__(self) -> IO[AnyStr]:
        try:
            return self.fd.__enter__()
        except AttributeError:
            pass

    def __exit__(self, *args, **kwargs):
        try:
            self.fd.__exit__(*args, **kwargs)
        except AttributeError:
            return


class SequentialReader(ArchiveBaseIO):
    def seek(self, *args, **kwargs) -> int:
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def seekable(self) -> bool:
        return False

    def truncate(self, *args, **kwargs) -> int:
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def writable(self) -> bool:
        return False

    def write(self, *args, **kwargs) -> int:
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def writelines(self, *args, **kwargs) -> None:
        raise OSError("Illegal operation on Sequential Read-Only IO")


def get_reader(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a reader for multiple format.
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)


def get_writer(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a writer for multiple format.
    """
    if is_binary:
        mode = "wb"
    else:
        mode = "wt"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)


def get_appender(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get an appender for multiple format.
    """
    if is_binary:
        mode = "ab"
    else:
        mode = "at"
    return ArchiveBaseIO(filename, mode=mode, **kwargs)
