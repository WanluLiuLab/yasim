import io
from types import TracebackType
from typing import Type, Iterator, AnyStr, TextIO, Optional, List, Union, Iterable

from commonutils.io.safe_io import get_reader
from commonutils.text_utils.typing import BaseStream
from commonutils.typing import FDType, PathOrFDType



class Cat(BaseStream):

    _fds:List[FDType]

    _closed:bool

    @property
    def closed(self) -> bool:
        return self._closed

    def __init__(self, *pathnames_or_fds:PathOrFDType):
        self._closed=False
        for pathname_or_fd in pathnames_or_fds:
            self._fds.append(get_reader(pathname_or_fd)) # FIXME: Refactor



    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[BaseException]):
        pass

    def __enter__(self) -> TextIO:
        return self

    def close(self) -> None:
        if self._closed:
            return
        for fd in self._fds:
            try:
                fd.close()
            except (OSError, io.UnsupportedOperation):
                pass
        self._closed = True

    def read(self, __n: int = ...) -> AnyStr:
        pass

    def readable(self) -> bool:
        pass

    def readline(self, __limit: int = ...) -> AnyStr:
        pass

    def readlines(self, __hint: int = ...) -> List[AnyStr]:
        pass

    def __next__(self) -> AnyStr:
        pass

    def __iter__(self) -> Iterator[AnyStr]:
        pass



    def __repr__(self):
        pass
