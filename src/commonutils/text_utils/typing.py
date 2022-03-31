import io
import threading
import time
from abc import ABC, abstractmethod
from typing import IO, Iterator, Optional, Type, List, Union, TextIO, Callable

from commonutils.io.safe_io import get_reader
from commonutils.stdlib_helper.docstring_helper import copy_doc
from commonutils.typing import PathOrFDType


class BaseStreamType(TextIO, ABC):
    """
    Restricted IO stream.
    """

    @abstractmethod
    def __repr__(self):
        pass

    @copy_doc(io.RawIOBase.__exit__)
    @abstractmethod
    def __exit__(self,
                 exc_type: Optional[Type[BaseException]],
                 exc_val: Optional[BaseException],
                 exc_tb: Optional[BaseException]):
        pass

    @property
    def mode(self) -> str:
        return "rt"

    @property
    def name(self) -> str:
        return repr(self)

    def fileno(self) -> int:
        """Disabled; will get io error"""
        raise io.UnsupportedOperation("fileno")

    def flush(self) -> None:
        """Function disabled."""
        pass

    def isatty(self) -> bool:
        """False"""
        return False

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def seekable(self) -> bool:
        """False"""
        return False

    def tell(self) -> int:
        """Disabled; will be set to -1"""
        return -1

    def truncate(self, *args, **kwargs) -> int:
        """Disabled; will be set to -1"""
        return -1

    def writable(self) -> bool:
        """False"""
        return False

    def write(self, *args, **kwargs) -> int:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    def writelines(self, *args, **kwargs) -> None:
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError("Illegal operation on Sequential Read-Only IO")

    @property
    @abstractmethod
    def closed(self) -> bool:
        pass

    @copy_doc(io.RawIOBase.close)
    @abstractmethod
    def close(self) -> None:
        pass

    @copy_doc(io.RawIOBase.read)
    @abstractmethod
    def read(self, size: int = -1) -> Optional[str]:
        pass

    @copy_doc(io.RawIOBase.readable)
    @abstractmethod
    def readable(self) -> bool:
        pass

    @copy_doc(io.RawIOBase.readline)
    @abstractmethod
    def readline(self, limit: int = -1) -> Optional[str]:
        pass

    @copy_doc(io.RawIOBase.readlines)
    @abstractmethod
    def readlines(self, hint: int = -1) -> Optional[List[str]]:
        pass

    @copy_doc(io.RawIOBase.__next__)
    @abstractmethod
    def __next__(self) -> str:
        pass

    @copy_doc(io.RawIOBase.__iter__)
    @abstractmethod
    def __iter__(self) -> Iterator[str]:
        pass

    @copy_doc(io.RawIOBase.__enter__)
    @abstractmethod
    def __enter__(self) -> IO[str]:
        pass


class BaseStreamProxy(BaseStreamType):
    _fd: Union[IO, io.IOBase]
    """
    The underlying file descriptor.
    """

    def __init__(self, path_or_fd: PathOrFDType):
        """
        Open a file.

        :param path_or_fd: The filename to be opened.
        """
        if isinstance(path_or_fd, str):
            self._fd = get_reader(path_or_fd)
        elif isinstance(path_or_fd, IO) or isinstance(path_or_fd, io.IOBase):
            self._fd = path_or_fd
        else:
            raise TypeError(f"Type {type(path_or_fd)} not supported!")

    def __repr__(self):
        return f"BaseStreamProxy of {self._fd}"

    @property
    def closed(self) -> bool:
        return self._fd.closed

    @copy_doc(io.RawIOBase.close)
    def close(self) -> None:
        return self._fd.close()

    @copy_doc(io.RawIOBase.read)
    def read(self, size: int = -1) -> str:
        return self._fd.read(size)

    @copy_doc(io.RawIOBase.readable)
    def readable(self) -> bool:
        return self._fd.readable()

    @copy_doc(io.RawIOBase.readline)
    def readline(self, limit: int = -1) -> str:
        return self._fd.readline(limit)

    @copy_doc(io.RawIOBase.readlines)
    def readlines(self, hint: int = -1) -> List[str]:
        return self._fd.readlines(hint)

    @copy_doc(io.RawIOBase.__next__)
    def __next__(self) -> str:
        return self._fd.__next__()

    @copy_doc(io.RawIOBase.__iter__)
    def __iter__(self) -> Iterator[str]:
        return self._fd.__iter__()

    @copy_doc(io.RawIOBase.__enter__)
    def __enter__(self) -> BaseStreamType[str]:
        try:
            self._fd.__enter__()
        except AttributeError:
            pass
        return self

    def __exit__(self,
                 exc_type: Optional[Type[BaseException]],
                 exc_val: Optional[BaseException],
                 exc_tb: Optional[BaseException]):
        try:
            self._fd.__exit__(exc_type, exc_val, exc_tb)
        except AttributeError:
            return


class BufferController(threading.Thread):
    buffer_size: int
    is_terminated: bool
    buffer: List[str]
    interval: float
    readline_from_input_hook: Callable[[], Optional[str]]

    def __init__(
            self,
            readline_from_input_hook: Callable[[], Optional[str]],
            buffer_size: int = 10,
            interval: float = 0.01
    ):
        super().__init__()
        self.is_terminated = False
        self.readline_from_input_hook = readline_from_input_hook
        self.buffer_size = buffer_size
        self.interval = interval

    def read(self, n: int) -> Optional[str]:
        if self.is_terminated:
            return None
        rets = ""
        while not self.is_terminated and len(rets) < n:
            if len(self.buffer) != 0:
                retsb = self.buffer.pop(0)
                if len(retsb) + len(rets) == n:
                    return rets + retsb
                elif len(retsb) + len(rets) > n:
                    ldiff = len(retsb) - len(rets)
                    self.buffer.insert(0, retsb[ldiff:])
                    return rets + retsb[:ldiff]
            time.sleep(self.interval)
        return None

    def readline(self, n: int) -> Optional[str]:
        if self.is_terminated:
            return None
        while not self.is_terminated:
            if len(self.buffer) != 0:
                rets = self.buffer.pop(0)
                if len(rets) > n:
                    self.buffer.insert(0, rets[n:])
                    return rets[:n]
            time.sleep(self.interval)
        return None

    def run(self) -> None:
        """
        Fill the buffer.

        Will set ``self.is_terminated`` after input stream gets terminated.
        """
        while not self.is_terminated:
            if len(self.buffer) < self.buffer_size:
                get_line = self.readline_from_input_hook()
                if get_line is None:
                    self.is_terminated = True
                else:
                    self.buffer.append(get_line)
            time.sleep(self.interval)


class BaseLineBufferedStream(BaseStreamType, ABC):
    _at_eof: bool
    _input_buffer_controller: BufferController

    def _setup_buffer(self):
        self._at_eof = False
        self._input_buffer_controller = BufferController(self._read_into_input_buffer)

    def _start_buffer_controller(self):
        self._input_buffer_controller.start()

    def _stop_buffer_controller(self):
        self._input_buffer_controller.is_terminated = True
        self._input_buffer_controller.join()

    @abstractmethod
    def _read_into_input_buffer(self) -> Optional[str]:
        pass

    def read(self, n: int = -1) -> Optional[str]:
        # Get some chars from buffer
        rets = self._input_buffer_controller.read(n)
        if rets is None:
            self._stop_buffer_controller(n)
        else:
            return rets

    def readline(self, n: int = -1) -> Optional[str]:
        # Get one line from buffer
        rets = self._input_buffer_controller.readline()
        if rets is None:
            self._stop_buffer_controller()
        else:
            return rets

    def readlines(self, hint: int = -1) -> List[str]:
        if self._at_eof:
            return []
        elif hint is None or hint <= 0:
            return [line for line in self]

    def __iter__(self) -> Iterator[str]:
        return self

    def __next__(self) -> str:
        rets = self.readline()
        if rets is None:
            raise StopIteration
        return rets
