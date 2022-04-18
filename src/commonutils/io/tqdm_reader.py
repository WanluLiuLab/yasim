"""
tqdm_reader.py -- Reader with Progress Bar

Here are wrappings for basic IO classes & functions in :py:mod:`commonutils.io` with additional progress bar.
"""

from typing import Iterator, AnyStr, List, Type

from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io import get_reader
from commonutils.io._ioproxy import RuleBasedIOProxy, SequentialReader
from commonutils.stdlib_helper.docstring_helper import copy_doc

__all__ = (
    "TqdmReader",
    "TqdmLineReader",
    "get_tqdm_reader",
    "get_tqdm_line_reader"
)


class _BaseTqdmReader(SequentialReader):
    _tqdm: Type[tqdm]

    @copy_doc(RuleBasedIOProxy.__enter__)
    def __enter__(self, *args, **kwargs):
        self._tqdm.__enter__(*args, **kwargs)
        return self

    @copy_doc(RuleBasedIOProxy.__exit__)
    def __exit__(self, *args, **kwargs):
        return self._tqdm.__exit__(*args, **kwargs)


class TqdmReader(_BaseTqdmReader):
    """
    A very simple tqdm reader with :py:func:``read`` and :py:func:``readline`` functions.
    """

    @copy_doc(RuleBasedIOProxy.__init__)
    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self._tqdm = tqdm(
            desc=f"Reading {filename}", total=shell_utils.wc_c_io(self._fd), unit='B', unit_scale=True,
            unit_divisor=1024
        )

    @copy_doc(RuleBasedIOProxy.readline)
    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(*args, **kwargs)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    @copy_doc(RuleBasedIOProxy.read)
    def read(self, size: int = -1) -> AnyStr:
        update_bytes = super().read(size)
        self._tqdm.update(len(update_bytes))
        return update_bytes

    @copy_doc(RuleBasedIOProxy.readlines)
    def readlines(self, *args, **kwargs) -> List[AnyStr]:
        update_bytes_arr = super().readlines(*args, **kwargs)
        self._tqdm.update(sum(map(len, update_bytes_arr)))
        return update_bytes_arr


class TqdmLineReader(_BaseTqdmReader):
    """
    A very simple tqdm reader with only :py:func:`readline` functions.

    .. warning:: This class has different ``__iter__`` method!
    """

    @copy_doc(RuleBasedIOProxy.__init__)
    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self._tqdm = tqdm(desc=f"Reading {filename}", total=shell_utils.wc_l_io(self._fd), unit='L')

    def read(self, *args, **kwargs):
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError('Illegal operation read.')

    def readlines(self, *args, **kwargs):
        """This function is disabled, will raise :py:class:`OSError`"""
        raise OSError('Illegal operation read.')

    @copy_doc(RuleBasedIOProxy.readline)
    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(-1)  # Size limit canceled.
        self._tqdm.update(1)
        return update_bytes

    def __iter__(self) -> Iterator[str]:
        """Iterator over lines with removal of trailing line feed (``LF``, ``\n``)/carriage return (``CR``, ``\r``)"""
        while True:
            line = self.readline()
            if not line:
                break
            yield line.rstrip('\n\r')


@copy_doc(get_reader)
def get_tqdm_reader(filename: str, is_binary: bool = False, **kwargs) -> TqdmReader:
    """
    Get a reader for multiple format.
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmReader(filename, mode=mode, **kwargs)


def get_tqdm_line_reader(filename: str, is_binary: bool = False, **kwargs) -> TqdmLineReader:
    """
    :py:func:`get_reader`-like wrapper for :py:class:`TqdmLineReader`
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmLineReader(filename, mode=mode, **kwargs)
