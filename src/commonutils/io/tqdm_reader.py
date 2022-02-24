from typing import IO, Iterator, AnyStr, List

from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io import SequentialReader


class TqdmReader(SequentialReader):
    """
    A very simple tqdm reader with :py:func:``read`` and :py:func:``readline`` functions.
    """

    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self.tqdm = tqdm(
            desc=f"Reading {filename}", total=shell_utils.wc_c_io(self.fd), unit='B', unit_scale=True,
            unit_divisor=1024
        )

    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(*args, **kwargs)
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def read(self, size: int = -1) -> AnyStr:
        update_bytes = super().read(size)
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def readlines(self, *args, **kwargs) -> List[AnyStr]:
        update_bytes_arr = super().readlines(*args, **kwargs)
        self.tqdm.update(sum(map(len, update_bytes_arr)))
        return update_bytes_arr

    def __enter__(self):
        self.tqdm.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        return self.tqdm.__exit__(*args, **kwargs)


def get_tqdm_reader(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a reader for multiple format.
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmReader(filename, mode=mode, **kwargs)


class TqdmLineReader(TqdmReader):
    """
    A very simple tqdm reader with only :py:func:``readline`` functions.
    """

    def __init__(self, filename: str, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)
        self.tqdm = tqdm(desc=f"Reading {filename}", total=shell_utils.wc_l_io(self.fd), unit='L')

    def read(self, *args, **kwargs):
        raise OSError('Illegal operation read.')

    def readlines(self, *args, **kwargs):
        update_bytes_arr = super().readlines(*args, **kwargs)
        self.tqdm.update(len(update_bytes_arr))
        return update_bytes_arr

    def readline(self, *args, **kwargs) -> AnyStr:
        update_bytes = super().readline(*args, **kwargs)
        self.tqdm.update(1)
        return update_bytes

    def __iter__(self) -> Iterator[str]:
        while True:
            line = self.readline()
            if not line:
                break
            yield line.rstrip('\n\r')


def get_tqdm_line_reader(filename: str, is_binary: bool = False, **kwargs) -> IO:
    """
    Get a reader for multiple format.
    """
    if is_binary:
        mode = "rb"
    else:
        mode = "rt"
    return TqdmLineReader(filename, mode=mode, **kwargs)
