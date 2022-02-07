from typing import IO

from commonutils import ioctl
from commonutils.tqdm_importer import tqdm


class tqdm_reader:
    """
    A very simple tqdm reader with :py:func:``read`` and :py:func:``readline`` functions.
    """
    fd: IO

    def __init__(self, path: str, underlying_opener=None, is_binary: bool = False, **kwargs):
        if underlying_opener == None:
            self.fd = ioctl.get_reader(path, is_binary=is_binary)
        else:
            if is_binary:
                self.fd = underlying_opener(path, 'rb', **kwargs)
            else:
                self.fd = underlying_opener(path, 'rt', **kwargs)
        self.tqdm = tqdm(desc=f"Reading {path}", total=ioctl.wc_c(path), unit='B', unit_scale=True, unit_divisor=1024)

    def readline(self):
        update_bytes = self.fd.readline()
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def read(self, size: int = -1):
        update_bytes = self.fd.read(size)
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def __enter__(self):
        self.tqdm.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        return self.tqdm.__exit__(*args, **kwargs)


class tqdm_line_reader:
    """
    A very simple tqdm reader with only :py:func:``readline`` functions.
    """
    fd: IO

    def __init__(self, path: str, underlying_opener=None, is_binary: bool = False, **kwargs):
        if underlying_opener == None:
            self.fd = ioctl.get_reader(path, is_binary=is_binary)
        else:
            if is_binary:
                self.fd = underlying_opener(path, 'rb', **kwargs)
            else:
                self.fd = underlying_opener(path, 'rt', **kwargs)
        self.tqdm = tqdm(desc=f"Reading {path}", total=ioctl.wc_l(path), unit='L')

    def readline(self):
        update_bytes = self.fd.readline()
        self.tqdm.update(1)
        return update_bytes

    def tell(self):
        return self.fd.tell()

    def __enter__(self):
        self.tqdm.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        return self.tqdm.__exit__(*args, **kwargs)
