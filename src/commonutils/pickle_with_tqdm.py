# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: pickle_with_tqdm.py -- Unpickle with progress bar
#
#  VERSION HISTORY:
#  2021-08-30 0.1  : Purposed and added by YU Zhejian
#
# ==============================================================================
"""
pickle_with_tqdm.py -- Unpickle with progress bar
"""

from pickle import Unpickler

from commonutils import ioctl
from commonutils.tqdm_importer import tqdm

__version__ = 0.1


class _UnpickleWithTQDM:
    """
    The underlying class with standard functions defined.
    """

    def __init__(self, path, **kwargs):
        super().__init__()
        self.fd = open(path, 'rb')
        self.tqdm = tqdm(desc=f"Unpickling {path}", total=ioctl.wc_c(path), **kwargs)

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


# noinspection all
def unpickle_with_tqdm(filename: str):
    """
    Unpickle a file with tqdm.

    :return: Picked object.
    """
    with _UnpickleWithTQDM(filename) as pbfd:
        up = Unpickler(pbfd)
        obj = up.load()
    return obj
