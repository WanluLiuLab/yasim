"""
_silent_tqdm.py -- A silent tqdm that does not pollutes stderr
"""

import sys
from typing import Optional, Iterable


class tqdm:
    """
    A silent tqdm that does not pollute stderr.

    This module will be imported by :py:mod:`commonutils.importer.tqdm_importer` if stderr is not a TTY,
    which usually happens when the stderr is connected to a logger.

    This module does not support rich functions as TQDM.
    It only supports 4 quarters (i.e., it will change only 4 times)

    Examples:

    Using tqdm as an iterator:

    >>> list(tqdm(iterable=range(10)))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    Using tqdm as context manager:

    >>> out_list = []
    >>> with tqdm(total=10) as pbar:
    ...    for i in range(10):
    ...        out_list.append(i)
    ...        pbar.update(1)
    >>> out_list
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    __slots__ = ['iterable', 'total', 'desc', '_n', '_quarters']

    iterable: Optional[Iterable]
    """
    The iterable to be iterated when using ``tqdm`` in ``for`` loops.
    """

    desc: Optional[str]
    """
    Description show at the front of the progress bar.
    """

    total: Optional[int]
    """
    The maximum number of items to be iterated.
    
    Set this variable if ``iterable`` is not sizeable or is ``None``.
    """

    _n: int
    _quarters: float

    def __init__(self,
                 iterable: Optional[Iterable] = None,
                 desc: Optional[str] = None,
                 total: Optional[int] = None,
                 **kwargs):
        _ = kwargs  # Stop PyCharm from complaining
        self.desc = desc
        self._n = 0
        self._quarters = 0
        if total is None and iterable is not None:
            try:
                self.total = len(iterable)
            except (TypeError, AttributeError):
                self.total = None
        else:
            self.total = total
        if total == float("inf"):
            self.total = None
        self.iterable = iterable

    def __iter__(self):
        for item in self.iterable:
            yield item
            self.update(1)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return

    def update(self, i: int = 1):
        self._n += i
        if self.total:
            percent = round(self._n / self.total, 2)
            if percent > self._quarters:
                total_len = 100
                print(
                    f"{self.desc}: {int(percent * 100)}% [{'=' * int(self._quarters * total_len)}|{' ' * int((1 - self._quarters) * total_len)}]",
                    file=sys.stderr)
                self._quarters += 0.25
