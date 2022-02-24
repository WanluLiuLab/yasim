# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: tqdm_importer.py -- Import tqdm without messing up stderr.
#
#  VERSION HISTORY:
#  2021-09-01 0.1  : Purposed and added by YU Zhejian.
#
# ==============================================================================
"""
tqdm_importer.py -- Import tqdm without messing up stderr.
"""

import sys

import tqdm as _external_tqdm

if sys.stderr.isatty():

    def tqdm(**kwargs):
        return _external_tqdm.tqdm(**kwargs)
else:
    class tqdm:

        __slots__ = ['iterable', 'total', 'desc', 'n', 'quaters']

        def __init__(self, iterable=None, desc="", total=None, **kwargs):
            self.desc = desc
            self.n = 0
            self.quaters = 0
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

        def __exit__(self, exc_type, exc_value, traceback):
            """
            """
            pass

        def update(self, i: int = 1):
            self.n += i
            if self.total:
                percent = round(self.n / self.total, 2)
                if percent > self.quaters:
                    total_len = 100 - len(self.desc)
                    print(
                        f"{self.desc}: {int(percent * 100)}% [{'=' * int(self.quaters * total_len)}|{' ' * int((1 - self.quaters) * total_len)}]",
                        file=sys.stderr)
                    self.quaters += 0.25
