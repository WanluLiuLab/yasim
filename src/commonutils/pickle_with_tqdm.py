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
import pickle
from pickle import Unpickler
from typing import Any

__version__ = 0.1

from commonutils.io.safe_io import get_writer

from commonutils.io.tqdm_reader import get_tqdm_reader


# noinspection all
def load(filename: str):
    """
    Unpickle a file with tqdm.

    :return: Picked object.
    """
    with get_tqdm_reader(filename, is_binary=True) as pbfd:
        up = Unpickler(pbfd)
        obj = up.load()
    return obj


def dump(obj: Any, filename: str):
    with get_writer(filename, is_binary=True) as writer:
        pickle.dump(obj, writer)
