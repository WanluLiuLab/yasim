"""
pickle_helper.py -- Pickle helper with compression and progress-bar.

Example:

We may firstly generate some random sequence:

>>> import random
>>> import tempfile
>>> from commonutils.shell_utils import rm_rf
>>> random_arr = [random.random() for _ in range(1000)]
>>> test_path = tempfile.mkdtemp()

Test with a normal file:

>>> pickle_fn = f'{test_path}/rd.pickle'
>>> dump(random_arr, pickle_fn)
>>> unpickle_obj = load(pickle_fn)
>>> assert unpickle_obj == random_arr

This module can also handle compressed pickle:

>>> pickle_fn = f'{test_path}/rd.pickle.xz'
>>> dump(random_arr, pickle_fn)
>>> unpickle_obj = load(pickle_fn)
>>> assert unpickle_obj == random_arr

Without a progress bar:
>>> unpickle_obj = load(pickle_fn, with_tqdm=False)
>>> assert unpickle_obj == random_arr

Clean up the environment.

>>> rm_rf(test_path)
"""
import pickle
from pickle import Unpickler
from typing import Any

from commonutils.io.safe_io import get_writer, get_reader
from commonutils.io.tqdm_reader import get_tqdm_reader


def load(filename: str, with_tqdm: bool = True) -> Any:
    """
    Unpickle a file with tqdm.

    :param filename: Filename to be load from
    :param with_tqdm: Whether to display a progress bar.
    :return: Picked object.
    """
    if with_tqdm:
        fd = get_tqdm_reader(filename, is_binary=True)
    else:
        fd = get_reader(filename, is_binary=True)
    up = Unpickler(fd)
    obj = up.load()
    fd.close()
    return obj


def dump(obj: Any, filename: str):
    """
    Pickle an object into a file.

    :param obj: The object to be pickled.
    :param filename: The filename to be written to.
    If the filename have compressed suffixes like ``xz``, it will be compressed.
    """
    with get_writer(filename, is_binary=True) as writer:
        pickle.dump(obj, writer)
