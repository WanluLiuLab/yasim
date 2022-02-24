"""
safe_io -- A Safe Wrapper for :py:mod:`commonutils.io`
======================================================

This is a "safe" IO. It does follow things:

On reader, it ensures existence of file being read by throwing errors.

On writer or appender, it ensures existence of file by :py:func:`touch`-ing it first.
"""
from typing import IO

from commonutils import shell_utils
from commonutils.io import get_appender as _get_appender
from commonutils.io import get_reader as _get_reader, file_system
from commonutils.io import get_writer as _get_writer
from commonutils.stdlib_helper.docstring_helper import copy_doc


@copy_doc(_get_reader)
def get_reader(filename: str, **kwargs) -> IO:
    if file_system.file_exists(filename, allow_special_paths=True):
        return _get_reader(filename, **kwargs)
    else:
        raise FileNotFoundError(f"File {filename} not found!")


@copy_doc(_get_writer)
def get_writer(filename: str, **kwargs) -> IO:
    if not file_system.file_exists(filename, allow_special_paths=True):
        shell_utils.touch(filename)
    return _get_writer(filename, **kwargs)


@copy_doc(_get_appender)
def get_appender(filename: str, **kwargs) -> IO:
    if file_system.file_exists(filename, allow_special_paths=True):
        return _get_appender(filename, **kwargs)
    else:
        shell_utils.touch(filename)
