"""
commonutils.io -- Enhanced Python IO Functions/Classes
======================================================

This module includes some enhanced IO functions,
like IO that automatically creates missing intermediate directories or IO with tqdm progress bar.
It also supports Python standard archive types like ``bz2``, ``gzip`` or ``lzma``.

For those who wish to get a ``gzip`` command-like utilities,
please visit :py:mod:`commonutils.shell_utils` for more information
"""

import os
from typing import IO

from commonutils.io._ioproxy import RuleBasedIOProxy, IOProxy

__all__ = (
    "determine_line_endings",
    "determine_file_line_endings",
    "get_reader",
    "get_writer",
    "get_appender"
)

from commonutils.io.typing import IOProxyType

from commonutils.typing import PathOrFDType


def determine_line_endings(fd: IO):
    """
    Determine line endings. If failed, will return OS default.

    """
    find_cr = False
    find_lf = False
    while True:
        c = fd.read(1)
        if c is None:
            break
        elif c == '\r':
            find_cr = True
            if find_lf:
                return '\n\r'
        elif c == '\n':
            find_lf = True
            if find_cr:
                return '\r\n'
        else:
            if find_cr:
                return '\r'
            if find_lf:
                return '\n'
            find_cr = False
            find_lf = False
    return os.linesep


def determine_file_line_endings(filename: str):
    return determine_line_endings(RuleBasedIOProxy(filename))


def get_reader(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxyType:
    """
    Get a reader for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.

    :param path_or_fd: Filename to be opened or IO that was opened..
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.

    .. warning::
        Do NOT specify ``mode`` keyword arguments!
    """
    if isinstance(path_or_fd, str):
        if is_binary:
            mode = "rb"
        else:
            mode = "rt"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)
    else:
        return IOProxy(path_or_fd)


def get_writer(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxyType:
    """
    Get a writer for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.
    
    :param path_or_fd: Filename to be opened or IO that was opened..
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.
    
    :: warning..
        Do NOT specify ``mode`` keyword arguments!
    """
    if isinstance(path_or_fd, str):
        if is_binary:
            mode = "wb"
        else:
            mode = "wt"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)
    else:
        return IOProxy(path_or_fd)


def get_appender(path_or_fd: PathOrFDType, is_binary: bool = False, **kwargs) -> IOProxyType:
    """
    Get an appender for multiple format.
    
    This function is for newbies or others who does not wish to have full control over what they opened.
    The IO wrapper given by this function may satisfy 95% of the needs.
    
    :param path_or_fd: Filename to be opened or IO that was opened.
    :param is_binary: Whether to read as binary.
    :param kwargs: Other arguments passed to underlying opener.
    
    :: warning..
        Do NOT specify ``mode`` keyword arguments!
    """
    if isinstance(path_or_fd, str):
        if is_binary:
            mode = "ab"
        else:
            mode = "at"
        return RuleBasedIOProxy(path_or_fd, mode=mode, **kwargs)
    else:
        return IOProxy(path_or_fd)
