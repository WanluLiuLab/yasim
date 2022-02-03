# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: sysctl.py -- System Control library.
#
#  VERSION HISTORY:
#  2021-08-23 0.1  : Migrated from ioctl.py
#
# ==============================================================================
"""
sysctl.py -- A rewritten of various sysadmin commands.

Due to the issues caused by default os.path modules, ww rewrote some of them to make it behave as is used in GNU Core
Utils. Reasons are provided.

This module is intended to work on Microsoft Windows.
"""

import ctypes
import os
import sys

_valid_os_name = (
    'freebsd',
    'win32',
    'cygwin',  # MSYS2, CYGWin
    'darwin',  # MacOS
    'unknown')
_valid_os_feature = ('nt', 'posix')
_this_platform = sys.platform


def get_os_name() -> str:
    """
    Fet a valid os name for current operating system. e.g. 'freebsd'
    """
    for name in _valid_os_name:
        if _this_platform.startswith(name):
            return name
    return 'unknown'


def is_windows() -> bool:
    """
    This piece of script determines whether the operating system in Microsoft Windows.
    """
    return get_os_name() == 'win32'


def is_user_admin() -> int:
    """
    To detect whether the user is root (admin) or not.
    A part of this function comes from
    <https://stackoverflow.com/questions/19672352/how-to-run-script-with-elevated-privilege-on-windows>

    :return: 0=no, 1=yes, -1=error
    """
    if os.name == 'nt':
        # WARNING: requires Windows XP SP2 or higher!
        try:
            return ctypes.windll.shell32.IsUserAnAdmin()
        except AttributeError:
            return -1
    elif os.name == 'posix':
        if os.getuid() == 0:  # root
            return 1
        else:
            return 0
    else:
        return -1


def is_64bit() -> bool:
    """
    Whether this is a 64bit system.
    """
    return sys.maxsize > 2 ** 32


def is_little_endian() -> bool:
    """
    Whether the system is little-endian.
    """
    return sys.byteorder == 'little'
