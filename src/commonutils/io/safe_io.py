from typing import IO

from commonutils import shell_utils
from commonutils.io import get_appender as _get_appender
from commonutils.io import get_reader as _get_reader, file_system
from commonutils.io import get_writer as _get_writer


def get_reader(filename: str, **kwargs) -> IO:
    """
    Get a reader for multiple format.
    """
    if file_system.file_exists(filename, allow_special_paths=True):
        return _get_reader(filename, **kwargs)
    else:
        raise FileNotFoundError(f"File {filename} not found!")


def get_writer(filename: str, **kwargs) -> IO:
    """
    Get a writer for multiple format.
    """
    if not file_system.file_exists(filename, allow_special_paths=True):
        shell_utils.touch(filename)
    return _get_writer(filename, **kwargs)


def get_appender(filename: str, **kwargs) -> IO:
    """
    Get an appender for multiple format.
    """
    if file_system.file_exists(filename, allow_special_paths=True):
        return _get_appender(filename, **kwargs)
    else:
        shell_utils.touch(filename)
