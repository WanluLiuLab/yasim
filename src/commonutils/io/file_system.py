import os
import stat

from commonutils.logger import chronolog


@chronolog(display_time=True)
def get_abspath(path: str) -> str:
    """
    Get the absolute path of a given relative path. Can be a non-existent file or directory.

    The default one provided by os.path.abspath() or os.path.realpath()
    cannot perform tide expand (That is, expand ~ to your HOME directory).

    :param path: The relative path
    :return: The absolute path
    """
    if path == '':
        return path
    path = os.path.abspath(os.path.expanduser(path))
    return path


def file_exists(path: str, allow_special_paths: bool = True) -> bool:
    """
    Check the existence of a file. Can check regular and special files.

    :param path: The path you wish to examine.
    :param allow_special_paths: Whether this is a special file, may be block device, etc.
    """
    if allow_special_paths:
        return os.path.exists(path) and not os.path.isdir(path)
    else:
        return os.path.isfile(path)


@chronolog(display_time=True)
def is_soft_link(path: str) -> bool:
    return stat.S_ISLNK(os.stat(path, follow_symlinks=False)[0])
