"""
This module defines pkgutil.resolve_name for Python >= 3.7, < 3.9

From https://pypi.org/project/pkgutil_resolve_name/
"""
import importlib
from typing import Any


def resolve_name(name: str) -> Any:
    """
    A Naive implementation that does not support wildcards.

    :param name: Name of an object.
    :return: Resolved object
    """
    parts = name.split('.')
    modname = parts.pop(0)
    # first part *must* be a module/package.
    mod = importlib.import_module(modname)
    while parts:
        p = parts[0]
        s = f'{modname}.{p}'
        try:
            mod = importlib.import_module(s)
            parts.pop(0)
            modname = s
        except ImportError:
            break
    # if we reach this point, mod is the module, already imported, and
    # parts is the list of parts in the object hierarchy to be traversed, or
    # an empty list if just the module is wanted.
    result = mod
    for p in parts:
        result = getattr(result, p)
    return result


try:
    from pkgutil import resolve_name
except ImportError:
    pass
