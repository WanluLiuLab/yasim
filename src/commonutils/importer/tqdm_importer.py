"""
tqdm_importer.py -- Import `tqdm` without messing up stderr

This module imports `tqdm`, the progress bar implementation on Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.
"""
import sys

__all__ = ("tqdm", "IMPORTED_TQDM_TYPE")

try:
    import tqdm as _external_tqdm
except ImportError:
    _external_tqdm = None
from commonutils.importer._silent_tqdm import tqdm as _silent_tqdm

IMPORTED_TQDM_TYPE = ""

if sys.stderr.isatty() and _external_tqdm is not None:
    IMPORTED_TQDM_TYPE = "official"
    tqdm = _external_tqdm.tqdm
else:
    IMPORTED_TQDM_TYPE = "silent"
    tqdm = _silent_tqdm
