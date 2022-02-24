"""
tqdm_importer.py -- Import `tqdm` without messing up stderr

This module imports `tqdm`, the progress bar implementation on Python.

If import is failed or stderr is not a Pseudo Terminal,
will use a home-made fallback which is more silent.
"""

import sys

try:
    import tqdm as _external_tqdm
except ImportError:
    _external_tqdm = None
from commonutils.importer._silent_tqdm import tqdm as _silent_tqdm

IMPORTED_TQDM_TYPE = ""

if sys.stderr.isatty() and _external_tqdm is not None:
    IMPORTED_TQDM_TYPE = "official"
else:
    IMPORTED_TQDM_TYPE = "silent"


def tqdm(*args, **kwargs):
    global IMPORTED_TQDM_TYPE
    if IMPORTED_TQDM_TYPE == "official":
        return _external_tqdm.tqdm(*args, **kwargs)
    else:
        return _silent_tqdm(*args, **kwargs)
