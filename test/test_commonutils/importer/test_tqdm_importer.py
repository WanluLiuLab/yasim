import sys
from importlib import reload

import pytest

from commonutils.importer import tqdm_importer

try:
    import pty
except ImportError:
    pytest.skip("pty module not found", allow_module_level=True)
    pty = None

try:
    open("/dev/null", "w").close()
except FileNotFoundError:
    pytest.skip("/dev/null not found", allow_module_level=True)
    pty = None


# test if tqdm is installed


def test_import_tty():
    """
    TODO: No idea how ro attach a TTY to ``sys.stdout``.
    """
    pass


def test_import_no_tty():
    sys.stderr = open("/dev/null", "w")
    reload(tqdm_importer)
    assert tqdm_importer.IMPORTED_TQDM_TYPE == "silent"
