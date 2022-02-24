import sys
from importlib import reload

import pytest

from commonutils import sysctl

pytest.mark.skipif(sysctl.is_windows(), "Not tested on Windows")

try:
    import pty
except ImportError:
    pytest.mark.skip("pty module not found")
    pty = None

# test if tqdm is installed

from commonutils.importer import tqdm_importer


def test_import_tty():
    """
    TODO: No idea how ro attach a TTY to ``sys.stdout``.
    """
    pass


def test_import_no_tty():
    sys.stderr = open("/dev/null", "w")
    reload(tqdm_importer)
    assert tqdm_importer.IMPORTED_TQDM_TYPE == "silent"
