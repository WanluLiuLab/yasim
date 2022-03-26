import os
import pathlib

import pytest

from commonutils import sysctl
from commonutils.io import file_system

if not sysctl.is_windows():
    pytest.skip("System is not widows", allow_module_level=True)


def test_get_abspath():
    assert file_system.get_abspath('') == ''
    home = str(pathlib.Path.home())
    assert file_system.get_abspath('~') == home
    assert file_system.get_abspath('~\\..\\') == os.path.dirname(home)
    assert file_system.get_abspath('~\\~') == home + '\\~'


def test_file_exists():
    pass  # TODO


def test_directory_exists():
    pass  # TODO


def test_is_soft_link():
    pass  # TODO
