import os

import pytest

from commonutils import sysctl
from commonutils.io import file_system

if sysctl.is_windows():
    pytest.skip("System is widows", allow_module_level=True)


def test_get_abspath():
    assert file_system.get_abspath('') == ''
    home_directory = os.environ.get("HOME")  # FIXME: In freebsd, /home -> /usr/home
    current_working_directory = os.environ.get("PWD")
    assert file_system.get_abspath('/a//b/c/') == '/a/b/c'
    assert file_system.get_abspath('/a//b/.././..////c/') == '/c'
    if current_working_directory is not None:
        assert file_system.get_abspath('.') == current_working_directory
    if home_directory is not None:
        assert file_system.get_abspath('~') == home_directory
        assert file_system.get_abspath('~/../') == os.path.dirname(home_directory)
        assert file_system.get_abspath('~/../~') == os.path.dirname(home_directory) + '/~'


def test_file_exists():
    pass  # TODO


def test_directory_exists():
    pass  # TODO


def test_is_soft_link():
    pass  # TODO
