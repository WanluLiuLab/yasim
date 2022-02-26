import os
import pathlib

import commonutils
import test_tetgs
from commonutils import sysctl
from commonutils.io import file_system


test_path = test_tetgs.initialize(__name__)




def test_get_abspath():
    assert file_system.get_abspath('') == ''
    if not sysctl.is_windows():
        HOME = os.environ['HOME']  # FIXME: In freebsd, /home -> /usr/home
        assert file_system.get_abspath('.') == os.environ['PWD']
        assert file_system.get_abspath('~') == HOME
        assert file_system.get_abspath('~/../') == os.path.dirname(HOME)
        assert file_system.get_abspath('~/../~') == os.path.dirname(HOME) + '/~'
        assert file_system.get_abspath('/a//b/c/') == '/a/b/c'
        assert file_system.get_abspath('/a//b/.././..////c/') == '/c'
    else:
        home = str(pathlib.Path.home())
        assert file_system.get_abspath('~') == home
        assert file_system.get_abspath('~/../') == os.path.dirname(home)
        assert file_system.get_abspath('~/~') == home + '\\~'

def test_file_exists():
    pass # TODO

def test_directory_exists():
    pass # TODO

def test_is_soft_link():
    pass # TODO
