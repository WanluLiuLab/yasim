# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_ioctl.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-10 0.1  : Added by YU Zhejian.
#
# ==============================================================================
"""
test_ioctl.py -- Unit test of corresponding module.
"""

import os
import pathlib

import pytest

import commonutils
import commonutils.io.file_system
import commonutils.shutil
import test_tetgs
from commonutils import sysctl
from commonutils.io.safe_io import get_writer, get_reader
from commonutils.sysctl import is_windows

test_path = test_tetgs.initialize(__name__)


def test_get_abspath():
    assert commonutils.io.file_system.get_abspath('') == ''
    if not commonutils.sysctl.is_windows():
        HOME = os.environ['HOME']  # FIXME: In freebsd, /home -> /usr/home
        assert commonutils.io.file_system.get_abspath('.') == os.environ['PWD']
        assert commonutils.io.file_system.get_abspath('~') == HOME
        assert commonutils.io.file_system.get_abspath('~/../') == os.path.dirname(HOME)
        assert commonutils.io.file_system.get_abspath('~/../~') == os.path.dirname(HOME) + '/~'
        assert commonutils.io.file_system.get_abspath('/a//b/c/') == '/a/b/c'
        assert commonutils.io.file_system.get_abspath('/a//b/.././..////c/') == '/c'
    else:
        home = str(pathlib.Path.home())
        assert commonutils.io.file_system.get_abspath('~') == home
        assert commonutils.io.file_system.get_abspath('~/../') == os.path.dirname(home)
        assert commonutils.io.file_system.get_abspath('~/~') == home + '\\~'


def test_mkdir_p_and_rm_f():
    if commonutils.sysctl.is_user_admin() == 0:
        if not commonutils.sysctl.is_windows():
            with pytest.raises(PermissionError):
                commonutils.shutil.mkdir_p('/__')  # FIXME: Error in Fedora
            with pytest.raises(PermissionError):
                commonutils.shutil.touch('/__')
        else:
            with pytest.raises(PermissionError):
                commonutils.shutil.mkdir_p('C:\\Windows\\__')
            with pytest.raises(PermissionError):
                commonutils.shutil.touch('C:\\Windows\\__')
    if not commonutils.sysctl.is_windows():
        assert commonutils.io.file_system.file_exists('/dev/null', allow_special_paths=True)
        assert not commonutils.io.file_system.file_exists('/dev/null', allow_special_paths=False)
    commonutils.shutil.rm_rf(test_path)
    aa = f'{test_path}/aa'
    commonutils.shutil.touch(aa)
    with pytest.raises(IsADirectoryError):
        commonutils.shutil.touch(test_path)
    os.path.isdir(test_path)
    os.path.isfile(aa)
    commonutils.shutil.rm_rf(aa)
    commonutils.shutil.touch(aa)
    commonutils.shutil.rm_rf(test_path)
    assert not os.path.isdir(test_path)
    assert not os.path.isfile(aa)


def test_wc_c():
    if sysctl.is_windows():
        pass
    aa = f'{test_path}/aa'
    commonutils.shutil.touch(aa)
    assert commonutils.shutil.wc_c(aa) == 0
    if not sysctl.is_windows():
        assert commonutils.shutil.wc_c('/dev/null') == 0
        assert commonutils.shutil.wc_c('/dev/zero') == 0
        assert commonutils.shutil.wc_c('/dev/stdin') == 0
    PROFILE = '~/.profile'
    if commonutils.io.file_system.file_exists(PROFILE):
        assert commonutils.shutil.wc_c(PROFILE) == os.path.getsize(PROFILE)
    commonutils.shutil.rm_rf(aa)


def test_readlink_f():
    assert commonutils.shutil.readlink_f('') == ''
    if not is_windows():
        commonutils.shutil.rm_rf(test_path)
        commonutils.shutil.touch(f'{test_path}/aa')
        os.symlink(f'{test_path}/aa', f'{test_path}/ab')
        assert commonutils.shutil.readlink_f(f'{test_path}/ab') == commonutils.io.file_system.get_abspath(
            f'{test_path}/aa')
        assert commonutils.shutil.readlink_f(
            commonutils.io.file_system.get_abspath(f'{test_path}/ab')) == commonutils.io.file_system.get_abspath(
            f'{test_path}/aa')
        commonutils.shutil.rm_rf(test_path)


contents = "A line.\nAnother line.\t\nSomething."
len_contents = len(contents)


def get_opener_family_assertions(suffix):
    filename = f"{test_path}/1.{suffix}"
    commonutils.shutil.rm_rf(filename)
    with get_writer(filename) as writer:
        writer.write(contents)
    with get_reader(filename) as reader:
        assert reader.read(len_contents) == contents
    commonutils.shutil.rm_rf(filename)


def test_txt():
    get_opener_family_assertions("txt")


def test_bz2():
    get_opener_family_assertions("bz2")


def test_gz():
    get_opener_family_assertions("gz")


def test_xz():
    get_opener_family_assertions("xz")


def test_lzma():
    get_opener_family_assertions("lzma")
