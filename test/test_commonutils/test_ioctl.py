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
import test_tetgs
from commonutils import ioctl, sysctl
from commonutils.sysctl import is_windows

test_path = test_tetgs.initialize(__name__)

def test_get_abspath():
    assert ioctl.get_abspath('') == ''
    if not commonutils.sysctl.is_windows():
        HOME = os.environ['HOME']  # FIXME: In freebsd, /home -> /usr/home
        assert ioctl.get_abspath('.') == os.environ['PWD']
        assert ioctl.get_abspath('~') == HOME
        assert ioctl.get_abspath('~/../') == os.path.dirname(HOME)
        assert ioctl.get_abspath('~/../~') == os.path.dirname(HOME) + '/~'
        assert ioctl.get_abspath('/a//b/c/') == '/a/b/c'
        assert ioctl.get_abspath('/a//b/.././..////c/') == '/c'
    else:
        home = str(pathlib.Path.home())
        assert ioctl.get_abspath('~') == home
        assert ioctl.get_abspath('~/../') == os.path.dirname(home)
        assert ioctl.get_abspath('~/~') == home + '\\~'


def test_mkdir_p_and_rm_f():
    if commonutils.sysctl.is_user_admin() == 0:
        if not commonutils.sysctl.is_windows():
            with pytest.raises(PermissionError):
                ioctl.mkdir_p('/__')  # FIXME: Error in Fedora
            with pytest.raises(PermissionError):
                ioctl.touch('/__')
        else:
            with pytest.raises(PermissionError):
                ioctl.mkdir_p('C:\\Windows\\__')
            with pytest.raises(PermissionError):
                ioctl.touch('C:\\Windows\\__')
    if not commonutils.sysctl.is_windows():
        assert ioctl.file_exists('/dev/null', allow_special_paths=True)
        assert not ioctl.file_exists('/dev/null', allow_special_paths=False)
    ioctl.rm_rf(test_path)
    aa = f'{test_path}/aa'
    ioctl.touch(aa)
    with pytest.raises(IsADirectoryError):
        ioctl.touch(test_path)
    os.path.isdir(test_path)
    os.path.isfile(aa)
    ioctl.rm_rf(aa)
    ioctl.touch(aa)
    ioctl.rm_rf(test_path)
    assert not os.path.isdir(test_path)
    assert not os.path.isfile(aa)


def test_wc_c():
    if sysctl.is_windows():
        pass
    aa = f'{test_path}/aa'
    ioctl.touch(aa)
    assert ioctl.wc_c(aa) == 0
    if not sysctl.is_windows():
        assert ioctl.wc_c('/dev/null') == 0
        assert ioctl.wc_c('/dev/zero') == 0
        assert ioctl.wc_c('/dev/stdin') == 0
    PROFILE = '~/.profile'
    if ioctl.file_exists(PROFILE):
        assert ioctl.wc_c(PROFILE) == os.path.getsize(PROFILE)
    ioctl.rm_rf(aa)


def test_readlink_f():
    assert ioctl.readlink_f('') == ''
    if not is_windows():
        ioctl.rm_rf(test_path)
        ioctl.touch(f'{test_path}/aa')
        os.symlink(f'{test_path}/aa', f'{test_path}/ab')
        assert ioctl.readlink_f(f'{test_path}/ab') == ioctl.get_abspath(f'{test_path}/aa')
        assert ioctl.readlink_f(ioctl.get_abspath(f'{test_path}/ab')) == ioctl.get_abspath(f'{test_path}/aa')
        ioctl.rm_rf(test_path)


contents = "A line.\nAnother line.\t\nSomething."
len_contents = len(contents)


def get_opener_family_assertions(suffix):
    filename = f"{test_path}/1.{suffix}"
    ioctl.rm_rf(filename)
    with ioctl.get_writer(filename) as writer:
        writer.write(contents)
    with ioctl.get_reader(filename) as reader:
        assert reader.read(len_contents) == contents
    ioctl.rm_rf(filename)


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
