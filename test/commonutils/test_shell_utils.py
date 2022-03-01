import os

import pytest

import commonutils
import test_tetgs
from commonutils import shell_utils, sysctl
from commonutils.io.file_system import file_exists, get_abspath
from commonutils.sysctl import is_windows

test_path = test_tetgs.initialize(__name__)


def test_mkdir_p_and_rm_f():
    if commonutils.sysctl.is_user_admin() == 0:
        if not commonutils.sysctl.is_windows():
            with pytest.raises(PermissionError):
                shell_utils.mkdir_p('/__')  # FIXME: Error in Fedora
            with pytest.raises(PermissionError):
                shell_utils.touch('/__')
        else:
            with pytest.raises(PermissionError):
                shell_utils.mkdir_p('C:\\Windows\\__')
            with pytest.raises(PermissionError):
                shell_utils.touch('C:\\Windows\\__')
    if not commonutils.sysctl.is_windows():
        assert file_exists('/dev/null', allow_special_paths=True)
        assert not file_exists('/dev/null', allow_special_paths=False)
    shell_utils.rm_rf(test_path)
    aa = f'{test_path}/aa'
    shell_utils.touch(aa)
    with pytest.raises(IsADirectoryError):
        shell_utils.touch(test_path)
    os.path.isdir(test_path)
    os.path.isfile(aa)
    shell_utils.rm_rf(aa)
    shell_utils.touch(aa)
    shell_utils.rm_rf(test_path)
    assert not os.path.isdir(test_path)
    assert not os.path.isfile(aa)


def test_wc_c():
    if sysctl.is_windows():
        pass
    aa = f'{test_path}/aa'
    shell_utils.touch(aa)
    assert shell_utils.wc_c(aa) == 0
    if not sysctl.is_windows():
        assert shell_utils.wc_c('/dev/null') == 0
        assert shell_utils.wc_c('/dev/zero') == 0
        assert shell_utils.wc_c('/dev/stdin') == 0
    PROFILE = '~/.profile'
    if file_exists(PROFILE):
        assert shell_utils.wc_c(PROFILE) == os.path.getsize(PROFILE)
    shell_utils.rm_rf(aa)


def test_readlink_f():
    assert shell_utils.readlink_f('') == ''
    if not is_windows():
        shell_utils.rm_rf(test_path)
        shell_utils.touch(f'{test_path}/aa')
        os.symlink(f'{test_path}/aa', f'{test_path}/ab')
        assert shell_utils.readlink_f(f'{test_path}/ab') == get_abspath(
            f'{test_path}/aa')
        assert shell_utils.readlink_f(
            get_abspath(f'{test_path}/ab')) == get_abspath(
            f'{test_path}/aa')
        shell_utils.rm_rf(test_path)
