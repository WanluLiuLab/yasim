import os

import pytest

import conftest
from commonutils import shell_utils, sysctl
from commonutils.io.file_system import file_exists, get_abspath
from commonutils.sysctl import is_windows


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> conftest.ModuleTestInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = conftest.ModuleTestInfo(session_test_info.base_test_dir, __name__)
    yield module_test_info
    module_test_info.teardown()


def test_mkdir_p_and_rm_f(initialize_module):
    if sysctl.is_user_admin() == 0:
        if not sysctl.is_windows():
            with pytest.raises(PermissionError):
                shell_utils.mkdir_p('/root/__')
            with pytest.raises(PermissionError):
                shell_utils.touch('/root/__')
        else:
            with pytest.raises(PermissionError):
                shell_utils.mkdir_p('C:\\Windows\\__')
            with pytest.raises(PermissionError):
                shell_utils.touch('C:\\Windows\\__')
    if not sysctl.is_windows():
        assert file_exists('/dev/null', allow_special_paths=True)
        assert not file_exists('/dev/null', allow_special_paths=False)
    aa = os.path.join(initialize_module.path, "aa")
    shell_utils.touch(aa)
    with pytest.raises(IsADirectoryError):
        shell_utils.touch(initialize_module.path)
    assert os.path.isdir(initialize_module.path)
    assert os.path.isfile(aa)
    shell_utils.rm_rf(aa)
    shell_utils.touch(aa)
    shell_utils.rm_rf(initialize_module.path)
    assert not os.path.isdir(initialize_module.path)
    assert not os.path.isfile(aa)


def test_wc_c(initialize_module):
    if sysctl.is_windows():
        return
    aa = os.path.join(initialize_module.path, "aa")
    shell_utils.touch(aa)
    assert shell_utils.wc_c(aa) == 0
    if not sysctl.is_windows():
        assert shell_utils.wc_c('/dev/null') == 0
        assert shell_utils.wc_c('/dev/zero') == 0
        assert shell_utils.wc_c('/dev/stdin') == 0
    dot_profile = '~/.profile'
    if file_exists(dot_profile):
        assert shell_utils.wc_c(dot_profile) == os.path.getsize(dot_profile)
    shell_utils.rm_rf(aa)


def test_readlink_f(initialize_module):
    assert shell_utils.readlink_f('') == ''
    test_path = initialize_module.path
    if not is_windows():
        shell_utils.rm_rf(test_path)
        shell_utils.touch(os.path.join(test_path, 'aa'))
        os.symlink(os.path.join(test_path, 'aa'), os.path.join(test_path, 'ab'))
        os.symlink(os.path.join(test_path, 'ab'), os.path.join(test_path, 'ac'))
        assert shell_utils.readlink_f(os.path.join(test_path, 'ab')) == get_abspath(
            os.path.join(test_path, 'aa'))
        assert shell_utils.readlink_f(
            get_abspath(os.path.join(test_path, 'ab'))) == get_abspath(
            os.path.join(test_path, 'aa'))
        assert shell_utils.readlink_f(os.path.join(test_path, 'ac')) == get_abspath(
            os.path.join(test_path, 'aa'))
        assert shell_utils.readlink_f(
            get_abspath(os.path.join(test_path, 'ac'))) == get_abspath(
            os.path.join(test_path, 'aa'))
        shell_utils.rm_rf(test_path)
