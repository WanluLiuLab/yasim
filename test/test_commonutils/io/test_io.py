import io
import os
import random
import string

import pytest

import conftest
from commonutils import shell_utils
from commonutils.io import ArchiveBaseIO
from commonutils.io.safe_io import get_writer, get_reader
from commonutils.io.tqdm_reader import get_tqdm_reader, get_tqdm_line_reader


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> conftest.ModuleTestInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = conftest.ModuleTestInfo(session_test_info.base_test_dir, __name__)
    yield module_test_info
    module_test_info.teardown()


available_chars = string.printable

contents = "".join((random.choice(available_chars) for _ in range(1024 * 1024)))
len_contents = len(contents)

contents_list = contents.splitlines()


def assess_archive_io(filename: str):
    # FIXME: Bugs here
    shell_utils.rm_rf(filename)
    with get_writer(filename, encoding="UTF-8", newline='') as writer:
        writer.write(contents)
    with get_reader(filename, encoding="UTF-8", newline='') as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename, encoding="UTF-8", newline='') as reader:
        assert reader.read(len_contents) == contents
        assert reader._tqdm.total == 1 + len(contents)
    with get_tqdm_line_reader(filename, encoding="UTF-8", newline='') as reader:
        i = 0
        assert reader._tqdm.total == 1 + len(contents_list)
        for line in reader:
            assert contents_list[i] == line
            i += 1
            assert reader._tqdm._n == i


extensions = (
    "txt", "xz", "bz2", "lzma", "gz"
)


@pytest.mark.parametrize(
    argnames="ext",
    argvalues=extensions,
    ids=["test_" + ext for ext in extensions]
)
def test_ext(initialize_module, ext: str):
    module_test_info = initialize_module
    filename = os.path.join(module_test_info.path, f"1.{ext}")
    assess_archive_io(filename)
    shell_utils.rm_rf(filename)


def test_string_io():
    # FIXME
    sio = io.StringIO(contents)
    bare_archive_io = ArchiveBaseIO(sio)
    assert bare_archive_io.read(len_contents) == contents
