"""
test_tetgs.py -- Unit test of corresponding module.
"""
import os
import tempfile

import pytest

from commonutils import shell_utils
from commonutils.stdlib_helper import logger_helper


class SessionTestInfo:
    __slots__ = (
        "base_test_dir"
    )
    base_test_dir: str

    def __init__(self):
        self.base_test_dir = tempfile.mkdtemp()
        lh = logger_helper.get_logger("TEST_TETGS")
        lh.info(f"Test dir {self.base_test_dir}")

    def teardown(self):
        shell_utils.rm_rf(self.base_test_dir)


class ModuleTestInfo:
    __slots__ = (
        "name",
        "path"
    )
    path: str

    def __init__(self, base_test_dir: str, name: str):
        self.name = name
        self.path = os.path.join(base_test_dir, name)
        shell_utils.mkdir_p(self.path)

    def teardown(self):
        shell_utils.rm_rf(self.path)


@pytest.fixture(scope="session")
def initialize_session():
    logger_helper.set_level(logger_helper.TRACE)
    session_test_info = SessionTestInfo()
    yield session_test_info
    session_test_info.teardown()
