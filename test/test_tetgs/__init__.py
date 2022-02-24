# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_tetgs.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-23 0.1  : Purposed and added by YU Zhejian.
#
# ==============================================================================
"""
test_tetgs.py -- Unit test of corresponding module.
"""
import logging
import tempfile

from commonutils import logger, shell_utils

TEST_DIR = tempfile.mkdtemp()
lh = logging.getLogger()
lh.info(f"Test dir {TEST_DIR}")


def initialize(name: str = __name__) -> str:
    test_path = TEST_DIR + name
    shell_utils.mkdir_p(test_path)
    logger.set_level(logger.TRACE)
    return test_path
