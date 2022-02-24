# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_logger.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-15 0.1  : Added by YU Zhejian.
#
# ==============================================================================
"""
test_logger.py -- Unit test of corresponding module.
"""
import logging

import pytest


def test_logger():
    lh = logging.getLogger()
    with pytest.raises(ValueError):
        logger.set_level('NON_EXISTENT_LEVEL')

    assert logging.getLevelName(logger.TRACE) == 'TRACE'
    logger.set_level(logging.DEBUG, quiet=False)
    assert lh.getEffectiveLevel() == logging.DEBUG
    assert lh.isEnabledFor(logging.DEBUG)
    assert not lh.isEnabledFor(logger.TRACE)
    logger.set_level('TRACE')
    assert lh.getEffectiveLevel() == logger.TRACE
    assert lh.isEnabledFor(logger.TRACE)
    lh.trace("TRACE")
    lh.debug("DEBUG")
    lh.info("INFO")
    lh.warning("WARNING")
    lh.error("ERROR")
    logger.set_level(logging.DEBUG)
