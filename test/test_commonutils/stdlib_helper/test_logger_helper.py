import logging

import pytest

from commonutils.stdlib_helper import logger_helper


def test_logger():
    lh = logging.getLogger()
    with pytest.raises(ValueError):
        logger_helper.set_level('NON_EXISTENT_LEVEL')

    assert logging.getLevelName(logger_helper.TRACE) == 'TRACE'
    logger_helper.set_level(logging.DEBUG, quiet=False)
    assert lh.getEffectiveLevel() == logging.DEBUG
    assert lh.isEnabledFor(logging.DEBUG)
    assert not lh.isEnabledFor(logger_helper.TRACE)
    logger_helper.set_level('TRACE')
    assert lh.getEffectiveLevel() == logger_helper.TRACE
    assert lh.isEnabledFor(logger_helper.TRACE)
    lh.trace("TRACE")
    lh.debug("DEBUG")
    lh.info("INFO")
    lh.warning("WARNING")
    lh.error("ERROR")
    logger_helper.set_level(logging.DEBUG)
