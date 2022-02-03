# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_cbisect.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-10-22 0.1  : Added by YU Zhejian.
#
# ==============================================================================

from commonutils.cbisect import test_cbisect_left as _test_cbisect_left
from commonutils.cbisect import test_cbisect_right as _test_cbisect_right


def test_cbisect_left():
    _test_cbisect_left()


def test_cbisect_right():
    _test_cbisect_right()
