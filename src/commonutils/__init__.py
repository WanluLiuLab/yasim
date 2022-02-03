# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: __init__.py -- Package indicator for util.
#
#  VERSION HISTORY:
#  2021-08-10 0.1  : Purposed and added by YU Zhejian.
#
# ==============================================================================
"""
Utils -- Common sysadmin utilities.
"""

import pyximport

pyximport.install(language_level=3, inplace=True, load_py_module_on_import_failure=True)
