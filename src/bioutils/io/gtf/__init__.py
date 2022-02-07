# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: new_gtf.py -- General-purposed GTF reader.
#
#  VERSION HISTORY:
#  VERSION HISTORY:
#  2021-08-16 0.1  : Purposed and added by YU Zhejian.
#  2021-08-24 0.1  : IntervalTree added.
#  2021-08-24 0.2  : Refactored as is designed by YUAN Ruihong.
#  2021-08-30 0.3  : Refactored as is designed by YUAN Ruihong.
#  2021-09-12 0.4  : Cythonized.
#
# ==============================================================================

"""
gtf.py -- General-purposed GTF reader.

Just as Fasta, this general-purposed reader should support multiple backends.
But currently it supports tetgs only.

.. warning::
    This file uses 1-based [] indexing!
"""

__version__ = 0.4

try:
    from bioutils.io.gtf._c_gtf_record import GtfRecord
except ImportError:
    from bioutils.io.gtf._py_gtf_record import GtfRecord
from bioutils.io.gtf._gtf_view import SimpleGtfView, GtfIterator

__all__ = [
    'GtfRecord',
    'SimpleGtfView',
    'GtfIterator',
    '__version__'
]

