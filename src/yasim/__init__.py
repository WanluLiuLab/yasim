"""
yasim -- Yet Another SIMulator for Alternative Splicing and Realistic Gene Expression Profile
"""

__version__ = "3.1.4"
__author__ = ",".join((
    "YU Zhejian"
))

description = __doc__.splitlines()[1]

try:
    import labw_utils

    _labw_utils_version = labw_utils.__version__

except ImportError as e:
    raise e

from labw_utils import UnmetDependenciesError

try:

    import numpy
    import numpy.typing as npt
except ImportError as e:
    raise UnmetDependenciesError("numpy") from e
try:

    import pandas
except ImportError as e:
    raise UnmetDependenciesError("pandas") from e

try:
    import joblib
except ImportError as e:
    raise UnmetDependenciesError("joblib") from e

try:
    import scipy
except ImportError as e:
    raise UnmetDependenciesError("scipy") from e
