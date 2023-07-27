"""
plot_utils.py -- Utility functions in GMM fitting process.
    
.. versionadded:: 3.1.5
"""

import numpy as np
import numpy.typing as npt

from labw_utils import UnmetDependenciesError
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import Optional, Any

try:

    import matplotlib
    from matplotlib import pyplot as plt
except ImportError as e:
    raise UnmetDependenciesError("matplotlib") from e

# Try to use Qt. If not available, will fall to defaults.
try:
    matplotlib.use('qtagg')
except ImportError as e:
    pass

__all__ = ("plot",)


def plot(data: npt.NDArray, model: Optional[Any] = None, title: str = ""):
    """
    Plot how fitted GMM model performs over actual data

    .. versionadded:: 3.1.5
    """

    def normalize(_simulated_data: npt.NDArray) -> npt.NDArray:
        return _simulated_data[_simulated_data >= 0][:len(data)]

    data = np.array(data)
    n_bins = np.arange(0, 10, 0.1)
    if model is not None:
        print(str(model))

        simulated_data = model.rvs(2 * len(data))
        simulated_data_n = normalize(simulated_data)

        _, axes = plt.subplots(1, 1)
        axes.hist(data, bins=n_bins, density=True, color="blue", alpha=0.2)
        axes.hist(simulated_data_n, bins=n_bins, density=True, color="red", alpha=0.2)
        print("data:                        " + describe(data))
        print("simulated_data:              " + describe(simulated_data))
        print("simulated_data (normalized): " + describe(simulated_data_n))
        plt.title(f"{title} Red: Simulated, Blue: Real")
        plt.show()
    else:
        plt.hist(data, bins=n_bins, density=True, color="blue", alpha=0.2)
        plt.title(f"{title} Blue: Real")
        plt.show()
