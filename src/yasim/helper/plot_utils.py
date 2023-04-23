"""
plot_utils.py -- Utility functions in GMM fitting process.
"""

import numpy as np
import numpy.typing as npt

from labw_utils import UnmetDependenciesError

try:

    import matplotlib
    from matplotlib import pyplot as plt
except ImportError as e:
    raise UnmetDependenciesError("matplotlib") from e


matplotlib.use('qtagg')

__all__ = ("plot",)


def plot(data: npt.ArrayLike, model, title: str = ""):
    """
    Plot how fitted GMM model performs over actual data
    """

    def normalize(_simulated_data: npt.ArrayLike) -> npt.ArrayLike:
        return _simulated_data[_simulated_data >= 0][:len(data)]

    print(str(model))
    # print(list(model.export()))
    n_bins = np.arange(0, 10, 0.1)

    simulated_data = model.rvs(2 * len(data))
    data_n = data
    simulated_data_n = normalize(simulated_data)

    fig, axes = plt.subplots(1, 1)
    axes.hist(data_n, bins=n_bins, density=True, color="blue", alpha=0.2)
    axes.hist(simulated_data_n, bins=n_bins, density=True, color="red", alpha=0.2)
    print(np.mean(data), np.mean(simulated_data))
    plt.title(f"{title} Red: Simulated, Blue: Real")
    plt.show()
