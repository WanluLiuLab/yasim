import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from numpy import typing as npt

matplotlib.use('qtagg')


def plot(data: npt.ArrayLike, model):
    def normalize(_simulated_data: npt.ArrayLike) -> npt.ArrayLike:
        return _simulated_data[_simulated_data > 1][:len(data)]  # list(map(int,

    print(str(model))
    print(list(model.export()))
    n_bins = np.arange(0, 10, 0.1)

    simulated_data = model.rvs(2 * len(data))
    data_n = data
    simulated_data_n = normalize(simulated_data)

    fig, axes = plt.subplots(1, 1)
    axes.hist(data_n, bins=n_bins, density=True, color="blue", alpha=0.2)
    axes.hist(simulated_data_n, bins=n_bins, density=True, color="red", alpha=0.2)
    print(np.mean(data), np.mean(simulated_data))
    plt.show()