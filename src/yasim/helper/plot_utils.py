import multiprocessing

import matplotlib
import numpy as np
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from numpy import typing as npt

matplotlib.use('qtagg')


def plot(data: npt.ArrayLike, model):
    def normalize(_simulated_data: npt.ArrayLike) -> npt.ArrayLike:
        return np.log(list(map(int, _simulated_data[_simulated_data > 1][:len(data)])))

    def generated_sim1(_model, _len: int) -> npt.ArrayLike:
        reta = np.array(Parallel(n_jobs=multiprocessing.cpu_count())(delayed(_model.rvs)(size=1) for _ in range(_len)))
        return reta

    print(str(model))
    n_bins = np.arange(0, 10, 0.1)

    simulated_data = model.rvs(2 * len(data))
    simulated_data1 = generated_sim1(model, 2 * len(data))

    data_n = np.log(data)
    simulated_data_n = normalize(simulated_data)
    simulated_data1_n = normalize(simulated_data1)

    fig, axes = plt.subplots(1, 1)
    axes.hist(data_n, bins=n_bins, density=True, color="blue", alpha=0.2)
    axes.hist(simulated_data_n, bins=n_bins, density=True, color="red", alpha=0.2)
    axes.hist(simulated_data1_n, bins=n_bins, density=True, color="black", alpha=0.2)
    print(np.mean(data), np.mean(simulated_data), np.mean(simulated_data1))
    plt.show()
