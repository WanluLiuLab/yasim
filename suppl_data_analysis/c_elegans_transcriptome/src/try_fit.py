import math
import multiprocessing
import queue
import warnings
from typing import Callable, Tuple, Optional, List

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import stats as ss
from scipy.stats import rv_continuous

from commonutils.stdlib_helper.parallel_helper import ParallelJobQueue, TimeOutKiller


class FitResult:
    dist_name: str
    _dist: rv_continuous
    dist_fit_param: Tuple
    dist_loc_scale: Tuple[float, float]
    _dist_aic: Optional[float]
    data: npt.ArrayLike
    _log_likelihood: Optional[float]

    def __init__(self, data: npt.ArrayLike, dist_name: str, dist_fit_param: Tuple, dist_loc_scale: Tuple[float, float]):
        self.data = data
        self.dist_name = dist_name
        self._dist = getattr(ss, self.dist_name)
        self.dist_fit_param = dist_fit_param
        self.dist_loc_scale = dist_loc_scale
        self._dist_aic = None
        self._log_likelihood = None

    @property
    def dist_aic(self):
        if self._dist_aic is None:
            self._dist_aic = 2 * len(self.dist_fit_param) - 2 * self.log_likelihood
        return self._dist_aic

    @property
    def log_likelihood(self):
        if self._log_likelihood is None:
            self._log_likelihood = np.sum(
                self._dist.logpdf(self.data, *self.dist_fit_param)
            )
        return self._log_likelihood

    def plot_histogram(self):
        warnings.filterwarnings(action="ignore")
        plt.figure(figsize=(10, 10))
        plt.hist(x=self.data, bins=150)
        data_max = int(np.max(self.data))
        x = range(0, data_max)
        y = self._dist.pdf(x, *self.dist_fit_param) * self.dist_loc_scale[1] + self.dist_loc_scale[0]
        plt.plot(x, y, lw=3)
        plt.yscale('log')
        plt.title(f"{self.dist_name}{self.dist_fit_param}, AIC={self.dist_aic}")
        plt.grid(True)
        plt.savefig(f"{self.dist_name}.png")
        plt.close()

    def __repr__(self):
        return f"Fitting {self.dist_name}, with AIC={self.dist_aic}"

    def __str__(self):
        return repr(self)


class AutoFitProcess(multiprocessing.Process):
    data: npt.ArrayLike
    dist_name: str
    _out_queue: multiprocessing.Queue
    _dist: rv_continuous
    _tok: TimeOutKiller

    def __init__(self, data: npt.ArrayLike, dist_name: str, out_queue: multiprocessing.Queue):
        super().__init__()
        self.data = data
        self.dist_name = dist_name
        self._dist = getattr(ss, self.dist_name)
        print(self._dist)
        self._out_queue = out_queue

    def run(self):
        self._tok = TimeOutKiller(process_or_pid=self, timeout=40)
        self._tok.start()
        warnings.filterwarnings(action="ignore")
        for fit_data in self.data, np.array(self.data, dtype=int):
            for method in ("MM", "MLE"):
                try:
                    dist_fit_param = self._dist.fit(fit_data, method=method)
                    dist_loc_scale = self._dist.fit_loc_scale(fit_data, *dist_fit_param)
                    fit_result = FitResult(fit_data, self.dist_name, dist_fit_param, dist_loc_scale)
                    fit_result.plot_histogram()
                    self._out_queue.put(fit_result)
                    print(f"PASSED {self.dist_name}")
                    return
                except Exception:
                    pass
        print(f"FAILED {self.dist_name}")


if __name__ == "__main__":
    fitable_functions: List[str] = []
    for spec in dir(ss):
        if hasattr(getattr(ss, spec), "fit_loc_scale") and \
                isinstance(getattr(getattr(ss, spec), "fit_loc_scale"), Callable) and \
                hasattr(getattr(ss, spec), "logpdf") and \
                isinstance(getattr(getattr(ss, spec), "logpdf"), Callable):
            fitable_functions.append(spec)
    print(fitable_functions)
    all_data = pd.read_table("../all_data.tsv")
    nanopore_coverage = all_data.loc[:, 'NANOPORE_AVG_DEPTH'].to_numpy(dtype=float)
    nanopore_coverage = nanopore_coverage[np.where(nanopore_coverage > 0)]
    fs: List[FitResult] = []
    fit_job_queue = ParallelJobQueue(pool_name="Fitting", pool_size=math.inf)
    sync_manager = multiprocessing.Manager()
    out_queue = sync_manager.Queue()
    for spec in fitable_functions:
        fit_job_queue.append(
            AutoFitProcess(nanopore_coverage, spec, out_queue)
        )
    # sync_manager.start()
    fit_job_queue.start()
    while fit_job_queue.is_alive():
        try:
            get_item = out_queue.get(timeout=0.1)
            fs.append(get_item)
        except (TimeoutError, queue.Empty):
            pass
    fs.sort(key=lambda x: x.dist_aic)
    print(fs)
    fit_job_queue.join()
    sync_manager.shutdown()
