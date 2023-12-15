"""
gmm.py -- Gaussian Mixture Model bu YUAN Ruihong

.. versionadded:: 3.1.5
"""

__all__ = ("GaussianMixture1D",)

import itertools
import multiprocessing
import secrets

import numpy as np
import numpy.typing as npt
import tqdm
from scipy.integrate import quad
from scipy.stats import norm

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import parallel_map
from labw_utils.typing_importer import (
    Optional,
    Union,
    Iterable,
    Tuple,
    List,
    Sequence,
    SequenceProxy,
    Any,
)

_lh = get_logger(__name__)


class GaussianMixture1D:
    """
    TODO docs

    .. versionadded:: 3.1.5
    """

    _n_components: int
    _n_iter: int
    _precision: float
    _weights: npt.NDArray
    _mu: npt.NDArray
    _sigma: npt.NDArray

    def __init__(
        self,
        n_components: int = 1,
        n_iter: int = 100,
        init_weights: Optional[npt.ArrayLike] = None,
        init_mus: Optional[npt.ArrayLike] = None,
        init_sigmas: Optional[npt.ArrayLike] = None,
        precision: float = 1e-4,
    ) -> None:
        self._n_components = n_components
        self._n_iter = n_iter
        self._precision = precision
        self._weights = init_weights or np.repeat(1 / n_components, n_components)
        self._mu = init_mus or np.zeros(n_components)
        self._sigma = init_sigmas or np.ones(n_components)

    def fit(self, data: npt.ArrayLike, show_progress_bar: bool = True):
        n = len(data)
        rdg = secrets.SystemRandom()

        # Use k-means to cluster points
        centers = np.array(rdg.choices(data, k=self._n_components))
        categories = np.repeat(0, n)
        bags: List[List[Union[int, float]]] = [[]] * self._n_components
        knn_range = range(self._n_iter)
        if show_progress_bar:
            knn_range = tqdm.tqdm(iterable=knn_range, desc="KNN")
        for _ in knn_range:
            converged = True
            for j in range(self._n_components):
                bags[j] = []

            for i in range(n):
                x = data[i]
                category = categories[i]
                least_dist = abs(x - centers[category])
                for j in range(self._n_components):
                    center = centers[j]
                    dist = abs(x - center)
                    if dist < least_dist:
                        converged = False
                        category = j  # assign to the nearest center
                        least_dist = dist
                categories[i] = category
                bags[category].append(x)

            for j in range(self._n_components):
                bag = bags[j]
                if len(bag) != 0:
                    centers[j] = np.mean(bag)

            if converged:
                break

        self._mu = centers
        for j in range(self._n_components):
            bag = bags[j]
            if len(bag) != 0:
                self._sigma[j] = np.std(bag)
                self._weights[j] = len(bag)
            else:
                self._sigma[j] = 0
                self._weights[j] = 0
        self._weights /= n

        last_ll = np.NINF
        densities = np.zeros((self._n_components, n))
        em_range = range(self._n_iter)
        if show_progress_bar:
            em_range = tqdm.tqdm(iterable=em_range, desc="EM")
        _last_range = self.export()
        for _ in em_range:
            # E step: compute ownership weights
            for j in range(self._n_components):
                model = norm(loc=self._mu[j], scale=self._sigma[j])
                densities[j, :] = model.pdf(data)

            a = np.transpose(np.multiply(np.transpose(densities), self._weights))
            b = np.transpose(self._weights) @ densities
            gamma = a / b

            # M step: compute means, variances and weights
            denom = np.sum(gamma, axis=1)  # row sum
            self._mu = (gamma @ data) / denom

            for j in range(self._n_components):
                for i in range(n):
                    gamma[j, i] *= (data[i] - self._mu[j]) ** 2

            self._sigma = np.sqrt(np.sum(gamma, axis=1) / denom)
            self._weights = denom / n

            ll = self.log_likelihood(data)
            if 0 <= ll - last_ll < self._precision:
                break
            last_ll = ll
            current_export = np.array(list(itertools.chain(*self.export())), dtype=float)
            if np.sum(np.bitwise_or(np.isnan(current_export), np.isinf(current_export))) != 0:
                _lh.warning("GMM: NaN observed, will rerun last successful result")
                return self.import_model(_last_range)
            else:
                _last_range = list(self.export())
        return self

    def log_likelihood(self, data):
        return np.sum(self.logpdf(data))

    def pdf(self, x: Union[int, float, npt.ArrayLike]) -> float:
        if isinstance(x, (int, float)):
            result = 0.0
        else:
            result = np.zeros(len(x))

        for j in range(self._n_components):
            model = norm(loc=self._mu[j], scale=self._sigma[j])
            result += self._weights[j] * model.pdf(x)
        return result

    def logpdf(self, x: Union[int, float, npt.ArrayLike]) -> float:
        if isinstance(x, (int, float)):
            result = 0.0
        else:
            result = np.zeros(len(x))

        for j in range(self._n_components):
            model = norm(loc=self._mu[j], scale=self._sigma[j])
            result += self._weights[j] * model.logpdf(x)
        return result

    def lintrans(self, a: float, b: float = 0):
        self._mu *= a
        self._mu += b
        self._sigma *= a
        return self

    def positive_mean(self, start: float = 0):
        models = []
        for j in range(self._n_components):
            models.append(norm(self._mu[j], self._sigma[j]))

        def f(x):
            result = 0.0
            for _j in range(self._n_components):
                result += self._weights[_j] * models[_j].pdf(x)
            return result * x

        return quad(f, start, np.inf)[0]

    def mean(self) -> float:
        return np.average(self._mu, weights=self._weights)

    def rvs(self, size: int = 1) -> npt.ArrayLike:
        indices = list(range(self._n_components))

        def _rvs(_: Any):
            rdg = secrets.SystemRandom()
            k = rdg.choices(indices, weights=self._weights)
            return rdg.normalvariate(mu=self._mu[k], sigma=self._sigma[k])

        result = np.array(
            list(parallel_map(_rvs, range(size), n_jobs=multiprocessing.cpu_count())),
            dtype=float,
        )
        return result

    def export(self) -> Sequence[Tuple[float, float, float]]:
        retl = []
        for j in range(self._n_components):
            retl.append((self._weights[j], self._mu[j], self._sigma[j]))
        return SequenceProxy(retl)

    @classmethod
    def import_model(cls, model: Iterable[Tuple[float, float, float]]):
        model = list(model)
        new_instance = cls(n_components=len(model))
        for j in range(len(model)):
            (
                new_instance._weights[j],
                new_instance._mu[j],
                new_instance._sigma[j],
            ) = model[j]
        return new_instance

    def __repr__(self):
        model_desc = "Model ="
        for j in self.export():
            model_desc += f" + {j[0]}*N({j[1]}, {j[2]})"
        return model_desc

    def __str__(self):
        return repr(self)
