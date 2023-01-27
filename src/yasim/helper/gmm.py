"""Gaussian Mixture Model bu YUAN Ruihong"""
import multiprocessing
from random import choices
from typing import Optional, Union, Iterable, Tuple

import numpy as np
from joblib import Parallel, delayed
from numpy.random import choice
from numpy.typing import ArrayLike
from scipy.integrate import quad
from scipy.stats import norm


class GaussianMixture1D:
    _n_components: int
    _n_iter: int
    _precision: float
    _weights: ArrayLike
    _mu: ArrayLike
    _sigma: ArrayLike

    def __init__(
            self, n_components: int = 1,
            n_iter: int = 100,
            init_weights: Optional[ArrayLike] = None,
            init_mus: Optional[ArrayLike] = None,
            init_sigmas: Optional[ArrayLike] = None,
            precision: float = 1E-4,
    ) -> None:
        self._n_components = n_components
        self._n_iter = n_iter
        self._precision = precision
        self._weights = init_weights or np.repeat(1 / n_components, n_components)
        self._mu = init_mus or np.zeros(n_components)
        self._sigma = init_sigmas or np.ones(n_components)

    def fit(self, data: ArrayLike, show_progress_bar: bool = True):
        N = len(data)

        # Use k-means to cluster points
        centers = np.random.choice(data, size=self._n_components)
        categories = np.repeat(0, N)
        bags = [[]] * self._n_components
        for _ in range(self._n_iter):
            converged = True
            for j in range(self._n_components):
                bags[j].clear()

            for i in range(N):
                x = data[i]
                category = categories[i]
                least_dist = abs(x - centers[category])
                for j in range(0, self._n_components):
                    center = centers[j]
                    dist = abs(x - center)
                    if dist < least_dist:
                        converged = False
                        category = j  # assign to nearest center
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
        self._weights /= N
        return self

    def log_likelihood(self, data):
        return np.sum(self.logpdf(data))

    def pdf(self, x: Union[int, float, ArrayLike]) -> float:
        if isinstance(x, (int, float)):
            result = 0.0
        else:
            result = np.zeros(len(x))

        for j in range(self._n_components):
            model = norm(loc=self._mu[j], scale=self._sigma[j])
            result += self._weights[j] * model.pdf(x)
        return result

    def logpdf(self, x: Union[int, float, ArrayLike]) -> float:
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
            for j in range(self._n_components):
                result += self._weights[j] * models[j].pdf(x)
            return result * x

        return quad(f, start, np.inf)[0]

    def mean(self) -> float:
        return np.average(self._mu, weights=self._weights)

    def rvs(self, size: int = 1) -> ArrayLike:
        indices = list(range(self._n_components))

        def _rvs():
            k = choices(indices, weights=self._weights)
            model = norm(self._mu[k], self._sigma[k])
            return model.rvs()

        result = np.array(Parallel(n_jobs=multiprocessing.cpu_count())(delayed(_rvs)() for _ in range(size)))
        return result

    def export(self) -> Iterable[Tuple[float, float, float]]:
        for j in range(self._n_components):
            yield self._weights[j], self._mu[j], self._sigma[j]

    @classmethod
    def import_model(cls, model: Iterable[Tuple[float, float, float]]):
        model = list(model)
        new_instance = cls(n_components=len(model))
        for j in range(len(model)):
            new_instance._weights[j], new_instance._mu[j], new_instance._sigma[j] = model[j]
        return new_instance

    def __repr__(self):
        model_desc = "Model ="
        for j in self.export():
            model_desc += f" + {j[0]}*N({j[1]}, {j[2]})"
        return model_desc

    def __str__(self):
        return repr(self)
