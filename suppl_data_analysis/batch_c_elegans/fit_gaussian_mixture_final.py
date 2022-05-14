"""
Fit Gaussian Mixture by YUAN Ruihong
"""

import argparse
import sys
from typing import Optional, Union
from random import random

import numpy as np
from numpy.typing import ArrayLike
from numpy.random import choice
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import pandas as pd
from scipy.stats import norm
from sklearn.preprocessing import scale

matplotlib.use('qtagg')

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
        self._weights = init_weights or np.repeat(1/n_components, n_components)
        self._mu = init_mus or np.zeros(n_components)
        self._sigma = init_sigmas or np.ones(n_components)

    def fit(self, data: ArrayLike):
        N = len(data)

        # Use k-means to cluster points
        centers = choice(data, size=self._n_components)
        categories = np.repeat(0, N)
        bags = [[]] * self._n_components
        for _ in range(self._n_iter):
            converged = True
            for j in range(self._n_components):
                bags[j] = []

            for i in range(N):
                x = data[i]
                category = categories[i]
                least_dist = abs(x - centers[category])
                for j in range(self._n_components):
                    center = centers[j]
                    dist = abs(x - center)
                    if dist < least_dist:
                        converged = False
                        category = j  # assign to nearest center
                        least_dist = dist
                categories[i] = category
                bags[category].append(x)

            for j in range(self._n_components):
                bag = np.array(bags[j])
                bags[j] = bag
                centers[j] = np.mean(bag)

            if converged:
                break

        self._mu = centers
        for j in range(self._n_components):
            self._sigma[j] = np.std(bags[j])
            self._weights[j] = len(bags[j])
        self._weights /= N

        last_ll = np.NINF
        densities = np.zeros((self._n_components, N))
        for epoch in range(self._n_iter):
            # E step: compute ownership weights
            for j in range(self._n_components):
                model = norm(loc=self._mu[j], scale=self._sigma[j])
                densities[j, :] = model.pdf(data)

            a = np.transpose(np.multiply(np.transpose(densities), self._weights))
            b = (np.transpose(self._weights) @ densities)
            gamma = a / b

            # M step: compute means, variances and weights
            denom = np.sum(gamma, axis=1)  # row sum
            self._mu = (gamma @ data) / denom

            for j in range(self._n_components):
                for i in range(N):
                    gamma[j, i] *= (data[i] - self._mu[j])**2

            self._sigma = np.sqrt(np.sum(gamma, axis=1) / denom)
            self._weights = denom / N

            ll = self.log_likelihood(data)
            if 0 <= ll - last_ll < self._precision:
                break
            last_ll = ll

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

def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data')
    parser.add_argument('-c', '--numComponents', type=int, required=False, default=2)
    parser.add_argument('--columnIndex', type=int, required=False, default=-1)
    parser.add_argument('--header', type=int, required=False)
    parser.add_argument('--transform', nargs='+', required=False,
                        choices=('scale', 'log'), default=[])
    parser.add_argument('--numIters', type=int, required=False, default=100)

    args = parser.parse_args()
    if args.data == 'dummy':
        mu1, sigma1 = random() * 100, random() * 10
        mu2, sigma2 = random() * 100, random() * 10
        p = random()
        model1 = norm(mu1, sigma1)
        model2 = norm(mu2, sigma2)
        data = np.append(model1.rvs(size=int(p * 200)), model2.rvs(size=int((1-p) * 200)))
        print(f"Ground truth is {p}*N({mu1}, {sigma1}) + {1-p}*N({mu2}, {sigma2})")
    else:
        df = pd.read_table(args.data, header=args.header or 'infer')
        df = df[df.iloc[:, args.columnIndex] != 0]  # filter zeros
        data = np.asarray(df.iloc[:, args.columnIndex])

    for transform_method in args.transform:
        if transform_method == 'scale':
            scale(data, copy=False)
        elif transform_method == 'log':
            data = np.log(data + 1)

    model = GaussianMixture1D(n_components=args.numComponents,
                              n_iter=args.numIters).fit(data)
    print(f"Fitted to Gaussian mixture model ({args.numComponents} components)")
    model_desc = "Model ="
    for j in range(model._n_components):
        model_desc += f" + {model._weights[j]}*N({model._mu[j]}, {model._sigma[j]})"
    print(model_desc)

    fig, axes = plt.subplots(1, 2)
    axes[0].yaxis.set_major_formatter(PercentFormatter(xmax=1))

    x_max = np.max(data)
    n_bins = min(40, int(x_max))
    xticks = np.linspace(np.min(data), x_max, 100)

    axes[0].hist(data, bins=n_bins, density=True)
    axes[1].plot(xticks, model.pdf(xticks))
    plt.show()

if __name__ == '__main__':
    sys.exit(_main())
