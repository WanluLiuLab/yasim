"""
Fit Gaussian Mixture by YUAN Ruihong
"""

import argparse
import sys
from random import random

import numpy as np
import pandas as pd
from scipy.stats import norm

from yasim.helper.gmm import GaussianMixture1D
from yasim.helper.plot_utils import plot


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data')
    parser.add_argument('-c', '--numComponents', type=int, required=False, default=2)
    parser.add_argument('--transform', nargs='+', required=False,
                        choices=('log',), default=[])
    parser.add_argument('--numIters', type=int, required=False, default=100)

    args = parser.parse_args()
    if args.data == 'dummy':
        mu1, sigma1 = random() * 100, random() * 10
        mu2, sigma2 = random() * 100, random() * 10
        p = random()
        model1 = norm(mu1, sigma1)
        model2 = norm(mu2, sigma2)
        data = np.append(model1.rvs(size=int(p * 200)), model2.rvs(size=int((1 - p) * 200)))
        print(f"Ground truth is {p}*N({mu1}, {sigma1}) + {1 - p}*N({mu2}, {sigma2})")
    else:
        df = pd.read_parquet(args.data)
        data = np.asarray(df.iloc[:, 1:]).flatten()

    for transform_method in args.transform:
        if transform_method == 'log':
            data = np.log10(data + 1) + 1

    model = GaussianMixture1D(
        n_components=args.numComponents,
        n_iter=args.numIters
    ).fit(data)
    print(f"Fitted to Gaussian mixture model ({args.numComponents} components)")
    export_str = '\n'.join(map(str,(model.export())))
    plot(data, model, f"GMM {args.numComponents} components\n{export_str}\n")
    print(list(model.export()))


if __name__ == '__main__':
    sys.exit(_main())
