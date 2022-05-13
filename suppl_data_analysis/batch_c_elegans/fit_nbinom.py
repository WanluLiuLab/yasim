"""
Fit negative binomial distribution by YUAN Ruihong
"""
import argparse
import sys
from random import random, randint
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib.ticker import PercentFormatter
from tqdm import tqdm

EPS = 0.01


# Fit under the assumption that maximum likelihood function
# is convex (increase, then decrease)
def nbinom_mle(data, n_epochs: int):
    n = 2
    N = len(data)
    b = np.sum(data)
    mll = np.NINF
    best_mll = np.NINF
    best_n: int = 0
    best_p: float = 0.0
    diff: float = 0.0

    learning_rates = np.linspace(1 / EPS, 1, n_epochs)
    cache = set()  # avoid oscillation
    first_part = np.sum(np.log(data + 1))  # for n=2
    progress_bar = tqdm(learning_rates)
    for r in progress_bar:
        a = n * N
        p = a / (a + b)

        last_mll = mll
        mll = first_part + a * np.log(p) + b * np.log(1 - p)
        diff = mll - last_mll
        if n not in cache and mll - best_mll > 0:
            best_mll = mll
            best_n = n
            best_p = p
            if diff < EPS:
                progress_bar.update(len(learning_rates))
                progress_bar.close()
                print("Converged")
                break

        if n in cache:
            k = int(np.sign(diff))  # enforce move by 1
            if (n + k) in cache:  # oscillation
                progress_bar.update(len(learning_rates))
                progress_bar.close()
                print("Broke out of oscillation")
                break
        else:
            k = int(np.clip(diff * r, -5, 5))
            if k == 0:
                k = int(np.sign(diff))
            if (n + k) < 2:
                k = 2 - n

        if k > 0:
            for j in range(n, n + k):
                for i in range(0, N):
                    first_part += np.log(data[i] / j + 1)
        else:  # k < 0
            for j in range(n + k, n):
                for i in range(0, N):
                    first_part -= np.log(data[i] / j + 1)
        cache.add(n)
        n += k

    return best_mll, best_n, best_p


def mll2aic(mll: float, n_params: int, n_data: Optional[int] = None) -> float:
    result = -2 * mll + 2 * n_params
    if n_data is not None:
        result += 2 * n_params * (n_params + 1) / (n_data - n_params - 1)
    return result


def mll2bic(mll: float, n_params: int, n_data: int) -> float:
    return -2 * mll + n_params * np.log(n_data)


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--depth", required=False,
                        help="Path to training dataset. If not specified, uses a random negative binomial model")
    parser.add_argument("-e", "--numEpochs", type=int,
                        required=False, default=100,
                        help="Number of training epochs (default 100)")
    parser.add_argument("-c", "--columnNo", type=int,
                        required=False, default=-1,
                        help="0-based index of column of expression level")
    parser.add_argument("--dataDisplayCutoff", type=int,
                        required=False, default=None,
                        help="Upper cut-off of data (default none)")
    parser.add_argument("--header", type=int,
                        required=False, default=None,
                        help="Header specification (default auto-infer)")
    parser.add_argument("--transform", required=False,
                        choices=('log',),
                        help="Transformation to the data (default identity)")

    args = parser.parse_args()

    if args.depth is None:
        n, p = randint(2, 100), random()
        print(f"Ground truth is n={n} and p={p}")
        data = scipy.stats.nbinom.rvs(n, p, size=20000)
    else:
        depth_df = pd.read_table(args.depth, header=args.header or "infer")
        if args.columnNo < 0:
            args.columnNo = depth_df.shape[1] - 1  # last column
        depth_df = depth_df[depth_df.iloc[:, args.columnNo] != 0]
        data = np.array(depth_df.iloc[:, args.columnNo])
    if args.dataDisplayCutoff is None:
        args.dataDisplayCutoff = np.max(data)
    if args.transform:
        if args.transform == 'log':
            data = np.log(data + 1)
            args.dataDisplayCutoff = int(np.log(args.dataDisplayCutoff + 1))
        else:
            raise ValueError("unrecognized transformation")

    mll, n, p = nbinom_mle(data, args.numEpochs)
    print(
        f"Fitted to nbinom(n={n}, p={p}) with LL={mll}, AICc={mll2aic(mll, 2, len(data))}, BIC={mll2bic(mll, 2, len(data))}")

    fitted = scipy.stats.nbinom(n, p)
    print(f"MEAN={fitted.mean()}")
    conf_a, conf_b = fitted.interval(0.95)
    n_outlier = len(data[data < conf_a]) + len(data[data > conf_b])
    print(f"Number of outliers (at p<0.05): {n_outlier}/{len(data)}" + "({:.2%})".format(n_outlier / len(data)))
    n_bins = min(args.dataDisplayCutoff, 40)
    fig, axes = plt.subplots(1, 2)
    rounded_step = int(args.dataDisplayCutoff / n_bins)
    xticks = np.arange(0, rounded_step * n_bins, rounded_step)
    axes[0].yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axes[1].yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axes[0].hist(data[data < args.dataDisplayCutoff], bins=n_bins, density=True)
    axes[1].vlines(xticks, 0, fitted.pmf(xticks))
    plt.show()
    plt.savefig("gep_fit_nb.png")


# n, p = 5, 0.5
# ground_truth = scipy.stats.nbinom(n, p)
# sample = ground_truth.rvs(size=1000)

if __name__ == "__main__":
    sys.exit(_main())
