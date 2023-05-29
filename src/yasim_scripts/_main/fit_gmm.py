"""
Fit Gaussian Mixture by YUAN Ruihong
"""

import argparse
import re
from random import random

import numpy as np
import pandas as pd
from scipy.stats import norm

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import List
from yasim.helper.gmm import GaussianMixture1D
from yasim.helper.plot_utils import plot

try:
    import pytest

    _ = pytest.importorskip("pyarrow")
except ImportError:
    pytest = None
    try:
        import pyarrow
    except ImportError as e:
        raise UnmetDependenciesError("pyarrow") from e

_lh = get_logger()


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_scripts fit_gmm",
        description=__doc__.splitlines()[1]
    )
    parser.add_argument(
        'data',
        help="Dataset in Apache Parquet Format"
    )
    parser.add_argument(
        '-c', '--num_components', type=int, required=False, default=2
    )
    parser.add_argument(
        '--col_regex', type=str, required=False, default="SRR.*"
    )
    parser.add_argument(
        '--filter_zero', action="store_true"
    )
    parser.add_argument(
        '--scale', action="store_true"
    )
    parser.add_argument(
        '--num_iters', type=int, required=False, default=100
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    _lh.info("GMM: Processing data...")
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
        col_names = df.columns[list(map(lambda x: re.match(args.col_regex, x) is not None, df.columns))]
        df: pd.DataFrame = df[col_names]
        df = df.fillna(0)
        if args.scale:
            df = df.apply(
                lambda x: np.log10(x / max(x.max(), 1E-9) * 5000 + 1),
                axis="columns"
            )

        data = np.asarray(df.to_numpy()).ravel()
        if args.filter_zero:
            data = data[data != 0]
    print(describe(data))
    # plot(data)

    _lh.info("GMM: Fitting START")
    # auto_select_num_components(data, 2, 8, args.num_iters)
    model = GaussianMixture1D(
        n_components=args.num_components,
        n_iter=args.num_iters
    ).fit(data)
    _lh.info("GMM: Fitting FIN")

    plot(data, model)
    print(model.export())
