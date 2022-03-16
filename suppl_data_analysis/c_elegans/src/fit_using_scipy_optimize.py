from math import e

import numpy as np
import pandas as pd
import scipy.optimize as opt


def d_flux_sim(x, k, Y_max, a, b):
    return Y_max * x ** k * e** (x / a * (x / b) ** 2)


if __name__ == "__main__":
    all_data = pd.read_table("all_data.tsv")
    nanopore_coverage = all_data.loc[:, 'NANOPORE_AVG_DEPTH'].to_numpy(dtype=float)
    nanopore_coverage = sorted(nanopore_coverage[np.where(nanopore_coverage > 0)], reverse=True)
    maxnc = max(nanopore_coverage)
    lnc = len(nanopore_coverage)
    f = opt.curve_fit(d_flux_sim, range(lnc), nanopore_coverage, bounds=([-0.9, maxnc-0.01, 1/10*lnc, 1/10*lnc], [-0.5, maxnc, 10*lnc, 10*lnc]))

#     upper=c(-0.5, 10*ldata, 10*ldata),
#     lower=c(-0.9, 1/10*ldata, 1/10*ldata)
