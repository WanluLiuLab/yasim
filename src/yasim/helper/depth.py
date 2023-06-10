"""
depth.py -- GEP Datastructure and Utils

Here contains generators for Gene Expression Profile (GEP).
"""

__all__ = (
    "GenerationFailureException",
    "simulate_depth_gmm_v2",
    "simulate_gene_level_depth_gmm",
    "simulate_isoform_variance_inside_a_gene",
    "generate_depth_replicates_uniform"
)

import functools
import math
from random import uniform

import numpy as np

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewType
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import List
from yasim.helper.depth_io import DepthType
from yasim.helper.gmm import GaussianMixture1D

_lh = get_logger(__name__)


class GenerationFailureException(RuntimeError):
    """Raised while generation of data was failed"""
    ...


def simulate_gene_level_depth_gmm(
        gv: GeneViewType,
        mu: float,
        low_cutoff: float,
        high_cutoff_ratio: float
) -> DepthType:
    """
    Simulate DGE using Gaussian mixture model. Used in YASIM 3.0
    The two ``cutoff`` parameters are for data generated by the model.
    After cutoff, the data would be scaled to make its mean equal to ``mu``.

    The data would be re-generated if it is all-zero or contains ``np.nan``.

    :param gv: GeneView of the GTF.
    :param mu: Targeted mean.
    :param low_cutoff: Depth lower than this value would be preserved to this value.
    :param high_cutoff_ratio: Depth higher than mu * high_cutoff_ratio would be preserved to this value.
    :return: Generated abundance.
    :raise GenerationFailureException: If the data was re-generated 20 times.
    """
    _lh.info("GEN GENE DEPTH: Loading GMM model...")
    gmm_model = GaussianMixture1D.import_model([
        # (0.14427447754677641, 1.5527662658235803, 0.5349602445403809),
        # (0.06838373058223621, 3.6990568545476674, 4.440892098500626e-16),
        # (0.3139879527820852, 2.3714543839859354, 0.3438604637707616),
        # (0.47335383908890194, 3.126802802022012, 0.3396120976402885)
        (0.1427655350146609, 1.5746800285632137, 0.5388586661428069),
        (0.06837661567650871, 3.699056854547668, 0.0),
        (0.31792430333258503, 2.369675925231799, 0.345735276372599),
        (0.4709335459762455, 3.1287227827640636, 0.33812788057291315)
    ])
    n_gene_ids = gv.number_of_genes
    depth = {}
    for i in range(20):
        _lh.info("GEN GENE DEPTH: Attempt %d: GMM...", i)
        data = np.power(10, gmm_model.rvs(size=2 * n_gene_ids) - 1) - 1
        _lh.info("GEN GENE DEPTH: Attempt %d: 1/4: Scaling...", i)
        data = data / np.mean(data) * mu  # Scale to similar mean; should have a ~10% error
        _lh.info("GEN GENE DEPTH: Attempt %d: 2/4: Filtering...", i)
        data = data[functools.reduce(
            np.logical_and,
            (
                data > 0,
                np.isfinite(data)
            )
        )]  # Filter data
        if len(data) < n_gene_ids:
            _lh.warning(
                "GEN GENE DEPTH: filtered data length (%d) smaller than required (%d); would regenerate",
                len(data), n_gene_ids
            )
            continue

        data = data[: n_gene_ids]
        _lh.info("GEN GENE DEPTH: Attempt %d: 3/4 Scaling...", i)
        data = data / np.mean(data) * mu  # Rescale to real mean
        _lh.info("GEN GENE DEPTH: Attempt %d: 4/4 Filtering...", i)
        data[mu * high_cutoff_ratio < data] = mu * high_cutoff_ratio
        data[data < low_cutoff] = low_cutoff
        break
    else:
        raise GenerationFailureException()
    _lh.info("GEN GENE DEPTH: Final distribution: %s", describe(data))
    for i, gene_id in enumerate(gv.iter_gene_ids()):
        depth[gene_id] = data[i]
    return depth


def simulate_isoform_variance_inside_a_gene(
        n: int,
        mu: float,
        alpha: int,
        low_cutoff: float,
        high_cutoff_ratio: float
) -> List[float]:
    """
    Generate isoform variance inside a gene using Zipf's Distribution

    :param n: Number of isoforms
    :param mu: Mean of expression data
    :param alpha: Zipf's Coefficient
    :param low_cutoff: Depth lower than this value would be preserved to this value.
    :param high_cutoff_ratio: Depth higher than mu * high_cutoff_ratio would be preserved to this value.
    :return: Generated abundance.
    :raise GenerationFailureException: If the data was re-generated 20 times.
    """
    if n == 1:
        return [mu]
    for _ in range(20):
        if alpha == 1:
            data = np.array([(rn ** (-1)) for rn in range(1, n + 1)])
        else:
            data = np.array([(alpha - 1) * (rn ** (-alpha)) for rn in range(1, n + 1)])
        if np.sum(np.isnan(data)) == 0 and np.sum(data) != 0:
            break
        else:
            _lh.warning("GEN ISOFORM DEPTH: NAN found in data; would regenerate")
    else:
        raise GenerationFailureException()
    data = data / np.mean(data) * mu
    data[data < low_cutoff] = low_cutoff
    data[data > mu * high_cutoff_ratio] = mu * high_cutoff_ratio
    np.random.shuffle(data)
    return data


def simulate_depth_gmm_v2(
        gv: GeneViewType,
        mu: float
) -> DepthType:
    """
    Simulate Isoform-level GEP using Gaussian mixture model using version 2 algorithm.

    :param gv: GeneView of the GTF.
    :param mu: Targeted mean.
    :return: Generated abundance.
    """
    gmm_model = GaussianMixture1D.import_model(
        [
            (0.2757332390925274, 1.115903346554058, 0.28730736063701984),
            (0.7242667609074726, 3.588018272411806, 1.5658127268321427)
        ]
    )
    transcript_ids = gv.iter_transcript_ids()
    n_transcript_ids = gv.number_of_transcripts
    depth = {}
    gmm_model.lintrans((math.log(mu) - 1) / gmm_model.positive_mean())
    data = np.exp(gmm_model.rvs(size=2 * n_transcript_ids) - 1)
    data = data[13 * mu >= data]
    data = data[data >= 0][:n_transcript_ids]
    i = 0
    for transcript_id in tqdm(iterable=transcript_ids, desc="Simulating..."):
        if data[i] != 0:
            depth[transcript_id] = data[i]
        i += 1
    return depth


def generate_depth_replicates_uniform(
        input_depth: DepthType,
        _range: float = 0.001
) -> DepthType:
    """
    Generate replications of depth in uniform distribution. Would be O = I + uniform(-_range, +_range) * I

    This is a naive method that modulates technical replicate.

    :param input_depth: Input depth.
    :param _range: Degree of variation to be added.
    :return: Mutated depth information.
    """
    return {
        k: v + uniform(-_range, +_range) * v
        for k, v in input_depth.items()
    }
