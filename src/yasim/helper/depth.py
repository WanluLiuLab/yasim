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

import math
from random import uniform

import numpy as np

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewType
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
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
    :param low_cutoff: Depth lower than this value would be 0.
    :param high_cutoff_ratio: Depth higher than mu * high_cutoff_ratio would be 0.
    :return: Generated abundance.
    :raise GenerationFailureException: If the data was re-generated 20 times.
    """
    gmm_model = GaussianMixture1D.import_model([
        (0.3043757984608804, 2.3107772552803634, 0.3956119888112459),
        (0.3106118627088962, 1.3815767710834788, 0.19646588205588317),
        (0.2905529799446133, 1.000000000000003, 3.108624468950436e-15),
        (0.09445935888561038, 3.04499166208647, 0.6200436778100588)
    ])
    n_gene_ids = gv.number_of_genes
    depth = {}
    for _ in range(20):
        data = np.power(10, gmm_model.rvs(size=n_gene_ids) - 1) - 1
        data[data < low_cutoff] = 0
        data[data > mu * high_cutoff_ratio] = 0
        data = data / np.mean(data) * mu
        if np.sum(np.isnan(data)) == 0 and np.sum(data) != 0:
            break
        else:
            _lh.warning("NAN/all zero found in data; would regenerate")
    else:
        raise GenerationFailureException()
    for i, gene_id in enumerate(tqdm(iterable=list(gv.iter_gene_ids()), desc="Simulating...")):
        depth[gene_id] = data[i]
    return depth


def simulate_isoform_variance_inside_a_gene(
        n: int,
        mu: float,
        alpha: int,
        low_cutoff: float
) -> List[float]:
    """
    Generate isoform variance inside a gene using Zipf's Distribution

    :param n: Number of isoforms
    :param mu: Mean of expression data
    :param alpha: Zipf's Coefficient
    :param low_cutoff: Depth lower than this value would be 0.
    :return: Generated abundance.
    :raise GenerationFailureException: If the data was re-generated 20 times.
    """
    if n == 1:
        return [mu]
    for _ in range(20):
        if alpha == 1:
            generated_abundance = np.array([(rn ** (-1)) for rn in range(1, n + 1)])
        else:
            generated_abundance = np.array([(alpha - 1) * (rn ** (-alpha)) for rn in range(1, n + 1)])
        if np.sum(np.isnan(generated_abundance)) == 0 and np.sum(generated_abundance) != 0:
            break
        else:
            _lh.warning("NAN found in data; would regenerate")
    else:
        raise GenerationFailureException()
    generated_abundance[generated_abundance < low_cutoff] = 0
    generated_abundance = generated_abundance / np.mean(generated_abundance) * mu
    np.random.shuffle(generated_abundance)
    return generated_abundance


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
