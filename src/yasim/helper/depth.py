"""
depth.py -- DGE Datastructure and Utils
"""
import math
from random import uniform
from typing import Dict, List

import numpy as np
from labw_utils.bioutils.datastructure.gene_view import GeneViewType
from labw_utils.commonutils import shell_utils
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer, get_reader

from yasim.helper.gmm import GaussianMixture1D

DepthType = Dict[str, float]
"""DGE type, is transcript_id -> coverage"""


def simulate_gene_level_depth_gmm(
        gv: GeneViewType,
        mu: float
):
    """
    Simulate DGE using Gaussian mixture model. Used in YASIM 3.0
    """
    gmm_model = GaussianMixture1D.import_model([
        (0.3043757984608804, 2.3107772552803634, 0.3956119888112459),
        (0.3106118627088962, 1.3815767710834788, 0.19646588205588317),
        (0.2905529799446133, 1.000000000000003, 3.108624468950436e-15),
        (0.09445935888561038, 3.04499166208647, 0.6200436778100588)
    ])
    n_gene_ids = gv.number_of_genes
    depth = {}
    data = np.power(10, gmm_model.rvs(size=n_gene_ids) - 1) - 1
    data[data < 0.001] = 0
    data = data / np.mean(data) * mu
    for i, gene_id in enumerate(tqdm(iterable=gv.iter_gene_ids(), desc="Simulating...")):
        depth[gene_id] = data[i]
    return depth


def simulate_isoform_variance_inside_a_gene(
        n: int,
        mu: float,
        alpha: int = 10
) -> List[float]:
    """
    Generate isoform variance inside a gene using Zipf's Distribution

    :param n: Number of isoforms
    :param mu: Mean of expression data
    :param alpha: Zipf's Coefficient
    :return: Generated abundance
    """
    if n == 1:
        return [mu]
    generated_abundance = np.array([(alpha - 1) * (rn ** (-alpha)) for rn in range(1, n + 1)])
    generated_abundance[generated_abundance < 0.01] = 0
    generated_abundance = generated_abundance / np.mean(generated_abundance) * mu
    np.random.shuffle(generated_abundance)
    return generated_abundance


def simulate_depth_gmm(
        gv: GeneViewType,
        mu: float
) -> DepthType:
    """
    Simulate DGE using Gaussian mixture model.
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


def write_depth(
        dge_data: DepthType,
        output_tsv: str,
        feature_name: str = "TRANSCRIPT_ID"
):
    """
    Write Depth information to file
    """
    with get_writer(output_tsv) as writer:
        writer.write(f"{feature_name}\tDEPTH\n")
        for transcript_id, d in dge_data.items():
            writer.write(f"{transcript_id}\t{d}\n")


def read_depth(input_tsv: str) -> DepthType:
    retd = {}
    total = shell_utils.wc_l(input_tsv)
    with get_reader(input_tsv) as reader:
        reader.readline()  # Skip line 1
        for line in tqdm(iterable=reader.readlines(), desc="Reading depth file...", total=total - 1):
            line = line.strip()
            lkv = line.split("\t")
            retd[lkv[0]] = float(lkv[1])
    return retd


def generate_depth_replicates_uniform(
        input_depth: DepthType,
        _range: float = 0.001
) -> DepthType:
    return {
        k: v + uniform(-_range, +_range) * v
        for k, v in input_depth.items()
    }
