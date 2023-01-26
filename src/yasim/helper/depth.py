"""
depth.py -- DGE Datastructure and Utils
"""
import math
from typing import Dict

import numpy as np
from labw_utils.bioutils.datastructure.gene_view import GeneViewType
from labw_utils.commonutils import shell_utils
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer, get_reader

from yasim.helper.gmm import GaussianMixture1D

DepthType = Dict[str, float]
"""DGE type, is transcript_id -> coverage"""


def simulate_gene_level_dge_gmm(
        gv: GeneViewType,
        mu: float
):
    """
    Simulate DGE using Gaussian mixture model. Used in YASIM 3.0
    """
    gmm_model = GaussianMixture1D.import_model(
        [
            (0.12827805183096117, 1.3597462971461431, 0.20564573367723052),
            (0.15548306830287628, 0.8132604801292427, 0.1932997409289587),
            (0.1182911269313193, 1.8663843368510689, 0.23483302069369835),
            (0.033586632398567316, 2.617052700717236, 0.5027565062645347),
            (0.16480393079998523, 0.3188344745000003, 0.12415694780744323),
            (0.3995571897362907, 1.0000000000000086e-06, 8.682087709356578e-21)
        ]
    )
    n_gene_ids = gv.number_of_genes
    depth = {}
    gmm_model.lintrans((math.log(mu) - 1) / gmm_model.positive_mean())
    data = np.power(10, gmm_model.rvs(size=2 * n_gene_ids) - -1E-6) - 1
    data = data[data >= 0][:n_gene_ids]
    for i, gene_id in enumerate(tqdm(iterable=gv.iter_genes(), desc="Simulating...")):
        depth[gene_id] = data[i]
    return depth


def simulate_dge_gmm(
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


def write_dge(
        dge_data: DepthType,
        output_tsv: str,
        feature_name: str = "TRANSCRIPT_ID"
):
    """
    Write DGE information to file
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
