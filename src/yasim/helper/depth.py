"""
depth.py -- DGE Datastructure and Utils
"""
from typing import Dict

from scipy.stats import nbinom, uniform

from bioutils.datastructure.gene_view import GeneViewType
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer, get_reader
from yasim.helper.gmm import GaussianMixture1D

DepthType = Dict[str, int]
"""DGE type, is transcript_id -> coverage"""


def simulate_dge_uniform(
        gv: GeneViewType,
        mu: int
) -> DepthType:
    """
    Simulate DGE using a uniform distribution.

    :param gv: Input GeneViewType. Read only.
    :param mu: Mean of sequencing depth
    """
    transcript_ids = list(gv.iter_transcript_ids())
    depth = {}

    generated_data = list(map(int, uniform.rvs(scale=mu * 2, size=len(transcript_ids)) + 1))

    for i in range(len(transcript_ids)):
        depth[transcript_ids[i]] = generated_data[i]
    return depth


nb_generator = nbinom(n=2, p=0.016291702363523883)
nb_generator_mean = nb_generator.mean()


def simulate_dge_nb(
        gv: GeneViewType,
        mu: int
) -> DepthType:
    """
    Simulate DGE using a negative binomial distribution.
    """
    transcript_ids = gv.iter_transcript_ids()
    depth = {}

    for transcript_id in tqdm(iterable=transcript_ids, desc="Simulating..."):
        d = 0
        while not d > 0:
            d = int(nb_generator.rvs() / nb_generator_mean * mu)
        depth[transcript_id] = d
    return depth


gmm_model = GaussianMixture1D.import_model(
    (
        (0.7533708725243431, 13.50517241736218, 15.935138269740976),
        (0.24662912747565696, 448.39584422522495, 997.3253918151006)
    )
)


def simulate_dge_gmm(
        gv: GeneViewType,
        mu: int
) -> DepthType:
    """
    Simulate DGE using Gaussian mixture model.
    """
    transcript_ids = gv.iter_transcript_ids()
    depth = {}
    gmm_model.lintrans(mu / gmm_model.mean())
    for transcript_id in tqdm(iterable=transcript_ids, desc="Simulating..."):
        d = 0
        while not d > 0:
            d = int(gmm_model.rvs())
        depth[transcript_id] = d
    return depth


def write_dge(dge_data: DepthType, output_tsv: str):
    """
    Write DGE information to file
    """
    with get_writer(output_tsv) as writer:
        writer.write(f"TRANSCRIPT_ID\tDEPTH\n")
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
            retd[lkv[0]] = int(lkv[1])
    return retd
