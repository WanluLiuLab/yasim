"""
depth.py -- DGE Datastructure and Utils
"""

import random
from typing import Dict, Union, Iterable, List

import numpy as np
from scipy.stats import nbinom, uniform, gengamma
import numpy.typing as npt

from bioutils.datastructure.gene_view import GeneView
from commonutils import shell_utils
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer, get_reader

DepthType = Dict[str, int]
"""DGE type, is transcript_id -> coverage"""


def simulate_dge_uniform(
        gv: GeneView,
        mu: int
) -> DepthType:
    """
    Simulate DGE using a uniform distribution.
    """
    transcript_ids = list(gv.transcripts.keys())
    depth = {}

    generated_data = list(map(int, uniform.rvs(scale=mu*2, size=len(transcript_ids)) + 1))

    for i in range(len(transcript_ids)):
        depth[transcript_ids[i]] = generated_data[i]
    return depth

def simulate_dge_nb(
        gv: GeneView,
        max_depth: int,
        levels: int = 100
) -> DepthType:
    """
    Simulate DGE using a negative binomial distribution.
    """
    transcript_ids = gv.transcripts.keys()
    depth = {}

    for transcript_id in tqdm(iterable=transcript_ids, desc="Simulating..."):
        d = int(nbinom.rvs(1, max_depth)) * levels // max_depth * int(max_depth / levels) + 1
        depth[transcript_id] = d
    return depth


def write_dge(dge_data: DepthType, output_tsv:str):
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
