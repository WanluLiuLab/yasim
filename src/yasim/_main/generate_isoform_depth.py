"""
generate_isoform_depth.py -- Generate Isoform-Level Depth using YASIM V3 API.

.. versionadded:: 3.1.5
"""

__all__ = ("main", "create_parser")

import argparse

import numpy as np

import yasim.helper.depth_io
from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import List
from yasim.helper import depth
from yasim.helper.frontend import patch_frontend_argument_parser

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim generate_isoform_depth",
        description=__doc__.splitlines()[1],
    )
    parser = patch_frontend_argument_parser(parser, "-g")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to output Isoform-Level Depth TSV. Can be compressed.",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-d",
        "--depth",
        required=True,
        help="Path to input Gene-Level Depth TSV. Can be compressed.",
        nargs="?",
        type=str,
        action="store",
    )
    parser = patch_frontend_argument_parser(parser, "--low_cutoff")
    parser = patch_frontend_argument_parser(parser, "--high_cutoff_ratio")
    parser.add_argument(
        "--alpha",
        required=False,
        help="Zipf's Coefficient, larger for larger differences",
        nargs="?",
        type=int,
        action="store",
        default=depth.DEFAULT_ALPHA,
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    gv = DiploidGeneTree.from_gtf_file(args.gtf, gene_implementation=DumbGene)
    gene_level_depth = yasim.helper.depth_io.read_depth(args.depth)
    transcript_level_depth = {}
    for gene in gv.gene_values:
        if gene.gene_id not in gene_level_depth:
            _lh.warning(
                "GEN ISOFORM DEPTH: Gene %s defined in GTF but not gene-level depth",
                gene.gene_id,
            )
        if gene_level_depth[gene.gene_id] == 0:
            for transcript in gene.transcript_values:
                transcript_level_depth[transcript.transcript_id] = 0
            continue
        try:
            this_transcript_level_depth = depth.simulate_isoform_variance_inside_a_gene(
                n=gene.number_of_transcripts,
                mu=gene_level_depth[gene.gene_id],
                low_cutoff=args.low_cutoff,
                alpha=args.alpha,
                high_cutoff_ratio=args.high_cutoff_ratio,
            )
        except depth.GenerationFailureException:
            _lh.error(
                "GEN ISOFORM DEPTH: Generation failed for gene %s -- SKIPPED",
                gene.gene_id,
            )
            continue
        for i, transcript in enumerate(gene.transcript_values):
            transcript_level_depth[transcript.transcript_id] = this_transcript_level_depth[i]
    _lh.info(
        "GEN ISOFORM DEPTH: Generation of isoform-level depth: Final distribution: %s",
        describe(np.array(list(transcript_level_depth.values()))),
    )
    yasim.helper.depth_io.write_depth(transcript_level_depth, args.out, "TRANSCRIPT_ID")
