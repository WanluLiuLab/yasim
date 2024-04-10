"""
generate_isoform_depth.py -- Generate Isoform-Level Depth using YASIM V3 API.

.. versionadded:: 3.1.5
"""

__all__ = ("main", "create_parser")

import argparse

import numpy as np

from labw_utils.bioutils.datastructure.gene_tree import DiploidGeneTree
from labw_utils.bioutils.datastructure.gv.gene import DumbGene
from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List
from yasim.helper import depth, depth_io, isoform_depth
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
    gene_level_depth = depth_io.read_depth(args.depth)
    transcript_level_depth = isoform_depth.generate_isoform_depth(
        gv.gene_id_transcript_ids_map(),
        gene_level_depth,
        args.alpha,
        args.low_cutoff,
        args.high_cutoff_ratio,
    )
    depth_io.write_depth(transcript_level_depth, args.out, "TRANSCRIPT_ID")
